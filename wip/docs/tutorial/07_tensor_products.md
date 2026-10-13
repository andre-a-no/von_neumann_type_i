**English** · [Español](es/07_tensor_products.md) · [Français](fr/07_tensor_products.md) · [Deutsch](de/07_tensor_products.md) · [中文](zh/07_tensor_products.md) · [日本語](ja/07_tensor_products.md) · [Русский](ru/07_tensor_products.md)

# 7. Tensor products and entanglement

Two quantum systems together are described by the tensor product of their algebras. With
sectors there are two products: the plain one, in which each part keeps its own charge, and the
*fused* one, in which only the total charge is conserved. This chapter proves the structure of
both, defines the partial trace as the dual of an embedding, and studies entanglement.

**You will learn** how sectors multiply and how the library orders them, why the partial trace is
a channel and, normalised, a conditional expectation of finite index, what the reduced state of a charge-conserving state is, how entanglement splits into a
number part and a configurational part, and the average purity and entropy of random states.

Throughout, $\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)$ acts on
$H_1=\bigoplus_c\mathbb C^{k_c}$ with minimal central projections $P_c$, and
$\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)$ acts on $H_2=\bigoplus_d\mathbb C^{l_d}$
with minimal central projections $Q_d$. $\operatorname{Tr}$ is the trace of the space on which an
operator acts; on $\mathcal M_i$ it is the blunt trace of Chapter 3. For $X$ on
$\mathbb C^k\otimes\mathbb C^l$ we write $X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle$,
and we identify $\mathbb C^k\otimes\mathbb C^l$ with $\mathbb C^{kl}$ by $e_a\otimes e_b\mapsto e_{(a-1)l+b}$
(the Kronecker convention), so that $(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}$.

## 7.1 The tensor product of two algebras

For $*$-algebras $\mathcal M_1$ on $H_1$ and $\mathcal M_2$ on $H_2$, the **tensor product**
$\mathcal M_1\otimes\mathcal M_2$ is the linear span of the operators $A\otimes B$ on
$H_1\otimes H_2$, $A\in\mathcal M_1$, $B\in\mathcal M_2$. Since
$(A\otimes B)(A'\otimes B')=AA'\otimes BB'$, $(A\otimes B)^*=A^*\otimes B^*$ and $1\otimes1=1$, it
is a $*$-algebra, hence a von Neumann algebra (Theorem 2.2).

**Theorem 7.1.** (a) $H_1\otimes H_2=\bigoplus_{c,d}\mathbb C^{k_c}\otimes\mathbb C^{l_d}$, an
orthogonal sum of the ranges of the projections $P_c\otimes Q_d$, and
$$
\mathcal M_1\otimes\mathcal M_2=\bigoplus_{c,d}M_{k_c}\otimes M_{l_d}\cong\bigoplus_{c,d}M_{k_cl_d}(\mathbb C):
$$
it consists of all operators $X=\sum_{c,d}(P_c\otimes Q_d)X(P_c\otimes Q_d)$, with arbitrary blocks
on the $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$. (b) Its minimal central projections are the $C_1C_2$
projections $P_c\otimes Q_d$, and $\dim(\mathcal M_1\otimes\mathcal M_2)=\dim\mathcal M_1\dim\mathcal M_2$.
(c) For traces $\operatorname{tr}_w$ on $\mathcal M_1$ and $\operatorname{tr}_v$ on $\mathcal M_2$ (Theorem 3.1),
$\operatorname{tr}_w(A)\operatorname{tr}_v(B)=\operatorname{tr}_{w\otimes v}(A\otimes B)$ with
$(w\otimes v)_{cd}=w_cv_d$. In particular `Tr_blunt`, `Tr_norm` and `tau_vN` of the product are the
products of the same traces of the factors.

*Proof.* (a) The $P_c\otimes Q_d$ are mutually orthogonal projections with sum $1\otimes1$, and the
range of $P_c\otimes Q_d$ is $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$. As $P_c$, $Q_d$ are central,
$(P_c\otimes Q_d)(A\otimes B)=AP_c\otimes BQ_d=(A\otimes B)(P_c\otimes Q_d)$, so every element of
$\mathcal M_1\otimes\mathcal M_2$ is block diagonal. Conversely, with the matrix units
$E^c_{aa'}\in\mathcal M_1$ and $F^d_{bb'}\in\mathcal M_2$ of the blocks $c$ and $d$,
$E^c_{aa'}\otimes F^d_{bb'}=|e_a\otimes e_b\rangle\langle e_{a'}\otimes e_{b'}|$ inside
$\mathbb C^{k_c}\otimes\mathbb C^{l_d}$; these are the matrix units of $M_{k_cl_d}(\mathbb C)$, so every
operator supported in the block $(c,d)$ lies in $\mathcal M_1\otimes\mathcal M_2$. (b) The centre of a
direct sum of full matrix algebras is spanned by the projections onto the summands (Chapter 2,
Exercise 1), and $\sum_{c,d}k_c^2l_d^2=\sum_ck_c^2\sum_dl_d^2$. (c)
$\operatorname{Tr}(X\otimes Y)=\sum_{a,b}X_{aa}Y_{bb}=\operatorname{Tr}X\operatorname{Tr}Y$, so
$\operatorname{tr}_{w\otimes v}(A\otimes B)=\sum_{c,d}w_cv_d\operatorname{Tr}A_c\operatorname{Tr}B_d=\operatorname{tr}_w(A)\operatorname{tr}_v(B)$.
The weights $1$, $1/k_c$, $1/(C_1k_c)$ of the three traces (Section 3.1) multiply to $1$,
$1/(k_cl_d)$, $1/(C_1C_2k_cl_d)$, the weights of the same traces on the product. $\square$

`tensor_product(alg1, alg2)` builds this algebra. Sector $(c,d)$ has index $s=cC_2+d$ (indices
from $0$ in the code) and charge $(q_c,r_d)$, and inside it $e_a\otimes e_b$ has index $a\,l_d+b$.
`kron(A, B, alg12)` computes $A\otimes B$. The sector order differs from the Kronecker order of
$H_1\otimes H_2$ by a permutation of the basis, which `kron_rows` below builds explicitly; with it
we compare with Kronecker products of full matrices.

```python
import math
import numpy as np
import torch
from torch_vn_algebra import (TypeIAlgebra, DensityMatrix, states, composite,
                              tensor_product, fused_tensor_product, kron, partial_trace)

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def full(X):
    """An operator of (+)_c M_{k_c} as a block-diagonal matrix on (+)_c C^{k_c}, shape (batch, D, D)."""
    return torch.stack([torch.block_diag(*[x[c, :k, :k] for c, k in enumerate(X.algebra.k_factors)])
                        for x in X.matrix])

def kron_rows(alg12):
    """For every sector of a product algebra, the positions of its basis vectors in H1 (x) H2 (Kronecker order)."""
    k1, k2 = alg12.tensor_factors[0].k_factors, alg12.tensor_factors[1].k_factors
    return [torch.tensor([(sum(k1[:c]) + a) * sum(k2) + sum(k2[:d]) + b
                          for c, d, off in pairs for a in range(k1[c]) for b in range(k2[d])])
            for pairs in alg12.product_blocks]              # pair blocks (c, d, offset), by offset

def full12(X):
    """An operator of a (plain or fused) product algebra as a matrix on H1 (x) H2 in Kronecker order."""
    f1, f2 = X.algebra.tensor_factors
    N = sum(f1.k_factors) * sum(f2.k_factors)
    out = torch.zeros(X.matrix.shape[0], N, N, dtype=X.matrix.dtype)
    for s, r in enumerate(kron_rows(X.algebra)):
        out[:, r[:, None], r[None, :]] = X.matrix[:, s, :len(r), :len(r)]
    return out

def random_operator(alg):
    return alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])

alg1 = TypeIAlgebra([1, 2], [1, 2], precision='double')            # C (+) M_2
alg2 = TypeIAlgebra([2, 3], [2, 3], precision='double')            # M_2 (+) M_3
alg12 = tensor_product(alg1, alg2)
print(alg12.k_factors, alg12.charges, sum(k * k for k in alg12.k_factors))     # 4 sectors, dimension 65
assert sum(k * k for k in alg12.k_factors) == 5 * 13                           # Theorem 7.1(b)

A, A2, B, B2 = random_operator(alg1), random_operator(alg1), random_operator(alg2), random_operator(alg2)
AB = kron(A, B, alg12)
assert close(full12(AB)[0], torch.kron(full(A)[0], full(B)[0]))               # Theorem 7.1(a)
assert close(AB @ kron(A2, B2, alg12), kron(A @ A2, B @ B2, alg12))
for tr in ('Tr_blunt', 'Tr_norm', 'tau_vN'):                                   # Theorem 7.1(c)
    assert (getattr(AB, tr)() - getattr(A, tr)() * getattr(B, tr)()).abs().max() < 1e-10
```

> **Physical meaning.** If each of two systems obeys its own superselection rule — two boxes
> that cannot exchange particles — the joint system has a sector for every *pair* of charges, and
> no observable connects different pairs. A joint sector counts with the product of the weights.

## 7.2 The partial trace

**Definition.** Let $\rho\in\mathcal M_1\otimes\mathcal M_2$. The **partial trace**
$\operatorname{Tr}_2\rho\in\mathcal M_1$ is the element with
$$
\operatorname{Tr}\big((\operatorname{Tr}_2\rho)\,A\big)=\operatorname{Tr}\big(\rho\,(A\otimes1)\big)\qquad\text{for all }A\in\mathcal M_1 .
$$
It exists and is unique by Theorem 3.2, applied to the linear functional
$A\mapsto\operatorname{Tr}(\rho(A\otimes1))$ and the faithful trace $\operatorname{Tr}$ of $\mathcal M_1$: the
partial trace is the dual, with respect to the blunt traces, of the embedding $A\mapsto A\otimes1$.
$\operatorname{Tr}_1$ is defined in the same way. For the full algebras $\mathcal M_i=B(H_i)$ this is the
**Hilbert-space partial trace** $\operatorname{Tr}_{H_2}:B(H_1\otimes H_2)\to B(H_1)$.

**Proposition 7.2.** (a) For $\rho=\bigoplus_{c,d}\rho_{cd}$,
$$
(\operatorname{Tr}_2\rho)_c=\sum_d\operatorname{tr}_{l_d}\rho_{cd},\qquad
(\operatorname{tr}_lX)_{aa'}=\sum_{b=1}^lX_{(a,b),(a',b)}\quad(X\text{ on }\mathbb C^k\otimes\mathbb C^l),
$$
and $\operatorname{Tr}_2\rho=\operatorname{Tr}_{H_2}\rho$, where $(\operatorname{Tr}_{H_2}X)_{ii'}=\sum_jX_{(i,j),(i',j)}$.
(b) $\operatorname{Tr}_2(A\otimes B)=(\operatorname{Tr}B)\,A$.
(c) $\operatorname{Tr}_2\big((A\otimes1)X(A'\otimes1)\big)=A\,(\operatorname{Tr}_2X)\,A'$ for $A,A'\in\mathcal M_1$.
(d) $\operatorname{Tr}(\operatorname{Tr}_2X)=\operatorname{Tr}X$, and $\operatorname{Tr}_2$ maps states to states.

*Proof.* (a) With $A=\bigoplus_cA_c$, $A\otimes1=\bigoplus_{c,d}A_c\otimes1_{l_d}$ and
$$
\operatorname{Tr}\big(\rho(A\otimes1)\big)=\sum_{c,d}\sum_{a,a',b}(\rho_{cd})_{(a,b),(a',b)}(A_c)_{a'a}
=\sum_c\operatorname{Tr}\Big(\Big(\sum_d\operatorname{tr}_{l_d}\rho_{cd}\Big)A_c\Big),
$$
so the element on the right satisfies the definition and is $\operatorname{Tr}_2\rho$. The same
computation for $B(H_i)$ gives the formula for $\operatorname{Tr}_{H_2}$; for block-diagonal $\rho$ the
entries $\rho_{(i,j),(i',j)}$ with $i,i'$ in different blocks vanish, so $\operatorname{Tr}_{H_2}\rho$ has
the blocks of (a). (b) $\operatorname{tr}_l(X\otimes Y)=(\operatorname{Tr}Y)X$. (c) For $B\in\mathcal M_1$,
$\operatorname{Tr}\big(\operatorname{Tr}_2((A\otimes1)X(A'\otimes1))B\big)=\operatorname{Tr}\big(X(A'BA\otimes1)\big)
=\operatorname{Tr}\big(A(\operatorname{Tr}_2X)A'B\big)$; use uniqueness. (d) Take $A=1$. If $\rho$ is a state,
$A\mapsto\operatorname{Tr}(\rho(A\otimes1))$ is positive, as $A^*A\otimes1=(A\otimes1)^*(A\otimes1)$, and
is $1$ at $A=1$, so its density $\operatorname{Tr}_2\rho$ is a state (Theorem 3.2). $\square$

Channels are the subject of Chapter 6; here we need only this definition: a linear map
$\Phi$ is **completely positive** if $\mathrm{id}_n\otimes\Phi$ is positive on
$M_n(\mathbb C)\otimes(\cdot)$ for every $n$, and a **channel** if it also preserves the trace.

**Theorem 7.3.** Let $j$ run over the standard basis of $H_2$ and $K_j=1\otimes\langle j|:H_1\otimes H_2\to H_1$,
$K_j(x\otimes y)=\langle j,y\rangle x$. Then
$$
\operatorname{Tr}_{H_2}X=\sum_jK_jXK_j^*,\qquad\sum_jK_j^*K_j=1,\qquad\sum_jK_j^*AK_j=A\otimes1 .
$$
Hence $\operatorname{Tr}_2$ is a channel from $\mathcal M_1\otimes\mathcal M_2$ to $\mathcal M_1$, and its
dual with respect to the blunt traces is $A\mapsto A\otimes1$.

*Proof.* $K_j^*x=x\otimes e_j$, so $\langle e_a,K_jXK_j^*e_{a'}\rangle=X_{(a,j),(a',j)}$; summing over $j$
gives Proposition 7.2(a). Further $\sum_jK_j^*K_j(x\otimes y)=\sum_j\langle j,y\rangle x\otimes e_j=x\otimes y$
and $\sum_jK_j^*AK_j=\sum_jA\otimes|j\rangle\langle j|=A\otimes1$. For $Y\ge0$,
$(\mathrm{id}_n\otimes\operatorname{Tr}_{H_2})(Y)=\sum_j(1_n\otimes K_j)Y(1_n\otimes K_j)^*\ge0$, because
$\langle x,TYT^*x\rangle=\langle T^*x,YT^*x\rangle\ge0$; this persists on the subalgebra
$\mathcal M_1\otimes\mathcal M_2$. Finally $\operatorname{Tr}\sum_jK_jXK_j^*A=\operatorname{Tr}X\sum_jK_j^*AK_j$,
which for $A=1$ is trace preservation and in general the duality. $\square$

`partial_trace(X, keep=1)` reshapes each block $\rho_{cd}$ to $(k_c,l_d,k_c,l_d)$ and traces;
`composite.partial_trace_channel(alg12, keep=1)` is the channel of Theorem 7.3, with Kraus operators
$1_{k_c}\otimes\langle j|$ from sector $(c,d)$ to sector $c$. We check both against the definition,
an independent `numpy.einsum`, and the duality:

```python
rho = states.random_density_matrix(alg12, batch_size=3)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
assert isinstance(rho1, DensityMatrix)                                      # states give states
R = full12(rho).numpy().reshape(3, 3, 5, 3, 5)                              # indices (batch, i, j, i', j')
assert np.abs(full(rho1).numpy() - np.einsum('xijkj->xik', R)).max() < 1e-12
assert np.abs(full(rho2).numpy() - np.einsum('xijil->xjl', R)).max() < 1e-12
one2 = alg2.identity()
assert (states.expectation(rho, kron(A, one2, alg12)) - states.expectation(rho1, A)).abs().max() < 1e-12
assert close(partial_trace(AB, keep=1), A * B.trace.item())                 # Proposition 7.2(b)

Ptr = composite.partial_trace_channel(alg12, keep=1)
print(Ptr)
assert Ptr.is_trace_preserving()
assert close(Ptr(rho), rho1)
assert close(Ptr.adjoint('blunt')(A), kron(A, one2, alg12))                 # Theorem 7.3: the dual
```

> **Physical meaning.** $\operatorname{Tr}_2\rho$ is the **reduced state**: all that an observer who can
> only measure the first system learns — the definition says it reproduces every expectation of an
> observable $A\otimes1$. Joint states with equal reduced states differ in their correlations.

**The partial trace as a conditional expectation.** Let $\mathcal M_2=M_m(\mathbb C)$, so that the
sectors of $\mathcal M_1\otimes M_m$ are those of $\mathcal M_1$ (Theorem 7.1), let
$\iota(A)=A\otimes1_m$ and $\mathcal E=\frac1m\operatorname{Tr}_2:\mathcal M_1\otimes M_m\to\mathcal M_1$.
Let $e_{ab}$ be the matrix units of $M_m$.

**Proposition 7.4.** (a) $\iota\circ\mathcal E$ is a conditional expectation of $\mathcal M_1\otimes M_m$
onto $\mathcal M_1\otimes1$ (Section 2.6), and for every faithful trace $\operatorname{tr}_w$ of
$\mathcal M_1\otimes M_m$ it is the only one that preserves $\operatorname{tr}_w$. (b) The $m^2$ elements
$u_{ab}=\sqrt m\,(1\otimes e_{ab})$ satisfy $X=\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*X))$ for all $X$,
and $\sum_{a,b}u_{ab}u_{ab}^*=m^2\,1$. (c) Every family $(v_k)$ with
$X=\sum_kv_k\,\iota(\mathcal E(v_k^*X))$ for all $X$ has $\sum_kv_kv_k^*=m^2\,1$.

*Proof.* (a) $\mathcal E$ is positive (Theorem 7.3) and so is $\iota$
($T^*T\otimes1=(T\otimes1)^*(T\otimes1)$); $\mathcal E(A\otimes1)=A$ by Proposition 7.2(b), and the
bimodule property is Proposition 7.2(c). Write $\operatorname{tr}_w$ also for the trace with the same
weights $w_c$ on $\mathcal M_1$; then $\operatorname{tr}_w(Y\otimes1)=m\operatorname{tr}_w(Y)$. The definition of
$\operatorname{Tr}_2$ applied to $AP_c$, multiplied by $w_c$ and summed over $c$, gives
$\operatorname{tr}_w(X(A\otimes1))=\operatorname{tr}_w((\operatorname{Tr}_2X)A)$, and $A=1$ gives
$\operatorname{tr}_w(\iota\mathcal E(X))=\operatorname{tr}_w(X)$. If $F=\iota\circ F_0$ is another conditional
expectation onto $\mathcal M_1\otimes1$ preserving $\operatorname{tr}_w$, then
$\operatorname{tr}_w(X(A\otimes1))=\operatorname{tr}_w(F(X)(A\otimes1))=m\operatorname{tr}_w(F_0(X)A)$, so
$\operatorname{tr}_w((\operatorname{Tr}_2X-mF_0(X))A)=0$ for all $A\in\mathcal M_1$, and $F_0=\mathcal E$ by faithfulness
(take $A=(\operatorname{Tr}_2X-mF_0(X))^*$). (b) Write $X=\sum_{a,b}X_{ab}\otimes e_{ab}$ with
$X_{ab}\in\mathcal M_1$ (Theorem 7.1(a)). By Proposition 7.2(b),
$\mathcal E(u_{ab}^*X)=\frac{\sqrt m}m\sum_{b'}\operatorname{Tr}_2(X_{ab'}\otimes e_{bb'})=X_{ab}/\sqrt m$, so
$u_{ab}\,\iota(\mathcal E(u_{ab}^*X))=X_{ab}\otimes e_{ab}$; and
$\sum_{a,b}u_{ab}u_{ab}^*=m\sum_{a,b}1\otimes e_{aa}=m^21$. (c) $\operatorname{Tr}_2(Y^*)=(\operatorname{Tr}_2Y)^*$ by
Proposition 7.2(a), so the adjoint of the expansion of $Y=u_{ab}$ is
$u_{ab}^*=\sum_k\iota(\mathcal E(u_{ab}^*v_k))\,v_k^*$, and by (b)
$\sum_{a,b}u_{ab}u_{ab}^*=\sum_k\big(\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*v_k))\big)v_k^*=\sum_kv_kv_k^*$. $\square$

A family as in (b) is a *quasi-basis* of $\mathcal E$, and $\sum_kv_kv_k^*$, which by (c) does not
depend on it, is the **Watatani index** $\operatorname{Ind}\mathcal E=m^2$ (Y. Watatani, Mem. Amer. Math.
Soc. 83, no. 424, 1990). For a factor, $M_k\otimes1\subset M_k\otimes M_m$ contains $\mathbb C^k$ with
multiplicity $m$ in $\mathbb C^{km}$: its inclusion matrix is $\Lambda=(m)$, and
$m^2=\|\Lambda\|^2=\dim(M_k\otimes M_m)/\dim M_k$ is its Jones index $[\mathcal M:\mathcal N]$ (V. F. R. Jones,
Invent. Math. 72, 1, 1983; F. M. Goodman, P. de la Harpe and V. F. R. Jones, *Coxeter graphs and
towers of algebras*, Springer, 1989). For general $\mathcal M_1$ the inclusion is the direct sum of the
$C_1$ inclusions $M_{k_c}\otimes1\subset M_{k_c}\otimes M_m$, $\Lambda=m1_{C_1}$. The module
`torch_vn_algebra.inclusions` implements inclusions, their conditional expectations and indices;
`from_tensor_product(alg1, m)` builds this one:

```python
from torch_vn_algebra import inclusions
m = 3
inc = inclusions.from_tensor_product(alg1, m)                # alg1 (x) 1 c alg1 (x) M_3
X = random_operator(inc.M)
for trace in ('Tr_blunt', 'Tr_norm', 'tau_vN'):
    E = inc.conditional_expectation(trace.replace('Tr_', ''))
    assert close(E(X), partial_trace(X, keep=1) * (1 / m))   # Proposition 7.4(a): the same E for every trace
    assert abs(getattr(inc.embed(E(X)), trace)() - getattr(X, trace)()).max() < 1e-10
    assert close(E.index(), inc.M.identity() * m**2)         # Proposition 7.4(c)

def unit(a, b):
    e = torch.zeros(m, m, dtype=dt)
    e[a, b] = 1
    return math.sqrt(m) * kron(alg1.identity(), inc.M.tensor_factors[1].from_blocks([e]), inc.M)

u = [unit(a, b) for a in range(m) for b in range(m)]         # the quasi-basis of Proposition 7.4(b)
assert close(sum((v @ inc.embed(E(v.adjoint() @ X))).matrix for v in u), X)
assert close(sum((v @ v.adjoint()).matrix for v in u), inc.M.identity() * m**2)
print(inc.matrix.tolist(), inc.component_indices())           # Lambda = 3 * 1_2, index 9 per component
```

## 7.3 Fusion of charges

Suppose now that sector $c$ of $\mathcal M_1$ carries a charge $q_c$, sector $d$ of $\mathcal M_2$ a
charge $r_d$, and that a **fusion rule** $f$ gives the total charge $f(q_c,r_d)$ — addition for
particle numbers, addition modulo $2$ for fermion parity, or any other map. Put
$$
H_q=\bigoplus_{f(q_c,r_d)=q}\mathbb C^{k_c}\otimes\mathbb C^{l_d},\qquad
\Pi_q=\sum_{f(q_c,r_d)=q}P_c\otimes Q_d,\qquad n_q=\sum_{f(q_c,r_d)=q}k_cl_d .
$$
The **fused product** $\mathcal M_1\otimes_f\mathcal M_2=\{\Pi_q\}'$ consists of the operators on
$H_1\otimes H_2$ that conserve the total charge.

**Proposition 7.5.** (a) $\mathcal M_1\otimes_f\mathcal M_2=\bigoplus_qB(H_q)\cong\bigoplus_qM_{n_q}(\mathbb C)$,
with minimal central projections $\Pi_q$. (b) $\mathcal M_1\otimes\mathcal M_2\subseteq\mathcal M_1\otimes_f\mathcal M_2$,
with equality iff every $H_q$ contains only one pair $(c,d)$. (c) If the charges are real numbers
and $f$ is addition, then $\mathcal M_1\otimes_f\mathcal M_2=\{N\}'$ with
$N=N_1\otimes1+1\otimes N_2$, $N_1=\sum_cq_cP_c$, $N_2=\sum_dr_dQ_d$.

*Proof.* (a) Every pair lies in exactly one $H_q$, so the $\Pi_q$ are mutually orthogonal
projections with sum $1$, and $X$ commutes with all of them iff $X=\sum_q\Pi_qX\Pi_q$, with
arbitrary blocks in $B(H_q)$; the centre is as in Theorem 7.1(b). (b) Every element of
$\mathcal M_1\otimes\mathcal M_2$ commutes with each $P_c\otimes Q_d$, hence with their sums $\Pi_q$.
If each $H_q$ is one $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$, both algebras are
$\bigoplus_{c,d}B(\mathbb C^{k_c}\otimes\mathbb C^{l_d})$; if $H_q$ contains two pairs, $B(H_q)$ contains
operators mapping one pair block into the other, which are not block diagonal. (c)
$N=\sum_{c,d}(q_c+r_d)P_c\otimes Q_d=\sum_qq\,\Pi_q$, so the nonzero $\Pi_q$ are the spectral
projections of $N$. They are polynomials in $N$ (proof of Theorem 1.5(a)) and $N$ is a combination
of them, so $X$ commutes with $N$ iff it commutes with every $\Pi_q$. $\square$

`fused_tensor_product(alg1, alg2, fuse)` (default `fuse`: addition) orders its sectors by the fused
charge and the pairs $(c,d)$ inside a sector lexicographically, each with the basis $e_a\otimes e_b$.
Example: hard-core bosons. The number-conserving operators of two sites form
$M_1\oplus M_2\oplus M_1$ (charges $0,1,2$), of three sites $M_1\oplus M_3\oplus M_3\oplus M_1$; the fused
product is that of five sites, up to the order of the basis: $n_q=\sum_c\binom2c\binom3{q-c}=\binom5q$.

```python
bos1 = TypeIAlgebra([1, 2, 1], [1, 2, 1], precision='double', charges=[0, 1, 2])
bos2 = TypeIAlgebra([1, 3, 3, 1], [1, 3, 3, 1], precision='double', charges=[0, 1, 2, 3])
F = fused_tensor_product(bos1, bos2)                         # U(1): total particle number
P = tensor_product(bos1, bos2)                               # particle numbers conserved separately
print(F.charges, F.k_factors, sum(k * k for k in F.k_factors), sum(k * k for k in P.k_factors))
assert F.k_factors == [math.comb(5, q) for q in range(6)]

N1, N2 = full(bos1.central([0.0, 1.0, 2.0]))[0], full(bos2.central([0.0, 1.0, 2.0, 3.0]))[0]
N = torch.kron(N1, torch.eye(8, dtype=dt)) + torch.kron(torch.eye(4, dtype=dt), N2)
X = full12(F.operator(torch.randn(1, F.C, F.k_max, F.k_max, dtype=dt)))[0]
assert close(X @ N, N @ X)                                   # Proposition 7.5(c)
A, B = random_operator(bos1), random_operator(bos2)
assert close(full12(kron(A, B, F)), full12(kron(A, B, P)))   # Proposition 7.5(b)
```

The fused algebra has dimension $252$, the plain product $120$; the difference consists of the
operators that move particles between the two parts, such as hopping terms.

A state of $\mathcal M_1\otimes_f\mathcal M_2$ can have coherences between pair blocks, so Proposition
7.2 does not apply. Let $E_1(X)=\sum_cP_cXP_c$ on $B(H_1)$ (the conditional expectation of Section
2.6, Physical meaning); it maps onto $\mathcal M_1$, and $\operatorname{Tr}(E_1(X)A)=\operatorname{Tr}(XA)$ for
$A\in\mathcal M_1$, because $A=\sum_cP_cAP_c$.

**Theorem 7.6.** Let $\rho$ be the density of a state $\omega$ on $\mathcal M_1\otimes_f\mathcal M_2$. The
density of its restriction $A\mapsto\omega(A\otimes1)$ to $\mathcal M_1$ is
$$
\rho_1=E_1\big(\operatorname{Tr}_{H_2}\rho\big).
$$
The block $(c,c')$ of $\operatorname{Tr}_{H_2}\rho$ is $\sum_d\operatorname{tr}_{l_d}\rho_{(c,d),(c',d)}$, summed over
the $d$ with $f(q_c,r_d)=f(q_{c'},r_d)$, where $\rho_{(c,d),(c',d)}$ is the block of $\rho$ from the
pair $(c',d)$ to the pair $(c,d)$. Hence $E_1$ acts trivially on all $\operatorname{Tr}_{H_2}\rho$ iff there
are no $c\neq c'$ and $d$ with $f(q_c,r_d)=f(q_{c'},r_d)$; for instance if $f$ is addition and the
$q_c$ are pairwise distinct.

*Proof.* $A\otimes1\in\mathcal M_1\otimes_f\mathcal M_2$ (Proposition 7.5(b)) and
$\omega(A\otimes1)=\operatorname{Tr}(\rho(A\otimes1))=\operatorname{Tr}\big((\operatorname{Tr}_{H_2}\rho)A\big)
=\operatorname{Tr}\big(E_1(\operatorname{Tr}_{H_2}\rho)A\big)$. As $E_1(\operatorname{Tr}_{H_2}\rho)\in\mathcal M_1$, it is
the density by the uniqueness in Theorem 3.2. The block formula is Proposition 7.2(a) for
$\operatorname{Tr}_{H_2}$, with $i$ in block $c$, $i'$ in block $c'$ and $j$ in block $d$; the block of
$\rho$ between $(c',d)$ and $(c,d)$ vanishes unless both pairs lie in the same $H_q$. If $c\neq c'$
and $d$ with $f(q_c,r_d)=f(q_{c'},r_d)=q$ exist, take unit vectors $x\in\mathbb C^{k_c}$,
$x'\in\mathbb C^{k_{c'}}$, $y\in\mathbb C^{l_d}$ and $\psi=(x+x')\otimes y/\sqrt2\in H_q$. Then
$|\psi\rangle\langle\psi|\in B(H_q)$, and $\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|=\frac12|x+x'\rangle\langle x+x'|$
has the nonzero block $\frac12|x\rangle\langle x'|$, which $E_1$ removes. $\square$

Both library implementations compute this map, for any fusion rule: `partial_trace` traces only
the diagonal pair blocks $\rho_{(c,d),(c,d)}$, and `partial_trace_channel` uses separate Kraus
indices for pairs of one sector with the same $c$; its dual is again $A\mapsto A\otimes1$. With
addition $E_1$ does nothing in our example; with the parity rule $f(q,r)=q+r\bmod2$ the charges $0$
and $2$ of `bos1` fuse identically with every $r$, and $E_1$ removes genuine coherences:

```python
G = fused_tensor_product(bos1, bos2, fuse=lambda q, r: (q + r) % 2)       # only the parity is conserved
print(G.k_factors)

def block_mask(alg):
    """Ones in the diagonal blocks of (+)_c C^{k_c}: X * block_mask(alg) = sum_c P_c X P_c."""
    return full(alg.operator(torch.ones(1, alg.C, alg.k_max, alg.k_max, dtype=dt)))

for alg, coherent in ((F, False), (G, True)):
    rho = states.random_density_matrix(alg, batch_size=2)
    R = full12(rho).reshape(2, 4, 8, 4, 8)
    for keep, X_H, alg_k in ((1, torch.einsum('xijkj->xik', R), bos1), (2, torch.einsum('xijil->xjl', R), bos2)):
        E_X = X_H * block_mask(alg_k)                                     # E(Tr_H rho), Theorem 7.6
        assert close(full(partial_trace(rho, keep=keep)), E_X)
        assert close(full(composite.partial_trace_channel(alg, keep=keep)(rho)), E_X)
        lost = (X_H - E_X).abs().max().item()                             # coherences removed by E
        print(f"keep={keep}: largest coherence removed: {lost:.3e}")
        assert (lost > 1e-3) if coherent else (lost == 0)
    assert close(composite.partial_trace_channel(alg, keep=1).adjoint('blunt')(A), kron(A, bos2.identity(), alg))
```

> **Physical meaning.** In the BCS description of a superconductor only the parity of the electron
> number is conserved, so the global state may superpose the particle numbers $0$ and $2$ of a
> region. An observer whose observables conserve the particle number of the region cannot see
> this coherence: for them the reduced state is the pinched one of Theorem 7.6, a mixture of particle numbers.

## 7.4 Entanglement of pure states

**Lemma 7.7 (Schmidt decomposition).** Every unit vector $\psi\in\mathbb C^a\otimes\mathbb C^b$ can be
written $\psi=\sum_{i=1}^rs_i\,u_i\otimes v_i$ with $s_i>0$, $\sum_is_i^2=1$ and orthonormal
$u_i\in\mathbb C^a$, $v_i\in\mathbb C^b$. Then $\operatorname{Tr}_2|\psi\rangle\langle\psi|=\sum_is_i^2|u_i\rangle\langle u_i|$
and $\operatorname{Tr}_1|\psi\rangle\langle\psi|=\sum_is_i^2|v_i\rangle\langle v_i|$: the two reduced states
have the same nonzero eigenvalues $s_i^2$.

*Proof.* Let $\Psi$ be the $a\times b$ matrix $\Psi_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,\psi\rangle$. By the
spectral theorem $\Psi\Psi^*=\sum_{i=1}^a\lambda_i|u_i\rangle\langle u_i|$ with an orthonormal basis
$u_i$ and $\lambda_i=\|\Psi^*u_i\|^2\ge0$, ordered so that $\lambda_i>0$ exactly for $i\le r$. Put
$s_i=\sqrt{\lambda_i}$ and $w_i=\Psi^*u_i/s_i$ for $i\le r$; then
$\langle w_i,w_j\rangle=\langle u_i,\Psi\Psi^*u_j\rangle/(s_is_j)=\delta_{ij}$, and $\Psi^*u_i=0$ for $i>r$.
Hence $\Psi=\sum_i|u_i\rangle\langle u_i|\Psi=\sum_{i\le r}|u_i\rangle\langle\Psi^*u_i|=\sum_{i\le r}s_i|u_i\rangle\langle w_i|$,
i.e. $\Psi_{\alpha\beta}=\sum_is_i(u_i)_\alpha\overline{(w_i)_\beta}$: the claim with $v_i=\bar w_i$, and
$\sum_is_i^2=\operatorname{Tr}\Psi\Psi^*=\|\psi\|^2=1$. Finally
$|\psi\rangle\langle\psi|=\sum_{i,j}s_is_j|u_i\rangle\langle u_j|\otimes|v_i\rangle\langle v_j|$, and
Proposition 7.2(b) gives the reduced states. $\square$

The **entanglement entropy** of $\psi$ is
$E(\psi)=S(\operatorname{Tr}_2|\psi\rangle\langle\psi|)=S(\operatorname{Tr}_1|\psi\rangle\langle\psi|)=-\sum_is_i^2\log s_i^2$.
It vanishes iff $r=1$, i.e. iff $\psi=u\otimes v$, and $E(\psi)\le\log\min(a,b)$ with equality iff the
smaller reduced state is maximally mixed (Corollary 3.6).

With charges, a pure state lies in one sector $H_q$ of the fused product; write
$\psi=\sum_{(c,d)}\psi_{cd}$ with $\psi_{cd}\in\mathbb C^{k_c}\otimes\mathbb C^{l_d}$.

**Proposition 7.8 (symmetry-resolved entanglement).** Let $f$ be addition and let the $q_c$ be
pairwise distinct and the $r_d$ pairwise distinct, so that for each $c$ at most one $d=d(c)$ has
$q_c+r_d=q$. Put $\psi_c=\psi_{c\,d(c)}$ (or $0$), $p_c=\|\psi_c\|^2$ and, for $p_c>0$,
$\hat\psi_c=\psi_c/\sqrt{p_c}$. Then $\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|$ lies in $\mathcal M_1$ (it
commutes with $N_1$) and equals $\bigoplus_cp_c\operatorname{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|$; the
reduced state on $\mathcal M_2$ has the sector probability $p_c$ at $d(c)$; and
$$
E(\psi)=H(p)+\sum_cp_c\,E(\hat\psi_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*Proof.* By the block formula of Theorem 7.6, the block $(c,c')$ of $\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|$
is the sum of $\operatorname{tr}_{l_d}|\psi_{cd}\rangle\langle\psi_{c'd}|$ over the $d$ with
$q_c+r_d=q=q_{c'}+r_d$, which forces $c=c'$ and $d=d(c)$. So the reduced state is
$\bigoplus_c\operatorname{Tr}_2|\psi_c\rangle\langle\psi_c|$, with block traces $p_c$; exchanging the roles of
the factors gives the statement for $\mathcal M_2$. By Lemma 7.7, $E(\psi)$ is the entropy of this
state, and Theorem 3.3 with $\sigma_c=\operatorname{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|$, $S(\sigma_c)=E(\hat\psi_c)$,
gives the formula. $\square$

$H(p)$ is the **number entanglement** and $\sum_cp_cE(\hat\psi_c)$ the **configurational
entanglement**. For a random state of $2$ particles on the five sites we compute $E(\psi)$ from the
partial traces on both sides and from the singular values of the $4\times8$ matrix $\Psi$:

```python
s = F.charges.index(2)                                      # the sector of two particles, n_q = 10
n = F.k_factors[s]
psi = torch.randn(n, dtype=dt)
psi = psi / psi.norm()
mat = torch.zeros(1, F.C, F.k_max, F.k_max, dtype=dt)
mat[0, s, :n, :n] = torch.outer(psi, psi.conj())
rho = DensityMatrix(F, matrix=mat)
r1, r2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)

psi_H = torch.zeros(4 * 8, dtype=dt)
psi_H[kron_rows(F)[s]] = psi                                # psi as a vector of H1 (x) H2
lam = torch.linalg.svdvals(psi_H.reshape(4, 8)) ** 2        # squared Schmidt coefficients
lam = lam[lam > 1e-14]
E_psi = -(lam * lam.log()).sum().item()
assert abs(r1.von_neumann_entropy().item() - E_psi) < 1e-10
assert abs(r2.von_neumann_entropy().item() - E_psi) < 1e-10  # Lemma 7.7

def entropy(sigma):
    w = torch.linalg.eigvalsh(sigma).clamp(min=1e-300)
    return -(w * w.log()).sum().item()

p = r1.sector_probabilities()[0]                            # c = 0, 1, 2 particles in the first part
assert torch.allclose(r2.sector_probabilities()[0, [2, 1, 0]], p)   # d(c) = 2 - c
H_num = -(p * p.log()).sum().item()
S_conf = sum(p[c].item() * entropy(r1.matrix[0, c, :k, :k] / p[c]) for c, k in enumerate(bos1.k_factors))
print(f"E = {E_psi:.4f} = H(p) {H_num:.4f} + configurational {S_conf:.4f}")
assert abs(E_psi - H_num - S_conf) < 1e-10                  # Proposition 7.8
```

> **Physical meaning.** Parties bound by the particle-number superselection rule can use only the
> configurational part as a resource: the *entanglement of particles* of H. M. Wiseman and J. A.
> Vaccaro (Phys. Rev. Lett. 91, 097902, 2003). $H(p)$ is measured in cold atoms as the *number
> entropy* (A. Lukin et al., Science 364, 256, 2019).

## 7.5 Random states and Page's formula

Let $g\in\mathbb C^N$ have independent coordinates with independent $N(0,\tfrac12)$ real and imaginary
parts, so $\mathbb E|g_i|^2=1$. A **Haar-random pure state** is $\psi=g/\|g\|$; its law is unitarily
invariant, as the density $\pi^{-N}e^{-\|g\|^2}$ and Lebesgue measure are (it is the only invariant
law on the sphere, a fact we do not use).

**Lemma 7.9.** $\|g\|$ and $\psi=g/\|g\|$ are independent.

*Proof.* Let $h\ge0$ and $\phi\ge0$ be measurable on $(0,\infty)$ and on the unit sphere $S^{2N-1}$.
Integration in polar coordinates in $\mathbb R^{2N}$, with $\sigma$ the surface measure, gives
$$
\mathbb E\big[h(\|g\|)\phi(\psi)\big]=\pi^{-N}\int_0^\infty h(t)e^{-t^2}t^{2N-1}dt\int_{S^{2N-1}}\phi\,d\sigma=:\alpha(h)\beta(\phi).
$$
Taking $\phi=1$, then $h=1$, then both, gives $\mathbb Eh(\|g\|)=\alpha(h)\beta(1)$,
$\mathbb E\phi(\psi)=\alpha(1)\beta(\phi)$ and $\alpha(1)\beta(1)=1$, hence
$\mathbb E[h(\|g\|)\phi(\psi)]=\mathbb Eh(\|g\|)\,\mathbb E\phi(\psi)$. $\square$

**Lemma 7.10.** $\mathbb E\|g\|^4=N(N+1)$. If $G$ is an $a\times b$ matrix with independent entries distributed
as the $g_i$ (a *Ginibre* matrix), then $\mathbb E\operatorname{Tr}(GG^*)^2=ab(a+b)$.

*Proof.* With $g_i=x+iy$, $\mathbb E|g_i|^4=\mathbb Ex^4+2\mathbb Ex^2\mathbb Ey^2+\mathbb Ey^4=\tfrac34+\tfrac12+\tfrac34=2$,
so $\mathbb E\|g\|^4=\sum_i\mathbb E|g_i|^4+\sum_{i\ne j}\mathbb E|g_i|^2\mathbb E|g_j|^2=2N+N(N-1)$. Next,
$\operatorname{Tr}(GG^*)^2=\sum_{i,i',j,j'}G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}$. The entries and their
conjugates are complex linear combinations of independent real Gaussians with
$\mathbb E[G_{ij}G_{kl}]=0$ and $\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}$, so Isserlis' formula
(proof of Proposition 1.10) gives $\mathbb E[G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}]=\delta_{ii'}+\delta_{jj'}$.
Summing gives $ab^2+a^2b$. $\square$

**Theorem 7.11 (E. Lubkin, J. Math. Phys. 19, 1028, 1978).** For a Haar-random
$\psi\in\mathbb C^a\otimes\mathbb C^b$ and $\rho_A=\operatorname{Tr}_2|\psi\rangle\langle\psi|$,
$$
\mathbb E\operatorname{Tr}\rho_A^2=\frac{a+b}{ab+1}.
$$

*Proof.* Let $G_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,g\rangle$, a Ginibre matrix; the matrix of $\psi$ is
$\Psi=G/\|g\|$, and by the proof of Lemma 7.7, $\rho_A=\Psi\Psi^*=GG^*/\|g\|^2$. Hence
$\operatorname{Tr}(GG^*)^2=\|g\|^4\operatorname{Tr}\rho_A^2$, where $\operatorname{Tr}\rho_A^2\in[0,1]$ is a function of
$\psi$. By Lemma 7.9, $\mathbb E\operatorname{Tr}(GG^*)^2=\mathbb E\|g\|^4\,\mathbb E\operatorname{Tr}\rho_A^2$, and by Lemma 7.10
$\mathbb E\operatorname{Tr}\rho_A^2=ab(a+b)/(ab(ab+1))$. $\square$

Chapter 6 obtains the same mean from the fourth moments of a Haar vector, as the output purity of a
random channel.

**Corollary 7.12.** $\log\frac{ab+1}{a+b}\le\mathbb E\,S(\rho_A)\le\log\min(a,b)$.

*Proof.* By Corollary 3.6, $S(\rho_A)\le\log a$, and by Lemma 7.7,
$S(\rho_A)=S(\operatorname{Tr}_1|\psi\rangle\langle\psi|)\le\log b$. Jensen's inequality for the concave $\log$ with the
weights $\lambda_i>0$, the nonzero eigenvalues of $\rho_A$, gives
$\sum_i\lambda_i\log\lambda_i\le\log\sum_i\lambda_i^2$, i.e. $S(\rho_A)\ge-\log\operatorname{Tr}\rho_A^2$; Jensen's
inequality for the convex $-\log$ and Theorem 7.11 give
$\mathbb E[-\log\operatorname{Tr}\rho_A^2]\ge-\log\mathbb E\operatorname{Tr}\rho_A^2$. $\square$

For fixed $a$ and $b\to\infty$ the bounds meet: a small part of a random pure state is almost
maximally mixed. The exact mean is Page's formula, which we state without proof.

**Theorem 7.13 (Page's formula).** For $a\le b$,
$$
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
$$
It was conjectured by D. N. Page (Phys. Rev. Lett. 71, 1291, 1993) and proved by S. K. Foong and
S. Kanno (Phys. Rev. Lett. 72, 1148, 1994), J. Sánchez-Ruiz (Phys. Rev. E 52, 5653, 1995) and
S. Sen (Phys. Rev. Lett. 77, 1, 1996). For $1\ll a\le b$ it is $\approx\log a-\frac a{2b}$; for two
qubits it is $\frac13$.

The law of $\rho_A$ is the **induced measure** (K. Życzkowski and H.-J. Sommers, J. Phys. A 34,
7111, 2001). `states.random_density_matrix(alg, rank=r)` forms $G_cG_c^*$ with $k_c\times r$ Ginibre
blocks $G_c$ and divides by $\sum_c\operatorname{Tr}G_cG_c^*$. For one block this is the law of $\rho_A$
above with $b=r$. For several blocks it is the restriction to $\mathcal M_1$ of a Haar-random pure state
of $H_1\otimes\mathbb C^r$: apply Theorem 7.6 with $\mathcal M_2=M_r$ and a constant fusion rule, so that
$\mathcal M_1\otimes_f\mathcal M_2=B(H_1\otimes\mathbb C^r)$, with $G_c$ the rows of $G$ in block $c$.

**Proposition 7.14.** For $\rho$ drawn by `random_density_matrix(alg, rank=r)` on $\bigoplus_cM_{k_c}$,
$D=\sum_ck_c$,
$$
\mathbb E\operatorname{Tr}\rho^2=\frac{\sum_ck_c(k_c+r)}{D(Dr+1)} .
$$

*Proof.* $\operatorname{Tr}\rho^2=\sum_c\operatorname{Tr}(G_cG_c^*)^2/\|g\|^4$ with $g\in\mathbb C^{Dr}$ the collection of
all entries; as in Theorem 7.11, Lemmas 7.9 and 7.10 give $\sum_ck_cr(k_c+r)/(Dr(Dr+1))$. $\square$

We check Theorems 7.11 and 7.13 on $10^4$ pure states of $\mathbb C^3\otimes\mathbb C^6$ reduced with
`partial_trace` and on the induced measure, and Proposition 7.14 for $M_1\oplus M_2\oplus M_3$, $r=2$:

```python
a, b, n_samples = 3, 6, 10_000
g = torch.randn(n_samples, a * b, dtype=dt)                 # complex Gaussian: E|g_i|^2 = 1
psi = g / g.norm(dim=-1, keepdim=True)
AB_alg = tensor_product(TypeIAlgebra([a], [a], precision='double'), TypeIAlgebra([b], [b], precision='double'))
pure = DensityMatrix(AB_alg, matrix=(psi[:, :, None] * psi.conj()[:, None, :])[:, None],
                     validate=False)                        # pure by construction
rho_A = partial_trace(pure, keep=1)
purity_exact = (a + b) / (a * b + 1)
page = sum(1 / k for k in range(b + 1, a * b + 1)) - (a - 1) / (2 * b)

def mean_se(x):
    return x.mean().item(), x.std().item() / math.sqrt(x.numel())

for name, rho_s in (('reduced pure states', rho_A),
                    ('induced measure', states.random_density_matrix(rho_A.algebra, 40_000, rank=b))):
    (m_p, e_p), (m_S, e_S) = mean_se(rho_s.purity()), mean_se(rho_s.von_neumann_entropy())
    print(f"{name}: E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {purity_exact:.4f}), "
          f"E S = {m_S:.4f} +- {e_S:.4f} (Page {page:.4f})")
    assert abs(m_p - purity_exact) < 5 * e_p and abs(m_S - page) < 5 * e_S     # Theorems 7.11, 7.13
    assert math.log((a * b + 1) / (a + b)) < m_S < math.log(a)                # Corollary 7.12

alg3 = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
m_p, e_p = mean_se(states.random_density_matrix(alg3, 40_000, rank=2).purity())
print(f"E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {26 / 78:.4f})")
assert abs(m_p - 26 / 78) < 5 * e_p                                          # Proposition 7.14
```

> **Physical meaning.** A typical pure state is nearly maximally entangled between any part of at
> most half the system and the rest. As the part grows, the entanglement rises like $\log a$ and
> falls beyond the middle: the *Page curve*, which describes how information leaves a black hole
> in the Hawking radiation, and why the parts of a closed system look thermal.

## 7.6 What to compare with

On one Hilbert space, QuTiP (`Qobj.ptrace`, `entropy_vn`, `rand_ket`, `rand_dm` with the Ginibre
distribution and a rank) and `numpy` (a reshape and an `einsum`, as above) compute partial traces,
entropies and induced states, and agree with the library on plain products. Tensor-network
libraries with conserved charges — ITensor (quantum-number block-sparse tensors), TeNPy (charged
legs, entanglement spectra by charge) and, for non-abelian symmetries, TensorKit.jl (fusion trees) —
fuse sectors in the same way for large systems and give symmetry-resolved entanglement. The
library adds the algebraic layer: products are von Neumann algebras with sectors and traces, the
fusion rule is any map of the labels, the partial trace is a channel with its dual for each trace
and the correct restriction for fused products (Theorems 7.3, 7.6), and all of it is batched. It has
no tensor networks, and `fused_tensor_product` takes a rule with one outcome per pair of charges, as
for abelian symmetries; a non-abelian rule such as $j_1\otimes j_2=\bigoplus_jj$ for $SU(2)$, with its
Clebsch–Gordan coefficients, is not a fused product here (the library treats $SU(2)$ for spin-$\frac12$
chains separately, Chapter 4).

## Exercises

1. Compute the dual of $\operatorname{Tr}_2$ with respect to $\tau_{\mathrm{vN}}$ on both
   $\mathcal M_1\otimes\mathcal M_2$ and $\mathcal M_1$, and check it with
   `partial_trace_channel(alg12).adjoint('tau_vN')`.
2. (Subadditivity.) Let $\rho$ be a state on $\mathcal M_1\otimes\mathcal M_2$ with invertible reduced
   states $\rho_1,\rho_2$. Show that $S(\rho)\le S(\rho_1)+S(\rho_2)$, with equality iff
   $\rho=\rho_1\otimes\rho_2$.
3. Let `site` be $\mathbb C\oplus\mathbb C$ with charges $0,1$ (an empty or occupied site), fused with
   itself by addition. Compute $E(\psi)$ and its two parts for $\psi=(|01\rangle+|10\rangle)/\sqrt2$, one
   particle shared by two sites, and describe the restriction of the vector state of
   $(|00\rangle+|11\rangle)/\sqrt2$ to the fused algebra.

<details>
<summary>Solutions</summary>

1. The weights are $1/(C_1k_c)$ on $\mathcal M_1$ and $1/(C_1C_2k_cl_d)$ on the product. For
   $A\in\mathcal M_1$, by Proposition 7.2(a),
   $\tau_{\mathrm{vN}}((\operatorname{Tr}_2X)A)=\sum_c\frac1{C_1k_c}\sum_d\operatorname{Tr}(X_{cd}(A_c\otimes1))
   =\sum_{c,d}\frac1{C_1C_2k_cl_d}\operatorname{Tr}\big(X_{cd}\,C_2l_d(A_c\otimes1)\big)$, so the dual is
   $A\mapsto\bigoplus_{c,d}C_2l_d\,A_c\otimes1_{l_d}=A\otimes Z$ with the central element
   $Z=\bigoplus_dC_2l_d1_{l_d}$ of $\mathcal M_2$. It is not unital: of the three traces of the
   library, only the blunt ones make $\operatorname{Tr}_2$ the dual of $A\mapsto A\otimes1$ (in general, $\operatorname{tr}_w$ on $\mathcal M_1$ and the weights
   $w_{cd}=w_c$ on the product do, by the same computation).

```python
dual_tau = composite.partial_trace_channel(alg12, keep=1).adjoint('tau_vN')
A = random_operator(alg1)
Z = alg2.central([alg2.C * l for l in alg2.k_factors])
assert close(dual_tau(A), kron(A, Z, alg12))
```

2. $\rho_1\otimes\rho_2$ is a state (Theorem 7.1(c)), invertible, with eigenvectors $u_i\otimes v_j$ and
   eigenvalues $\lambda_i\mu_j$, so $\log(\rho_1\otimes\rho_2)=\log\rho_1\otimes1+1\otimes\log\rho_2$
   (Theorem 1.5(a)). By the definition of the partial trace
   $\operatorname{Tr}\rho(\log\rho_1\otimes1)=\operatorname{Tr}\rho_1\log\rho_1$, and likewise for the second
   factor. Hence $S(\rho\|\rho_1\otimes\rho_2)=S(\rho_1)+S(\rho_2)-S(\rho)$, which is $\ge0$ with
   equality iff $\rho=\rho_1\otimes\rho_2$ (Corollary 3.6). This difference is the mutual information.

```python
rho = states.random_density_matrix(alg12, batch_size=4)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
mutual = rho1.von_neumann_entropy() + rho2.von_neumann_entropy() - rho.von_neumann_entropy()
assert (mutual > 0).all()
assert torch.allclose(mutual, states.relative_entropy(rho, kron(rho1, rho2, alg12)))
```

3. The sectors are $N=0,1,2$, of sizes $1,2,1$. $\psi$ lies in $N=1$, with $\psi_0=|01\rangle/\sqrt2$,
   $\psi_1=|10\rangle/\sqrt2$, so $p=(\frac12,\frac12)$ and, the blocks being one-dimensional,
   $E(\hat\psi_c)=0$: by Proposition 7.8, $E(\psi)=H(p)=\log2$ is all number entanglement, and the
   usable configurational part is $0$. The vector $(|00\rangle+|11\rangle)/\sqrt2$ superposes $N=0$ and
   $N=2$; on the block-diagonal fused algebra its vector state is the mixture
   $\frac12|00\rangle\langle00|+\frac12|11\rangle\langle11|$ of product states (Chapter 3, Exercise 3).

```python
site = TypeIAlgebra([1, 1], [1, 1], precision='double', charges=[0, 1])
S2 = fused_tensor_product(site, site)
assert S2.k_factors == [1, 2, 1]
mat = torch.zeros(1, 3, 2, 2, dtype=dt)
mat[0, 1] = 0.5                                       # |psi><psi| with psi = (|01> + |10>)/sqrt 2
r1 = partial_trace(DensityMatrix(S2, matrix=mat), keep=1)
assert torch.allclose(r1.sector_probabilities()[0], torch.tensor([0.5, 0.5], dtype=torch.float64))
assert abs(r1.von_neumann_entropy().item() - math.log(2)) < 1e-12
```

</details>

**Next:** Chapter 8 (in preparation) treats dynamics: the Schrödinger and Lindblad equations in an
algebra with sectors, jumps between sectors and post-selection.
