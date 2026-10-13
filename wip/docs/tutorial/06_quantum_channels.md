**English** · [Español](es/06_quantum_channels.md) · [Français](fr/06_quantum_channels.md) · [Deutsch](de/06_quantum_channels.md) · [中文](zh/06_quantum_channels.md) · [日本語](ja/06_quantum_channels.md) · [Русский](ru/06_quantum_channels.md)

# 6. Quantum channels

A quantum channel is the most general way in which a state can change: a unitary evolution,
a measurement whose result is forgotten, the interaction with an environment, the loss of a
particle. Mathematically it is a completely positive, trace-preserving linear map. This chapter
proves the standard descriptions of such maps — Choi matrices, Kraus operators, the Stinespring
dilation — for maps between finite-dimensional von Neumann algebras, and then studies what is
special about direct sums: maps between sectors, duals with respect to the traces of Chapter 3,
and the classical part of a channel. Chapters 2 and 3 already used the library's channels; here
is their theory.

**You will learn** why positivity is not enough, how to pass between Choi matrices, Kraus
operators, superoperators and isometries, why "trace preserving" depends on the trace once a map
changes the sector, which observables survive decoherence, and what a random channel does to a
pure state on average.

Throughout, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ and
$\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)$ with minimal central projections $P_c$ and
$Q_d$; $A=\bigoplus_cA_c$; $M_{k'\times k}$ is the space of complex $k'\times k$ matrices;
$E_{ij}=|i\rangle\langle j|$ are matrix units; $\operatorname{Tr}$ is the blunt trace
$\sum_c\operatorname{Tr}A_c$, and $\operatorname{tr}_w=\sum_cw_c\operatorname{Tr}A_c$ a faithful
trace (Theorem 3.1). Block matrices act on $\mathbb C^n\otimes\mathbb C^k$, block index first.

## 6.1 Completely positive maps and Choi's theorem

A linear map $\Phi:\mathcal M\to\mathcal N$ is **positive** if $A\ge0$ implies $\Phi(A)\ge0$.
The $n\times n$ matrices $[A_{ij}]$ with entries in $\mathcal M$ form again a direct sum of matrix
algebras, $M_n(\mathcal M)\cong\bigoplus_cM_n(M_{k_c})$ via $[A_{ij}]\mapsto([(A_{ij})_c]_{ij})_c$,
so $[A_{ij}]\ge0$ means that every block matrix $[(A_{ij})_c]\in M_{nk_c}$ is positive. $\Phi$ is
**$n$-positive** if $\Phi_n:[A_{ij}]\mapsto[\Phi(A_{ij})]$ is positive, and **completely positive**
(CP) if it is $n$-positive for every $n$. It is **trace preserving** (TP) if
$\operatorname{Tr}\Phi(A)=\operatorname{Tr}A$, **unital** if $\Phi(1)=1$, and a **channel** if it is
CP and TP. We use repeatedly that $X\ge0$ implies $TXT^*\ge0$ for every rectangular $T$, since
$\langle y,TXT^*y\rangle=\langle T^*y,XT^*y\rangle$.

**Proposition 6.1.** For $k\ge2$, the transpose $\theta(A)=A^T$ on $M_k$ is positive but not
$2$-positive.

*Proof.* For $A\ge0$, $\langle x,A^Tx\rangle=\sum_{ij}\bar x_iA_{ji}x_j=\langle\bar x,A\bar x\rangle\ge0$.
With $\Omega=e_1\otimes e_1+e_2\otimes e_2\in\mathbb C^2\otimes\mathbb C^k$,
$|\Omega\rangle\langle\Omega|=[E_{ij}]_{i,j=1}^2\ge0$, and $\theta_2$ maps it to
$F=[E_{ji}]_{i,j=1}^2$, which sends $e_1\otimes e_2-e_2\otimes e_1$ to its negative. $\square$

**Theorem 6.2 (Choi).** For a linear map $\Phi:M_k\to M_{k'}$ the following are equivalent:
(a) $\Phi$ is CP; (b) $\Phi$ is $k$-positive; (c) the **Choi matrix**
$$
J(\Phi)=\sum_{i,j=1}^kE_{ij}\otimes\Phi(E_{ij})=[\Phi(E_{ij})]_{i,j=1}^k\in M_k(M_{k'})
$$
is positive; (d) there are $K_1,\dots,K_r\in M_{k'\times k}$ (**Kraus operators**) with
$\Phi(A)=\sum_lK_lAK_l^*$, and one can take $r=\operatorname{rank}J(\Phi)\le kk'$.

*Proof.* (a)$\Rightarrow$(b) is trivial. (b)$\Rightarrow$(c): for
$\Omega=\sum_ie_i\otimes e_i\in\mathbb C^k\otimes\mathbb C^k$,
$|\Omega\rangle\langle\Omega|=[E_{ij}]\ge0$ and $J(\Phi)=\Phi_k([E_{ij}])$. (c)$\Rightarrow$(d): by
Theorem 1.3, $J(\Phi)=\sum_{l=1}^rv_lv_l^*$ with $r=\operatorname{rank}J(\Phi)$ (the eigenvectors of
the nonzero eigenvalues, scaled by their square roots). Every $v\in\mathbb C^k\otimes\mathbb C^{k'}$
is $v=\sum_ie_i\otimes Ke_i$ for exactly one $K\in M_{k'\times k}$ ($Ke_i$ is the $i$-th block of
$v$), and the $(i,j)$ block of $vv^*$ is $Ke_ie_j^*K^*=KE_{ij}K^*$. Writing $v_l$ in this way with
$K_l$ and comparing blocks gives $\Phi(E_{ij})=\sum_lK_lE_{ij}K_l^*$, hence (d) by linearity.
(d)$\Rightarrow$(a): $\Phi_n([A_{ij}])=\sum_l(1_n\otimes K_l)[A_{ij}](1_n\otimes K_l)^*\ge0$. $\square$

**Corollary 6.3.** If $\Phi(A)=\sum_{l=1}^sL_lAL_l^*$, then $s\ge\operatorname{rank}J(\Phi)$. So the
minimal number of Kraus operators, the **Kraus rank**, is $\operatorname{rank}J(\Phi)$.

*Proof.* Reading the last computation backwards, $J(\Phi)=\sum_{l=1}^sw_lw_l^*$ with
$w_l=\sum_ie_i\otimes L_le_i$, a sum of $s$ matrices of rank $\le1$. $\square$

A map between direct sums is a matrix of maps between blocks. Let $\iota_c:M_{k_c}\to\mathcal M$
put a matrix into block $c$, $\pi_d:\mathcal N\to M_{k'_d}$ take block $d$, and
$\Phi_{dc}=\pi_d\circ\Phi\circ\iota_c:M_{k_c}\to M_{k'_d}$. Since $A=\sum_c\iota_c(A_c)$,
$\Phi(A)_d=\sum_c\Phi_{dc}(A_c)$.

**Proposition 6.4.** $\Phi:\mathcal M\to\mathcal N$ is $n$-positive iff every $\Phi_{dc}$ is. Hence
$\Phi$ is CP iff there are $K^{dc}_l\in M_{k'_d\times k_c}$, $l=1,\dots,r$, with
$$
\Phi(A)_d=\sum_c\sum_{l=1}^rK^{dc}_l\,A_c\,K^{dc*}_l ,
$$
and $r=\max_{d,c}\operatorname{rank}J(\Phi_{dc})$ suffices. $\Phi$ is TP iff
$\sum_{d,l}K^{dc*}_lK^{dc}_l=1_{k_c}$ for every $c$.

*Proof.* If all $\Phi_{dc}$ are $n$-positive and $[A_{ij}]\ge0$, then each $X^c=[(A_{ij})_c]\ge0$, and
block $d$ of $\Phi_n([A_{ij}])$ is $\sum_c(\Phi_{dc})_n(X^c)\ge0$. Conversely, for $[X_{ij}]\ge0$ in
$M_n(M_{k_c})$, $[\iota_c(X_{ij})]\ge0$ in $M_n(\mathcal M)$ (its block $c$ is $[X_{ij}]$, the others
vanish), and block $d$ of its image is $[\Phi_{dc}(X_{ij})]$. The Kraus form follows from Theorem 6.2
for each pair $(d,c)$, padded with zero operators to a common $r$. Finally, by cyclicity,
$\operatorname{Tr}\Phi(A)=\sum_c\operatorname{Tr}\big(A_c\sum_{d,l}K^{dc*}_lK^{dc}_l\big)$, which equals
$\sum_c\operatorname{Tr}A_c$ for all $A$ iff each sum is $1_{k_c}$ (take matrix units:
$\operatorname{Tr}(E_{ij}X)=X_{ji}$). $\square$

The library has two classes. `Channel` stores **sector-preserving** maps ($C'=C$, $k'=k$,
$\Phi_{dc}=0$ for $d\ne c$), whose Kraus operators $K_l=\bigoplus_cK^{cc}_l$ lie in $\mathcal M$, as a
tensor of shape `(batch, r, C, k_max, k_max)`; `InterSectorChannel` stores the general map of
Proposition 6.4 with shape `(batch, r, C', C, k'_max, k_max)`. `Channel.choi()` returns
$J(\Phi_{cc})$ for every sector and `Channel.superoperator()` the matrix $S_c$ of $\Phi_{cc}$ on
row-major vectors, $\operatorname{vec}(X)_{(a,b)}=X_{ab}$. Both contain the numbers
$\Phi(E_{ij})_{ab}=J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$ in a different order. `Channel.from_superoperator`
diagonalises the Choi matrices as in the proof of Theorem 6.2 and raises `ValueError` if a Choi
matrix has an eigenvalue below $-10^{-10}\max(1,\|J\|)$ (in double precision): by Theorem 6.2 it
refuses the maps that are not CP, up to this rounding tolerance.

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, states, channels
from torch_vn_algebra.channels import Channel, InterSectorChannel

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def unit(i, j, k):
    X = torch.zeros(k, k, dtype=dt)
    X[i, j] = 1
    return X

def choi(f, k):
    """Choi matrix sum_ij E_ij (x) f(E_ij) of a linear map f on k x k matrices."""
    return sum(torch.kron(unit(i, j, k), f(unit(i, j, k)).contiguous()) for i in range(k) for j in range(k))

def to_superoperator(J, k):
    """S[(a,b),(i,j)] = J[(i,a),(j,b)], with a leading batch axis."""
    return J.reshape(k, k, k, k).permute(1, 3, 0, 2).reshape(1, k * k, k * k)

ev = torch.linalg.eigvalsh(choi(lambda X: X.T, 2))
print(ev)                                                 # the Choi matrix of the transpose is the swap
assert abs(ev[0] + 1) < 1e-12                             # Proposition 6.1

alg = TypeIAlgebra([2, 3], [2, 3], precision='double')    # M_2 (+) M_3
Phi = channels.random_channel(alg, kraus_rank=2)          # Kraus tensor of shape (1, 2, 2, 3, 3)

def block(Phi, c):
    """The block Phi_cc of a sector-preserving map, as a function on k_c x k_c matrices."""
    ks = Phi.algebra.k_factors
    def f(X):
        blocks = [torch.zeros(k, k, dtype=dt) for k in ks]
        blocks[c] = X
        return Phi(Phi.algebra.from_blocks(blocks)).matrix[0, c, :ks[c], :ks[c]]
    return f

J, S = Phi.choi(), Phi.superoperator()                    # one matrix per sector, (1, k_c^2, k_c^2)
rho = states.random_density_matrix(alg)
out = Phi(rho)
for c, k in enumerate(alg.k_factors):
    assert close(J[c][0], choi(block(Phi, c), k))         # the convention of Theorem 6.2
    w = torch.linalg.eigvalsh(J[c][0])
    assert w.min() > -1e-12 and int((w > 1e-10).sum()) == 2      # J >= 0, Kraus rank 2 (Corollary 6.3)
    assert close(S[c][0], to_superoperator(J[c][0], k)[0])
    assert close(S[c][0] @ rho.matrix[0, c, :k, :k].reshape(-1), out.matrix[0, c, :k, :k].reshape(-1))
assert Phi.is_trace_preserving()
assert close(Channel.from_superoperator(alg, S)(rho), out)       # back from the superoperator
try:
    Channel.from_superoperator(alg, [to_superoperator(choi(lambda X: X.T, k), k) for k in (2, 3)])
    raise AssertionError("the transpose was accepted")
except ValueError as e:
    print(e)                                              # superoperator is not completely positive
```

> **Physical meaning** (physical interpretation, not a theorem of this course). A system is never alone: it may be correlated with an untouched ancilla
> with $n$ levels. A physical process acts as $\Phi_n=\mathrm{id}_n\otimes\Phi$ on the pair and must
> map states to states, so $\Phi$ must be completely positive. The transpose is positive, but it maps
> the entangled state $|\Omega\rangle\langle\Omega|/2$ to an operator with eigenvalue $-1/2$. The
> Choi matrix is the image of the maximally entangled state: one state decides complete positivity.

## 6.2 Stinespring dilation and the freedom in the Kraus operators

For $\Phi(A)=\sum_{l=1}^rK_lAK_l^*$ from $M_k$ to $M_{k'}$ put $\Phi^*(B)=\sum_lK_l^*BK_l$; by
cyclicity $\operatorname{Tr}(\Phi(A)B)=\operatorname{Tr}(A\Phi^*(B))$, so $\Phi^*$ is the dual of
$\Phi$ for the trace (Section 6.3) and does not depend on the choice of the $K_l$ (if
$\operatorname{Tr}(AX)=\operatorname{Tr}(AX')$ for all $A$, then $X=X'$: take $A=E_{ji}$). On
$\mathbb C^r\otimes\mathbb C^{k'}$ (the **environment** $\mathbb C^r$ first) let
$\operatorname{Tr}_EX=\sum_l(e_l^*\otimes1)X(e_l\otimes1)$.

**Theorem 6.5 (Stinespring, finite dimensions).** Let
$V=\sum_le_l\otimes K_l:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}$. Then
$$
\Phi(A)=\operatorname{Tr}_E(VAV^*),\qquad\Phi^*(B)=V^*(1_r\otimes B)V,
$$
and $\Phi$ is TP iff $V$ is an isometry. Conversely, every $V$ defines a CP map by the first formula.

*Proof.* $(e_l^*\otimes1)V=K_l$, so $\operatorname{Tr}_E(VAV^*)=\sum_lK_lAK_l^*$; conversely
$\operatorname{Tr}_E(VAV^*)$ is a Kraus form with $K_l=(e_l^*\otimes1)V$. Next
$V^*(1\otimes B)V=\sum_{l,l'}\langle e_l,e_{l'}\rangle K_l^*BK_{l'}=\Phi^*(B)$, and
$V^*V=\sum_lK_l^*K_l$, which is $1$ iff $\Phi$ is TP (Proposition 6.4). $\square$

The formula $\Phi^*(B)=V^*\pi(B)V$ with the representation $\pi(B)=1_r\otimes B$ is the form in
which W. F. Stinespring (Proc. AMS 6 (1955) 211) stated the theorem, in a far more general setting
than the finite-dimensional one proved here (historical attribution; further reading: V. Paulsen,
*Completely Bounded Maps and Operator Algebras*, Cambridge 2002, Ch. 4). For a map between direct
sums, let $\iota'_d$ be the inclusion of the $d$-th summand of $\bigoplus_d\mathbb C^{k'_d}$ and
$V_c=\sum_{d,l}e_l\otimes\iota'_dK^{dc}_l:\mathbb C^{k_c}\to\mathbb C^r\otimes\bigoplus_d\mathbb C^{k'_d}$.
Then $(e_l^*\otimes1)V_c=\sum_d\iota'_dK^{dc}_l$, so the $(d,d)$ block of
$\sum_c\operatorname{Tr}_E(V_cA_cV_c^*)=\sum_{c,l}\big(\sum_d\iota'_dK^{dc}_l\big)A_c\big(\sum_{d'}\iota'_{d'}K^{d'c}_l\big)^*$
is $\sum_{c,l}K^{dc}_lA_cK^{dc*}_l=\Phi(A)_d$; and, as $\iota'^*_d\iota'_{d'}=\delta_{dd'}1$,
$V_c^*V_c=\sum_{d,l}K^{dc*}_lK^{dc}_l$, so $\Phi$ is TP iff every $V_c$ is an isometry
(Proposition 6.4). This is how `random_channel` and `random_inter_sector_channel` draw channels
(Section 6.5).

**Lemma 6.6.** If $X,Y\in M_{p\times n}$ and $XX^*=YY^*$, then $Y=XW$ for a unitary $W\in M_n$.

*Proof.* $\|X^*z\|^2=\langle z,XX^*z\rangle=\|Y^*z\|^2$ for all $z$, so $X^*z\mapsto Y^*z$ is a
well-defined linear isometry from the range of $X^*$ onto the range of $Y^*$. As in the proof of
Theorem 1.8 it extends to a unitary $U$ with $UX^*=Y^*$; take $W=U^*$. $\square$

**Theorem 6.7 (unitary freedom).** Let $K_1,\dots,K_n$ and $L_1,\dots,L_n$ be in $M_{k'\times k}$
(pad the shorter family with zeros). Then $\sum_lK_lAK_l^*=\sum_jL_jAL_j^*$ for all $A$ iff
$L_j=\sum_lu_{jl}K_l$ for a unitary $u\in M_n$.

*Proof.* If $L_j=\sum_lu_{jl}K_l$, then
$\sum_jL_jAL_j^*=\sum_{l,l'}(u^*u)_{l'l}K_lAK_{l'}^*=\sum_lK_lAK_l^*$. Conversely, equal maps have
equal Choi matrices, $\sum_lv_lv_l^*=\sum_jw_jw_j^*$ with $v_l=\sum_ie_i\otimes K_le_i$ and
$w_j=\sum_ie_i\otimes L_je_i$. Lemma 6.6, applied to the matrices with columns $v_l$ and $w_j$, gives
$w_j=\sum_lW_{lj}v_l$; as $K\mapsto\sum_ie_i\otimes Ke_i$ is linear and injective,
$L_j=\sum_lW_{lj}K_l$, and $u=W^T$ is unitary. $\square$

In the code the environment index of $V$ is the first index of the Kraus tensor:

```python
K = Phi.kraus[0, :, 1, :3, :3]                            # the two Kraus operators of block c = 1
V = K.reshape(2 * 3, 3)                                   # V = sum_l e_l (x) K_l : C^3 -> C^2 (x) C^3
assert close(V.conj().T @ V, torch.eye(3, dtype=dt))      # TP <=> V is an isometry
r1 = rho.matrix[0, 1, :3, :3]
TrE = torch.einsum('lalb->ab', (V @ r1 @ V.conj().T).reshape(2, 3, 2, 3))
assert close(TrE, out.matrix[0, 1, :3, :3])               # Phi(rho) = Tr_E(V rho V^*)
B1 = torch.randn(3, 3, dtype=dt)
Bop = alg.from_blocks([torch.zeros(2, 2, dtype=dt), B1])
assert close(V.conj().T @ torch.kron(torch.eye(2, dtype=dt), B1) @ V,
             Phi.adjoint()(Bop).matrix[0, 1, :3, :3])    # Phi^*(B) = V^*(1 (x) B)V

u = alg.random_unitary(3)                                 # Theorem 6.7 with n = 3: pad with a zero operator
Kpad = torch.cat([Phi.kraus, torch.zeros_like(Phi.kraus[:, :1])], dim=1)
PhiL = Channel(alg, torch.einsum('jl,blcxy->bjcxy', u, Kpad))
assert close(PhiL(rho), out) and all(close(a, b) for a, b in zip(PhiL.choi(), J))
```

> **Physical meaning** (physical interpretation, not a theorem of this course). Every channel is an
> isometric interaction of the system with an environment
> in a fixed state, followed by ignoring the environment. Measuring the environment in the basis
> $e_l$ gives outcome $l$ with probability $\operatorname{Tr}K_l\rho K_l^*$ and leaves the state
> $K_l\rho K_l^*/\operatorname{Tr}(K_l\rho K_l^*)$: the Kraus operators are the possible "quantum jumps".
> Another basis gives other jumps but the same channel (Theorem 6.7); the smallest environment has
> dimension $\operatorname{rank}J$ (Theorem 6.5 and Corollary 6.3).

## 6.3 Duals with respect to a trace

**Definition.** Let $\operatorname{tr}_w$, $\operatorname{tr}_v$ be faithful traces on $\mathcal M$,
$\mathcal N$. The **dual** of a linear map $\Phi:\mathcal M\to\mathcal N$ with respect to $(w,v)$ is
the linear map $\Phi^\sharp:\mathcal N\to\mathcal M$ with
$$
\operatorname{tr}_v\big(\Phi(x)\,y\big)=\operatorname{tr}_w\big(x\,\Phi^\sharp(y)\big)\qquad(x\in\mathcal M,\ y\in\mathcal N).
$$
It exists and is unique: for fixed $y$, $x\mapsto\operatorname{tr}_v(\Phi(x)y)$ is a linear functional,
so by Theorem 3.2 it equals $\operatorname{tr}_w(\rho x)=\operatorname{tr}_w(x\rho)$ for exactly one
$\rho=:\Phi^\sharp(y)$, and uniqueness makes $\Phi^\sharp$ linear. For blunt traces we write $\Phi^*$.
Let $z_w=\sum_cw_cP_c$ and $z_v=\sum_dv_dQ_d$ (central and invertible).

**Proposition 6.8.** If $\Phi$ has Kraus operators $K^{dc}_l$ (Proposition 6.4), then
$$
\Phi^\sharp(y)_c=\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_l\,y_d\,K^{dc}_l,\qquad\text{i.e.}\qquad
\Phi^\sharp(y)=z_w^{-1}\,\Phi^*(z_v\,y).
$$
In particular $\Phi^\sharp$ is CP, with Kraus operators $(v_d/w_c)^{1/2}K^{dc*}_l$.

*Proof.* By cyclicity of the trace in each block,
$$
\operatorname{tr}_v(\Phi(x)y)=\sum_dv_d\sum_{c,l}\operatorname{Tr}\big(K^{dc}_lx_cK^{dc*}_ly_d\big)
=\sum_cw_c\operatorname{Tr}\Big(x_c\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_ly_dK^{dc}_l\Big),
$$
and uniqueness gives the formula. For $v=w=1$ it reads $\Phi^*(y)_c=\sum_{d,l}K^{dc*}_ly_dK^{dc}_l$;
the factors $v_d$ and $1/w_c$ are multiplication by $z_v$ and $z_w^{-1}$. $\square$

$\Phi$ is **TP with respect to $(w,v)$** if $\operatorname{tr}_v\Phi(x)=\operatorname{tr}_w(x)$ for all $x$.

**Corollary 6.9.** (a) $\Phi$ is TP with respect to $(w,v)$ iff $\Phi^\sharp(1)=1$, i.e.
$\sum_{d,l}(v_d/w_c)K^{dc*}_lK^{dc}_l=1_{k_c}$ for every $c$. (b) $\Phi$ is unital iff $\Phi^\sharp$ is
TP with respect to $(v,w)$. (c) The dual of $\Phi^\sharp$ with respect to $(v,w)$ is $\Phi$. (d) If
$\Phi:\mathcal M\to\mathcal M$ is sector preserving and $v=w$, then $\Phi^\sharp=\Phi^*$ for every $w$,
and $\Phi$ is TP for one trace iff it is TP for all.

*Proof.* (a) $\operatorname{tr}_v\Phi(x)=\operatorname{tr}_w(x\Phi^\sharp(1))$, and
$\operatorname{tr}_w(x(\Phi^\sharp(1)-1))=0$ for all $x$ iff $\Phi^\sharp(1)=1$ (take
$x=(\Phi^\sharp(1)-1)^*$; $\operatorname{tr}_w$ is faithful); then use Proposition 6.8. (b)
$\operatorname{tr}_w\Phi^\sharp(y)=\operatorname{tr}_v(\Phi(1)y)$ equals $\operatorname{tr}_v(y)$ for all
$y$ iff $\Phi(1)=1$. (c) $\operatorname{tr}_w(\Phi^\sharp(y)x)=\operatorname{tr}_v(y\Phi(x))$ by the trace
property; use uniqueness. (d) Only the pairs $d=c$ occur, where $v_d/w_c=1$. $\square$

For maps that change the sector, (d) fails. Let $\mathcal M=\mathcal N=M_1\oplus M_2$ (sector $0$: no
particle; sector $1$: one particle on two sites) and let $\Phi$ be **particle loss**,
$\Phi(x)=(x_0+\operatorname{Tr}x_1)\oplus0$, with Kraus operators $K^{00}=1$, $K^{01}_l=\langle l|$.
It is TP for the blunt traces. For `Tr_norm`, $w=v=(1,\tfrac12)$, $\Phi^\sharp(1)=1\oplus2\cdot1_2$
and $\operatorname{tr}_v\Phi(x)=x_0+\operatorname{Tr}x_1\ne x_0+\tfrac12\operatorname{Tr}x_1=\operatorname{tr}_w(x)$.
In the library, `Channel.adjoint()` is $\Phi^*$, and `InterSectorChannel.adjoint(trace=...)` is
$\Phi^\sharp$ for the same named trace on both algebras (`'blunt'`, `'norm'`, `'tau_vN'`; the weights
$1/(Ck_c)$ and $1/(C'k'_d)$ of `tau_vN` also see the numbers of sectors). We check Proposition 6.8 for
a random map $M_2\oplus M_3\to M_1\oplus M_2\oplus M_2$ and the three traces:

```python
alg_out = TypeIAlgebra([1, 2, 2], [1, 2, 2], precision='double')
Psi = channels.random_inter_sector_channel(alg, alg_out, kraus_rank=2)     # Kraus tensor (1, 2, 3, 2, 2, 3)
x = alg.operator(torch.randn(1, 2, 3, 3, dtype=dt))
y = alg_out.operator(torch.randn(1, 3, 2, 2, dtype=dt))
weights = {'blunt': lambda k, C: 1.0, 'norm': lambda k, C: 1 / k, 'tau_vN': lambda k, C: 1 / (C * k)}
trace = {'blunt': lambda X: X.Tr_blunt(), 'norm': lambda X: X.Tr_norm(), 'tau_vN': lambda X: X.tau_vN()}

def z(a, t):                                              # the central element sum_c w_c P_c of trace t
    return a.central(torch.tensor([weights[t](k, a.C) for k in a.k_factors], dtype=torch.float64))

for t in trace:
    dual = Psi.adjoint(trace=t)
    assert (trace[t](Psi(x) @ y) - trace[t](x @ dual(y))).abs().max() < 1e-10   # the definition
    assert close(dual(y), z(alg, t).inverse() @ Psi.adjoint()(z(alg_out, t) @ y))   # Proposition 6.8
dev = [(Psi.adjoint(trace=t)(alg_out.identity()) - alg.identity()).matrix.abs().max().item() for t in trace]
print(dev)                         # TP only for the blunt traces: the other duals are not unital
assert Psi.is_trace_preserving() and dev[0] < 1e-12 and min(dev[1:]) > 0.1

alg01 = TypeIAlgebra([1, 2], [1, 2], precision='double')
loss = InterSectorChannel.from_blocks(alg01, alg01, {(0, 0): [torch.ones(1, 1)],
                                                     (0, 1): [torch.tensor([[1.0, 0.0]]), torch.tensor([[0.0, 1.0]])]})
assert loss.is_trace_preserving()
assert close(loss.adjoint(trace='norm')(alg01.identity()), alg01.central([1.0, 2.0]))
s = states.random_density_matrix(alg01)
assert abs(loss(s).Tr_norm() - (s @ alg01.central([1.0, 2.0])).Tr_norm()) < 1e-12   # Corollary 6.9(a)
print(s.Tr_norm().real.item(), loss(s).Tr_norm().real.item())      # Tr_norm is not conserved
```

> **Physical meaning** (physical interpretation, not a theorem of this course). The trace-independent object is the Heisenberg map
> $\Psi:\mathcal N\to\mathcal M$ on observables: a state $\omega$ before the process becomes
> $\omega\circ\Psi$ after it, and $\Psi$ is CP and unital. Densities depend on the trace
> (Theorem 3.2); the map sending the density of $\omega$ to that of $\omega\circ\Psi$ is the dual of
> $\Psi$ for the chosen traces, which is TP for them by Corollary 6.9(b). If the algebras act with
> multiplicities, the physical trace has weights $w_c=m_c$ (Chapter 2). Particle loss above conserves
> probability for the blunt trace, which counts states, not for `Tr_norm`, which weighs sectors equally.

## 6.4 Decoherence

The library's standard channels act inside each block ($p,\gamma\in[0,1]$):

| channel | $\Phi(\rho)_c$ | Kraus operators in block $c$ |
|---|---|---|
| `dephasing_channel` | $(1-p)\rho_c+p\operatorname{diag}\rho_c$ | $\sqrt{1-p}\,1$, $\sqrt p\,E_{ii}$ |
| `depolarizing_channel` | $(1-p)\rho_c+p\operatorname{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1-p}\,1$, $\sqrt{p/k_c}\,E_{ij}$ |
| `center_expectation` | $\operatorname{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1/k_c}\,E_{ij}$ |
| `amplitude_damping_channel` | levels $i>0$ decay to $0$ | $E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}$, $\sqrt\gamma\,E_{0i}$ |
| `random_mixed_unitary_channel` | $\sum_jp_jU_j\rho U_j^*$ | $\sqrt{p_j}\,U_j$ |

The Kraus forms are checked by direct computation, e.g.
$\sum_{ij}\frac1kE_{ij}\rho E_{ji}=\frac1k\sum_{ij}\rho_{jj}E_{ii}$. All five are TP, and all but
amplitude damping with $\gamma>0$ are unital (Proposition 6.4: compute $\sum_lK_l^*K_l$ and
$\sum_lK_lK_l^*$; for amplitude damping in a block of size $k$ they are $1$ and
$1+(k-1)\gamma E_{00}-\gamma\sum_{i>0}E_{ii}$, and $\sum_jp_j=1$ for the mixed unitary channel).
Explicitly, amplitude damping maps $\rho_{00}\mapsto\rho_{00}+\gamma\sum_{i>0}\rho_{ii}$,
$\rho_{ii}\mapsto(1-\gamma)\rho_{ii}$ and $\rho_{ij}\mapsto(1-\gamma)\rho_{ij}$ for $i\ne j$ both
positive, and $\rho_{0j}\mapsto\sqrt{1-\gamma}\,\rho_{0j}$, $\rho_{j0}\mapsto\sqrt{1-\gamma}\,\rho_{j0}$ for $j>0$. The conditional expectation $E$ onto the centre
(Proposition 2.8) is depolarization with $p=1$. Unlike the measurement of the sector,
$X\mapsto\sum_cP_cXP_c$, which is the identity on $\mathcal M$ (Section 3.3), it keeps only $p_c$.

**Theorem 6.10 (fixed points of unital channels).** Let $\Phi(X)=\sum_lK_lXK_l^*$ with $K_l\in\mathcal M$
be TP and unital. Then $\Phi(X)=X$ iff $X$ commutes with every $K_l$ and $K_l^*$. In particular the
fixed points form a von Neumann subalgebra of $\mathcal M$.

*Proof.* If $X$ commutes with all $K_l$, then $\Phi(X)=X\sum_lK_lK_l^*=X$. Conversely, $\Phi(X^*)=\Phi(X)^*$,
so the self-adjoint parts $(X+X^*)/2$ and $(X-X^*)/2i$ of a fixed point are fixed points; let
$X=X^*$. With $Y_l=K_lX-XK_l$,
$$
\sum_lY_lY_l^*=\sum_l\big(K_lX^2K_l^*-K_lXK_l^*X-XK_lXK_l^*+XK_lK_l^*X\big)=\Phi(X^2)-X^2,
$$
whose trace is $0$ because $\Phi$ is TP. As $\operatorname{Tr}Y_lY_l^*=\|Y_l\|_2^2\ge0$, every $Y_l=0$;
taking adjoints, $X$ commutes with $K_l^*$ too. These $X$ form the commutant of a self-adjoint set,
intersected with $\mathcal M$: a $*$-algebra containing $1$, hence a von Neumann algebra
(Theorem 2.2). $\square$

So dephasing ($p>0$) keeps exactly the diagonal of every block, and depolarization ($p>0$) and $E$
keep exactly the centre. Without unitality the theorem fails: amplitude damping fixes $E_{00}$,
which does not commute with the jump $E_{01}$. For $0<\gamma\le1$ and $k\ge2$ its fixed points in a block
and the commutant of its Kraus operators both have dimension $1$, but they are $\mathbb CE_{00}$ and
$\mathbb C1$: by the explicit form above, a fixed point has $\rho_{ij}=0$ for $i\ne j$ and
$\rho_{ii}=0$ for $i>0$ (these entries are multiplied by $\sqrt{1-\gamma}\ne1$ or $1-\gamma\ne1$), and every
multiple of $E_{00}$ is fixed; and anything commuting with all $E_{0i}$, $E_{i0}$ ($i>0$) commutes
with their products $E_{ij}=E_{i0}E_{0j}$ and $E_{00}=E_{01}E_{10}$, hence with $M_k$, so it is a multiple
of $1$ (Chapter 2, Exercise 1).

**Theorem 6.11 (unital channels increase entropy).** Let $\Phi:\mathcal M\to\mathcal N$ be positive,
TP and unital, and let $\rho$ and $\Phi(\rho)$ be states with eigenvalues $\lambda_i$ and $\mu_j$.
Then $\mu=\Delta\lambda$ for a doubly stochastic matrix $\Delta$, so $\sum_jf(\mu_j)\ge\sum_if(\lambda_i)$ for
every concave $f$ on $[0,1]$. In particular $S(\Phi(\rho))\ge S(\rho)$ and
$\operatorname{Tr}\Phi(\rho)^2\le\operatorname{Tr}\rho^2$.

*Proof.* Diagonalise each block: $\rho=\sum_i\lambda_i|e_i\rangle\langle e_i|$ and
$\Phi(\rho)=\sum_j\mu_j|f_j\rangle\langle f_j|$ with orthonormal bases of $\bigoplus_c\mathbb C^{k_c}$
and $\bigoplus_d\mathbb C^{k'_d}$ whose vectors lie in one summand each, so that
$|e_i\rangle\langle e_i|\in\mathcal M$. Then $\mu_j=\langle f_j,\Phi(\rho)f_j\rangle=\sum_i\Delta_{ji}\lambda_i$
with $\Delta_{ji}=\langle f_j,\Phi(|e_i\rangle\langle e_i|)f_j\rangle\ge0$,
$\sum_j\Delta_{ji}=\operatorname{Tr}\Phi(|e_i\rangle\langle e_i|)=1$ and
$\sum_i\Delta_{ji}=\langle f_j,\Phi(1)f_j\rangle=1$ (so $\Delta$ is square: sum all entries). By the
finite Jensen inequality $f(\sum_i\Delta_{ji}\lambda_i)\ge\sum_i\Delta_{ji}f(\lambda_i)$ (induction on
the number of terms from the definition of concavity, as in the solution of Exercise 2 in Chapter 3),
$\sum_jf(\mu_j)\ge\sum_j\sum_i\Delta_{ji}f(\lambda_i)=\sum_if(\lambda_i)$. Take $f(t)=-t\log t$ and
$f(t)=-t^2$. $\square$

The classical part of a map between sectors is its **transition matrix**
$T_{dc}(\rho)=\operatorname{Tr}\Phi_{dc}(\rho_c)/\operatorname{Tr}\rho_c$ ($\rho_c\ne0$); the library's
`transition_matrix()` is $T=T(1)$, for the maximally mixed state of every input sector.

**Proposition 6.12.** Let $\Phi:\mathcal M\to\mathcal N$ be positive and TP, $\rho$ a state with
$p_c=\operatorname{Tr}\rho_c>0$, and $q_d=\operatorname{Tr}\Phi(\rho)_d$. Then $T(\rho)$ is column
stochastic and $q=T(\rho)p$. With the conditional expectations $E_{\mathcal M}$, $E_{\mathcal N}$ onto
the centres,
$$
E_{\mathcal N}\Phi E_{\mathcal M}(\rho)=\bigoplus_d\frac{(Tp)_d}{k'_d}\,1_{k'_d},\qquad T=T(1).
$$

*Proof.* $T_{dc}\ge0$ as $\Phi_{dc}$ is positive, and
$\sum_dT_{dc}=\operatorname{Tr}\Phi(\iota_c\rho_c)/\operatorname{Tr}\rho_c=1$. Next
$q_d=\sum_c\operatorname{Tr}\Phi_{dc}(\rho_c)=\sum_cT_{dc}(\rho)p_c$. Finally
$E_{\mathcal M}\rho=\bigoplus_cp_c1_{k_c}/k_c$, so the trace of
$\Phi(E_{\mathcal M}\rho)_d=\sum_cp_c\Phi_{dc}(1_{k_c}/k_c)$ is $\sum_cT_{dc}p_c$, and $E_{\mathcal N}$
replaces block $d$ by its trace times $1/k'_d$. $\square$

```python
def fixed_dims(Phi):
    """dim {X in M_{k_c} : Phi_cc(X) = X} for every sector: the null space of S_c - 1."""
    return [int((torch.linalg.svdvals(S_c[0] - torch.eye(k * k, dtype=dt)) < 1e-9).sum())
            for S_c, k in zip(Phi.superoperator(), Phi.algebra.k_factors)]

def commutant_dims(Phi):
    """dim of the commutant of {K_l, K_l^*} in every block (row-major vec, as in Chapter 2)."""
    out = []
    for c, k in enumerate(Phi.algebra.k_factors):
        I, Ks = torch.eye(k, dtype=dt), Phi.kraus[0, :, c, :k, :k]
        ops = [A.contiguous() for A in torch.cat([Ks, Ks.mH])]
        L = torch.cat([torch.kron(A, I) - torch.kron(I, A.T.contiguous()) for A in ops])
        out.append(k * k - int((torch.linalg.svdvals(L) > 1e-9).sum()))
    return out

zoo = {'dephasing': channels.dephasing_channel(alg, 0.3),
       'depolarizing': channels.depolarizing_channel(alg, 0.3),
       'centre': channels.center_expectation(alg),
       'mixed unitary': channels.random_mixed_unitary_channel(alg, 2),
       'amplitude damping': channels.amplitude_damping_channel(alg, 0.4)}
for name, ch in zoo.items():
    print(f"{name:18s} unital {ch.is_unital()!s:5s} fixed {fixed_dims(ch)} commutant {commutant_dims(ch)}")
    assert ch.is_trace_preserving()
    if ch.is_unital():
        assert fixed_dims(ch) == commutant_dims(ch)                         # Theorem 6.10
assert fixed_dims(zoo['dephasing']) == [2, 3] and fixed_dims(zoo['depolarizing']) == [1, 1]
E00 = alg.from_blocks([unit(0, 0, 2), unit(0, 0, 3)])
E01 = alg.from_blocks([unit(0, 1, 2), unit(0, 1, 3)])
assert close(zoo['amplitude damping'](E00), E00) and not close(E00 @ E01, E01 @ E00)   # equal dims, different spaces
assert fixed_dims(zoo['amplitude damping']) == commutant_dims(zoo['amplitude damping']) == [1, 1]

rhos = states.random_density_matrix(alg, batch_size=100)
for name in ('dephasing', 'depolarizing', 'centre', 'mixed unitary'):                 # Theorem 6.11
    assert (zoo[name](rhos).von_neumann_entropy() >= rhos.von_neumann_entropy() - 1e-12).all()
    assert (zoo[name](rhos).purity() <= rhos.purity() + 1e-12).all()
mm = states.maximally_mixed_state(alg)
S_after = channels.amplitude_damping_channel(alg, 1.0)(mm).von_neumann_entropy().item()
print(mm.von_neumann_entropy().item(), S_after)          # not unital: the entropy drops
assert abs(S_after + 0.4 * math.log(0.4) + 0.6 * math.log(0.6)) < 1e-12

T = Psi.transition_matrix()[0]                            # (C', C) = (3, 2)
p, q = rho.sector_probabilities()[0], Psi(rho).sector_probabilities()[0]
assert torch.allclose(T.sum(0), torch.ones(2, dtype=T.dtype))                       # Proposition 6.12
assert torch.allclose(q, Psi.transition_matrix(rho)[0] @ p)
classical = channels.center_expectation(alg_out) @ Psi @ channels.center_expectation(alg)
Tp = T @ p
assert close(classical(rho), alg_out.central([Tp[d] / k for d, k in enumerate(alg_out.k_factors)]))
assert torch.allclose(loss.transition_matrix()[0], torch.tensor([[1.0, 1.0], [0.0, 0.0]], dtype=T.dtype))
```

Under complete amplitude damping the entropy of the maximally mixed state of $M_2\oplus M_3$ drops
from $\log5\approx1.609$ to $H(\tfrac25,\tfrac35)\approx0.673$: a non-unital channel can purify.

**Qubit coherences.** On $M_2$ let $D_p$ be dephasing and $A_\gamma$ amplitude damping. By the table
and the explicit form above, $D_p$ multiplies $\rho_{01}$ and $\rho_{10}$ by $1-p$ and fixes the
diagonal, while $A_\gamma$ multiplies $\rho_{01}$ and $\rho_{10}$ by $\sqrt{1-\gamma}$ and maps the
diagonal to $(\rho_{00}+\gamma\rho_{11},(1-\gamma)\rho_{11})$, independently of the coherences. Hence
$D_pA_\gamma=A_\gamma D_p$: both multiply the coherences by $(1-p)\sqrt{1-\gamma}$ and act on the diagonal
as $A_\gamma$. In particular, for $T_1,T_\varphi>0$, $t\ge0$, $1-p=e^{-t/T_\varphi}$ and $1-\gamma=e^{-t/T_1}$,
the composition multiplies the coherences by $e^{-t/T_2}$ with
$$
\frac1{T_2}=\frac1{2T_1}+\frac1{T_\varphi},
$$
and these families compose additively in $t$: $D_{p(t)}D_{p(s)}=D_{p(t+s)}$ and
$A_{\gamma(t)}A_{\gamma(s)}=A_{\gamma(t+s)}$, since the factors $e^{-t/T}$ multiply and the trace is preserved.

```python
q2 = TypeIAlgebra([2], [2], precision='double')
p_, g_ = 0.3, 0.4
Dp, Ag = channels.dephasing_channel(q2, p_), channels.amplitude_damping_channel(q2, g_)
r2 = states.random_density_matrix(q2)
assert close(Dp(Ag(r2)), Ag(Dp(r2)))                                         # D_p A_gamma = A_gamma D_p
coh = Dp(Ag(r2)).matrix[0, 0, 0, 1] / r2.matrix[0, 0, 0, 1]
assert abs(coh - (1 - p_) * math.sqrt(1 - g_)) < 1e-12                       # coherence factor
T1, Tphi, t1, t2 = 2.0, 3.0, 0.7, 1.1
gam = lambda t: 1 - math.exp(-t / T1)
assert close(channels.amplitude_damping_channel(q2, gam(t1))(channels.amplitude_damping_channel(q2, gam(t2))(r2)),
             channels.amplitude_damping_channel(q2, gam(t1 + t2))(r2))       # A_gamma(t) A_gamma(s) = A_gamma(t+s)
coh_t = math.exp(-t1 / Tphi) * math.sqrt(1 - gam(t1))
assert abs(coh_t - math.exp(-t1 * (1 / (2 * T1) + 1 / Tphi))) < 1e-15       # 1/T_2 = 1/(2T_1) + 1/T_phi
```

> **Physical meaning** (physical interpretation, not a theorem of this course). For a qubit,
> dephasing models pure dephasing with time $T_\varphi$ and amplitude damping the decay of the
> excited level with time $T_1$; the computation above then gives the coherence time $T_2$ of a qubit
> subject to both. Theorem 6.10 is the algebraic
> form of decoherence-free information: what commutes with all jump operators survives (Chapter 10).
> Proposition 6.12 says that, seen through the centre, a map between sectors is a classical Markov
> matrix on the charges.

## 6.5 Random channels

A **Haar-random isometry** $V:\mathbb C^k\to\mathbb C^N$ is, as in the library (`random_isometry`),
the factor $Q$ of $G=QR$, where $G$ is an $N\times k$ Ginibre matrix (independent entries with
density $\pi^{-1}e^{-|z|^2}$), $Q^*Q=1$, and $R$ is upper triangular with positive diagonal (on the
CPU the library takes a LAPACK factorisation $G=Q'R'$ and multiplies the $i$-th column of $Q'$ by
$R'_{ii}/|R'_{ii}|$, which gives this factorisation, as in Theorem 5.3; for large batches on a GPU it
uses Gram–Schmidt, which gives it directly).
`random_channel(alg, r)` draws such $V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}$
independently in every block and uses the Kraus operators of Theorem 6.5.

**Lemma 6.13.** $V$ has the law of the first $k$ columns of a Haar unitary of $U(N)$. In
particular, for every unitary $U\in M_N$, $UV$ has the same law as $V$.

*Proof.* Let $Z$ be an $N\times N$ Ginibre matrix with first $k$ columns $G$, and $E\in M_{N\times k}$
the first $k$ columns of $1_N$. By Theorem 5.3, $Z$ is invertible almost surely and $Z=Q_ZR_Z$ with
$Q_Z$ Haar distributed. As $R_Z$ is upper triangular, $G=ZE=(Q_ZE)(E^*R_ZE)$ is a factorisation of
the required type, and it is the only one: if $G=QR=Q'R'$, then $R^*R=G^*G=R'^*R'$, so $R'R^{-1}$
is unitary, upper triangular with positive diagonal, hence $1$ (proof of Lemma 5.2). So $V=Q_ZE$,
and $UV=(UQ_Z)E$ has the law of $V$ because $UQ_Z$ has that of $Q_Z$. $\square$

**Lemma 6.14.** Let $u$ be a random unit vector in $\mathbb C^N$ whose law is invariant under all
unitaries. Then $\mathbb E[u_a\bar u_b]=\delta_{ab}/N$ and
$$
\mathbb E\big[u_a\bar u_b\,u_c\bar u_d\big]=\frac{\delta_{ab}\delta_{cd}+\delta_{ad}\delta_{cb}}{N(N+1)} .
$$

*Proof.* Invariance under $\operatorname{diag}(e^{i\theta_1},\dots,e^{i\theta_N})$ multiplies the
moments by $e^{i(\theta_a-\theta_b)}$ and $e^{i(\theta_a-\theta_b+\theta_c-\theta_d)}$, so they vanish
unless $a=b$, resp. $\{a,c\}=\{b,d\}$ as multisets. By invariance under permutation matrices
$\mathbb E|u_a|^2=1/N$, and only $\alpha=\mathbb E|u_1|^4$ and $\beta=\mathbb E|u_1|^2|u_2|^2$ remain.
First, $\sum_c|u_c|^2=1$ gives $\alpha+(N-1)\beta=\mathbb E|u_1|^2=1/N$. Second, $(u_1+u_2)/\sqrt2$ is
the first coordinate of $Uu$ for a unitary $U$, so $\mathbb E|u_1+u_2|^4=4\alpha$. With
$s=|u_1|^2+|u_2|^2$ and $z=u_1\bar u_2$, $|u_1+u_2|^4=s^2+2s(z+\bar z)+(z+\bar z)^2$; the
invariances $u_1\mapsto-u_1$ and $u_1\mapsto iu_1$ give $\mathbb E[s(z+\bar z)]=0$ and
$\mathbb Ez^2=0$, so $4\alpha=(2\alpha+2\beta)+2\beta$, i.e. $\alpha=2\beta$. Hence
$\beta=1/(N(N+1))$ and $\alpha=2\beta$, which is the formula in all cases. (This used $N\ge2$; for
$N=1$, $|u_1|=1$ and both sides equal $1$.) $\square$

**Theorem 6.15 (output purity of a random channel).** Let $\psi$ be a pure state in a sector of size
$k$ and $\Phi$ drawn by `random_channel` with Kraus rank $r$. Then
$$
\mathbb E\operatorname{Tr}\Phi(\psi)^2=\frac{k+r}{kr+1},\qquad\mathbb E\,\Phi(\psi)=\frac1k1_k .
$$

*Proof.* In the sector of $\psi$, $\Phi(\psi)=\operatorname{Tr}_E|u\rangle\langle u|$ with $u=V\psi$, a
unit vector whose law is unitarily invariant by Lemma 6.13. Index $\mathbb C^r\otimes\mathbb C^k$ by
pairs $(e,x)$. Then $\Phi(\psi)_{xy}=\sum_eu_{ex}\bar u_{ey}$, whose mean is $r\delta_{xy}/(kr)$, and
$\operatorname{Tr}\Phi(\psi)^2=\sum_{x,y,e,f}u_{ex}\bar u_{ey}u_{fy}\bar u_{fx}$. By Lemma 6.14 with
$N=kr$, the first term forces $x=y$ ($kr^2$ index tuples) and the second $e=f$ ($k^2r$ tuples), so
the mean is $(kr^2+k^2r)/(kr(kr+1))$. $\square$

For $r=1$ the channel is unitary and the purity is $1$; as $r\to\infty$ it tends to $1/k$, the
purity of the maximally mixed state. The proof used only the unitary invariance of the law of $u$,
so it also gives Lubkin's mean purity of the reduced state of a random pure state $u=g/\|g\|$,
$g$ Gaussian (whose law is unitarily invariant by Theorem A.20(b)), which Chapter 7, Theorem 7.11
proves with Gaussian integrals; the numbers $\alpha,\beta$ of Lemma 6.14 are those of Corollary 5.11
for the entries of a Haar unitary, and Corollary A.21 recomputes them for $u=g/\|g\|$.

```python
B = 20_000
for r in (1, 2, 4):
    Phi_r = channels.random_channel(alg, kraus_rank=r, batch_size=B)
    for c, k in enumerate(alg.k_factors):
        psi = torch.randn(k, 1, dtype=dt)
        psi = psi / psi.norm()
        out_r = Phi_r(alg.from_blocks([psi @ psi.conj().T if d == c else torch.zeros(kk, kk, dtype=dt)
                                       for d, kk in enumerate(alg.k_factors)]))     # a batch of B outputs
        purity = (out_r @ out_r).trace.real
        mean, err = purity.mean().item(), purity.std().item() / math.sqrt(B)
        exact = (k + r) / (k * r + 1)
        print(f"k = {k}, r = {r}: E Tr Phi(psi)^2 = {mean:.4f} +- {err:.4f}   (exact {exact:.4f})")
        assert abs(mean - exact) < 5 * err + 1e-12
        avg = out_r.matrix.mean(0)[c, :k, :k]
        assert (avg - torch.eye(k, dtype=dt) / k).abs().max() < 0.02          # E Phi(psi) = 1/k
```

> **Physical meaning** (physical interpretation, not a theorem of this course). A random channel
> models generic, unstructured noise from an environment of
> dimension $r$. Already $r=4$ brings the purity of a pure qutrit state to $7/13\approx0.54$ on
> average (Theorem 6.15), against $1$ without noise and the minimum $1/3$ (by Cauchy–Schwarz,
> $\operatorname{Tr}\rho^2\ge(\operatorname{Tr}\rho)^2/k$). In a superselected system the noise acts in each sector
> separately, and the relevant dimension is that of the sector, not of the whole space.

## 6.6 What to compare with

QuTiP (`to_choi`, `to_kraus`, `to_super`, `to_stinespring`), Qiskit (`qiskit.quantum_info`: the
classes `Kraus`, `Choi`, `SuperOp`, `Stinespring` and `random_quantum_channel`; noise models such as
`amplitude_damping_error` in Qiskit Aer), toqito (`toqito.channel_ops`, `toqito.channel_props`) and
QETLAB (`ChoiMatrix`, `KrausOperators`, `ApplyMap`, `RandomSuperoperator`) convert between the
representations of Sections 6.1–6.2, test complete positivity and (most of them) draw random
channels, for maps $M_n\to M_m$ between full matrix algebras.
Their conventions differ (column- or row-stacking, the order of the factors of the Choi matrix), so
compare conventions before comparing numbers. A sector-preserving channel on $\mathcal M$ can be
treated there as a channel on $M_D$, $D=\sum_ck_c$, with superoperators of size $D^2\times D^2$
instead of $k_c^2\times k_c^2$ per sector. The library adds the block structure itself: Kraus
operators stored by sector and in batches, maps between different algebras (`InterSectorChannel`),
the duals for each of the traces of Chapter 3, and the transition matrix between sectors.

## Exercises

1. (Werner–Holevo) Show that $\Phi(\rho)=\big(\operatorname{Tr}(\rho)1-\rho^T\big)/(k-1)$ on $M_k$,
   $k\ge2$, is a unital channel with Kraus rank $k(k-1)/2$, although it is built from the
   transpose. Check it with `Channel.from_superoperator`.
2. Show that for $0<p\le1$ the dephasing channel on $M_k$ has Kraus rank $k$, one less than the
   number of Kraus operators the library stores.
3. Let $X$ be the cyclic shift on $\mathbb C^k$ and $Z=\operatorname{diag}(1,\omega,\dots,\omega^{k-1})$,
   $\omega=e^{2\pi i/k}$. Show that $\frac1{k^2}\sum_{a,b=0}^{k-1}X^aZ^b\rho Z^{-b}X^{-a}=\operatorname{Tr}(\rho)1/k$:
   the conditional expectation onto the centre of $M_k$ is a mixed-unitary channel.
4. Show that $(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp$ when the same trace is used on the
   middle algebra, and check it for the three traces of the library.

<details>
<summary>Solutions</summary>

1. With $F$ the swap, $J(\theta)=\sum_{ij}E_{ij}\otimes E_{ji}=F$ and $J(\rho\mapsto\operatorname{Tr}(\rho)1)=1\otimes1$,
   so $J(\Phi)=(1-F)/(k-1)=2P_-/(k-1)$, where $P_-=(1-F)/2$ is the projection onto the
   antisymmetric vectors ($F=F^*$, $F^2=1$), of rank $\operatorname{Tr}P_-=(k^2-\operatorname{Tr}F)/2$ with
   $\operatorname{Tr}F=\sum_{ij}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k$. So $J\ge0$ with rank $k(k-1)/2$ (Theorem 6.2,
   Corollary 6.3), and $\operatorname{Tr}\Phi(\rho)=\operatorname{Tr}\rho$, $\Phi(1)=1$ directly.

```python
k = 3
M3 = TypeIAlgebra([k], [k], precision='double')
J_WH = choi(lambda X: (X.trace() * torch.eye(k, dtype=dt) - X.T) / (k - 1), k)
w = torch.linalg.eigvalsh(J_WH)
assert w.min() > -1e-12 and int((w > 1e-10).sum()) == k * (k - 1) // 2
WH = Channel.from_superoperator(M3, [to_superoperator(J_WH, k)])
assert WH.is_trace_preserving() and WH.is_unital()
r3 = states.random_density_matrix(M3).matrix[0, 0]
assert close(WH(M3.operator(r3[None, None])).matrix[0, 0], (torch.eye(k, dtype=dt) - r3.T) / (k - 1))
```

2. $J=(1-p)|\Omega\rangle\langle\Omega|+p\sum_iE_{ii}\otimes E_{ii}$ vanishes outside the span of the
   $e_i\otimes e_i$, and on this $k$-dimensional space it is the matrix $(1-p)\mathbf 1\mathbf 1^T+p1$
   ($\mathbf 1$ the vector of ones), with eigenvalues $p$ and $(1-p)k+p$, all positive. So
   $\operatorname{rank}J=k$, while the library stores $1+k_{\max}$ operators.

```python
ranks = [int((torch.linalg.eigvalsh(J_c[0]) > 1e-10).sum()) for J_c in zoo['dephasing'].choi()]
print(ranks, zoo['dephasing'].kraus_rank)                 # [2, 3] versus 4 stored operators
assert ranks == alg.k_factors
```

3. $\operatorname{Tr}((X^aZ^b)^*X^{a'}Z^{b'})=k\,\delta_{aa'}\delta_{bb'}$: $Z^{-b}X^{a'-a}Z^{b'}$ has zero
   diagonal unless $a=a'$, and $\operatorname{Tr}Z^{b'-b}=\sum_j\omega^{j(b'-b)}=k\delta_{bb'}$. So the
   $B_{ab}=X^aZ^b/\sqrt k$ are an orthonormal basis of $M_k$ for $\operatorname{Tr}(A^*B)$, and any
   two orthonormal bases are related by a unitary; by Theorem 6.7, $\sum_mB_m\rho B_m^*$ is the
   same for all orthonormal bases. For the basis $E_{ab}$ it is
   $\sum_{ab}E_{ab}\rho E_{ba}=\operatorname{Tr}(\rho)1$. Divide by $k$.

```python
Xs = torch.roll(torch.eye(k, dtype=dt), 1, dims=0)
Zs = torch.diag(torch.exp(2j * math.pi * torch.arange(k, dtype=torch.float64) / k)).to(dt)
W = [torch.linalg.matrix_power(Xs, a) @ torch.linalg.matrix_power(Zs, b) / k for a in range(k) for b in range(k)]
twirl = channels.kraus_channel([M3.operator(Wab[None, None]) for Wab in W])
r3op = states.random_density_matrix(M3)
assert close(twirl(r3op), channels.center_expectation(M3)(r3op))
```

4. $\operatorname{tr}(\Psi(\Phi(x))z)=\operatorname{tr}(\Phi(x)\Psi^\sharp(z))=\operatorname{tr}(x\,\Phi^\sharp(\Psi^\sharp(z)))$
   for all $x$ and $z$; uniqueness of the dual gives the claim.

```python
Psi2 = channels.random_inter_sector_channel(alg_out, alg, kraus_rank=2)
for t in trace:
    assert close((Psi2 @ Psi).adjoint(trace=t)(x), Psi.adjoint(trace=t)(Psi2.adjoint(trace=t)(x)))
```

</details>

**Next:** [Chapter 7](07_tensor_products.md) treats tensor products and entanglement: the partial trace —
a channel between different algebras whose dual is the embedding $A\mapsto A\otimes1$ — Page's
formula, and the fusion of charges.
