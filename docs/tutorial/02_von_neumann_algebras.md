**English** · [Español](es/02_von_neumann_algebras.md) · [Français](fr/02_von_neumann_algebras.md) · [Deutsch](de/02_von_neumann_algebras.md) · [中文](zh/02_von_neumann_algebras.md) · [日本語](ja/02_von_neumann_algebras.md) · [Русский](ru/02_von_neumann_algebras.md)

# 2. Finite-dimensional von Neumann algebras

A von Neumann algebra is an algebra of operators that is closed under taking the
"commutant twice". In finite dimensions this is a purely algebraic condition, and
the main theorem of this chapter says that every such algebra is, up to a unitary
change of basis, a direct sum of full matrix algebras, each repeated with some
multiplicity. This is the structure the library is built on.

**You will learn** the double commutant theorem, the structure theorem with a
complete proof, how a symmetry group produces such an algebra, what the centre is,
and how the library represents all this.

Throughout, *$*$-algebra on $\mathbb C^n$* means a complex linear subspace
$\mathcal M\subseteq M_n(\mathbb C)$ closed under products and adjoints and containing
the identity $1$.

## 2.1 Commutants and the double commutant theorem

For a set $S\subseteq M_n(\mathbb C)$, the **commutant** is
$S'=\{X\in M_n(\mathbb C):XA=AX\text{ for all }A\in S\}$. It is always a unital algebra;
if $S$ is closed under adjoints ($S^*=S$), then $S'$ is a $*$-algebra, since
$XA=AX$ for all $A\in S$ implies $X^*A^*=A^*X^*$ for all $A\in S$. Clearly
$S\subseteq S''$ and $S'=S'''$.

**Definition.** A **von Neumann algebra** on $\mathbb C^n$ is a $*$-algebra
$\mathcal M$ with $\mathcal M=\mathcal M''$. Its **centre** is
$Z(\mathcal M)=\mathcal M\cap\mathcal M'$. A **factor** is a von Neumann algebra with
$Z(\mathcal M)=\mathbb C1$.

In infinite dimensions $\mathcal M=\mathcal M''$ is a topological condition (closure in
the weak operator topology). In finite dimensions it is automatic:

**Lemma 2.1.** Let $\mathcal M$ be a $*$-algebra on $\mathbb C^n$ and $K\subseteq\mathbb C^n$
a subspace with $\mathcal MK\subseteq K$. Then the orthogonal projection $P$ onto $K$ lies
in $\mathcal M'$.

*Proof.* For $A\in\mathcal M$, $AK\subseteq K$ and $A^*K\subseteq K$ (as $A^*\in\mathcal M$).
The second inclusion gives $AK^\perp\subseteq K^\perp$: for $y\in K^\perp$, $x\in K$,
$\langle x,Ay\rangle=\langle A^*x,y\rangle=0$. Hence $A$ maps $K$ into $K$ and $K^\perp$
into $K^\perp$, which means $AP=PAP=PA$. $\square$

**Theorem 2.2 (double commutant theorem).** Every $*$-algebra $\mathcal M$ on
$\mathbb C^n$ satisfies $\mathcal M=\mathcal M''$.

*Proof.* We have $\mathcal M\subseteq\mathcal M''$; let $T\in\mathcal M''$.

*Step 1: one vector.* For $x\in\mathbb C^n$, the subspace $K=\mathcal Mx$ is invariant
under $\mathcal M$ and contains $x=1x$. By Lemma 2.1 its projection $P$ is in $\mathcal M'$,
so $TP=PT$ and $Tx=TPx=PTx\in K$: there is $A\in\mathcal M$ with $Tx=Ax$.

*Step 2: all vectors at once.* On $\mathbb C^n\otimes\mathbb C^n=(\mathbb C^n)^{\oplus n}$
consider the $*$-algebra $\tilde{\mathcal M}=\{A\oplus\dots\oplus A:A\in\mathcal M\}$. A block
matrix $[X_{ij}]_{i,j=1}^n$ commutes with every $A\oplus\dots\oplus A$ iff $X_{ij}A=AX_{ij}$
for all $i,j$, so $\tilde{\mathcal M}'=\{[X_{ij}]:X_{ij}\in\mathcal M'\}$. The operator
$\tilde T=T\oplus\dots\oplus T$ commutes with every such $[X_{ij}]$, because
$TX_{ij}=X_{ij}T$. Hence $\tilde T\in\tilde{\mathcal M}''$. Apply Step 1 to $\tilde{\mathcal M}$,
$\tilde T$ and the vector $\xi=(e_1,\dots,e_n)$: there is $A\in\mathcal M$ with
$(Te_1,\dots,Te_n)=(Ae_1,\dots,Ae_n)$, i.e. $T=A\in\mathcal M$. $\square$

So in finite dimensions "von Neumann algebra" and "unital $*$-algebra of matrices"
are the same thing. The double commutant theorem is nevertheless the right point of
view: physically, an algebra of observables is very often *given* as a commutant —
the operators that commute with a symmetry (Section 2.4).

## 2.2 The centre

**Lemma 2.3.** Let $\mathcal Z$ be a commutative $*$-algebra on $\mathbb C^n$. There are
projections $P_1,\dots,P_C\in\mathcal Z$, nonzero, pairwise orthogonal, with
$\sum_cP_c=1$, such that $\mathcal Z=\operatorname{span}\{P_1,\dots,P_C\}$. They are the
minimal nonzero projections of $\mathcal Z$ and are unique up to order.

*Proof.* Every $z\in\mathcal Z$ is normal ($z^*\in\mathcal Z$ commutes with $z$). Let
$z_1,\dots,z_r$ span $\mathcal Z$. We decompose $\mathbb C^n$ into joint eigenspaces. The
eigenspaces of $z_1$ are mutually orthogonal and span $\mathbb C^n$ (Theorem 1.3); they are
invariant under $z_2,\dots,z_r$ and their adjoints (which commute with $z_1$), so each
$z_j$ restricted to an eigenspace of $z_1$ is again normal, and we can decompose each
eigenspace of $z_1$ into eigenspaces of $z_2$, and so on. The result is an orthogonal
decomposition $\mathbb C^n=E_1\oplus\dots\oplus E_C$ into nonzero subspaces on which every
$z\in\mathcal Z$ acts as a scalar $\chi_c(z)$, where we have merged subspaces with the same
function $\chi_c:\mathcal Z\to\mathbb C$, so that $\chi_c\neq\chi_d$ for $c\neq d$. Let
$P_c$ be the projection onto $E_c$. Then $z=\sum_c\chi_c(z)P_c$ for every $z\in\mathcal Z$.

It remains to show $P_c\in\mathcal Z$. For $d\neq c$ choose $z_{cd}\in\mathcal Z$ with
$\chi_c(z_{cd})\neq\chi_d(z_{cd})$. The product
$\prod_{d\neq c}\big(z_{cd}-\chi_d(z_{cd})1\big)/\big(\chi_c(z_{cd})-\chi_d(z_{cd})\big)$
lies in $\mathcal Z$, acts as $1$ on $E_c$ and as $0$ on every $E_d$, $d\neq c$; it is
$P_c$. Minimality and uniqueness: a projection $\sum_c a_cP_c$ has $a_c\in\{0,1\}$, so
the nonzero projections of $\mathcal Z$ are the sums of nonempty sets of the $P_c$. $\square$

Applied to $\mathcal Z=Z(\mathcal M)$, the projections $P_c$ are the **minimal central
projections** of $\mathcal M$. Every $A\in\mathcal M$ decomposes as
$A=\sum_cAP_c$ with $AP_c=P_cAP_c$, so
$\mathcal M=\bigoplus_c\mathcal MP_c$ (an orthogonal direct sum of algebras acting on
$E_c=P_c\mathbb C^n$).

**Lemma 2.4.** Each $\mathcal M_c=\mathcal MP_c$, as an algebra on $E_c$ with unit $P_c$, is a
factor.

*Proof.* Let $y\in\mathcal M_c$ commute with $\mathcal M_c$. For $A\in\mathcal M$,
$Ay=AP_cy=(AP_c)y=y(AP_c)=yP_cA=yA$, using that $P_c$ is central and $y=yP_c=P_cy$.
So $y\in Z(\mathcal M)$, $y=\sum_d\chi_d(y)P_d$, and $y=yP_c=\chi_c(y)P_c$. $\square$

## 2.3 The structure theorem

**Lemma 2.5 (factors).** Let $\mathcal N$ be a factor on a space $E$. There are integers
$k,m\ge1$ and a unitary $V:\mathbb C^k\otimes\mathbb C^m\to E$ with
$V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$. Moreover $V^*\mathcal N'V=1_k\otimes M_m(\mathbb C)$.

*Proof.* *Step 1: a minimal projection.* $\mathcal N$ contains nonzero projections (e.g. $1$).
Let $e\in\mathcal N$ be a nonzero projection of minimal rank. We claim
$e\mathcal Ne=\mathbb Ce$. Let $x=x^*\in e\mathcal Ne$. Each spectral projection $q$ of $x$
for a nonzero eigenvalue is a polynomial $p(x)$ with $p(0)=0$ (Theorem 1.5), so
$q\in\mathcal N$ and $q=eqe$, i.e. $q\le e$. By minimality $q=e$ or $q=0$; hence $x$ has at
most one nonzero eigenvalue and its eigenvector space is the range of $e$, so
$x=\lambda e$. A general element of $e\mathcal Ne$ is $x_1+ix_2$ with
$x_{1,2}=x_{1,2}^*\in e\mathcal Ne$.

*Step 2: $e$ generates $E$.* Let $F=\operatorname{span}\mathcal Ne E$. It is invariant under
$\mathcal N$, so its projection $f$ lies in $\mathcal N'$ (Lemma 2.1). It is also invariant
under $\mathcal N'$: for $b\in\mathcal N'$, $bAe\xi=Aeb\xi\in F$. So $f\in\mathcal N''=\mathcal N$
(Theorem 2.2), $f\in Z(\mathcal N)=\mathbb C1$, and $f\neq0$, so $F=E$.

*Step 3: matrix units.* Call a $v\in\mathcal N$ with $v^*v=e$ an *$e$-isometry*; then $vv^*$
is a projection in $\mathcal N$ (as $(vv^*)^2=vev^*=vv^*$, since $ve=vv^*v=v$). Choose
$e$-isometries $v_1=e,v_2,\dots,v_k$ with mutually orthogonal ranges $e_i=v_iv_i^*$, $k$
maximal (the ranges are orthogonal subspaces of the same dimension $\operatorname{rank}e$,
so $k\le n$). We claim $\sum_ie_i=1$. Otherwise $g=1-\sum_ie_i\neq0$ is a projection in
$\mathcal N$. By Step 2, $gAe\neq0$ for some $A\in\mathcal N$ (if $gAe=0$ for all $A$, then
$g$ vanishes on $F=E$). Let $y=gAe$. Then $y^*y\in e\mathcal Ne$ is positive and nonzero,
so $y^*y=\lambda e$ with $\lambda>0$, and $v=y/\sqrt\lambda$ is an $e$-isometry with range
in $gE$, orthogonal to all $e_i$ — contradicting maximality.

For $i\neq j$, $v_i^*v_j=v_i^*e_ie_jv_j=0$. Set $E_{ij}=v_iv_j^*$. Then
$E_{ij}E_{kl}=v_iv_j^*v_kv_l^*=\delta_{jk}v_iev_l^*=\delta_{jk}E_{il}$,
$E_{ij}^*=E_{ji}$ and $\sum_iE_{ii}=1$. For $A\in\mathcal N$,
$A=\sum_{ij}e_iAe_j=\sum_{ij}v_i(v_i^*Av_j)v_j^*$ and $v_i^*Av_j\in e\mathcal Ne=\mathbb Ce$,
so $A=\sum_{ij}a_{ij}E_{ij}$ with $a_{ij}\in\mathbb C$. Thus $\mathcal N=\operatorname{span}\{E_{ij}\}$.

*Step 4: the unitary.* Let $m=\operatorname{rank}e$ and $f_1,\dots,f_m$ an orthonormal basis
of $eE$. The vectors $v_if_r$ are orthonormal:
$\langle v_if_r,v_jf_s\rangle=\langle f_r,v_i^*v_jf_s\rangle=\delta_{ij}\langle f_r,ef_s\rangle=\delta_{ij}\delta_{rs}$.
They span $E$, since $E=\bigoplus_ie_iE$ and $e_iE=v_ieE$. Define $V(e_i\otimes f_r)=v_if_r$.
Then $E_{ij}v_lf_r=\delta_{jl}v_if_r$, i.e. $V^*E_{ij}V=|i\rangle\langle j|\otimes1_m$, and
$V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$.

*Step 5: the commutant.* Write $X\in M_k\otimes M_m$ as a block matrix
$X=\sum_{ij}|i\rangle\langle j|\otimes X_{ij}$. For fixed $a,b$,
$(|a\rangle\langle b|\otimes1)X=\sum_j|a\rangle\langle j|\otimes X_{bj}$ and
$X(|a\rangle\langle b|\otimes1)=\sum_i|i\rangle\langle b|\otimes X_{ia}$. If they are equal,
comparing the blocks in position $(a,j)$ with $j\neq b$ gives $X_{bj}=0$, and in position
$(a,b)$ gives $X_{bb}=X_{aa}$. As $a,b$ are arbitrary, $X=1_k\otimes Y$ with $Y=X_{11}$;
conversely every $1_k\otimes Y$ commutes with $M_k\otimes1_m$. $\square$

**Theorem 2.6 (structure of finite-dimensional von Neumann algebras).** Let
$\mathcal M$ be a von Neumann algebra on $\mathbb C^n$. There are $C\ge1$, integers
$k_c,m_c\ge1$ with $\sum_ck_cm_c=n$, and a unitary
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$ such that
$$
V^*\mathcal MV=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},\qquad
V^*\mathcal M'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),
$$
and $Z(\mathcal M)=\operatorname{span}\{P_1,\dots,P_C\}$, $P_c$ the projection onto the $c$-th
summand. In particular $\mathcal M\cong\bigoplus_cM_{k_c}(\mathbb C)$ as $*$-algebras, and
$\dim\mathcal M=\sum_ck_c^2$, $\dim\mathcal M'=\sum_cm_c^2$.

*Proof.* By Lemma 2.3 and the remark after it, $\mathcal M=\bigoplus_c\mathcal M_c$ with
$\mathcal M_c$ acting on $E_c$, and each $\mathcal M_c$ is a factor (Lemma 2.4). Lemma 2.5
gives unitaries $V_c$ with $V_c^*\mathcal M_cV_c=M_{k_c}\otimes1_{m_c}$; put
$V=\bigoplus_cV_c$. For the commutant: $X\in\mathcal M'$ commutes with each $P_c$, so
$X=\sum_cP_cXP_c$, and $P_cXP_c$ commutes with $\mathcal M_c$; conversely any such sum
commutes with $\mathcal M$. Lemma 2.5 identifies $\mathcal M_c'$ on $E_c$. $\square$

The numbers $k_c$ describe the algebra itself. The **multiplicities** $m_c$ describe how
it sits in $M_n(\mathbb C)$: they are invisible to the algebraic structure but determine
the trace of $\mathbb C^n$ restricted to $\mathcal M$,
$$
\operatorname{Tr}_{\mathbb C^n}(A)=\sum_cm_c\operatorname{Tr}(A_c),\qquad A=\bigoplus_cA_c\otimes1_{m_c}.
$$
This is why the library works with a family of traces (Chapter 3).

> **Physical meaning.** The $P_c$ are **superselection sectors** or **symmetry
> sectors**: no observable in $\mathcal M$ connects different sectors, so a measurement
> of "which sector" is compatible with every other measurement. The centre is the
> *classical* part of the system. The multiplicity $m_c$ is a degeneracy: every
> eigenvalue of an observable in sector $c$ appears $m_c$ times in $\mathbb C^n$.

## 2.4 Example: a symmetry group and its commutant

Let a finite group $G$ act on $\mathbb C^n$ by a unitary representation $\pi$. The
operators that commute with the symmetry form $\pi(G)'$, a von Neumann algebra; its
commutant $\pi(G)''$ is the span of the $\pi(g)$ (Theorem 2.2 applied to that span,
which is a $*$-algebra). Decomposing $\pi$ into irreducible representations
$\pi\cong\bigoplus_\lambda\rho_\lambda\otimes1_{m_\lambda}$ shows (Schur's lemma, or
Theorem 2.6 applied to $\pi(G)''$) that
$$
\pi(G)''\cong\bigoplus_\lambda M_{\dim\rho_\lambda},\qquad
\pi(G)'\cong\bigoplus_\lambda M_{m_\lambda}.
$$
The roles of $k_c$ and $m_c$ are exchanged between the algebra and its commutant.

We take $G=S_3$ permuting the three factors of $(\mathbb C^2)^{\otimes3}$, i.e. three
qubits or three spins $\tfrac12$. Representation theory predicts: the symmetric subspace
(total spin $\tfrac32$, dimension $4$) carries the trivial representation with
multiplicity $4$; its complement (dimension $4$) carries the two-dimensional standard
representation with multiplicity $2$; the sign representation does not occur. So
$\pi(S_3)''\cong M_1\oplus M_2$ (dimension $5$) and $\pi(S_3)'\cong M_4\oplus M_2$
(dimension $20$). We verify this by linear algebra: the commutant of a set $S$ is the
null space of the linear map $X\mapsto(AX-XA)_{A\in S}$.

```python
import itertools
import torch

torch.manual_seed(0)
dt = torch.complex128

def permutation_operator(perm, n_qubits=3):
    """Unitary on (C^2)^{(x)3} that permutes the tensor factors."""
    dim = 2 ** n_qubits
    P = torch.zeros(dim, dim, dtype=dt)
    for bits in itertools.product([0, 1], repeat=n_qubits):
        src = int(''.join(map(str, bits)), 2)
        dst = int(''.join(str(bits[perm[i]]) for i in range(n_qubits)), 2)
        P[dst, src] = 1
    return P

group = [permutation_operator(p) for p in itertools.permutations(range(3))]

def commutant_basis(ops, tol=1e-10):
    """Orthonormal basis (Hilbert-Schmidt) of {X : X A = A X for all A in ops}."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    # vec(AX - XA) = (A (x) 1 - 1 (x) A^T) vec(X) for row-major vec
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    _, s, Vh = torch.linalg.svd(L)
    null = Vh[(s > tol).sum():].conj().resolve_conj()   # rows span the null space
    return [v.reshape(d, d) for v in null]

def span_dim(ops, tol=1e-10):
    return int((torch.linalg.svdvals(torch.stack([A.reshape(-1) for A in ops])) > tol).sum())

comm = commutant_basis(group)
print(len(comm), span_dim(group), len(commutant_basis(comm)))
assert len(comm) == 20                                   # pi(S3)'  = M_4 (+) M_2
assert span_dim(group) == 5                              # pi(S3)'' = M_1 (+) M_2
assert len(commutant_basis(comm)) == 5                   # double commutant theorem
```

The multiplicities are visible in the spectrum. Averaging over the group,
$X\mapsto\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$, maps any operator into the commutant
(Exercise 3). A random self-adjoint element of $\pi(S_3)'$ therefore has $4$ simple
eigenvalues on the symmetric subspace and $2$ doubly degenerate ones on its complement:
the factor $1_2$ in $M_2\otimes1_2$.

```python
def group_average(X):
    return sum(g @ X @ g.conj().T for g in group) / len(group)

H0 = torch.randn(8, 8, dtype=dt)
Hsym = group_average(H0 + H0.conj().T)                  # a random symmetric Hamiltonian
assert all(torch.allclose(g @ Hsym, Hsym @ g) for g in group)

P_sym = sum(group) / len(group)                         # projection onto the symmetric subspace
w, Q = torch.linalg.eigh(P_sym)
Q_rest = Q[:, w.real < 0.5]                             # orthonormal basis of the complement (dim 4)
e_rest = torch.linalg.eigvalsh(Q_rest.conj().T @ Hsym @ Q_rest)
print(e_rest)
assert Q_rest.shape[1] == 4
assert torch.allclose(e_rest[0], e_rest[1]) and torch.allclose(e_rest[2], e_rest[3])
```

## 2.5 The algebra in the library

The library stores the abstract algebra $\bigoplus_cM_{k_c}(\mathbb C)$ of Theorem 2.6 on
$\bigoplus_c\mathbb C^{k_c}$, without the multiplicities: an operator is the list of its
blocks $A_c$, padded to a common size. The commutant of the example is
`TypeIAlgebra([4, 2], [4, 2])`.

```python
from torch_vn_algebra import TypeIAlgebra, channels

alg = TypeIAlgebra([4, 2], [4, 2], precision='double')    # M_4 (+) M_2
print(alg.k_factors, alg.C)

A = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
print(A.matrix.shape)              # (1, 2, 4, 4): block 1 occupies the top-left 2 x 2 corner
assert A.matrix[0, 1, 2:, :].abs().max() == 0 and A.matrix[0, 1, :, 2:].abs().max() == 0
```

The centre consists of the elements $\sum_cz_cP_c$; `alg.central` builds them, and
they commute with everything:

```python
def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

P0, P1 = alg.central([1.0, 0.0]), alg.central([0.0, 1.0])   # minimal central projections
Z = alg.central([2.0, -1.0])
assert P0.is_projection and P1.is_projection
assert close(P0 + P1, alg.identity())
assert close(Z @ A, A @ Z)
assert close(A @ P0 + A @ P1, A)                         # A = sum_c A P_c
```

## 2.6 Conditional expectations

A **conditional expectation** from $\mathcal M$ onto a von Neumann subalgebra
$\mathcal N\subseteq\mathcal M$ is a linear map $E:\mathcal M\to\mathcal N$ that is positive
($A\ge0\Rightarrow E(A)\ge0$), satisfies $E(B)=B$ for $B\in\mathcal N$, and is an
$\mathcal N$-bimodule map: $E(B_1AB_2)=B_1E(A)B_2$ for $B_1,B_2\in\mathcal N$.

**Proposition 2.7.** For $\mathcal M=\bigoplus_cM_{k_c}$ and $\mathcal N=Z(\mathcal M)$,
$$
E(A)=\sum_c\frac{\operatorname{Tr}A_c}{k_c}\,P_c
$$
is a conditional expectation. It preserves every trace of the form
$\sum_cw_c\operatorname{Tr}A_c$.

*Proof.* Linearity is clear; $A\ge0$ implies $\operatorname{Tr}A_c\ge0$, so $E(A)\ge0$. For
$B=\sum_cb_cP_c$, $E(B)=\sum_cb_ck_cP_c/k_c=B$. For $B_i=\sum_cb^{(i)}_cP_c$,
$\operatorname{Tr}(B_1AB_2)_c=b^{(1)}_cb^{(2)}_c\operatorname{Tr}A_c$, giving the bimodule
property. Finally $\operatorname{Tr}(E(A)_c)=\operatorname{Tr}A_c$ for every $c$. $\square$

In the library this map is the channel `channels.center_expectation(alg)`. Applying it
to operators gives the classical (sector-averaged) part:

```python
E = channels.center_expectation(alg)
EA = E(A)
assert close(E(EA), EA)                                  # idempotent
assert close(E(Z @ A @ P1), Z @ EA @ P1)                 # bimodule property
assert close(EA @ A, A @ EA)                             # E(A) is central
for c, k in enumerate(alg.k_factors):                    # Tr E(A)_c = Tr A_c
    assert abs(EA.matrix[0, c].diagonal().sum() - A.matrix[0, c].diagonal().sum()) < 1e-10
```

The group average of Section 2.4 is the conditional expectation onto the commutant
$\pi(G)'$ inside $M_n(\mathbb C)$ (Exercise 3). Both are instances of the same idea:
projecting onto the operators that "do not see" some structure.

## 2.7 What to compare with

The examples of Section 2.4 can be checked with any linear-algebra package;
representation-theory software (GAP, SageMath) gives the decomposition
$\pi\cong\bigoplus\rho_\lambda\otimes1_{m_\lambda}$ symbolically. Exact-diagonalisation
codes in physics (QuSpin, the symmetry modules of QuTiP) use the same block structure to
diagonalise Hamiltonians sector by sector. The library takes the structure as the
starting point: once the blocks are known, all further computations — functions,
traces, states, channels, dynamics — stay inside the algebra.

## Exercises

1. Show that the centre of $M_k(\mathbb C)$ is $\mathbb C1$ (use the matrix units
   $|i\rangle\langle j|$), and that $Z(\bigoplus_cM_{k_c})$ has dimension $C$.
2. Let $S$ be the cyclic shift on $\mathbb C^3$, $Se_i=e_{i+1\bmod3}$. Determine
   $\{S\}'$ and $\{S\}''$ and check the dimensions numerically.
3. For a unitary representation $\pi$ of a finite group, show that
   $E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$ is a conditional expectation of
   $M_n(\mathbb C)$ onto $\pi(G)'$.

<details>
<summary>Solutions</summary>

1. If $X$ commutes with all $|i\rangle\langle j|$, then
   $X|i\rangle\langle j|=|i\rangle\langle j|X$; applying both sides to $e_j$ and comparing
   gives $Xe_i=X_{jj}e_i$ for all $i,j$, so $X=X_{11}1$. In a direct sum, an element
   $\bigoplus X_c$ is central iff every $X_c$ is, so the centre is
   $\{\bigoplus_cz_c1_{k_c}\}$, of dimension $C$.

2. $S$ is unitary with eigenvalues $1,\omega,\omega^2$, $\omega=e^{2\pi i/3}$, each simple.
   An operator commuting with $S$ preserves its eigenspaces, so in the eigenbasis it is
   diagonal: $\{S\}'\cong\mathbb C^3$, the circulant matrices. $\{S\}''$ is then also
   $\cong\mathbb C^3$ (the commutant of a maximal abelian algebra is itself), spanned by
   $1,S,S^2$.

```python
S = torch.roll(torch.eye(3, dtype=dt), 1, dims=0)
assert len(commutant_basis([S])) == 3
assert len(commutant_basis(commutant_basis([S]))) == 3
assert span_dim([torch.eye(3, dtype=dt), S, S @ S]) == 3
```

3. $E$ is linear and positive (a sum of maps $X\mapsto UXU^*$). For $h\in G$,
   $\pi(h)E(X)\pi(h)^*=\frac1{|G|}\sum_g\pi(hg)X\pi(hg)^*=E(X)$, so $E(X)\in\pi(G)'$. If
   $X\in\pi(G)'$, every term equals $X$, so $E(X)=X$. For $B_1,B_2\in\pi(G)'$,
   $\pi(g)B_1XB_2\pi(g)^*=B_1\pi(g)X\pi(g)^*B_2$, which gives the bimodule property.

</details>

**Next:** [Chapter 3](03_traces_and_states.md) classifies the traces on
$\bigoplus_cM_{k_c}$, introduces states and densities, and proves the entropy and
trace inequalities used later.
