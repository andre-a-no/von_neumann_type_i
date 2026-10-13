**English** · [Español](es/B_schur_weyl.md) · [Français](fr/B_schur_weyl.md) · [Deutsch](de/B_schur_weyl.md) · [中文](zh/B_schur_weyl.md) · [日本語](ja/B_schur_weyl.md) · [Русский](ru/B_schur_weyl.md)

# Appendix B. Schur–Weyl duality and the Weingarten calculus

The unitary group $U(d)$ acts on $(\mathbb C^d)^{\otimes n}$ by $V\mapsto V^{\otimes n}$, and the symmetric group
$S_n$ by permuting the tensor factors. Each of the two families of operators is the commutant of the other
(**Schur–Weyl duality**). This appendix proves the duality completely and derives from it the Weingarten formula
for the moments of a Haar unitary of every order. Along the way it proves the facts that Chapters 4, 5, 9 and 11
take from here: Theorem 4.14(b) for all $L$, the irreducibility of the $S_L$-representations on the multiplicity
spaces of the total spin, the Catalan numbers of Section 4.5, Theorem 5.12, the Wigner–Eckart theorem for $SU(2)$
used in Section 9.3, and the Haar measure of the orthogonal group used in Chapter 11.

The route is elementary: a polarisation identity identifies the commutant of the permutations with the span of all
$A^{\otimes n}$, $A\in M_d(\mathbb C)$; a Cayley-transform argument shows that the unitaries $V^{\otimes n}$ span the same
space; the double commutant theorem (Theorem 2.2) turns this into the duality. The Weingarten formula then follows
because the Haar average is the orthogonal projection onto the commutant.

Throughout, $d,n\ge1$, $\mathcal H_n=(\mathbb C^d)^{\otimes n}$ with the orthonormal basis
$e_{\mathbf i}=e_{i_1}\otimes\dots\otimes e_{i_n}$, where a multi-index $\mathbf i=(i_1,\dots,i_n)\in[d]^n$ is read as a map
$[n]\to[d]$, $m\mapsto i_m$ ($[d]=\{1,\dots,d\}$). $\operatorname{End}\mathcal H_n=M_{d^n}(\mathbb C)$ carries the
Hilbert–Schmidt inner product $\langle X,Y\rangle=\operatorname{Tr}(X^*Y)$. For $X\in M_d$ and $1\le m\le n$,
$X^{(m)}=1\otimes\dots\otimes X\otimes\dots\otimes1$ ($X$ in position $m$) and $d\Pi(X)=\sum_{m=1}^nX^{(m)}$. We write
$$
\mathcal T_n=\operatorname{span}\{A^{\otimes n}:A\in M_d(\mathbb C)\},\qquad
\mathcal U_n=\operatorname{span}\{V^{\otimes n}:V\in U(d)\}.
$$
Products of Kronecker products are computed factorwise, $(A_1\otimes\dots\otimes A_n)(B_1\otimes\dots\otimes B_n)=A_1B_1\otimes\dots\otimes A_nB_n$,
and $(A_1\otimes\dots\otimes A_n)^*=A_1^*\otimes\dots\otimes A_n^*$. Permutations are composed as maps,
$(\sigma\tau)(m)=\sigma(\tau(m))$, and $\#\text{cycles}(\sigma)$ counts the cycles of $\sigma$, fixed points included.

## B.1 Permutations of tensor factors

For $\sigma\in S_n$ let $\pi(\sigma)e_{\mathbf i}=e_{\mathbf i\circ\sigma^{-1}}$, i.e. the basis vector in slot $m$ moves to slot
$\sigma(m)$. This is the convention of `decompose.permutation_representation` and of Chapters 2 and 4.

**Lemma B.1.** (a) $\pi$ is a unitary representation of $S_n$, and
$\pi(\sigma)(x_1\otimes\dots\otimes x_n)=x_{\sigma^{-1}(1)}\otimes\dots\otimes x_{\sigma^{-1}(n)}$ for all $x_m\in\mathbb C^d$.
(b) $\pi(\sigma)(A_1\otimes\dots\otimes A_n)\pi(\sigma)^*=A_{\sigma^{-1}(1)}\otimes\dots\otimes A_{\sigma^{-1}(n)}$ for all
$A_m\in M_d$. In particular $\pi(\sigma)$ commutes with every $A^{\otimes n}$ and every $d\Pi(X)$, and
$\pi(\sigma)X^{(m)}\pi(\sigma)^*=X^{(\sigma(m))}$.
(c) $\operatorname{Tr}\pi(\sigma)=d^{\#\text{cycles}(\sigma)}$ and $\langle\pi(\sigma),\pi(\tau)\rangle=d^{\#\text{cycles}(\sigma^{-1}\tau)}$.
The number of cycles is invariant under conjugation and inversion, and two permutations are conjugate iff they
have the same cycle type.
(d) The $n!$ operators $\pi(\sigma)$ are linearly independent iff $n\le d$.

*Proof.* (a) $\mathbf i\mapsto\mathbf i\circ\sigma^{-1}$ is a bijection of $[d]^n$, so $\pi(\sigma)$ permutes an orthonormal basis
and is unitary. $\pi(\sigma)\pi(\tau)e_{\mathbf i}=\pi(\sigma)e_{\mathbf i\circ\tau^{-1}}=e_{\mathbf i\circ\tau^{-1}\circ\sigma^{-1}}=\pi(\sigma\tau)e_{\mathbf i}$.
Expanding $x_m=\sum_ax_m(a)e_a$, $x_1\otimes\dots\otimes x_n=\sum_{\mathbf i}\prod_mx_m(i_m)\,e_{\mathbf i}$; applying $\pi(\sigma)$ and
substituting $\mathbf j=\mathbf i\circ\sigma^{-1}$ gives $\sum_{\mathbf j}\prod_mx_m(j_{\sigma(m)})\,e_{\mathbf j}=\sum_{\mathbf j}\prod_{p}x_{\sigma^{-1}(p)}(j_p)\,e_{\mathbf j}$,
which is the claimed product vector.
(b) By (a), on a product vector: $\pi(\sigma)^*=\pi(\sigma^{-1})$ puts $x_{\sigma(m)}$ into slot $m$, the operator applies
$A_m$ there, and $\pi(\sigma)$ moves slot $m$ to slot $\sigma(m)$; so slot $p$ ends up holding $A_{\sigma^{-1}(p)}x_p$. Product
vectors span $\mathcal H_n$. For $A^{(m)}$ the factor $A$ moves from slot $m$ to slot $\sigma(m)$, so $d\Pi(X)$ is fixed.
(c) $\pi(\sigma)$ is a permutation matrix, and $e_{\mathbf i}$ is fixed iff $i_{\sigma^{-1}(p)}=i_p$ for all $p$, i.e. iff
$\mathbf i$ is constant on the cycles of $\sigma$: there are $d^{\#\text{cycles}(\sigma)}$ such $\mathbf i$. Then
$\langle\pi(\sigma),\pi(\tau)\rangle=\operatorname{Tr}(\pi(\sigma^{-1})\pi(\tau))=\operatorname{Tr}\pi(\sigma^{-1}\tau)$. The cycles of
$\kappa\rho\kappa^{-1}$ are the images under $\kappa$ of the cycles of $\rho$, and $\rho^{-1}$ runs through the same cycles
backwards. Conversely, if $\rho,\rho'$ have the same cycle type, write their cycles one below the other, ordered by
length, and let $\kappa$ map each entry of the first row to the entry below it; then $\kappa\rho\kappa^{-1}=\rho'$.
(d) If $n\le d$, let $\mathbf u=(1,2,\dots,n)$. The vectors $\pi(\sigma)e_{\mathbf u}=e_{\mathbf u\circ\sigma^{-1}}$ are distinct basis
vectors for distinct $\sigma$, so $\sum_\sigma c_\sigma\pi(\sigma)=0$ applied to $e_{\mathbf u}$ gives all $c_\sigma=0$. If $n>d$,
let $\mathcal A=\sum_\sigma\operatorname{sgn}(\sigma)\pi(\sigma)$. Every $\mathbf i\in[d]^n$ has $i_a=i_b$ for some $a\neq b$; for the
transposition $t=(a\,b)$, $\pi(t)e_{\mathbf i}=e_{\mathbf i}$ and $\mathcal A\pi(t)=\sum_\sigma\operatorname{sgn}(\sigma)\pi(\sigma t)=-\mathcal A$.
Hence $\mathcal Ae_{\mathbf i}=\mathcal A\pi(t)e_{\mathbf i}=-\mathcal Ae_{\mathbf i}=0$ for all $\mathbf i$: $\mathcal A=0$ is a nontrivial
linear relation. $\square$

## B.2 Polarisation: the commutant of the permutations

**Lemma B.2 (polarisation).** For $A_1,\dots,A_n\in M_d(\mathbb C)$, with $|\varepsilon|=\sum_m\varepsilon_m$,
$$
\sum_{\sigma\in S_n}A_{\sigma(1)}\otimes\dots\otimes A_{\sigma(n)}=
\sum_{\varepsilon\in\{0,1\}^n}(-1)^{n-|\varepsilon|}\Big(\sum_{m=1}^n\varepsilon_mA_m\Big)^{\otimes n}.
$$

*Proof.* By multilinearity, $\big(\sum_m\varepsilon_mA_m\big)^{\otimes n}=\sum_f\varepsilon_{f(1)}\cdots\varepsilon_{f(n)}\,A_{f(1)}\otimes\dots\otimes A_{f(n)}$,
the sum over all maps $f:[n]\to[n]$. Since $\varepsilon_m^r=\varepsilon_m$ for $r\ge1$, the coefficient of the term $f$ on the
right side is
$$
\sum_{\varepsilon}\prod_{m=1}^n(-1)^{1-\varepsilon_m}\prod_{m\in f([n])}\varepsilon_m
=\prod_{m\notin f([n])}\Big(\sum_{\varepsilon_m=0,1}(-1)^{1-\varepsilon_m}\Big)\prod_{m\in f([n])}\Big(\sum_{\varepsilon_m=0,1}(-1)^{1-\varepsilon_m}\varepsilon_m\Big).
$$
The first kind of factor is $-1+1=0$, the second is $1$. So the coefficient is $1$ if $f$ is onto, i.e. a
permutation, and $0$ otherwise. $\square$

The **symmetriser** is $\mathcal S(X)=\frac1{n!}\sum_\sigma\pi(\sigma)X\pi(\sigma)^*$ on $\operatorname{End}\mathcal H_n$.

**Proposition B.3.** $\mathcal S$ is the orthogonal projection of $\operatorname{End}\mathcal H_n$ onto $\pi(S_n)'$, and
$$
\pi(S_n)'=\mathcal T_n,\qquad\dim\pi(S_n)'=\binom{d^2+n-1}n .
$$

*Proof.* $\mathcal S$ is the twirl of the finite group $S_n$ in the representation $\pi$, hence the orthogonal
projection onto $\pi(S_n)'$ (Theorem 12.5 with $\pi_1=\pi_2=\pi$; also Section 4.2). By Lemma B.1(b),
$\mathcal T_n\subseteq\pi(S_n)'$. Conversely, $\operatorname{End}\mathcal H_n$ is spanned by the products
$A_1\otimes\dots\otimes A_n$ of matrix units, and by Lemma B.1(b) and the substitution $\sigma\to\sigma^{-1}$,
$\mathcal S(A_1\otimes\dots\otimes A_n)=\frac1{n!}\sum_\sigma A_{\sigma(1)}\otimes\dots\otimes A_{\sigma(n)}$, which lies in $\mathcal T_n$ by
Lemma B.2. So $\pi(S_n)'=\mathcal S(\operatorname{End}\mathcal H_n)\subseteq\mathcal T_n$.

*Dimension.* Let $P=[d]\times[d]$ and $E_\alpha=E_{\alpha_1}\otimes\dots\otimes E_{\alpha_n}$ for $\alpha\in P^n$, with the matrix units
$E_{(a,b)}=|e_a\rangle\langle e_b|$: a basis of $\operatorname{End}\mathcal H_n$. By Lemma B.1(b),
$\mathcal S(E_\alpha)=\frac1{n!}\sum_\sigma E_{\alpha\circ\sigma^{-1}}$, which is $|O_\alpha|^{-1}\sum_{\beta\in O_\alpha}E_\beta$ for the orbit
$O_\alpha=\{\alpha\circ\sigma\}$ (each element of the orbit occurs for the same number $n!/|O_\alpha|$ of $\sigma$). Sums over
different orbits involve disjoint sets of basis elements, so the range of $\mathcal S$ has a basis indexed by the
orbits. Two tuples lie in the same orbit iff every element of $P$ occurs equally often in both: the orbits are the
multisets of size $n$ from the $N=d^2$ elements of $P$. Ordering $P$ as $1,\dots,N$, such a multiset is a sequence
$1\le c_1\le\dots\le c_n\le N$, and $c_m\mapsto c_m+m-1$ is a bijection onto the strictly increasing sequences in
$[N+n-1]$, i.e. the $n$-element subsets; there are $\binom{N+n-1}n$ of them. $\square$

## B.3 Schur–Weyl duality

**Lemma B.4 (polynomials that vanish on the unitary group).** Let $p:M_d(\mathbb C)\to\mathbb C$ be a polynomial in the
$d^2$ entries, homogeneous of degree $n\ge1$. If $p(V)=0$ for all $V\in U(d)$, then $p=0$.

*Proof.* *Step 1 (Cayley transform).* For $A\in M_d$ let $q(A)=p\big((1+iA)\operatorname{adj}(1-iA)\big)$, where
$\operatorname{adj}$ is the adjugate matrix; $q$ is a polynomial in the entries of $A$. Let $H$ be Hermitian. Its
eigenvalues are real (Corollary 1.4), so $1-iH$ is invertible, $\operatorname{adj}(1-iH)=\det(1-iH)(1-iH)^{-1}$ (Cramer's
rule), and by homogeneity $q(H)=\det(1-iH)^n\,p(C)$ with $C=(1+iH)(1-iH)^{-1}$. The matrices $1\pm iH$ and their inverses
are functions of $H$ and commute (Theorem 1.5), and $(1+iH)^*=1-iH$, so
$C^*C=(1+iH)^{-1}(1-iH)(1+iH)(1-iH)^{-1}=1$. Hence $q(H)=0$.

*Step 2 (from Hermitian to all matrices).* Choose a basis $B_1,\dots,B_{d^2}$ of the real vector space of Hermitian
matrices. It is also a basis of $M_d(\mathbb C)$ over $\mathbb C$, since every $A$ is $\frac{A+A^*}2+i\frac{A-A^*}{2i}$ and the
dimensions agree. $r(t)=q(\sum_kt_kB_k)$ is a polynomial in $t\in\mathbb C^{d^2}$ that vanishes on $\mathbb R^{d^2}$.
Such a polynomial has all coefficients $0$: for one variable it has infinitely many roots; for $N$ variables
write $r=\sum_jr_j(t')t_N^j$ — for each real $t'$ the polynomial in $t_N$ vanishes on $\mathbb R$, so every
$r_j$ vanishes on $\mathbb R^{N-1}$, and induction applies. Hence $q(A)=r(t)=0$ for every $A=\sum_kt_kB_k\in M_d$.

*Step 3 (inverting the transform).* Let $B\in M_d$ with $1+B$ invertible, and $A=-i(1+B)^{-1}(B-1)$. Then
$1-iA=1-(1+B)^{-1}(B-1)=2(1+B)^{-1}$ and $1+iA=(1+B)^{-1}2B$; as $B$ commutes with $1+B$,
$(1+iA)(1-iA)^{-1}=(1+B)^{-1}B(1+B)=B$. So $0=q(A)=\det(1-iA)^np(B)$ with $\det(1-iA)\neq0$: $p(B)=0$.

*Step 4.* For any $B\neq0$ and $0<t<1/\|B\|$, $-1$ is not an eigenvalue of $tB$ (if $tBv=-v$ then
$\|v\|=t\|Bv\|\le t\|B\|\|v\|<\|v\|$ unless $v=0$), so $t^np(B)=p(tB)=0$. Finally $p(0)=0$ by homogeneity. $\square$

**Lemma B.5 (Newton's identities).** Let $x_1,\dots,x_n$ be pairwise commuting elements of an algebra, $x_S=\prod_{m\in S}x_m$,
$e_k=\sum_{|S|=k}x_S$ ($e_0=1$) and $p_j=\sum_mx_m^j$. Then $k\,e_k=\sum_{j=1}^k(-1)^{j-1}e_{k-j}\,p_j$ for $1\le k\le n$.

*Proof.* For $1\le j\le k$ let $a_j=\sum x_Sx_m^j$, summed over the pairs $(S,m)$ with $|S|=k-j$ and $m\notin S$, and
$a_{k+1}=0$. In $e_{k-j}p_j=\sum_{|S|=k-j}\sum_mx_Sx_m^j$, the terms with $m\notin S$ give $a_j$, and the terms with $m\in S$,
written with $S'=S\setminus\{m\}$, are $x_{S'}x_m^{j+1}$ and give $a_{j+1}$. So the alternating sum telescopes to $a_1$, and
$a_1=\sum_{|T|=k}|T|\,x_T=k\,e_k$ (each $k$-set $T$ arises from the $k$ choices of $m\in T$). $\square$

**Theorem B.6 (Schur–Weyl duality).** Let $d,n\ge1$.
(a) $\pi(S_n)'=\mathcal T_n=\mathcal U_n=\operatorname{span}\{V^{\otimes n}:V\in SU(d)\}$. This space is a $*$-algebra, so it is the
algebra generated by the $V^{\otimes n}$, $V\in U(d)$.
(b) For every set $\mathcal X\subseteq M_d$ whose span is $M_d$, the unital algebra generated by
$\{d\Pi(X):X\in\mathcal X\}$ is $\mathcal T_n$.
(c) $\{V^{\otimes n}:V\in U(d)\}'=\{V^{\otimes n}:V\in SU(d)\}'=\{d\Pi(X):X\in M_d\}'=\operatorname{span}\pi(S_n)$.

*Proof.* (a) $\pi(S_n)'=\mathcal T_n$ is Proposition B.3, and $\mathcal U_n\subseteq\mathcal T_n$. Suppose $\mathcal U_n\neq\mathcal T_n$. Then
some $A_0^{\otimes n}\notin\mathcal U_n$; let $Y\neq0$ be the component of $A_0^{\otimes n}$ orthogonal to $\mathcal U_n$. The function
$p(A)=\langle Y,A^{\otimes n}\rangle$ is a homogeneous polynomial of degree $n$ in the entries of $A$ (each entry of
$A^{\otimes n}$ is a product of $n$ entries of $A$), it vanishes on $U(d)$, and $p(A_0)=\|Y\|^2\neq0$, contradicting Lemma B.4.
For $V\in U(d)$ choose $\theta$ with $\det V=e^{i\theta}$; then $W=e^{-i\theta/d}V\in SU(d)$ and $V^{\otimes n}=e^{in\theta/d}W^{\otimes n}$,
so the two spans agree. $\mathcal U_n$ contains $1=1^{\otimes n}$, $V^{\otimes n}W^{\otimes n}=(VW)^{\otimes n}$ and
$(V^{\otimes n})^*=(V^*)^{\otimes n}$, so it is a $*$-algebra, and every algebra containing the $V^{\otimes n}$ contains their span.

(b) Let $\mathcal A$ be that algebra. $d\Pi$ is linear, so $\mathcal A$ contains $d\Pi(X)$ for all $X\in M_d$. By Lemma B.1(b)
every $d\Pi(X)$ lies in the unital algebra $\pi(S_n)'=\mathcal T_n$, so $\mathcal A\subseteq\mathcal T_n$. Conversely, fix $A\in M_d$ and
apply Lemma B.5 to the pairwise commuting $x_m=A^{(m)}$: then $e_n=A^{\otimes n}$ and $p_j=d\Pi(A^j)\in\mathcal A$, and by
induction on $k$ every $e_k=\frac1k\sum_j(-1)^{j-1}e_{k-j}p_j$ lies in $\mathcal A$. So $A^{\otimes n}\in\mathcal A$ and
$\mathcal T_n\subseteq\mathcal A$.

(c) The commutant of a set equals the commutant of the algebra it generates. By (a) and (b) the three sets
generate $\mathcal T_n=\pi(S_n)'$, whose commutant is $\pi(S_n)''$. The span of $\pi(S_n)$ is a $*$-algebra
($\pi(\sigma)\pi(\tau)=\pi(\sigma\tau)$, $\pi(\sigma)^*=\pi(\sigma^{-1})$, $\pi(e)=1$), so $\pi(S_n)''=\operatorname{span}\pi(S_n)$ by
Theorem 2.2. $\square$

**Corollary B.7 (Theorem 4.14(b)).** For the total spin $S_a=\sum_is_a^{(i)}$ of $L$ spins $\frac12$, the algebra generated
by $S_x,S_y,S_z$ is $\{S_a\}''=\pi(S_L)'$, and $\{S_a\}'=\operatorname{span}\pi(S_L)$.

*Proof.* $M_2=\operatorname{span}\{1,s_x,s_y,s_z\}$ (the Pauli matrices and $1$ are a basis), $d\Pi(1)=L\cdot1$ and $d\Pi(s_a)=S_a$.
By Theorem B.6(b) the unital algebra $\mathcal A$ generated by $S_x,S_y,S_z$ is $\mathcal T_L=\pi(S_L)'$. Since the $S_a$ are
Hermitian, $\mathcal A$ is a $*$-algebra (the adjoint of a product of generators is the reversed product), so
$\mathcal A=\mathcal A''\supseteq\{S_a\}''$ by Theorem 2.2, and $\{S_a\}''\supseteq\mathcal A$ because $\{S_a\}''$ is a unital algebra
containing the $S_a$. Finally $\{S_a\}'=\{S_a\}'''=\pi(S_L)''=\operatorname{span}\pi(S_L)$ (Section 2.1 and Theorem B.6(c)). $\square$

By Theorem B.6(a) for $d=2$, the same algebra is spanned by the rotations $V^{\otimes L}$, $V\in SU(2)$; this does not
need the exponential parametrisation of $SU(2)$ mentioned in Section 4.5.

The code compares, for small $d$ and $n$, four dimensions: of $\operatorname{span}\pi(S_n)$, of the commutant of two random
$V^{\otimes n}$ (which contains $\{V^{\otimes n}\}'=\operatorname{span}\pi(S_n)$, and coincides with it in all cases printed), and of
the spans of many random $A^{\otimes n}$ and $V^{\otimes n}$ against $\binom{d^2+n-1}n$. It also checks that the library's
`permutation_representation` follows the convention of Lemma B.1(a), Lemma B.2, and the Newton recursion of the proof of
Theorem B.6(b).

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

## B.4 The decomposition of $(\mathbb C^d)^{\otimes n}$

**Theorem B.8.** There are $C\ge1$, integers $k_c,m_c\ge1$ with $\sum_ck_cm_c=d^n$, and a unitary
$W:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathcal H_n$ such that
$$
W^*\pi(\sigma)W=\bigoplus_{c=1}^C\rho_c(\sigma)\otimes1_{m_c},\qquad W^*V^{\otimes n}W=\bigoplus_{c=1}^C1_{k_c}\otimes\tau_c(V)
\qquad(\sigma\in S_n,\ V\in U(d)),
$$
where the $\rho_c$ are pairwise inequivalent irreducible unitary representations of $S_n$ on $\mathbb C^{k_c}$, and the
$\tau_c$ pairwise inequivalent irreducible unitary representations of $U(d)$ on $\mathbb C^{m_c}$. Moreover
$W^*\operatorname{span}\pi(S_n)W=\bigoplus_cM_{k_c}\otimes1_{m_c}$, $W^*\mathcal U_nW=\bigoplus_c1_{k_c}\otimes M_{m_c}$,
$$
\sum_ck_c^2=\dim\operatorname{span}\pi(S_n),\qquad\sum_cm_c^2=\binom{d^2+n-1}n,
$$
and $\sum_ck_c^2=n!$ if $n\le d$.

*Proof.* Proposition 2.7 for the finite group $S_n$ and Theorem 2.6 for $\mathcal M=\operatorname{span}\pi(S_n)=\pi(S_n)''$ give
$W$, the $\rho_c$ with the stated properties, and $W^*\pi(S_n)'W=\bigoplus_c1_{k_c}\otimes M_{m_c}$. By Theorem B.6(a),
$\pi(S_n)'=\mathcal U_n$. So $W^*V^{\otimes n}W=\bigoplus_c1\otimes\tau_c(V)$ for some $\tau_c(V)\in M_{m_c}$, and $V\mapsto\tau_c(V)$ is
a unitary representation because $V\mapsto V^{\otimes n}$ is one. Since the $V^{\otimes n}$ span $\mathcal U_n$, the linear map
$\sum_la_lV_l^{\otimes n}\mapsto(\sum_la_l\tau_c(V_l))_c$ is onto $\bigoplus_cM_{m_c}$; irreducibility and inequivalence of the
$\tau_c$ follow word for word as in the proof of Proposition 2.7. The dimensions are Theorem 2.6, Proposition B.3
and Lemma B.1(d). $\square$

For $d=2$ the sectors are identified by the total spin:

**Corollary B.9 ($L$ spins $\frac12$).** Let $V$ be a unitary with $V^*S_aV=\bigoplus_SS^{(S)}_a\otimes1_{m_S}$ for the total
spin of $L$ spins $\frac12$ (Corollary 4.12, Proposition 4.13). Then
$V^*\pi(\sigma)V=\bigoplus_S1_{2S+1}\otimes\rho_S(\sigma)$, where the $\rho_S$ are pairwise inequivalent irreducible unitary
representations of $S_L$ of dimension
$$
m_S=\binom Lj-\binom L{j-1},\qquad j=\tfrac L2-S\in\{0,1,\dots,\lfloor L/2\rfloor\}.
$$
In particular $\dim\operatorname{span}\pi(S_L)=\sum_Sm_S^2$ and $\dim\pi(S_L)'=\sum_S(2S+1)^2=\binom{L+3}3$.

*Proof.* By Theorem 4.3, $V^*\{S_a\}'V=\bigoplus_S1_{2S+1}\otimes M_{m_S}$, and $\{S_a\}'=\operatorname{span}\pi(S_L)$ by Corollary B.7. So
$V^*\pi(\sigma)V=\bigoplus_S1\otimes\rho_S(\sigma)$ with a unitary representation $\rho_S$, and the map
$\sum_\sigma a_\sigma\pi(\sigma)\mapsto(\sum_\sigma a_\sigma\rho_S(\sigma))_S$ is onto $\bigoplus_SM_{m_S}$; irreducibility and inequivalence
follow as in the proof of Proposition 2.7. The values of $m_S$ are Proposition 4.13. By Theorem 4.3,
$\dim\{S_a\}'=\sum_Sm_S^2$ and $\dim\{S_a\}''=\sum_S(2S+1)^2$; by Corollary B.7, $\{S_a\}'=\operatorname{span}\pi(S_L)$ and
$\{S_a\}''=\pi(S_L)'$, whose dimension is $\binom{4+L-1}L=\binom{L+3}3$ by Proposition B.3 with $d=2$. $\square$

**Corollary B.10 (Catalan numbers).** For $L\ge1$, $\displaystyle\sum_Sm_S^2=\binom{2L}L\frac1{L+1}$; this is
$\dim\{S_a\}'=\dim\operatorname{span}\pi(S_L)$ on $(\mathbb C^2)^{\otimes L}$.

*Proof.* Write $N=2^L$ and let $\operatorname{vec}:M_N\to\mathbb C^N\otimes\mathbb C^N$ be the linear bijection
$|x\rangle\langle y|\mapsto x\otimes\bar y$. Then $\operatorname{vec}(AX-XA)=(A\otimes1-1\otimes A^T)\operatorname{vec}X$, since
$A|x\rangle\langle y|B=|Ax\rangle\langle B^*y|\mapsto Ax\otimes B^T\bar y$. So $\{S_a\}'$ is mapped onto the common kernel of
$S_a\otimes1-1\otimes S_a^T$, $a=x,y,z$, where $S_a^T=\sum_i(s_a^T)^{(i)}$. From $s_x^T=s_x$, $s_y^T=-s_y$, $s_z^T=s_z$ and
$\sigma_ys_x\sigma_y=-s_x$, $\sigma_ys_y\sigma_y=s_y$, $\sigma_ys_z\sigma_y=-s_z$ we get $\sigma_y(-s_a^T)\sigma_y=s_a$, so the unitary
$Y=1\otimes\sigma_y^{\otimes L}$ gives $Y(S_a\otimes1-1\otimes S_a^T)Y^*=S_a\otimes1+1\otimes S_a$, the total spin of $2L$ spins
$\frac12$ on $(\mathbb C^2)^{\otimes2L}=\mathbb C^N\otimes\mathbb C^N$. In the decomposition of Corollary 4.12 for $2L$ spins, a vector
annihilated by all three components is annihilated by $\mathbf S^2$, which is $S(S+1)\ne0$ on the summands of spin
$S\ge\frac12$; conversely the components vanish on the spin-$0$ summand. So the common kernel is that summand, of
dimension $m_0(2L)=\binom{2L}L-\binom{2L}{L-1}=\binom{2L}L\big(1-\frac L{L+1}\big)$ (Proposition 4.13). The first
equality of dimensions is Theorem 4.3, the second Corollary B.7. $\square$

The code decomposes $\pi(S_3)$ on $(\mathbb C^3)^{\otimes3}$ and checks the dimensions of Theorem B.8, then decomposes
$\pi(S_L)$ on $(\mathbb C^2)^{\otimes L}$ for $L=5,6$ and compares the pairs $(k_c,m_c)$ with $(m_S,2S+1)$ of Corollary B.9.

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

## B.5 The Weingarten calculus

**Lemma B.11 (the Haar twirl is the orthogonal projection onto the commutant).** Let $U$ be Haar distributed on
$U(d)$ and $\Phi(X)=\mathbb E[U^{\otimes n}X(U^{\otimes n})^*]$ for $X\in\operatorname{End}\mathcal H_n$. Then $\Phi$ is the orthogonal
projection of $\operatorname{End}\mathcal H_n$ onto $\{V^{\otimes n}:V\in U(d)\}'=\operatorname{span}\pi(S_n)$.

*Proof.* The entries of $U^{\otimes n}X(U^{\otimes n})^*$ are polynomials in the entries of $U$ and $\bar U$, which are bounded
by $1$, so $\Phi$ is a well-defined linear map. For $V\in U(d)$, $VU$ has the law of $U$, so
$V^{\otimes n}\Phi(X)(V^{\otimes n})^*=\mathbb E[(VU)^{\otimes n}X((VU)^{\otimes n})^*]=\Phi(X)$: the range lies in the commutant. If $Y$ is in
the commutant, $U^{\otimes n}Y(U^{\otimes n})^*=Y$ and $\Phi(Y)=Y$. So $\Phi$ is idempotent with range the commutant. By
Lemma 5.1, $U^*=U^{-1}$ has the law of $U$, hence
$\langle Y,\Phi(X)\rangle=\mathbb E\operatorname{Tr}\big(((U^*)^{\otimes n}YU^{\otimes n})^*X\big)=\langle\Phi(Y),X\rangle$: $\Phi$ is self-adjoint. A
self-adjoint idempotent is the orthogonal projection onto its range
($\langle X-\Phi(X),\Phi(Z)\rangle=\langle\Phi(X-\Phi(X)),Z\rangle=0$). The commutant is $\operatorname{span}\pi(S_n)$ by Theorem B.6(c). $\square$

To write this projection in terms of the spanning family $\pi(\sigma)$, which is linearly dependent for $n>d$, we use the
**pseudo-inverse** of a Hermitian matrix $G$: $G^+=f(G)$ with $f(\lambda)=1/\lambda$ for $\lambda\neq0$ and $f(0)=0$
(functional calculus, Theorem 1.5). It is Hermitian, $G^+=G^{-1}$ if $G$ is invertible, $GG^+=\mathbf 1_{\neq0}(G)$ is
the projection onto $(\ker G)^\perp$, $G^+v=v/\lambda$ if $Gv=\lambda v$ with $\lambda\neq0$, and
$(R^*GR)^+=R^*G^+R$ for unitary $R$ (if $G=Q\Lambda Q^*$, then $R^*GR=(R^*Q)\Lambda(R^*Q)^*$).

**Lemma B.12.** Let $v_1,\dots,v_N$ be vectors in a finite-dimensional inner product space $\mathcal K$, $G_{ab}=\langle v_a,v_b\rangle$
their Gram matrix and $W=G^+$. Then the orthogonal projection onto $\operatorname{span}\{v_a\}$ is
$x\mapsto\sum_{a,b}v_a\,W_{ab}\,\langle v_b,x\rangle$.

*Proof.* Let $T:\mathbb C^N\to\mathcal K$, $Tc=\sum_ac_av_a$. Then $T^*x=(\langle v_b,x\rangle)_b$, $G=T^*T$, and the map in question is
$Q=TG^+T^*$. As $\langle c,Gc\rangle=\|Tc\|^2$, $\ker G=\ker T$, so $TG^+G=T(1-\Pi)=T$ for the projection $\Pi$ onto
$\ker T$. Hence $QT=T$: $Q$ is the identity on $\operatorname{ran}T$. $Q$ vanishes on $(\operatorname{ran}T)^\perp=\ker T^*$, and $Q^*=Q$. So
$Q$ is the orthogonal projection onto $\operatorname{ran}T$. $\square$

**Theorem B.13 (Weingarten formula of every order).** Let $n,d\ge1$, $U$ Haar distributed on $U(d)$, $G$ the
$n!\times n!$ matrix $G_{\sigma\tau}=d^{\#\text{cycles}(\sigma^{-1}\tau)}$ ($\sigma,\tau\in S_n$) and $W=G^+$.
(a) $W_{\sigma\tau}=\mathrm{Wg}(\sigma^{-1}\tau)$ for the function $\mathrm{Wg}(\rho)=\mathrm{Wg}(\rho,d)=W_{e\rho}$ on $S_n$, which depends only
on the cycle type of $\rho$; in particular $\mathrm{Wg}(\rho^{-1})=\mathrm{Wg}(\rho)$.
(b) $G$ is invertible iff $n\le d$, and then $W=G^{-1}$.
(c) For $X\in\operatorname{End}\mathcal H_n$,
$\mathbb E[U^{\otimes n}X(U^{\otimes n})^*]=\sum_{\sigma,\tau}\mathrm{Wg}(\sigma^{-1}\tau)\operatorname{Tr}(\pi(\tau)^*X)\,\pi(\sigma)$.
(d) For all $\mathbf i,\mathbf j,\mathbf i',\mathbf j'\in[d]^n$,
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\tau\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\tau(m)}}\;\mathrm{Wg}(\sigma\tau^{-1}),
$$
i.e. the sum runs over the pairs with $\mathbf i=\mathbf i'\circ\sigma$ and $\mathbf j=\mathbf j'\circ\tau$.
(e) Moments with different numbers $p\neq q$ of factors $U$ and $\bar U$ vanish.

This is Theorem 5.12 (with $k=d$), now proved.

*Proof.* (a) Write $g(\rho)=d^{\#\text{cycles}(\rho)}$, so $G_{\sigma\tau}=g(\sigma^{-1}\tau)$, and let $\Lambda_\kappa,R_\kappa$ be the
permutation matrices on $\mathbb C^{S_n}$ with $\Lambda_\kappa e_\sigma=e_{\kappa\sigma}$, $R_\kappa e_\sigma=e_{\sigma\kappa}$. Then
$(\Lambda_\kappa^*G\Lambda_\kappa)_{\sigma\tau}=g(\sigma^{-1}\kappa^{-1}\kappa\tau)=G_{\sigma\tau}$ and
$(R_\kappa^*GR_\kappa)_{\sigma\tau}=g(\kappa^{-1}\sigma^{-1}\tau\kappa)=G_{\sigma\tau}$ (Lemma B.1(c)). By the properties of the
pseudo-inverse, $W$ has the same two invariances: $W_{\kappa\sigma,\kappa\tau}=W_{\sigma\tau}$ and
$W_{\sigma\kappa,\tau\kappa}=W_{\sigma\tau}$. The first with $\kappa=\sigma^{-1}$ gives $W_{\sigma\tau}=\mathrm{Wg}(\sigma^{-1}\tau)$; the second
gives $\mathrm{Wg}(\kappa^{-1}\rho\kappa)=\mathrm{Wg}(\rho)$, so $\mathrm{Wg}$ is constant on conjugacy classes, which are the cycle types
(Lemma B.1(c)); $\rho^{-1}$ has the cycle type of $\rho$.
(b) $G$ is the Gram matrix of the $\pi(\sigma)$ (Lemma B.1(c)), and with $T$ as in the proof of Lemma B.12, $G=T^*T$ is
invertible iff $T$ is injective, i.e. iff the $\pi(\sigma)$ are linearly independent: iff $n\le d$ by Lemma B.1(d).
(c) Lemma B.11 and Lemma B.12 with $\mathcal K=\operatorname{End}\mathcal H_n$ and $v_\sigma=\pi(\sigma)$.
(d) The entries of a Kronecker product are $\langle e_{\mathbf i},U^{\otimes n}e_{\mathbf j}\rangle=\prod_mU_{i_mj_m}$, so the moment is
$\langle e_{\mathbf i},\Phi(X)e_{\mathbf i'}\rangle$ for $X=|e_{\mathbf j}\rangle\langle e_{\mathbf j'}|$. Here
$\operatorname{Tr}(\pi(\tau)^*X)=\langle\pi(\tau)e_{\mathbf j'},e_{\mathbf j}\rangle=\delta_{\mathbf j,\mathbf j'\circ\tau^{-1}}$ and
$\langle e_{\mathbf i},\pi(\sigma)e_{\mathbf i'}\rangle=\delta_{\mathbf i,\mathbf i'\circ\sigma^{-1}}$. By (c) the moment is
$\sum_{\sigma,\tau}\delta_{\mathbf i,\mathbf i'\circ\sigma^{-1}}\delta_{\mathbf j,\mathbf j'\circ\tau^{-1}}\mathrm{Wg}(\sigma^{-1}\tau)$; substituting
$\sigma\to\sigma^{-1}$, $\tau\to\tau^{-1}$ gives the formula.
(e) For real $\theta$, $e^{i\theta}1\in U(d)$, so $e^{i\theta}U$ has the law of $U$, and the moment equals $e^{i(p-q)\theta}$ times
itself. Take $\theta$ with $e^{i(p-q)\theta}\neq1$. $\square$

For $n=2\le d$, $G=\begin{pmatrix}d^2&d\\d&d^2\end{pmatrix}$ and $G^{-1}$ gives $\mathrm{Wg}(e)=\frac1{d^2-1}$,
$\mathrm{Wg}(\tau)=-\frac1{d(d^2-1)}$: Theorem 5.10. For $d=1$ every $G_{\sigma\tau}=1$, $W=\frac1{(n!)^2}G$, and the moment of
the phase $U=U_{11}$ is $\sum_{\sigma,\tau}\frac1{(n!)^2}=1=|U_{11}|^{2n}$, as it must be.

**Corollary B.14.** For all $n,d\ge1$, $\mathbb E|U_{11}|^{2n}=\dfrac{n!}{d(d+1)\cdots(d+n-1)}$.

*Proof.* With all indices $1$ in Theorem B.13(d), the moment is $\sum_{\sigma,\tau}\mathrm{Wg}(\sigma\tau^{-1})=\sum_{\sigma,\tau}W_{\sigma^{-1},\tau^{-1}}=\mathbf 1^TW\mathbf 1$ (part (a)). Every row
sum of $G$ is $c_n(d)=\sum_{\rho\in S_n}d^{\#\text{cycles}(\rho)}\ge d^n>0$, so $G\mathbf 1=c_n(d)\mathbf 1$, $W\mathbf 1=\mathbf 1/c_n(d)$
and the moment is $n!/c_n(d)$. Finally $c_n(x)=x(x+1)\cdots(x+n-1)$ (Chapter 5, Exercise 2, solution). $\square$

This extends Chapter 5, Exercise 2, to $n>d$, where $G$ is singular. The code computes $\mathrm{Wg}$ from the
pseudo-inverse, checks (a), (b) and Theorem 5.10, and compares the whole moment tensor of order $n=3$ for $d=2$ — a
case with singular $G$ — and the twirl of a random $X$ with Monte Carlo averages over Haar unitaries.

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

## B.6 Tensor operators and the Wigner–Eckart theorem for $SU(2)$

Section 9.3 builds spin correlations from the Wigner–Eckart theorem. We prove it for the spin representations of
Section 4.5. Spin $k$ is the representation $S^{(k)}_a$ of Theorem 4.11 on $\mathbb C^{2k+1}$ with basis $|k,q\rangle$,
$q=k,\dots,-k$; recall $S^{(k)}_\pm|k,q\rangle=\sqrt{k(k+1)-q(q\pm1)}\,|k,q\pm1\rangle$ (for $S_-=S_+^*$ this follows from the
formula for $S_+$).

**Definition.** Let $(S_a)$ be a spin representation on $\mathbb C^N$ and $k\in\{0,\frac12,1,\dots\}$. A **tensor operator of
rank $k$** is a family $(T_q)_{q=-k,\dots,k}$ of $N\times N$ matrices with
$$
[S_z,T_q]=q\,T_q,\qquad[S_\pm,T_q]=\sqrt{k(k+1)-q(q\pm1)}\;T_{q\pm1}
$$
(the coefficient vanishes for $q=\pm k$, where $T_{\pm(k+1)}$ is undefined). A **vector operator** is a triple
$(T_x,T_y,T_z)$ with $[S_a,T_b]=i\sum_c\varepsilon_{abc}T_c$, i.e. $[S_x,T_y]=iT_z$, $[S_x,T_z]=-iT_y$ and cyclically, and
$[S_a,T_a]=0$.

**Lemma B.15.** (a) If $(T_a)$ is a vector operator, its **spherical components** $T_{\pm1}=\mp\frac1{\sqrt2}(T_x\pm iT_y)$,
$T_0=T_z$ form a tensor operator of rank $1$, and for two vector operators
$\sum_aT_aT'_a=\sum_q(-1)^qT_qT'_{-q}$.
(b) For $L$ spins $\frac12$ and the total spin $S_a$, each site spin $(s^{(i)}_a)_a$ and $(S_a)_a$ itself are vector operators.
(c) $(T_q)$ is a tensor operator of rank $k$ iff the map $\Psi:\mathbb C^N\otimes\mathbb C^{2k+1}\to\mathbb C^N$,
$\Psi(v\otimes|k,q\rangle)=T_qv$, satisfies $\Psi\,(S_a\otimes1+1\otimes S^{(k)}_a)=S_a\Psi$ for $a=x,y,z$.

*Proof.* (a) From the relations, $[S_z,T_x\pm iT_y]=\pm(T_x\pm iT_y)$, $[S_\pm,T_z]=\mp(T_x\pm iT_y)$,
$[S_\pm,T_x\mp iT_y]=\pm2T_z$ and $[S_\pm,T_x\pm iT_y]=0$ (e.g. $[S_x+iS_y,T_x-iT_y]=-i[S_x,T_y]+i[S_y,T_x]=T_z+T_z$). In
spherical components: $[S_z,T_q]=qT_q$, $[S_\pm,T_0]=\sqrt2\,T_{\pm1}$, $[S_\pm,T_{\mp1}]=\sqrt2\,T_0$, $[S_\pm,T_{\pm1}]=0$,
which are the defining relations for $k=1$ ($\sqrt{2-q(q\pm1)}$ is $\sqrt2$ for $q=0$ and $q=\mp1$, and $0$ for $q=\pm1$).
For the product, $T_{+1}T'_{-1}+T_{-1}T'_{+1}=-\frac12\big((T_x+iT_y)(T'_x-iT'_y)+(T_x-iT_y)(T'_x+iT'_y)\big)=-(T_xT'_x+T_yT'_y)$.
(b) $[S_a,s^{(i)}_b]=[s^{(i)}_a,s^{(i)}_b]$ (other sites commute), and the spin-$\frac12$ matrices satisfy the spin relations;
sum over $i$ for $S_a$. (c) $\Psi(S_a\otimes1+1\otimes S^{(k)}_a)(v\otimes|k,q\rangle)=T_qS_av+\sum_{q'}(S^{(k)}_a)_{q'q}T_{q'}v$, so the
intertwining relation for $S_a$ says $[S_a,T_q]=\sum_{q'}(S^{(k)}_a)_{q'q}T_{q'}$ for all $q$. For $S_z,S_\pm$ these are the
defining relations, and the relations for $S_x,S_y,S_z$ and for $S_z,S_\pm$ are linear combinations of each other. $\square$

**Theorem B.16 (Wigner–Eckart).** Let $(S_a)$ be a spin representation on $\mathbb C^N$, $V$ a unitary with
$V^*S_aV=\bigoplus_SS^{(S)}_a\otimes1_{m_S}$ (Corollary 4.12), $|S,M;b\rangle=V(|S,M\rangle\otimes f_b)$ the corresponding orthonormal
basis ($f_b$ the standard basis of $\mathbb C^{m_S}$), and $(T_q)$ a tensor operator of rank $k$. Let $S,S'$ occur ($m_S,m_{S'}\ge1$).
(a) If $S'\notin\{|S-k|,|S-k|+1,\dots,S+k\}$, then $\langle S',M';a|T_q|S,M;b\rangle=0$ for all $M',M,q,a,b$.
(b) Otherwise the intertwiners from spin $S\otimes$ spin $k$ to spin $S'$ form a one-dimensional space; for every nonzero
such $C=C^{S'}_{S,k}:\mathbb C^{2S+1}\otimes\mathbb C^{2k+1}\to\mathbb C^{2S'+1}$ there is a unique $m_{S'}\times m_S$ matrix
$R=\langle S'\|T\|S\rangle$, the **reduced matrix element**, with
$$
\langle S',M';a|\,T_q\,|S,M;b\rangle=\langle S',M'|\,C\,(|S,M\rangle\otimes|k,q\rangle)\rangle\;R_{ab}\qquad\text{for all }M',M,q,a,b .
$$
(c) The element vanishes unless $M'=M+q$. In case (b) there is a $q$ with $\langle S',S'|C(|S,S'-q\rangle\otimes|k,q\rangle)\rangle\neq0$,
so $R$, and with it the whole block of $T$ between the multiplets $S$ and $S'$, is determined by the matrix elements
$\langle S',S';a|T_q|S,S'-q;b\rangle$ between the highest-weight vectors of $S'$ and one vector of each copy of $S$.

*Proof.* Fix $a\le m_{S'}$, $b\le m_S$, and let $J:\mathbb C^{2S+1}\to\mathbb C^N$, $Jx=V(x\otimes f_b)$ (in the summand $S$), and
$K:\mathbb C^N\to\mathbb C^{2S'+1}$ the adjoint of $x\mapsto V(x\otimes f_a)$ (in the summand $S'$). Write $\alpha\in\{x,y,z\}$ for the
component index. From the form of $V^*S_\alpha V$, $S_\alpha J=JS^{(S)}_\alpha$ and, taking adjoints of the same relation for $S'$
(the spin matrices are Hermitian), $KS_\alpha=S^{(S')}_\alpha K$. Let $t=K\Psi(J\otimes1)$ with $\Psi$ of Lemma B.15(c); its matrix
elements are $\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle=\langle S',M';a|T_q|S,M;b\rangle$, and
$$
t\,(S^{(S)}_\alpha\otimes1+1\otimes S^{(k)}_\alpha)=K\Psi(S_\alpha\otimes1+1\otimes S^{(k)}_\alpha)(J\otimes1)=KS_\alpha\Psi(J\otimes1)=S^{(S')}_\alpha\,t,
$$
so $t$ is an intertwiner from $\rho=$ spin $S\otimes$ spin $k$ (the spin representation
$S^{(S)}_\alpha\otimes1+1\otimes S^{(k)}_\alpha$, Section 4.5) to $\sigma=$ spin $S'$. Taking adjoints is a linear bijection between
$\operatorname{Hom}(\rho,\sigma)$ and $\operatorname{Hom}(\sigma,\rho)$ (the adjoint of $t\rho(\alpha)=\sigma(\alpha)t$ is
$\rho(\alpha)t^*=t^*\sigma(\alpha)$, as all matrices involved are Hermitian). By Corollary 4.12 and Chapter 4, Exercise 2 (proved
in its solution from Corollary 4.12 alone), $\rho$ is unitarily equivalent to $\bigoplus_j(\text{spin }j)$, one copy of each
$j\in\{|S-k|,|S-k|+1,\dots,S+k\}$. These spins are irreducible and pairwise inequivalent (Theorem 4.11), so this is a
decomposition as in Theorem 4.3, and its last statement gives $\dim\operatorname{Hom}(\sigma,\rho)=1$ if $S'$ is one of these
$j$, and $\operatorname{Hom}(\sigma,\rho)=0$ otherwise (spin $S'$ is then equivalent to none of them, by Theorem 4.11).
This gives (a), and in case (b) $t=R_{ab}C$ for a unique number $R_{ab}$.
(c) The intertwining relation for $S_z$ gives $M'\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle=(M+q)\langle S',M'|t(|S,M\rangle\otimes|k,q\rangle)\rangle$.
The range of $C\neq0$ is a nonzero invariant subspace of the irreducible spin $S'$, hence everything, so
$\langle S',S'|C(|S,M\rangle\otimes|k,q\rangle)\rangle\neq0$ for some basis vector, and then $M=S'-q$ by the same relation for
the intertwiner $C$. Divide by this number. $\square$

By Lemma B.15(a), the block between multiplets $S$ and $S'$ of a scalar product such as
$\mathbf s^{(i)}\cdot\mathbf s^{(j)}=\sum_q(-1)^q s^{(i)}_qs^{(j)}_{-q}$ is
$\sum_{S''}\sum_q(-1)^q(\text{block }S''\to S'\text{ of }s^{(i)}_q)(\text{block }S\to S''\text{ of }s^{(j)}_{-q})$, and each factor is given by
Theorem B.16 through a reduced matrix element and a Clebsch–Gordan map. This is the principle behind
`su2.spin_correlation`; the correctness of its implementation is what Chapter 9 checks numerically against the
`'direct'` method. The code checks Lemma B.15 for the site spins of six spins $\frac12$ and Theorem B.16 in a basis-free
form: in the decomposition computed by the library (whose sector bases need not be the $|S,M\rangle$), the block of
$(T_q)_q$ between two multiplets, as a matrix with rows $(M',M,q)$ and columns $(a,b)$, has rank at most one, and it is zero
outside the range of (a).

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

## B.7 Haar measures used in the course

The course uses the Haar measures of the following groups: $U(k)$ (Theorem 5.3), $U(1)$ (the case $k=1$; Chapter 5,
Exercise 3), $U(\mathcal M)\cong\prod_cU(k_c)$ (Corollary 5.4), $SU(k)$ (Section 5.1, after Corollary 5.4), the orthogonal
group $O(k)$ (Section 5.1, Chapter 11, Exercise 2, and real algebras), and the diagonal groups $T$ of Lemma C.2. Uniqueness
for all of them is Lemma 5.1, since each is a closed subgroup of some $U(n)$; existence is an explicit construction in
each case (for $T$ in Lemma C.2), and the construction for $O(k)$ is the following. No general existence theorem for the
Haar measure of a compact group is used anywhere in the course.

**Theorem B.17 (Haar measure of the orthogonal group).** Let $Z$ be a $k\times k$ real matrix with independent
$N(0,1)$ entries. Then $Z$ is invertible with probability one, $Z=OR$ with $O$ orthogonal and $R$ real upper triangular
with positive diagonal, and $O$ is Haar distributed on $O(k)$. If $Z=O'R'$ is any real QR decomposition ($O'$ orthogonal,
$R'$ real upper triangular), then $O'\operatorname{diag}(\operatorname{sgn}R'_{ii})=O$. For Haar distributed $O$ and
$X\in M_k(\mathbb C)$, $\mathbb E[O]=0$ and $\mathbb E[OXO^T]=\frac{\operatorname{Tr}X}k1$; in particular
$\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k$.

*Proof.* $\det Z$ is a polynomial in the $k^2$ real entries, not identically zero ($\det1=1$), and $Z$ has a density, so
$\det Z\neq0$ almost surely (Corollary A.7(b), as in the proof of Theorem 5.3). Gram–Schmidt applied to real columns produces real $q_j$ and real coefficients, i.e. a
real decomposition $Z=OR$ with $R_{jj}>0$; by the uniqueness in Lemma 5.2 it is the decomposition of Lemma 5.2, so $O$
depends continuously on $Z$. The density $(2\pi)^{-k^2/2}e^{-\|Z\|_2^2/2}$ of $Z$ is invariant under $Z\mapsto VZ$ for
$V\in O(k)$, a linear map of $\mathbb R^{k^2}$ that preserves $\|Z\|_2^2=\operatorname{Tr}(Z^TZ)$, hence orthogonal with Jacobian $1$;
so $VZ$ has the law of $Z$ (Proposition A.14(a)). As $VZ=(VO)R$, uniqueness gives $O(VZ)=VO(Z)$, so the law of $O$ is left-invariant on the
compact group $O(k)$, and Lemma 5.1 applies. The sign correction is the end of the proof of Theorem 5.3, with
$\operatorname{sgn}R'_{ii}=R'_{ii}/|R'_{ii}|$ ($R'_{ii}\neq0$ as $Z$ is invertible). Since $-1\in O(k)$, $-O$ has the law of $O$, so
$\mathbb E[O]=-\mathbb E[O]=0$ (the entries are bounded by $1$). Finally $T=\mathbb E[OXO^T]$ satisfies $VTV^T=T$ for $V\in O(k)$, since $VO$ has the law of
$O$. With $V$ the diagonal matrix with $-1$ in position $i$ and $1$ elsewhere, $T_{ij}=-T_{ij}$ for $j\neq i$; with
permutation matrices, all $T_{ii}$ are equal (as in Lemma 5.6). So $T=\lambda1$, and $k\lambda=\operatorname{Tr}T=\operatorname{Tr}X$. For
$X=|e_j\rangle\langle e_m|$, $(OXO^T)_{il}=O_{ij}O_{lm}$. $\square$

The code checks the second moments, $\mathbb E\operatorname{Tr}O=0$ and $\mathbb E(\operatorname{Tr}O)^2=\sum_{i,l}\mathbb E[O_{ii}O_{ll}]=1$ by
Monte Carlo. It also records a *numerical observation* about the QR routine of `torch.linalg.qr` (a statement about that
routine, not a theorem): without the sign correction its factor $Q'$ has $\mathbb E\operatorname{Tr}Q'$ many standard errors
away from $0$, so it is not Haar distributed.

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

**Back to:** [Chapter 4](04_symmetries.md) (Theorem 4.14), [Chapter 5](05_random_matrices.md) (Theorem 5.12),
[Chapter 9](09_spin_chains.md) (correlations).
