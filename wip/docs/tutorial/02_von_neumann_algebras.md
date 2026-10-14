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
$XA=AX$ for all $A\in S$ implies $X^*A^*=A^*X^*$ for all $A\in S$. Every element of $S$
commutes with every element of $S'$, so $S\subseteq S''$; applied to $S'$ this gives
$S'\subseteq S'''$, and taking commutants in $S\subseteq S''$ (which reverses inclusions)
gives $S'''\subseteq S'$. Hence $S'=S'''$.

When an algebra acts on a subspace $E$ of $\mathbb C^n$, we identify $E$ with
$\mathbb C^{\dim E}$ by an orthonormal basis, so that all results stated on $\mathbb C^n$
apply to it.

**Definition.** A **von Neumann algebra** on $\mathbb C^n$ is a $*$-algebra
$\mathcal M$ with $\mathcal M=\mathcal M''$. Its **centre** is
$Z(\mathcal M)=\mathcal M\cap\mathcal M'$. A **factor** is a von Neumann algebra with
$Z(\mathcal M)=\mathbb C1$.

The definition is made for operators on infinite-dimensional Hilbert spaces, which this course
does not treat. In finite dimensions the condition $\mathcal M=\mathcal M''$ is automatic:

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
$P_c$ (for $C=1$, $P_1=1$: the empty product is $1$). Minimality and uniqueness: a projection $\sum_c a_cP_c$ has $a_c\in\{0,1\}$, so
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
for a nonzero eigenvalue $\lambda$ is $q=p(x)$, where $p$ is the Lagrange polynomial with
$p(\lambda)=1$ and $p(\mu)=0$ for $\mu\in(\operatorname{spec}x\cup\{0\})\setminus\{\lambda\}$
(Theorem 1.5(a),(b): $p$ agrees on $\operatorname{spec}x$ with the indicator function of $\{\lambda\}$).
As $p(0)=0$, $p(t)=t\,r(t)$ for a polynomial $r$, so $q=x\,r(x)=r(x)\,x\in\mathcal N$, and $x=ex=xe$ gives
$q=eq=qe$; so the range of $q$ lies in the range of $e$. If $q\neq0$, minimality gives
$\operatorname{rank}q\ge\operatorname{rank}e$, so the two ranges coincide and $q=e$. Since the spectral
projections for distinct eigenvalues are orthogonal, $x$ has at most one nonzero eigenvalue
$\lambda$, with spectral projection $e$, so $x=\lambda e$ (and $x=0=0\cdot e$ if there is none). A general element of $e\mathcal Ne$ is $x_1+ix_2$ with
$x_{1,2}=x_{1,2}^*\in e\mathcal Ne$.

*Step 2: $e$ generates $E$.* Let $F=\operatorname{span}\mathcal Ne E$. It is invariant under
$\mathcal N$, so its projection $p_F$ lies in $\mathcal N'$ (Lemma 2.1). It is also invariant
under the $*$-algebra $\mathcal N'$: for $b\in\mathcal N'$, $bAe\xi=Aeb\xi\in F$; so by Lemma 2.1
applied to $\mathcal N'$, $p_F\in\mathcal N''=\mathcal N$ (Theorem 2.2). Hence
$p_F\in Z(\mathcal N)=\mathbb C1$, and $p_F\neq0$ (as $eE\subseteq F$), so $F=E$.

*Step 3: matrix units.* Call a $v\in\mathcal N$ with $v^*v=e$ an *$e$-isometry*; then $vv^*$
is a projection in $\mathcal N$. (First $v=ve$: for every $x$,
$\|v(1-e)x\|^2=\langle(1-e)x,v^*v(1-e)x\rangle=\langle(1-e)x,e(1-e)x\rangle=0$. Hence
$vv^*v=ve=v$ and $(vv^*)^2=vev^*=vv^*$. For the $v_i$ below this gives $e_iv_i=v_i$ and
$v_i^*e_i=v_i^*$.) Choose
$e$-isometries $v_1=e,v_2,\dots,v_k$ with mutually orthogonal ranges $e_i=v_iv_i^*$, $k$
maximal (the ranges are orthogonal subspaces of the same dimension $\operatorname{rank}e$,
so $k\le\dim E/\operatorname{rank}e$). We claim $\sum_ie_i=1$. Otherwise $g=1-\sum_ie_i\neq0$ is a projection in
$\mathcal N$. By Step 2, $gAe\neq0$ for some $A\in\mathcal N$ (if $gAe=0$ for all $A$, then
$gAe\xi=0$ for all $A,\xi$, so $g$ vanishes on $F=E$). Let $y=gAe$. Then $y^*y\in e\mathcal Ne$ is positive and nonzero,
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
They span $E$, since $E=\bigoplus_ie_iE$ and $e_iE=v_ieE$. Define $V(|i\rangle\otimes f_r)=v_if_r$, where $|1\rangle,\dots,|k\rangle$ is
the standard basis of $\mathbb C^k$ and we identify $\mathbb C^m$ with $eE$ via the $f_r$.
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
and $V^*Z(\mathcal M)V=\operatorname{span}\{V^*P_1V,\dots,V^*P_CV\}$, where $P_1,\dots,P_C$ are the
minimal central projections of $\mathcal M$ (Lemma 2.3) and $V^*P_cV$ is the projection onto the
$c$-th summand. The pairs $(k_c,m_c)$ are unique up to order: for *any* unitary $V$ with the
first property, $V(\cdot)V^*$ maps the projections onto the summands to the $P_c$ (see the proof),
so $k_c^2=\dim\mathcal MP_c$ and $m_c=\operatorname{rank}P_c/k_c$ are determined by $\mathcal M$. In particular $\mathcal M\cong\bigoplus_cM_{k_c}(\mathbb C)$ as $*$-algebras, and
$\dim\mathcal M=\sum_ck_c^2$, $\dim\mathcal M'=\sum_cm_c^2$.

*Proof.* By Lemma 2.3 and the remark after it, $\mathcal M=\bigoplus_c\mathcal M_c$ with
$\mathcal M_c$ acting on $E_c$, and each $\mathcal M_c$ is a factor (Lemma 2.4). Lemma 2.5
gives unitaries $V_c$ with $V_c^*\mathcal M_cV_c=M_{k_c}\otimes1_{m_c}$; put
$V=\bigoplus_cV_c$. For the commutant: $X\in\mathcal M'$ commutes with each $P_c$, so
$X=\sum_cP_cXP_c$, and $P_cXP_c$ commutes with $\mathcal M_c$; conversely any such sum
commutes with $\mathcal M$. Lemma 2.5 identifies $\mathcal M_c'$ on $E_c$, and conjugating by $V$ gives
the formula for $V^*\mathcal M'V$.

*Uniqueness.* Let $V$ be any unitary with $V^*\mathcal MV=\mathcal B:=\bigoplus_{c=1}^{C'}M_{k_c}\otimes1_{m_c}$ and
let $Q_c$ be the projection onto the $c$-th summand. An element $\bigoplus_cX_c$ of $\mathcal B$ ($X_c\in M_{k_c}\otimes1_{m_c}$)
is central iff each $X_c$ commutes with $M_{k_c}\otimes1_{m_c}$, i.e. (Lemma 2.5, Step 5) $X_c\in(M_{k_c}\otimes1_{m_c})\cap(1_{k_c}\otimes M_{m_c})=\mathbb C1$
(an element $a\otimes1=1\otimes b$ has all diagonal blocks equal to $b$ and all off-diagonal blocks $0$, so $a$ is a scalar).
Hence $Z(\mathcal B)=\operatorname{span}\{Q_1,\dots,Q_{C'}\}$, and $Z(\mathcal M)=VZ(\mathcal B)V^*$ is spanned by the nonzero,
pairwise orthogonal projections $VQ_cV^*$ with sum $1$. Each of them is a nonzero projection of $Z(\mathcal M)$, hence (Lemma 2.3)
a sum of a nonempty set $S_c$ of the $P_d$; the $S_c$ are disjoint (orthogonality) and cover all $d$ (sum $1$), and as
$C'=\dim Z(\mathcal M)=C$ each $S_c$ is a single element. So, after renumbering, $VQ_cV^*=P_c$, whence
$\dim\mathcal MP_c=\dim(M_{k_c}\otimes1_{m_c})=k_c^2$ and $\operatorname{rank}P_c=k_cm_c$. $\square$

The numbers $k_c$ describe the algebra itself. The **multiplicities** $m_c$ describe how
it sits in $M_n(\mathbb C)$: they are invisible to the algebraic structure but determine
the trace of $\mathbb C^n$ restricted to $\mathcal M$,
$$
\operatorname{Tr}_{\mathbb C^n}(A)=\sum_cm_c\operatorname{Tr}(A_c),\qquad A=V\Big(\bigoplus_cA_c\otimes1_{m_c}\Big)V^*.
$$
This is why the library works with a family of traces (Chapter 3). Likewise, for self-adjoint $A$,
the restriction of $A$ to $V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})$ has the eigenvalues of the block $A_c$, each with
$m_c$ times its multiplicity in $A_c$, since $A_c\otimes1_{m_c}$ is $m_c$ copies of $A_c$.

> **Physical meaning** (physical interpretation, not a theorem of this course). The $P_c$ are **superselection sectors** or **symmetry
> sectors**: no observable in $\mathcal M$ connects different sectors, so a measurement
> of "which sector" is compatible with every other measurement. The centre is the
> *classical* part of the system. The multiplicity $m_c$ is a degeneracy of every observable in $\mathcal M$.

## 2.4 Example: a symmetry group and its commutant

Let a finite group $G$ act on $\mathbb C^n$ by a unitary representation $\pi$. The
operators that commute with the symmetry form $\pi(G)'$, a von Neumann algebra; its
commutant $\pi(G)''$ is the span of the $\pi(g)$ (Theorem 2.2 applied to that span, which
is a $*$-algebra: it contains $\pi(e)=1$, products $\pi(g)\pi(h)=\pi(gh)$ and adjoints
$\pi(g)^*=\pi(g^{-1})$).

**Proposition 2.7.** Apply Theorem 2.6 to $\mathcal M=\pi(G)''$. Then
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$, where each $\rho_c$ is an irreducible
unitary representation of $G$ on $\mathbb C^{k_c}$, and $\rho_c,\rho_d$ are inequivalent for
$c\neq d$. Consequently
$$
\pi\cong\bigoplus_c\rho_c\otimes1_{m_c},\qquad
\pi(G)''\cong\bigoplus_cM_{\dim\rho_c},\qquad
\pi(G)'\cong\bigoplus_cM_{m_c}:
$$
the blocks of the symmetric operators $\pi(G)'$ have the sizes of the multiplicities of the
irreducible representations, and the roles of $k_c$ and $m_c$ are exchanged between an
algebra and its commutant.

*Proof.* $\pi(g)\in\mathcal M$, so $V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$ for some
$\rho_c(g)\in M_{k_c}$, and $g\mapsto\rho_c(g)$ is a unitary representation because $\pi$ is.
Since the $\pi(g)$ span $\mathcal M$, the $\rho_c(g)$ span $M_{k_c}$, and more strongly the map
$\sum_ga_g\pi(g)\mapsto\big(\sum_ga_g\rho_c(g)\big)_c$ is onto $\bigoplus_cM_{k_c}$.
*Irreducible:* a subspace of $\mathbb C^{k_c}$ invariant under all $\rho_c(g)$ is invariant
under their span $M_{k_c}$, hence $0$ or $\mathbb C^{k_c}$. *Inequivalent:* if
$\rho_d(g)=T\rho_c(g)T^{-1}$ for all $g$ and some invertible $T$, $c\neq d$, then every element
$(x_{c'})_{c'}$ of the image satisfies $x_d=Tx_cT^{-1}$, which excludes the element with
$x_c=1$, $x_d=0$ — a contradiction to surjectivity. The last two isomorphisms are
Theorem 2.6. $\square$

The decomposition of a representation into irreducibles is unique up to equivalence, so that
the $m_c$ are *the* multiplicities in any decomposition. *Proof.* Let $\rho,\rho'$ be irreducible
and $T$ an intertwiner, $T\rho(g)=\rho'(g)T$ for all $g$. Kernel and range of $T$ are invariant,
so $T=0$ or $T$ is invertible (Schur's lemma); if $\rho'=\rho$, the same applies to $T-\lambda1$
for an eigenvalue $\lambda$ of $T$, so $T=\lambda1$. Hence the intertwiners from $\rho$ to $\rho'$
form a space of dimension $1$ if $\rho\cong\rho'$ (for an invertible intertwiner $S$, $T\mapsto S^{-1}T$
identifies them with the intertwiners from $\rho$ to itself, which are $\mathbb C1$) and $0$ otherwise (a nonzero
one would be invertible, an equivalence). If $\pi\cong\bigoplus_j\sigma_j$
with irreducible $\sigma_j$ (repetitions allowed), an intertwiner from an irreducible $\rho$ to
$\bigoplus_j\sigma_j$ is a tuple $(T_j)_j$ of intertwiners from $\rho$ to $\sigma_j$, and composing with the equivalence
does not change dimensions; so the intertwiners from $\rho$ to
$\pi$ form a space of dimension $\#\{j:\sigma_j\cong\rho\}$, a number that depends only
on $\pi$ and $\rho$. For the decomposition of Proposition 2.7, whose $\rho_c$ are pairwise
inequivalent, it is $m_c$ if $\rho\cong\rho_c$ and $0$ if $\rho$ is equivalent to no $\rho_c$. $\square$

We take $G=S_3$ permuting the three factors of $(\mathbb C^2)^{\otimes3}$, i.e. three
qubits or three spins $\tfrac12$: $\pi(\sigma)|x\rangle=|\sigma\cdot x\rangle$ for $x\in\{0,1\}^3$, where
$(\sigma\cdot x)_{\sigma(i)}=x_i$. Let $F=\{v:\pi(g)v=v\text{ for all }g\}$ be the **symmetric subspace**. The
decomposition is: $F$ (dimension $4$) carries the trivial representation with multiplicity $4$; its
orthogonal complement (dimension $4$) carries the two-dimensional **standard representation** (the
permutation representation of $S_3$ on the vectors of $\mathbb C^3$ with coordinate sum $0$) with
multiplicity $2$; the sign representation does not occur. So
$\pi(S_3)''\cong M_1\oplus M_2$ (dimension $5$) and $\pi(S_3)'\cong M_4\oplus M_2$
(dimension $20$).

*Proof.* Write $V^*\pi(S_3)''V=\bigoplus_cM_{k_c}\otimes1_{m_c}$ (Theorem 2.6), so that
$\sum_ck_cm_c=8$ and $\sum_cm_c^2=\dim\pi(S_3)'$.

*Step 1: $\dim\pi(S_3)'=20$.* Since $\pi(g)^*|g\cdot x\rangle=|x\rangle$, the matrix
$\pi(g)X\pi(g)^*$ has the entry $X_{xy}$ at position $(g\cdot x,g\cdot y)$. Hence $X\in\pi(S_3)'$, i.e.
$\pi(g)X\pi(g)^*=X$ for all $g$, iff $(x,y)\mapsto X_{xy}$ is constant on the orbits of $S_3$ acting
on pairs $(x,y)$ by $g\cdot(x,y)=(g\cdot x,g\cdot y)$; the indicator matrices of the orbits form a basis of
$\pi(S_3)'$, whose dimension is therefore the number of orbits. A pair $(x,y)$ is a triple
$\big((x_1,y_1),(x_2,y_2),(x_3,y_3)\big)$ of elements of the four-element set $\{0,1\}^2$, and $S_3$
permutes the three entries, so an orbit is a multiset of size $3$ from four elements: $4$ of the form
$\{a,a,a\}$, $4\cdot3=12$ of the form $\{a,a,b\}$ with $b\neq a$, and $4$ of the form $\{a,b,c\}$ with
distinct elements, $20$ in all.

*Step 2: $k=(1,2)$, $m=(4,2)$.* As $k_c\ge1$, $\sum_cm_c\le\sum_ck_cm_c=8$. Let $M=\max_cm_c$.
Then $20=\sum_cm_c^2\le M\sum_cm_c\le8M$ and $M^2\le20$, so $M\in\{3,4\}$. If $M=4$, the other
$m_c$ have sum of squares $4$, so they are $(2)$ or $(1,1,1,1)$. If $M=3$, the other $m_c$ have sum
of squares $11$ and sum at most $5$; if all of them were $\le2$, their sum of squares would be at most
$2\cdot5=10$, so one of them is $3$ and the remaining ones have sum of squares $2$: they are $(1,1)$.
Hence $m$ is, up to order, $(4,2)$, $(4,1,1,1,1)$ or $(3,3,1,1)$. In the last two cases
$\sum_cm_c=8=\sum_ck_cm_c$ forces all $k_c=1$, and then $\pi(S_3)''\cong\bigoplus_cM_1$ would be
commutative, whereas the transpositions $(12)$ and $(23)$ act non-commutingly:
$\pi((12))\pi((23))|001\rangle=|100\rangle$ but $\pi((23))\pi((12))|001\rangle=|010\rangle$. So
$m=(4,2)$ (numbering the summands so that $m_1=4$), and $4k_1+2k_2=8$ with $k_1,k_2\ge1$ gives $k=(1,2)$. In particular
$\dim\pi(S_3)''=\sum_ck_c^2=5$.

*Step 3: the trivial and the sign representation.* $\rho_2$ fixes no nonzero vector $u$: the line
$\mathbb Cu$ would be an invariant subspace of $\mathbb C^2$, against irreducibility
(Proposition 2.7). Write $V^*v=v_1\oplus v_2$ with $v_1\in\mathbb C^1\otimes\mathbb C^4$ and
$v_2=\sum_{r=1}^2u_r\otimes f_r\in\mathbb C^2\otimes\mathbb C^2$ ($f_1,f_2$ the standard basis). Then
$v\in F$ iff $\rho_1(g)v_1=v_1$ and $\rho_2(g)u_r=u_r$ for all $g$ and $r$, so every $v\in F$ has
$v_2=0$. The nonzero vector $|000\rangle$ lies in $F$, so the scalar $\rho_1(g)$ fixes a nonzero vector:
$\rho_1$ is trivial, $F=V(\mathbb C^1\otimes\mathbb C^4)$ has dimension $4$, and its orthogonal
complement is $V(\mathbb C^2\otimes\mathbb C^2)$. The sign representation is not equivalent to
$\rho_1$ (as $\operatorname{sgn}(12)=-1$) nor to $\rho_2$ (the dimensions differ), so by the uniqueness of
the multiplicities shown above it occurs $0$ times.

*Step 4: $\rho_2$ is the standard representation.* The vectors $e_1=|100\rangle$, $e_2=|010\rangle$,
$e_3=|001\rangle$ are permuted by the $\pi(\sigma)$ as $\sigma$ permutes the positions, so
$W=\{\sum_ix_ie_i:\sum_ix_i=0\}$ is invariant and $\pi$ acts on it by the standard representation.
$W\perp F$: for $v\in F$ and $w=\sum_ix_ie_i\in W$, $\langle v,w\rangle=\langle\pi(g)v,w\rangle=\langle v,\pi(g)^{-1}w\rangle$
for every $g$, so $\langle v,w\rangle=\langle v,\frac16\sum_g\pi(g)w\rangle$, and
$\frac16\sum_g\pi(g)w=\frac13\big(\sum_ix_i\big)(e_1+e_2+e_3)=0$ (each $e_i$ is mapped to each $e_j$ by
exactly two of the six permutations). Hence $V^*W\subseteq\mathbb C^2\otimes\mathbb C^2$. The maps
$T_r:\mathbb C^2\otimes\mathbb C^2\to\mathbb C^2$, $u\otimes f\mapsto\langle f_r,f\rangle u$, satisfy
$T_r(\rho_2(g)\otimes1)=\rho_2(g)T_r$. Choose $r$ with $T_rV^*|_W\neq0$ (if both vanished, $V^*W$ would
be $0$). Then $T_rV^*|_W$ is a nonzero intertwiner from $\pi|_W$ to $\rho_2$; its range is a nonzero
invariant subspace of $\mathbb C^2$, hence $\mathbb C^2$, and as $\dim W=2$ it is bijective. So $\rho_2$ is
equivalent to the standard representation (which is, in passing, irreducible). $\square$

The code below confirms the two dimensions numerically: the commutant of a set $S$ is the null
space of the linear map $X\mapsto(AX-XA)_{A\in S}$, whose dimension is computed here in floating
point, by counting the singular values above a tolerance. This is a check of the values just proved,
not part of the proof.

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
assert len(comm) == 20                                   # Step 1: pi(S3)'  = M_4 (+) M_2
assert span_dim(group) == 5                              # Step 2: pi(S3)'' = M_1 (+) M_2
assert len(commutant_basis(comm)) == 5                   # double commutant theorem
```

The multiplicities are visible in the spectrum. Every self-adjoint $K\in\pi(S_3)'$ has the form
$V^*KV=(1_1\otimes\kappa_1)\oplus(1_2\otimes\kappa_2)$ with self-adjoint $\kappa_1\in M_4$, $\kappa_2\in M_2$
(Theorem 2.6). On $F$ its eigenvalues are those of $\kappa_1$; on the complement of $F$ they are those
of $\kappa_2$, each repeated twice, since $1_2\otimes\kappa_2$ is two copies of $\kappa_2$. So *every* eigenvalue of
$K$ on the complement of $F$ has even multiplicity. For a random $K$ more is true. Averaging over the
group, $E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$, maps $M_8(\mathbb C)$ onto the commutant and is the
identity there (Exercise 3), so for a Ginibre matrix $Z$ (independent standard complex Gaussian
entries) $K=E(Z+Z^*)$ has a density on the real vector space $\pi(S_3)'_h$ of self-adjoint elements
(Theorem A.15 and the remark after it). By Proposition 4.7(a) (a forward reference; its proof uses
Theorem 2.6 and Appendix A, not this paragraph), with probability one the six
eigenvalues of $\kappa_1$ and $\kappa_2$ are pairwise distinct: $K$ then has $4$ simple eigenvalues on $F$ and
$2$ doubly degenerate ones on its complement. The code checks the degenerate pairs (which hold for every $K$)
and, for the sampled $K$, that the two pairs and the four eigenvalues on $F$ are distinct, using that
$P_F=\frac1{|G|}\sum_g\pi(g)$ is the orthogonal projection onto $F$ ($\pi(h)P_F=P_F$ for all $h$, so its
range lies in $F$; $P_Fv=v$ for $v\in F$; and $P_F^*=P_F$, by the substitution $g\to g^{-1}$).

```python
def group_average(X):
    return sum(g @ X @ g.conj().T for g in group) / len(group)

H0 = torch.randn(8, 8, dtype=dt)                        # Ginibre: Re, Im independent N(0, 1/2)
Hsym = group_average(H0 + H0.conj().T)                  # K = E(Z + Z^*), a random symmetric Hamiltonian
assert all(torch.allclose(g @ Hsym, Hsym @ g) for g in group)

P_sym = sum(group) / len(group)                         # P_F, the projection onto the symmetric subspace
w, Q = torch.linalg.eigh(P_sym)
Q_rest = Q[:, w.real < 0.5]                             # orthonormal basis of the complement (dim 4)
e_rest = torch.linalg.eigvalsh(Q_rest.conj().T @ Hsym @ Q_rest)
print(e_rest)
assert Q_rest.shape[1] == 4
assert torch.allclose(e_rest[0], e_rest[1]) and torch.allclose(e_rest[2], e_rest[3])   # even multiplicities
e_sym = torch.linalg.eigvalsh(Q[:, w.real > 0.5].conj().T @ Hsym @ Q[:, w.real > 0.5])   # eigenvalues on F
assert e_sym.shape[0] == 4 and e_rest[2] - e_rest[1] > 1e-6                           # this sample: generic
assert (e_sym[1:] - e_sym[:-1]).min() > 1e-6                                           # 4 simple eigenvalues on F
assert (e_sym[:, None] - e_rest[None, :]).abs().min() > 1e-6                           # all six distinct
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
$\mathcal N\subseteq\mathcal M$ (a $*$-algebra contained in $\mathcal M$, with the same unit) is a linear map $E:\mathcal M\to\mathcal N$ that is positive
($A\ge0\Rightarrow E(A)\ge0$), satisfies $E(B)=B$ for $B\in\mathcal N$, and is an
$\mathcal N$-bimodule map: $E(B_1AB_2)=B_1E(A)B_2$ for $B_1,B_2\in\mathcal N$.

**Proposition 2.8.** For $\mathcal M=\bigoplus_cM_{k_c}$ on $\bigoplus_c\mathbb C^{k_c}$, with $P_c$ the projection onto
the $c$-th block, and $\mathcal N=Z(\mathcal M)=\operatorname{span}\{P_1,\dots,P_C\}$ (Exercise 1),
$$
E(A)=\sum_c\frac{\operatorname{Tr}A_c}{k_c}\,P_c
$$
is a conditional expectation. It preserves every trace of the form
$\sum_cw_c\operatorname{Tr}A_c$.

*Proof.* Linearity is clear; $A\ge0$ implies $\operatorname{Tr}A_c=\sum_i\langle e_i,Ae_i\rangle\ge0$ (sum over the
standard basis vectors $e_i$ of the $c$-th block), so $E(A)\ge0$. For
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
A2 = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
assert close(EA @ A2, A2 @ EA)                           # E(A) commutes with a fresh element
for c, k in enumerate(alg.k_factors):                    # E(A)_c = (Tr A_c / k_c) 1, so Tr E(A)_c = Tr A_c
    trA = A.matrix[0, c].diagonal().sum()
    assert close(EA.matrix[0, c, :k, :k], trA / k * torch.eye(k, dtype=dt))
    assert abs(EA.matrix[0, c].diagonal().sum() - trA) < 1e-10
```

The group average of Section 2.4 is a conditional expectation onto the commutant
$\pi(G)'$ inside $M_n(\mathbb C)$, and it preserves $\operatorname{Tr}$ (Exercise 3). Similarly, on the full
matrix algebra $M_D(\mathbb C)$, $D=\sum_ck_c$, which contains $\mathcal M=\bigoplus_cM_{k_c}$ as its
block-diagonal matrices, the map $E_{\mathcal M}(X)=\sum_cP_cXP_c$ is a conditional expectation onto
$\mathcal M$ that preserves $\operatorname{Tr}$. *Proof.* $P_cXP_c$ vanishes outside the $c$-th diagonal block,
so $E_{\mathcal M}(X)\in\mathcal M$. For $X\ge0$, $\langle v,P_cXP_cv\rangle=\langle P_cv,XP_cv\rangle\ge0$. An element
$B\in\mathcal M$ is block diagonal, so $BP_c=P_cB$ and $B=\sum_cP_cBP_c$; hence $E_{\mathcal M}(B)=B$ and
$E_{\mathcal M}(B_1XB_2)=\sum_cB_1P_cXP_cB_2=B_1E_{\mathcal M}(X)B_2$ for $B_1,B_2\in\mathcal M$. Finally
$\sum_c\operatorname{Tr}(P_cXP_c)=\sum_c\operatorname{Tr}(P_cX)=\operatorname{Tr}X$. $\square$

For a **state** $\omega$ on $\mathcal M$ (a linear functional with $\omega(A^*A)\ge0$ for all $A$ and $\omega(1)=1$;
states are studied in Chapter 3), Proposition 2.8 and linearity give
$\omega(E(A))=\sum_c\omega(P_c)\operatorname{Tr}(A_c)/k_c$, and $\omega\circ E=\omega$ on the centre, since $E(z)=z$
there. (Conditional expectations are not unique in general: onto $\mathbb C1$, every state
$\varphi$ gives one, $A\mapsto\varphi(A)1$; on $M_2(\mathbb C)$ the states $A\mapsto A_{11}$ and $A\mapsto A_{22}$
give two different ones.) All are instances of the same idea:
projecting onto the operators that "do not see" some structure.

> **Physical meaning** (physical interpretation, not a theorem of this course). The map $E_{\mathcal M}$ erases the coherences between superselection
> sectors, which is what happens when the charge (the sector) is measured and the result is not
> recorded. The conditional expectation $E$ onto the centre goes further: it keeps only the
> classical information, the weight $\omega(P_c)$ of each sector, and replaces the state inside each
> sector by the maximally mixed one (the normalised trace $\operatorname{Tr}(A_c)/k_c$).

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
   $M_n(\mathbb C)$ onto $\pi(G)'$ that preserves $\operatorname{Tr}$.

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
   $\cong\mathbb C^3$, because in the eigenbasis of $S$ a matrix that commutes with all
   diagonal matrices is diagonal (compare entries of $XD$ and $DX$ for $D$ with distinct
   diagonal entries); it is spanned by $1,S,S^2$ (Vandermonde: the eigenvalues are distinct).

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
   Finally $\operatorname{Tr}(\pi(g)X\pi(g)^*)=\operatorname{Tr}X$ for every $g$, so
   $\operatorname{Tr}E(X)=\operatorname{Tr}X$.

</details>

**Next:** [Chapter 3](03_traces_and_states.md) classifies the traces on
$\bigoplus_cM_{k_c}$, introduces states and densities, and proves the entropy and
trace inequalities used later.
