**English** · [Español](es/09_spin_chains.md) · [Français](fr/09_spin_chains.md) · [Deutsch](de/09_spin_chains.md) · [中文](zh/09_spin_chains.md) · [日本語](ja/09_spin_chains.md) · [Русский](ru/09_spin_chains.md)

# 9. Spin chains

A chain of $`L`$ spins $`\frac12`$ lives on $`(\mathbb C^2)^{\otimes L}`$, and dense linear algebra in dimension $`2^L`$
becomes expensive beyond $`L\approx14`$ (a dense matrix has $`4^L`$ entries). Symmetries push the limit: the magnetisation splits the space into $`L+1`$ sectors
(Chapter 4), translations split each sector again, and $`SU(2)`$ reduces the blocks once more. This chapter
builds these sectors for the XXZ chain, proves what they say about ground states, checks the library against
free fermions and the Bethe ansatz, and goes beyond dense matrices (Krylov methods) and clean chains (disorder).

**You will learn** why the algebra of a chain is $`\bigoplus_NM_{\binom LN}`$ and why the tower of these
algebras has Pascal's triangle as Bratteli diagram and index $`4\cos^2\frac\pi{2L+4}`$, how momentum sectors
are counted, the Marshall sign rule and its consequences for the momentum and spin of ground states
(Lieb–Mattis), the Jordan–Wigner solution of the XX chain, the two-magnon Bethe equations, why Lanczos
converges monotonically, and how to measure level statistics and free-fermion entanglement. The numerical
checks of this chapter are in [Tutorial 9](../tutorial/09_spin_chains.md).

Throughout, the sites are $`0,\dots,L-1`$, the basis of $`(\mathbb C^2)^{\otimes L}`$ consists of the
configurations $`|b\rangle=|b_0\dots b_{L-1}\rangle`$ with $`b_i=1`$ for spin up (the library's convention: bit
$`i`$ of a label is site $`i`$), $`s^a_i`$ ($`a=x,y,z`$) are the spin matrices $`\frac12\sigma_a`$ on site $`i`$,
$`s^\pm_i=s^x_i\pm is^y_i`$, $`n_i=s^z_i+\frac12`$, $`S^a=\sum_is^a_i`$, and $`N=S^z+\frac L2`$ is the number of
up spins. A chain is **open** (bonds $`(i,i+1)`$, $`i<L-1`$) or **periodic** (also the bond $`(L-1,0)`$; then
$`L\ge3`$, indices mod $`L`$). $`J>0`$ (antiferromagnetic), $`\hbar=1`$, $`\log`$ is natural.

<a id="s-9-1"></a>

## 9.1 The XXZ chain and its $`U(1)`$ sectors

The **XXZ chain** with anisotropy $`\Delta\in\mathbb R`$ and fields $`h_i\in\mathbb R`$ is
```math
H=J\sum_{(i,j)}\Big[\tfrac12\big(s^+_is^-_j+s^-_is^+_j\big)+\Delta\,s^z_is^z_j\Big]+\sum_ih_is^z_i ,
```
the sum over the bonds; $`\frac12(s^+_is^-_j+s^-_is^+_j)=s^x_is^x_j+s^y_is^y_j`$ moves an up spin along the bond.
For $`\Delta=1`$ it is the Heisenberg chain $`J\sum\mathbf s_i\cdot\mathbf s_j`$, for $`\Delta=0`$ the XX chain.

<a id="r-9-1"></a>

**Proposition 9.1.** $`[H,S^z]=0`$. Hence $`H\in U(1)'\cong\bigoplus_{N=0}^LM_{\binom LN}(\mathbb C)`$, the
block $`N`$ acting on the span of the configurations with $`N`$ up spins.

*Proof.* The diagonal terms commute with the diagonal $`S^z`$. For $`i\neq j`$, $`s^+_is^-_j`$ commutes with
$`s^z_l`$, $`l\ne i,j`$, and $`[s^z_i+s^z_j,s^+_is^-_j]=[s^z_i,s^+_i]s^-_j+s^+_i[s^z_j,s^-_j]=s^+_is^-_j-s^+_is^-_j=0`$ by
$`[s^z,s^\pm]=\pm s^\pm`$ (Section 4.5); take adjoints for $`s^-_is^+_j`$. Hence $`H`$ commutes with
$`N=S^z+\frac L2`$, which has the integer eigenvalues $`0,\dots,L`$ with eigenspaces spanned by the configurations
with $`N`$ up spins, of dimension $`\binom LN`$; the rest is Proposition 4.10(2) with $`Q=N`$. $`\square`$

*Computation:* Proposition 9.1, and the blocks $`H_N`$ that `SpinChain.xxz` builds from the configurations of each sector (`chain.labels[N]`), against the full $`2^L\times2^L`$ matrix — [tutorial, Section 9.1](../tutorial/09_spin_chains.md#s-9-1).

**The tower of chains.** Let $`\mathcal A_L=\bigoplus_NM_{\binom LN}`$ act on $`L`$ sites. Adding a site,
$`x\mapsto x\otimes1_2`$, maps $`\mathcal A_L`$ into $`\mathcal A_{L+1}`$, because $`x\otimes1`$ commutes with
$`S^z\otimes1+1\otimes s^z`$. This gives a tower $`\mathcal A_0\subset\mathcal A_1\subset\mathcal A_2\subset\cdots`$.

<a id="r-9-2"></a>

**Proposition 9.2.** The inclusion $`\mathcal A_L\subset\mathcal A_{L+1}`$ has the $`(L+2)\times(L+1)`$
inclusion matrix $`\Lambda_L`$ with $`\lambda_{N,N}=\lambda_{N,N-1}=1`$ and all other entries $`0`$; the block
sizes obey Pascal's rule $`\binom{L+1}N=\binom LN+\binom L{N-1}`$, the Bratteli diagram of the tower is
Pascal's triangle, and the bipartite graph of $`\Lambda_L`$ is a path with $`2L+3`$ vertices.

*Proof.* The configurations with $`N`$ up spins on $`L+1`$ sites are those with $`N`$ up spins on the first
$`L`$ sites and the new spin down, and those with $`N-1`$ and the new spin up. On these two subspaces
$`x\otimes1`$ acts as $`x_N`$ and $`x_{N-1}`$, so block $`N`$ of $`\mathcal A_{L+1}`$ contains blocks $`N`$ and $`N-1`$ of
$`\mathcal A_L`$ once each; the dimensions add up to Pascal's rule. In the bipartite graph the sector $`N`$ of
$`\mathcal A_{L+1}`$ is joined to $`N-1`$ and $`N`$ of $`\mathcal A_L`$, which gives the path
$`0'-0-1'-1-\cdots-L'-L-(L+1)'`$ (primes for $`\mathcal A_{L+1}`$). $`\square`$

A path from the top of the Bratteli diagram to sector $`N`$ at level $`L`$ is a configuration $`b`$; there are
$`\binom LN`$ of them, the size of the block. For a connected inclusion with inclusion matrix $`\Lambda`$ we
**define** the **index** as $`\|\Lambda\|^2`$ (operator norm); for $`M_k\otimes1\subset M_k\otimes M_m`$ the inclusion
matrix is the $`1\times1`$ matrix $`(m)`$, so the index is the $`m^2`$ of Section 7.2. `Inclusion.index` computes it.
(The name is borrowed from V. Jones' index theory, *Invent. Math.* 72 (1983) 1; in this course the index of an
inclusion of finite-dimensional algebras is this definition and nothing more.)

*In the library:* the index is `Inclusion.index` — [tutorial, Section 9.1](../tutorial/09_spin_chains.md#s-9-1).

<a id="r-9-3"></a>

**Lemma 9.3.** The adjacency matrix $`A_n`$ of the path with vertices $`1,\dots,n`$ has the eigenvalues
$`2\cos\frac{\pi k}{n+1}`$, $`k=1,\dots,n`$, with eigenvectors $`v_k(j)=\sin\frac{\pi jk}{n+1}`$; hence
$`\|A_n\|=2\cos\frac\pi{n+1}`$.

*Proof.* Put $`v_k(0)=v_k(n+1)=0`$. Then $`(A_nv_k)(j)=v_k(j-1)+v_k(j+1)=2\cos\frac{\pi k}{n+1}\,v_k(j)`$ by
$`\sin(a-b)+\sin(a+b)=2\sin a\cos b`$, and $`v_k\ne0`$ since $`v_k(1)=\sin\frac{\pi k}{n+1}\neq0`$. The $`n`$ eigenvalues are
distinct ($`\cos`$ is strictly decreasing on $`(0,\pi)`$), so they are all of them, and the norm of the
symmetric $`A_n`$ is the largest modulus, attained at $`k=1`$ and $`k=n`$. $`\square`$

<a id="r-9-4"></a>

**Corollary 9.4.** $`[\mathcal A_{L+1}:\mathcal A_L]=\|\Lambda_L\|^2=4\cos^2\frac\pi{2L+4}`$, increasing to $`4`$.

*Proof.* Ordering the vertices along the path of Proposition 9.2, the adjacency matrix
$`B=\begin{pmatrix}0&\Lambda_L\\\Lambda_L^T&0\end{pmatrix}`$ becomes $`A_{2L+3}`$ (a permutation similarity). By Lemma 1.1,
$`\|B\|^2=\|B^2\|`$ and $`\|\Lambda_L\Lambda_L^T\|=\|\Lambda_L^T\Lambda_L\|=\|\Lambda_L\|^2`$, and $`B^2=\Lambda_L\Lambda_L^T\oplus\Lambda_L^T\Lambda_L`$, so
$`\|B\|=\|\Lambda_L\|`$; Lemma 9.3 with $`n=2L+3`$ gives the value. $`\square`$

*Computation:* Proposition 9.2, Lemma 9.3 and Corollary 9.4 for $`L\le8`$, with `inclusions.chain_tower` and `inclusions.bratteli_diagram` — [tutorial, Section 9.1](../tutorial/09_spin_chains.md#s-9-1).

> **Physical meaning** (physical interpretation, not a theorem of this course). $`S^z`$ is conserved whenever the interaction is invariant under rotations about the
> $`z`$ axis; $`N`$ counts up spins, or particles after the Jordan–Wigner map of Section 9.4. In the tower each
> charge sector of the longer chain is fed by two sectors of the shorter one.

<a id="s-9-2"></a>

## 9.2 Translations, momentum and the Marshall sign rule

On a periodic chain let $`T`$ be the translation $`T|b_0b_1\dots b_{L-1}\rangle=|b_{L-1}b_0\dots b_{L-2}\rangle`$,
which moves site $`i`$ to site $`i+1`$. $`T`$ is a permutation matrix, so unitary, $`T^L=1`$ and
$`Ts^a_iT^*=s^a_{i+1}`$. Hence $`TS^zT^*=S^z`$ and $`THT^*=H`$ for a uniform field $`h_i=h`$: $`T`$ maps bonds to bonds.
Its eigenvalues are therefore of the form $`e^{ik}`$ with **momenta** $`k=2\pi m/L`$, $`m=0,\dots,L-1`$.

<a id="r-9-5"></a>

**Proposition 9.5 (momentum sectors).** Fix $`N`$. For a configuration $`a`$ with $`N`$ up spins let $`R_a`$ be the
size of its orbit $`\{T^ja\}`$. (a) $`R_a`$ divides $`L`$; for $`k=2\pi m/L`$ with $`mR_a\equiv0\pmod L`$ the vector
$`|a,k\rangle=\sum_{j=0}^{R_a-1}e^{-ikj}T^j|a\rangle`$ is nonzero and $`T|a,k\rangle=e^{ik}|a,k\rangle`$, and these
$`R_a`$ vectors form an orthogonal basis of the span of the orbit. (b) The eigenvalue $`e^{2\pi im/L}`$ of $`T`$ on
sector $`N`$ has multiplicity $`d_{N,m}=\#\{\text{orbits }:mR_a\equiv0\pmod L\}`$, and every operator commuting with
$`S^z`$ and $`T`$ — e.g. the translation-invariant XXZ chain — is block diagonal, in $`\bigoplus_N\bigoplus_mM_{d_{N,m}}`$.

*Proof.* (a) The stabiliser of $`a`$ in the cyclic group $`\{T^j\}\cong\mathbb Z_L`$ is a subgroup, and $`R_a`$ is its
index, a divisor of $`L`$; the vectors $`T^j|a\rangle`$, $`0\le j<R_a`$, are orthonormal and $`T^{R_a}|a\rangle=|a\rangle`$. If
$`e^{ikR_a}=1`$, then $`T|a,k\rangle=e^{ik}\sum_{j=1}^{R_a}e^{-ikj}T^j|a\rangle=e^{ik}|a,k\rangle`$, since the term $`j=R_a`$ equals
the term $`j=0`$; the coefficient of $`|a\rangle`$ is $`1`$. The admissible $`m\in\{0,\dots,L-1\}`$ are the multiples of
$`L/R_a`$, $`R_a`$ of them, with distinct eigenvalues of the unitary $`T`$, so the vectors are orthogonal and span the
$`R_a`$-dimensional span of the orbit. (b) Sector $`N`$ is the orthogonal sum of the spans of its orbits, each
invariant under $`T`$, and by (a) the eigenvalue $`e^{2\pi im/L}`$ occurs once in the span of the orbit of $`a`$ if
$`mR_a\equiv0\pmod L`$ and not at all otherwise. An $`X`$
commuting with $`S^z`$ and $`T`$ commutes with their spectral projections $`P_N`$, $`Q_m`$ (Theorem 1.5(d)), so
$`X=\sum_{N,m}P_NQ_mXP_NQ_m`$. $`\square`$

*Computation:* Proposition 9.5(b) for $`L=12`$, and the momentum blocks built by `SpinChain.xxz_momentum` — [tutorial, Section 9.2](../tutorial/09_spin_chains.md#s-9-2).

<a id="r-9-6"></a>

**Lemma 9.6 (Perron–Frobenius, symmetric case).** Let $`M`$ be a real symmetric matrix with $`M_{ij}\ge0`$ for
$`i\ne j`$ whose graph ($`i\sim j`$ iff $`M_{ij}>0`$) is connected. Then its largest eigenvalue $`\lambda`$ is simple,
with an eigenvector whose entries are all $`>0`$.

*Proof.* $`\lambda=\max_{\|x\|=1}\langle x,Mx\rangle`$, and the maximum is attained exactly on the unit vectors
of the $`\lambda`$-eigenspace (expand in an orthonormal eigenbasis). For such $`x`$, the vector $`|x|`$ of moduli
satisfies $`\langle|x|,M|x|\rangle\ge\mathrm{Re}\langle x,Mx\rangle=\lambda`$, because the diagonal terms agree and
$`M_{ij}|x_i||x_j|\ge\mathrm{Re}\,M_{ij}\bar x_ix_j`$ for $`i\ne j`$; so $`|x|`$ is an eigenvector. If $`|x|_i=0`$, then
$`0=(\lambda-M_{ii})|x|_i=\sum_{j\ne i}M_{ij}|x|_j`$ forces $`|x|_j=0`$ for all neighbours $`j`$ of $`i`$, and by
connectedness $`x=0`$. So no nonzero eigenvector has a zero entry; but an eigenspace of dimension $`\ge2`$
contains a nonzero vector with $`x_1=0`$. Hence $`\lambda`$ is simple, and $`|x|>0`$ is an eigenvector. $`\square`$

We prove the sign rule for any bipartite arrangement of couplings; the chain is a special case.

<a id="r-9-7"></a>

**Theorem 9.7 (Marshall sign rule).** Let the sites $`\{0,\dots,L-1\}`$ be split into two sets $`A\sqcup B`$, and let
```math
H=\sum_{i<j}\Big[\tfrac{J_{ij}}2\big(s^+_is^-_j+s^-_is^+_j\big)+K_{ij}\,s^z_is^z_j\Big]+\sum_ih_is^z_i
```
with real $`K_{ij},h_i`$ and real $`J_{ij}`$ such that $`J_{ij}\ge0`$ if $`i,j`$ lie in different sets, $`J_{ij}\le0`$ if they lie
in the same set, and the graph $`G`$ on the sites with an edge $`i\sim j`$ iff $`J_{ij}\ne0`$ is connected. Let $`N_B(b)`$
be the number of up spins of $`b`$ on $`B`$. In every sector $`N`$ the lowest eigenvalue of $`H_N`$ is simple, and its
eigenvector can be chosen with $`(-1)^{N_B(b)}\psi(b)>0`$ for every configuration $`b`$ with $`N`$ up spins.

In particular this holds for the XXZ chain with $`J>0`$ and arbitrary $`\Delta,h_i`$ if the chain is open, or
periodic with $`L`$ even, with $`A`$ the even and $`B`$ the odd sites: then $`J_{ij}=J`$, $`K_{ij}=J\Delta`$ on the bonds and
$`0`$ otherwise, every bond joins an even and an odd site (for the ring because $`L`$ is even), and $`G`$ contains the
path $`0-1-\cdots-(L-1)`$.

*Proof.* Let $`V|b\rangle=(-1)^{N_B(b)}|b\rangle`$; $`V`$ is diagonal, $`V=V^*=V^{-1}`$. The operator $`s^+_is^-_j`$ moves an up
spin from $`j`$ to $`i`$; it changes $`N_B`$ by $`\pm1`$ if $`i,j`$ lie in different sets and leaves it unchanged otherwise,
so $`Vs^+_is^-_jV=\mp s^+_is^-_j`$ accordingly, and the same for $`s^-_is^+_j`$. Two distinct configurations $`b,b'`$ are
joined by such a term iff $`b'`$ arises from $`b`$ by exchanging the opposite spins at two sites $`i,j`$, and then $`i,j`$
are determined by $`b,b'`$ (the sites where they differ). Hence $`M=-VH_NV`$ has the off-diagonal entries
$`J_{ij}/2`$ (different sets) or $`-J_{ij}/2`$ (same set), i.e. $`|J_{ij}|/2\ge0`$ in both cases, positive iff $`i\sim j`$ in
$`G`$, and $`0`$ between configurations that differ elsewhere. The graph of $`M`$ is connected. Indeed, write
$`(i\,j)`$ for the exchange of the spins at sites $`i,j`$; for distinct $`a,b,c`$, $`(a\,c)=(a\,b)(b\,c)(a\,b)`$, so by induction
along a path in $`G`$ from $`i`$ to $`j`$ every $`(i\,j)`$ is a composition of exchanges along edges of $`G`$, and each such
exchange either fixes a configuration (equal spins) or moves along an edge of the graph of $`M`$. If $`b\ne b'`$ have
$`N`$ up spins, there are sites $`i`$, up in $`b`$ and down in $`b'`$, and $`j`$, down in $`b`$ and up in $`b'`$; $`(i\,j)b`$ differs
from $`b'`$ at two sites fewer than $`b`$, so finitely many such steps lead from $`b`$ to $`b'`$. By Lemma 9.6 the
largest eigenvalue of $`M`$ — minus the lowest of $`H_N`$ — is simple with an eigenvector $`\varphi>0`$; take
$`\psi=V\varphi`$. $`\square`$

<a id="r-9-8"></a>

**Corollary 9.8 (momentum of the ground state).** On a periodic chain with $`L`$ even, $`J>0`$ and a uniform field,
the ground state of sector $`N`$ has momentum $`k=0`$ if $`N`$ is even and $`k=\pi`$ if $`N`$ is odd.

*Proof.* $`T`$ commutes with $`H_N`$ and the ground state $`\psi=V\varphi`$ is simple, so $`T\psi=e^{ik}\psi`$. Now
$`(T\psi)(b)=\psi(T^{-1}b)`$, and $`T^{-1}b`$ has on site $`i`$ the spin $`b_{i+1}`$; as $`L`$ is even, $`i+1`$ is even for odd
$`i`$, so $`N_B(T^{-1}b)=N-N_B(b)`$ and $`(VT\psi)(b)=(-1)^{N_B(b)}(-1)^{N-N_B(b)}\varphi(T^{-1}b)=(-1)^N\varphi(T^{-1}b)`$. Also
$`VT\psi=e^{ik}\varphi`$. Both $`\varphi`$ and $`\varphi\circ T^{-1}`$ are positive, so $`e^{ik}=(-1)^N`$. $`\square`$

*Computation:* Theorem 9.7 for open and periodic chains with random fields, and Corollary 9.8 for $`L=12`$ — [tutorial, Section 9.2](../tutorial/09_spin_chains.md#s-9-2).

> **Physical meaning** (physical interpretation, not a theorem of this course). In an antiferromagnet the
> lattice fixes the signs of the ground state: moving a spin across a bond changes the sign. Physicists use
> this to run quantum Monte Carlo without a sign problem on bipartite antiferromagnets; for frustrated chains,
> where an antiferromagnetic coupling joins two sites of the same set, the hypothesis of Theorem 9.7 fails.

<a id="s-9-3"></a>

## 9.3 $`SU(2)`$ at $`\Delta=1`$

For $`\Delta=1`$ and $`h=0`$, $`H=J\sum\mathbf s_i\cdot\mathbf s_j`$ commutes with $`S^x,S^y,S^z`$ (Proposition 4.15), so
$`H\cong\bigoplus_S1_{2S+1}\otimes H_S`$ with blocks $`H_S\in M_{m_S}`$, $`m_S=\binom L{L/2-S}-\binom L{L/2-S-1}`$
(Proposition 4.13), and the module `su2` builds $`H_S`$ on the highest-weight vectors (Proposition 4.15(b)).
The same holds for every $`H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j`$ with real $`J_{ij}`$: by Proposition 4.15(a) each
$`\mathbf s_i\cdot\mathbf s_j`$ lies in $`\{S_a\}'`$. Let $`E_0(S)`$ be the lowest eigenvalue of $`H_S`$, for the spins $`S`$ that occur.

*In the library:* the blocks $`H_S`$ are built by `SU2Chain` — [tutorial, Section 9.3](../tutorial/09_spin_chains.md#s-9-3).

<a id="r-9-9"></a>

**Theorem 9.9 (Lieb–Mattis).** Let $`H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j`$ with $`A,B,J_{ij}`$ as in Theorem 9.7
(so $`K_{ij}=J_{ij}`$, $`h_i=0`$), $`n_A=|A|\ge n_B=|B|`$, and $`S_*=\frac{n_A-n_B}2`$.
(a) For $`S_*\le S\le\frac L2`$ (with $`S-\frac L2`$ an integer) the ground state of the sector $`S^z=S`$ has total spin
$`S`$; for $`0\le M\le S_*`$ (with $`M-\frac L2`$ an integer) the ground state of the sector $`S^z=M`$ has total spin $`S_*`$.
(b) Hence $`E_0(S)`$ is a simple eigenvalue of $`H_S`$ for $`S\ge S_*`$, $`E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)`$, and
$`E_0(S)>E_0(S_*)`$ for $`S<S_*`$: the ground space of $`H`$ is a single multiplet, of spin $`S_*`$.
For the Heisenberg chain with $`J>0`$, open or periodic with $`L`$ even, $`A`$ and $`B`$ are the even and odd sites,
$`S_*=0`$ for even $`L`$ and $`S_*=\frac12`$ for odd $`L`$ (open chain), and $`E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)`$.

*Proof.* (a) Let $`M\ge0`$ with $`M-\frac L2\in\mathbb Z`$, $`s=\max(M,S_*)`$, and let $`\psi`$ be the ground state of the sector
$`S^z=M`$, simple by Theorem 9.7: $`\psi=V\varphi`$ with $`\varphi>0`$. $`\mathbf S^2`$ commutes with $`H`$ and with $`S^z`$, so it
preserves the one-dimensional ground space of the sector: $`\psi`$ has a total spin $`S'`$. We build a vector $`\chi`$ of
the sector with $`\mathbf S^2\chi=s(s+1)\chi`$ and $`\langle\psi,\chi\rangle\ne0`$; as eigenvectors of $`\mathbf S^2`$ for distinct
eigenvalues are orthogonal, this gives $`S'=s`$. Put $`p=\frac L2-s`$; then $`0\le p\le n_B`$, because $`s\ge S_*`$ means
$`p\le\frac{n_A+n_B}2-\frac{n_A-n_B}2`$. Choose distinct $`a_1,\dots,a_p\in A`$ and $`b_1,\dots,b_p\in B`$, let $`R`$ be the
remaining $`2s`$ sites, and let
```math
\chi_0=\bigotimes_{l=1}^p\tfrac1{\sqrt2}\big(|{\uparrow}_{a_l}{\downarrow}_{b_l}\rangle-|{\downarrow}_{a_l}{\uparrow}_{b_l}\rangle\big)\otimes|{\uparrow\cdots\uparrow}\rangle_R ,
\qquad \chi=(S^-)^{s-M}\chi_0 .
```
Each singlet is annihilated by $`s^\pm_{a_l}+s^\pm_{b_l}`$ and the up spins of $`R`$ by $`s^+_i`$, so $`S^+\chi_0=0`$; with
$`S^z\chi_0=s\chi_0`$ this gives $`\mathbf S^2\chi_0=(S^-S^++(S^z)^2+S^z)\chi_0=s(s+1)\chi_0`$ (Section 4.5), and
$`\mathbf S^2\chi=s(s+1)\chi`$ since $`\mathbf S^2`$ commutes with $`S^-`$. Writing $`S^-=\sum_l(s^-_{a_l}+s^-_{b_l})+S^-_R`$, the
summands commute and the first ones annihilate the singlets, so $`\chi`$ is the product of the singlets with
$`(S^-_R)^{s-M}|{\uparrow\cdots\uparrow}\rangle_R=(s-M)!\sum_D|D\rangle`$, the sum over the configurations of $`R`$ with
exactly $`s-M`$ spins down (they exist, as $`s-M\le2s=|R|`$); thus $`\chi\ne0`$, $`S^z\chi=M\chi`$, and the coefficients of $`\chi`$ are real. Its sign: a
singlet contributes $`+`$ when its up spin is on $`A`$ and $`-`$ when it is on $`B`$, and the $`R`$-part is positive. If
$`M=s`$, the sites of $`R`$ are all up; if $`M<s`$, then $`s=S_*`$, $`p=n_B`$ and $`R\subseteq A`$. In both cases $`N_B(b)`$ is the
number of singlets with the up spin on $`B`$ plus a constant $`|R\cap B|`$ (if $`M=s`$) or $`0`$, so $`(-1)^{N_B(b)}\chi(b)`$ has
the same sign for all $`b`$ with $`\chi(b)\ne0`$. Hence $`\langle\psi,\chi\rangle=\langle\varphi,V\chi\rangle=\pm\sum_b\varphi(b)|\chi(b)|\ne0`$.

(b) By Proposition 4.15(a), Corollary 4.12 and Theorem 4.11, a multiplet of spin $`S'`$ meets the sector $`S^z=M`$
in one vector if $`S'\ge|M|`$ and not at all otherwise, so the spectrum of $`H`$ on the sector $`S^z=M\ge0`$ is the
union of the spectra of the $`H_{S'}`$, $`S'\ge M`$, with multiplicities. By (a) its lowest eigenvalue is simple and
belongs to spin $`s=\max(M,S_*)`$: so $`E_0(s)`$ is simple in $`H_s`$ and $`E_0(s)<E_0(S')`$ for all $`S'\ge M`$, $`S'\ne s`$.
For $`M=S\ge S_*`$ this gives simplicity and $`E_0(S)<E_0(S')`$ for $`S'>S`$; for $`M=S<S_*`$ it gives $`E_0(S_*)<E_0(S)`$.
So $`E_0(S_*)`$ is the lowest eigenvalue of $`H`$, attained only in $`H_{S_*}`$, where it is simple. For the chain the
even sites outnumber the odd ones by $`0`$ ($`L`$ even) or $`1`$ ($`L`$ odd). $`\square`$

Theorem 9.7 goes back to W. Marshall (*Proc. Roy. Soc. A* 232 (1955) 48) and Theorem 9.9 to E. Lieb and
D. Mattis (*J. Math. Phys.* 3 (1962) 749). The **singlet–triplet gap** $`E_0(1)-E_0(0)`$ of the even Heisenberg
ring is positive by Theorem 9.9.

*Computation:* Theorem 9.9 against the $`U(1)`$ sectors for $`L=12`$, the singlet–triplet gap, and Theorems 9.7 and 9.9 on a bipartite graph that is not a chain — [tutorial, Section 9.3](../tutorial/09_spin_chains.md#s-9-3).

**Correlations.** On every even Heisenberg ring the ground-state correlations $`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$
alternate in sign; Theorems 9.7 and 9.9 force this. Let $`\psi`$ be the real normalised
ground state, a singlet (Theorem 9.9), hence the ground state of the sector $`S^z=0`$. From
$`\mathbf S^2\psi=0`$, $`\|S^\pm\psi\|^2=\langle\psi,S^\mp S^\pm\psi\rangle=\langle\psi,(\mathbf S^2-(S^z)^2\mp S^z)\psi\rangle=0`$, so $`S^\pm\psi=0`$
and $`\langle\psi,[S^+,[S^-,X]]\psi\rangle=0`$ for every $`X`$ (each of the four terms has $`S^\pm`$ next to $`\psi`$). For
$`X=s^z_0s^z_r`$, $`r\ne0`$, the relations $`[s^\pm,s^z]=\mp s^\pm`$, $`[s^+,s^-]=2s^z`$ give $`[S^-,X]=s^-_0s^z_r+s^z_0s^-_r`$ and
$`[S^+,[S^-,X]]=4s^z_0s^z_r-s^+_0s^-_r-s^-_0s^+_r`$. Hence $`\langle s^z_0s^z_r\rangle=\frac14\langle s^+_0s^-_r+s^-_0s^+_r\rangle`$ and
```math
\langle\mathbf s_0\cdot\mathbf s_r\rangle=\langle s^z_0s^z_r\rangle+\tfrac12\langle s^+_0s^-_r+s^-_0s^+_r\rangle
=\tfrac34\sum_{b_0\ne b_r}\psi(b)\psi(b') ,
```
with $`b'`$ obtained from $`b`$ by exchanging the spins of sites $`0`$ and $`r`$. By Theorem 9.7, $`\psi(b)\psi(b')`$ has the
sign of $`(-1)^{N_B(b)+N_B(b')}`$, and $`N_B(b')-N_B(b)`$ is $`\pm1`$ for odd $`r`$ and $`0`$ for even $`r`$; the sum is not
empty, so $`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$ is nonzero with the sign $`(-1)^r`$ for $`0<r<L`$.

*Computation:* the correlations $`\langle\mathbf s_0\cdot\mathbf s_r\rangle`$ for $`L=12`$, with `SU2Chain.spin_correlation`, which uses the Wigner–Eckart theorem (Theorem B.16) — [tutorial, Section 9.3](../tutorial/09_spin_chains.md#s-9-3).

<a id="s-9-4"></a>

## 9.4 Exact checks: free fermions and the Bethe ansatz

**Jordan–Wigner.** On an open chain put $`c_j=\big(\prod_{l<j}(-1)^{n_l}\big)s^-_j`$.

<a id="r-9-10"></a>

**Theorem 9.10 (the XX chain is free).** (a) $`\{c_i,c_j^*\}=\delta_{ij}`$, $`\{c_i,c_j\}=0`$ and $`c_j^*c_j=n_j`$.
(b) The open XX chain is $`H=\frac J2\sum_{j<L-1}(s^+_js^-_{j+1}+s^-_js^+_{j+1})=\sum_{i,j}h_{ij}c_i^*c_j`$ with
$`h=\frac J2A_L`$. (c) Let $`u_k`$ be orthonormal real eigenvectors of $`h`$ with eigenvalues
$`\varepsilon_k=J\cos\frac{\pi k}{L+1}`$ (Lemma 9.3) and $`b_k=\sum_ju_k(j)c_j`$. For $`K=\{k_1<\dots<k_N\}`$ the vectors
$`b_K^*|0\rangle=b_{k_1}^*\cdots b_{k_N}^*|0\rangle`$ ($`|0\rangle`$ all spins down) form an orthonormal eigenbasis of sector
$`N`$, with energies $`\sum_{k\in K}\varepsilon_k`$.

*Proof.* (a) Operators on different sites commute, and on one site $`(-1)^n`$ anticommutes with $`s^\pm`$. For
$`i<j`$, $`c_i`$ and $`c_j^*`$ (or $`c_j`$) meet on site $`i`$ as $`s^-_i`$ and $`(-1)^{n_i}`$ — one anticommutation — and
otherwise commute, so they anticommute. For $`i=j`$ the strings square to $`1`$ and
$`c_jc_j^*+c_j^*c_j=s^-_js^+_j+s^+_js^-_j=1`$, $`c_j^*c_j=s^+_js^-_j=n_j`$, $`c_j^2=(s^-_j)^2=0`$; the case $`i>j`$ follows
by exchanging the roles of $`i,j`$ (and taking adjoints). (b) $`c_j^*c_{j+1}=s^+_j(-1)^{n_j}s^-_{j+1}=s^+_js^-_{j+1}`$,
because $`s^+_j\neq0`$ only on $`n_j=0`$; take adjoints. (c) $`\{b_k,b_l^*\}=\sum_ju_k(j)u_l(j)=\delta_{kl}`$ and
$`\{b_k,b_l\}=0`$, and $`H=\sum_k\varepsilon_kb_k^*b_k`$ by the spectral decomposition of $`h`$. From the CAR,
$`[b_l^*b_l,b_k^*]=\delta_{lk}b_k^*`$, so $`[H,b_k^*]=\varepsilon_kb_k^*`$; with $`H|0\rangle=0`$ this gives
$`Hb_K^*|0\rangle=(\sum_K\varepsilon_k)b_K^*|0\rangle`$. Since $`b_k|0\rangle=0`$ (each $`c_j`$ lowers a spin), moving annihilators
to the right with the CAR gives $`\langle b_{K'}^*0,b_K^*0\rangle=\delta_{KK'}`$ by induction on $`N`$. Each $`c_j^*`$ raises
$`N`$ by one, so $`b_K^*|0\rangle`$ lies in sector $`|K|`$: $`\binom LN`$ orthonormal vectors in a space of that dimension. $`\square`$

*Computation:* all $`1024`$ levels of the open XX chain with $`L=10`$, against Theorem 9.10 — [tutorial, Section 9.4](../tutorial/09_spin_chains.md#s-9-4).

**Magnons.** On the ring with $`h=0`$ the fully polarised state $`|0\rangle`$ has energy $`E_F=J\Delta L/4`$. One up
spin at $`x`$ breaks two bonds: $`H|x\rangle=(E_F-J\Delta)|x\rangle+\frac J2(|x-1\rangle+|x+1\rangle)`$ ($`L\ge3`$).

<a id="r-9-11"></a>

**Proposition 9.11 (one magnon).** In sector $`N=1`$ of the ring, $`|k\rangle=\sum_xe^{-ikx}|x\rangle`$,
$`k=2\pi m/L`$, satisfies $`T|k\rangle=e^{ik}|k\rangle`$ and $`H|k\rangle=\big(E_F+J(\cos k-\Delta)\big)|k\rangle`$.

*Proof.* $`\sum_xe^{-ikx}(|x-1\rangle+|x+1\rangle)=\sum_y(e^{-ik(y+1)}+e^{-ik(y-1)})|y\rangle=2\cos k\,|k\rangle`$, using
$`e^{ikL}=1`$; the same reindexing gives $`T|k\rangle=e^{ik}|k\rangle`$. $`\square`$

Two magnons scatter. Write $`|x_1,x_2\rangle`$, $`x_1<x_2`$, for the up spins at $`x_1,x_2`$; a pair of neighbours
breaks two bonds, a separated pair four.

<a id="r-9-12"></a>

**Proposition 9.12 (two magnons; Bethe equations).** Let $`k_1,k_2\in\mathbb C`$, $`K=k_1+k_2`$, and
$`(A_{12},A_{21})\neq0`$ with (i) $`A_{12}(1+e^{iK}-2\Delta e^{ik_2})+A_{21}(1+e^{iK}-2\Delta e^{ik_1})=0`$ and
(ii) $`A_{21}=e^{ik_2L}A_{12}`$, $`A_{12}=e^{ik_1L}A_{21}`$. Let
$`\psi(x_1,x_2)=A_{12}e^{i(k_1x_1+k_2x_2)}+A_{21}e^{i(k_2x_1+k_1x_2)}`$ on $`\mathbb Z^2`$ and
$`\Psi=\sum_{0\le x_1<x_2<L}\psi(x_1,x_2)|x_1,x_2\rangle`$. If $`\Psi\neq0`$, then
$`H\Psi=\big(E_F+J(\cos k_1+\cos k_2-2\Delta)\big)\Psi`$. If the denominators do not vanish, (i) and (ii) say
```math
e^{ik_1L}=-\frac{1+e^{iK}-2\Delta e^{ik_1}}{1+e^{iK}-2\Delta e^{ik_2}},\qquad
e^{ik_2L}=-\frac{1+e^{iK}-2\Delta e^{ik_2}}{1+e^{iK}-2\Delta e^{ik_1}} .
```

*Proof.* Let $`E=E_F+J(\cos k_1+\cos k_2-2\Delta)`$. (1) For all $`(x_1,x_2)\in\mathbb Z^2`$, each plane wave gives
$`\psi(x_1\pm1,x_2)`$ and $`\psi(x_1,x_2\pm1)`$ summing to $`2(\cos k_1+\cos k_2)\psi(x_1,x_2)`$. (2) Substituting,
$`\psi(x,x)+\psi(x+1,x+1)-2\Delta\psi(x,x+1)`$ is $`e^{iKx}`$ times the left side of (i), so it vanishes for all $`x`$.
(3) By (ii), $`\psi(x_2,x_1+L)=A_{12}e^{ik_2L}e^{i(k_2x_1+k_1x_2)}+A_{21}e^{ik_1L}e^{i(k_1x_1+k_2x_2)}=\psi(x_1,x_2)`$, and
applying this twice, $`\psi(x_1+L,x_2+L)=\psi(x_1,x_2)`$. The representatives of a configuration $`\{y_1<y_2\}`$ in
$`D=\{x_1<x_2<x_1+L\}`$ are $`(y_1+nL,y_2+nL)`$ and $`(y_2+nL,y_1+(n+1)L)`$, so $`\psi`$ is a function of configurations.
Let $`(x_1,x_2)\in D`$. If $`2\le x_2-x_1\le L-2`$, all four hops stay in $`D`$ and the diagonal is $`E_F-2J\Delta`$, so
$`(H\Psi)(x_1,x_2)=(E_F-2J\Delta)\psi+\frac J2\cdot2(\cos k_1+\cos k_2)\psi=E\psi`$ by (1). If $`x_2=x_1+1`$, the hops
$`(x_1+1,x_2)`$, $`(x_1,x_2-1)`$ are blocked and the diagonal is $`E_F-J\Delta`$; by (1) and (2),
$`(H\Psi)(x_1,x_2)=(E_F-J\Delta)\psi+\frac J2\big(2(\cos k_1+\cos k_2)\psi-\psi(x_1,x_1)-\psi(x_1+1,x_1+1)\big)=E\psi`$.
If $`x_2=x_1+L-1`$, use the representative $`(x_2,x_1+L)`$, which is of the previous kind ($`L\ge3`$). $`\square`$

Proposition 9.12 only says that solutions give eigenvectors; it does not say that all eigenvectors arise
this way. Write $`k_{1,2}=K/2\pm q`$ and $`c=\cos\frac K2`$. Then $`1+e^{iK}-2\Delta e^{ik_{1,2}}=2e^{iK/2}(c-\Delta e^{\pm iq})`$,
so $`A_{12}=c-\Delta e^{iq}`$, $`A_{21}=-(c-\Delta e^{-iq})`$ solve (i). By (ii) and $`(A_{12},A_{21})\ne0`$ both
coefficients are nonzero, so multiplying the two equations (ii) gives
$`e^{iKL}=1`$, $`K=2\pi n/L`$, and then (ii) is equivalent to $`F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}=0`$,
one equation in one variable. We look for its roots
$`q\in(0,\pi)`$, where $`e^{-iqL/2}F(q)`$ is real or imaginary because $`c-\Delta e^{\pm iq}`$ are complex conjugate
(scattering states), and on the lines $`q=i\kappa`$, $`\pi+i\kappa`$, $`\kappa>0`$, where $`F`$ is real (bound states).

For even $`L`$ and $`\Delta\ne0`$ the Bethe vectors are not complete. If $`\Delta\ne0`$,
$`\Phi=\sum_x(-1)^x|x,x+1\rangle`$ satisfies $`T\Phi=-\Phi`$ and $`H\Phi=(E_F-J\Delta)\Phi`$ (the hops to $`|x-1,x+1\rangle`$ and
$`|x,x+2\rangle`$ cancel in pairs). A Bethe vector has $`T\Psi=e^{-iK}\Psi`$, as $`\psi(x_1-1,x_2-1)=e^{-iK}\psi(x_1,x_2)`$. If
$`e^{-iK}\ne-1`$, $`\Psi\perp\Phi`$ because they are eigenvectors of the unitary $`T`$ for distinct eigenvalues; if
$`e^{iK}=-1`$, then $`\cos k_1+\cos k_2=2\cos\frac K2\cos\frac{k_1-k_2}2=0`$, the energy of $`\Psi`$ is $`E_F-2J\Delta\ne E_F-J\Delta`$,
and $`\Psi\perp\Phi`$ as eigenvectors of $`H`$ for distinct eigenvalues. So $`\Phi`$ is orthogonal to all Bethe vectors of
Proposition 9.12, which therefore never span the sector $`N=2`$ of an even ring. The two-magnon computation is the case $`N=2`$
of the ansatz of H. Bethe (*Z. Phys.* 71 (1931) 205); we use no other case.

*Computation:* the roots of $`F`$ for $`L=9`$, $`\Delta=0.6`$, which give all $`36`$ states of the sector $`N=2`$, and the vector $`\Phi`$ for $`L=8`$ — [tutorial, Section 9.4](../tutorial/09_spin_chains.md#s-9-4).

> **Physical meaning** (physical interpretation, not a theorem of this course). A magnon is a spin wave: one
> flipped spin delocalised over the ring with momentum $`k`$. Two magnons scatter with a phase shift (real $`q`$)
> or bind (complex $`q`$).

<a id="s-9-5"></a>

## 9.5 Lanczos and Krylov methods

The sector $`N=L/2`$ has dimension $`12870`$ for $`L=16`$ and $`184756`$ for $`L=20`$, but $`H`$ is sparse ($`O(L)`$ nonzero
entries per row, `xxz_sparse`), and Krylov methods need only products $`Hv`$. Let $`H`$ be Hermitian of size $`n`$ and
$`v\ne0`$. The **Krylov space** is $`\mathcal K_m=\mathrm{span}\{v,Hv,\dots,H^{m-1}v\}`$.

<a id="r-9-13"></a>

**Proposition 9.13 (Lanczos).** Suppose $`\dim\mathcal K_m=m`$. Let $`q_1=v/\|v\|`$ and, for $`j<m`$, let $`q_{j+1}`$ be
the normalised component of $`Hq_j`$ orthogonal to $`\mathcal K_j`$. Then $`q_1,\dots,q_m`$ is an orthonormal basis of
$`\mathcal K_m`$; with $`Q_m=(q_1\cdots q_m)`$, $`T_m=Q_m^*HQ_m`$ is real symmetric tridiagonal with diagonal
$`\alpha_j=\langle q_j,Hq_j\rangle`$ and off-diagonal $`\beta_j=\langle q_{j+1},Hq_j\rangle>0`$; and $`HQ_m=Q_mT_m+r_me_m^T`$, where
$`r_m\perp\mathcal K_m`$ vanishes iff $`\mathcal K_{m+1}=\mathcal K_m`$ (otherwise $`r_m=\beta_mq_{m+1}`$, $`q_{m+1}=r_m/\|r_m\|`$).

*Proof.* Suppose $`\mathrm{span}\{q_1,\dots,q_j\}=\mathcal K_j`$ (true for $`j=1`$). Then $`H\mathcal K_j\subseteq\mathcal K_j+\mathbb CHq_j`$,
as $`Hq_i\in\mathcal K_{i+1}`$; so $`\mathcal K_{j+1}=\mathcal K_j+\mathbb CH^jv\subseteq\mathcal K_j+\mathbb CHq_j\subseteq\mathcal K_{j+1}`$, and
$`\mathcal K_{j+1}=\mathcal K_j+\mathbb CHq_j`$. For $`j<m`$ the dimension grows, so the component of $`Hq_j`$ orthogonal to
$`\mathcal K_j`$ is nonzero and the induction proceeds; for $`j=m`$ that component, $`r_m`$, vanishes iff
$`\mathcal K_{m+1}=\mathcal K_m`$. For $`i>j+1`$, $`q_i\perp\mathcal K_{j+1}\ni Hq_j`$, so $`(T_m)_{ij}=0`$, and $`(T_m)_{ji}=0`$ as $`T_m`$ is
Hermitian. $`\alpha_j`$ is real, and $`\beta_j=\langle q_{j+1},Hq_j\rangle`$ is the norm of the component, $`>0`$. So
$`Hq_j=\beta_{j-1}q_{j-1}+\alpha_jq_j+\beta_jq_{j+1}`$ for $`j<m`$ and $`Hq_m=\beta_{m-1}q_{m-1}+\alpha_mq_m+r_m`$. $`\square`$

The eigenvalues $`\theta^{(m)}_1\le\dots\le\theta^{(m)}_m`$ of $`T_m`$ are the **Ritz values**.

<a id="r-9-14"></a>

**Theorem 9.14 (monotone convergence).** Let $`\lambda_1\le\dots\le\lambda_n`$ be the eigenvalues of $`H`$ and
$`\dim\mathcal K_m=m`$. Then $`\theta^{(m)}_j\ge\lambda_j`$, and $`\theta^{(m)}_j\ge\theta^{(m+1)}_j`$ if $`\dim\mathcal K_{m+1}=m+1`$
($`j\le m`$). If $`v`$ has a nonzero component in the eigenspace of $`\lambda_1`$ and $`\mathcal K_{m+1}=\mathcal K_m`$, then
$`\theta^{(m)}_1=\lambda_1`$.

*Proof.* By Corollary A.9 (Courant–Fischer for compressions, Theorem A.8) applied to $`Q=Q_m`$,
$`\theta^{(m)}_j=\min_{U}\max_{x\in U,\|x\|=1}\langle x,Hx\rangle`$ over the $`j`$-dimensional $`U\subseteq\mathcal K_m=\mathrm{ran}\,Q_m`$.
The same formula over $`U\subseteq\mathcal K_{m+1}`$ gives $`\theta^{(m+1)}_j`$ and over $`U\subseteq\mathbb C^n`$ gives $`\lambda_j`$
(Theorem A.8); minimising over more subspaces gives a smaller value. If
$`\mathcal K_{m+1}=\mathcal K_m`$, then $`\mathcal K_m`$ is $`H`$-invariant and contains $`p(H)v`$ for every polynomial $`p`$, in
particular $`P_1v\ne0`$ ($`P_1`$ the spectral projection of $`\lambda_1`$, a polynomial in $`H`$ by Lagrange
interpolation, as in Theorem 1.5); so $`\theta^{(m)}_1\le\langle x,Hx\rangle=\lambda_1`$ at $`x=P_1v/\|P_1v\|\in\mathcal K_m`$. $`\square`$

A small residual certifies an eigenvalue nearby: for Hermitian $`H`$, a number $`E`$ and $`\|\psi\|=1`$,
$`\|H\psi-E\psi\|^2=\sum_i|\lambda_i-E|^2|\langle u_i,\psi\rangle|^2\ge\min_i|\lambda_i-E|^2`$ ($`u_i`$ an orthonormal eigenbasis). It
does not certify that this eigenvalue is the lowest one. For $`L=4`$ the ground-state energy of the Heisenberg ring is
exact: with $`\mathbf S_A=\mathbf s_0+\mathbf s_2`$, $`\mathbf S_B=\mathbf s_1+\mathbf s_3`$ the
ring is $`H=J\,\mathbf S_A\cdot\mathbf S_B=\frac J2(\mathbf S^2-\mathbf S_A^2-\mathbf S_B^2)\ge\frac J2(0-2-2)=-2J`$, since $`\mathbf S^2\ge0`$ and
$`\mathbf S_A^2,\mathbf S_B^2`$ have the eigenvalues $`0,2`$ (Proposition 4.13 with $`L=2`$); equality holds on a spin-$`0`$ vector
of the product of the two triplets, which exists by Chapter 4, Exercise 2 ($`1\otimes1\cong2\oplus1\oplus0`$).

*Computation:* Theorem 9.14 for the open XX chain, with `krylov.lanczos` on the sparse sector of `SpinChain.xxz_sparse`, and ground states of the Heisenberg ring up to $`L=16`$ with `krylov.ground_state` and its residual — [tutorial, Section 9.5](../tutorial/09_spin_chains.md#s-9-5).

**Time evolution.** For $`e^{-i\tau H}v`$ one projects onto $`\mathcal K_m`$: $`w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1`$.

<a id="r-9-15"></a>

**Proposition 9.15.** Let $`\dim\mathcal K_m=m`$. For every $`\tau`$, $`\|w_\tau\|=\|v\|`$ and $`\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle`$. If
$`\mathcal K_m`$ is $`H`$-invariant, $`w_\tau=e^{-i\tau H}v`$.

*Proof.* $`Q_m`$ is an isometry and $`e^{-i\tau T_m}`$ is unitary ($`T_m`$ is real symmetric), which gives the norm.
Since $`Q_m^*HQ_m=T_m`$ commutes with $`e^{-i\tau T_m}`$,
$`\langle w_\tau,Hw_\tau\rangle=\|v\|^2\langle e^{-i\tau T_m}e_1,T_me^{-i\tau T_m}e_1\rangle=\|v\|^2\langle e_1,T_me_1\rangle=\langle v,Hv\rangle`$. If
$`\mathcal K_m`$ is invariant, $`HQ_m=Q_mT_m`$ (Proposition 9.13, $`r_m=0`$), so $`H^jQ_m=Q_mT_m^j`$ and, summing the exponential series (Proposition 8.2(1)),
$`e^{-i\tau H}Q_m=Q_me^{-i\tau T_m}`$; apply to
$`e_1`$, with $`Q_me_1=v/\|v\|`$. $`\square`$

*Computation:* Proposition 9.15 for `krylov.evolve`, which steps in Krylov spaces, on a quench of the Néel state, against the dense propagator — [tutorial, Section 9.5](../tutorial/09_spin_chains.md#s-9-5).

<a id="s-9-6"></a>

## 9.6 Disorder: level statistics and entanglement

The **random-field Heisenberg chain** has independent fields $`h_i`$ uniform on $`[-W,W]`$. It conserves $`S^z`$
(Proposition 9.1), so level statistics are computed per sector $`N`$ (Theorem 5.15 and Section 5.4: mixing sectors
fakes independent levels). If two fields differ, $`THT^*-H=\sum_i(h_{i-1}-h_i)s^z_i\ne0`$, so translations are no longer
a symmetry.

*Computation:* the mean ratio $`\langle r\rangle`$ of the random-field chain, compared with GOE matrices and with independent levels (numerical observations) — [tutorial, Section 9.6](../tutorial/09_spin_chains.md#s-9-6).

**Entanglement of free fermions.** For the eigenstates of the open XX chain the entanglement entropy of
the first $`n`$ sites is determined by an $`n\times n`$ matrix.

<a id="r-9-16"></a>

**Theorem 9.16 (Peschel's formula).** Let $`\Phi=b_{k_1}^*\cdots b_{k_N}^*|0\rangle`$ as in Theorem 9.10,
$`C_{ij}=\langle\Phi,c_i^*c_j\Phi\rangle`$ and $`C_A=(C_{ij})_{i,j<n}`$. Then $`0\le C_A\le1`$, and the entanglement entropy
of the sites $`A=\{0,\dots,n-1\}`$ is, with $`\nu`$ running over the eigenvalues of $`C_A`$ counted with multiplicity,
```math
S_A=\sum_\nu \eta(\nu),\qquad \eta(\nu)=-\nu\log\nu-(1-\nu)\log(1-\nu)\quad(0\log0=0).
```

*Proof.* For a real vector $`x`$ put $`b(x)=\sum_jx_jc_j`$; then $`b(x)^*b(y)^*=-b(y)^*b(x)^*`$ and
$`\{b(x),b(y)^*\}=\langle x,y\rangle`$. Let $`U\subseteq\mathbb R^L`$ be the span of the $`u_{k_a}`$ and $`P_A`$, $`P_B=1-P_A`$ the coordinate projections on
$`\{0,\dots,n-1\}`$ and its complement. Diagonalising the form $`\langle P_Ax,P_Ay\rangle`$ on $`U`$ gives an orthonormal basis
$`\phi_1,\dots,\phi_N`$ of $`U`$ with $`\langle P_A\phi_a,P_A\phi_b\rangle=\nu_a\delta_{ab}`$, where $`\nu_a=\|P_A\phi_a\|^2\in[0,1]`$; hence
$`\langle P_B\phi_a,P_B\phi_b\rangle=(1-\nu_a)\delta_{ab}`$. Since $`(x_1,\dots,x_N)\mapsto b(x_1)^*\cdots b(x_N)^*|0\rangle`$ is
multilinear and alternating, $`\Phi=\pm\prod_a\big(\sqrt{\nu_a}f_a^*+\sqrt{1-\nu_a}\,g_a^*\big)|0\rangle`$ with
$`f_a=b(P_A\phi_a/\sqrt{\nu_a})`$, $`g_a=b(P_B\phi_a/\sqrt{1-\nu_a})`$ (terms with a zero coefficient omitted).
Expanding, $`\Phi=\sum_I\pm\prod_{a\in I}\sqrt{\nu_a}\prod_{a\notin I}\sqrt{1-\nu_a}\,F_I^*G_{I^c}^*|0\rangle`$ over the
subsets $`I`$, with $`F_I^*`$, $`G_{I^c}^*`$ the ordered products. For $`j\ge n`$ the string of $`c_j^*`$ contains all of
$`A`$: $`c_j^*=(-1)^{N_A}\otimes\tilde c_j^*`$ with Jordan–Wigner operators $`\tilde c_j`$ of the sites $`\ge n`$, while for
$`i<n`$, $`c_i^*=\hat c_i^*\otimes1`$. As $`(-1)^{N_A}`$ is $`1`$ on $`|0_A\rangle`$,
$`F_I^*G_{I^c}^*|0\rangle=\xi_I\otimes\chi_{I^c}`$ with $`\xi_I=\hat F_I^*|0_A\rangle`$, $`\chi_{I^c}=\tilde G_{I^c}^*|0_B\rangle`$, which are
orthonormal families by the CAR computation of Theorem 9.10. This is a Schmidt decomposition (Lemma 7.7) with
weights $`p_I=\prod_{a\in I}\nu_a\prod_{a\notin I}(1-\nu_a)`$ (the signs go into $`\xi_I`$), the law of independent
Bernoulli variables, so, as $`\log p_I=\sum_{a\in I}\log\nu_a+\sum_{a\notin I}\log(1-\nu_a)`$ and $`\sum_{I\ni a}p_I=\nu_a`$,
$`S_A=-\sum_Ip_I\log p_I=\sum_a\eta(\nu_a)`$. Finally $`C_{ij}=\sum_a\phi_a(i)\phi_a(j)`$
(expand $`c_i^*c_j`$ in the $`b(\phi_a)`$), so $`C_A=\sum_a(P_A\phi_a)(P_A\phi_a)^T`$ has the eigenvalues $`\nu_a`$ on the
orthogonal vectors $`P_A\phi_a`$ and $`0`$ elsewhere, and $`\eta(0)=0`$. $`\square`$

Theorem 9.16 is due to I. Peschel (*J. Phys. A* 36 (2003) L205).

*Computation:* Theorem 9.16 for the ground states of the open XX chain with $`L=12`$, against `SpinChain.entanglement_entropy` — [tutorial, Section 9.6](../tutorial/09_spin_chains.md#s-9-6).

<a id="s-9-7"></a>

## 9.7 What to compare with

*Computation:* the comparison with other software (QuSpin, exact-diagonalisation codes, tensor networks) — [tutorial, Section 9.7](../tutorial/09_spin_chains.md#s-9-7).

<a id="x-9"></a>

## Exercises

1. (Momentum sectors by hand.) For $`L=4`$ and $`N=2`$ find the orbits of the six configurations under $`T`$
   and the dimensions $`d_{2,m}`$, $`m=0,\dots,3`$.
2. (Sum rules.) Let $`\psi`$ be the ground state of the Heisenberg ring with $`L`$ even. Show that
   $`E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle`$ and $`\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0`$.
3. (Too much symmetry.) Show that in sector $`N=L/2`$ the vector with all coefficients equal to $`1`$ is a
   multiple of $`(S^-)^{L/2}|{\uparrow\cdots\uparrow}\rangle`$ and an eigenvector of the Heisenberg ring with energy
   $`JL/4`$, and that every Krylov space $`\mathcal K_m(H,v)`$ lies in each $`H`$-invariant subspace containing $`v`$.
   What does `krylov.ground_state` return when started from this vector?

<details>
<summary>Solutions</summary>

1. $`0011\to1001\to1100\to0110\to0011`$ is one orbit with $`R=4`$, and $`0101\to1010\to0101`$ one with $`R=2`$. By
   Proposition 9.5 the first contributes to every $`m`$, the second to $`m`$ with $`2m\equiv0\pmod4`$, i.e. $`m=0,2`$:
   $`d_{2,m}=2,1,2,1`$.

2. The ground state is a simple singlet (Theorem 9.9) and $`T`$ commutes with $`H`$, so $`T\psi=e^{ik}\psi`$. As
   $`\mathbf s_i\cdot\mathbf s_{i+1}=T^i(\mathbf s_0\cdot\mathbf s_1)T^{*i}`$, $`\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle`$ is the same for all $`L`$ bonds, and
   $`E_0=J\sum_i\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle`$. As $`\mathbf S^2\psi=0`$,
   $`0=\langle\mathbf S^2\rangle=\sum_{i,j}\langle\mathbf s_i\cdot\mathbf s_j\rangle=L\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle`$, by translation
   invariance again.

3. $`S^-`$ lowers one of the up spins; applying it $`L/2`$ times to the fully polarised state reaches each
   configuration with $`L/2`$ up spins in $`(L/2)!`$ orders, all with coefficient $`1`$. $`|{\uparrow\cdots\uparrow}\rangle`$ has all
   $`L`$ bonds parallel, energy $`JL/4`$, and $`S^-`$ commutes with $`H`$, so the vector is an eigenvector with the same
   energy. If $`W`$ is invariant and $`v\in W`$, then $`H^jv\in W`$ for all $`j`$. Here $`\mathcal K_m=\mathbb Cv`$ for all $`m`$, and
   `ground_state` returns $`E=JL/4`$, the top of the spectrum (each $`\mathbf s_i\cdot\mathbf s_j=\frac12F_{ij}-\frac14\le\frac14`$ by
   Proposition 4.15(a), as $`F_{ij}`$ is a symmetric permutation matrix with eigenvalues $`\pm1`$). A random start vector
   avoids this trap: a Gaussian $`v`$ has a density, so it has a nonzero component $`P_1v`$ in the eigenspace of
   $`\lambda_1`$ almost surely (Corollary A.7(c); for a real $`v`$ apply it in $`W=\mathbb R^n`$ to the proper subspace
   $`\{v:P_1v=0\}`$, proper because the real matrix $`P_1\ne0`$ has a nonzero column), and then by Theorem 9.14 an
   $`H`$-invariant Krylov space $`\mathcal K_m(H,v)`$ has the lowest Ritz value $`\lambda_1`$, whereas here the invariant
   $`\mathcal K_1=\mathbb Cv`$ has the Ritz value $`JL/4`$.

</details>

*Computation:* Exercises 1–3 numerically — [tutorial, Exercises](../tutorial/09_spin_chains.md#x-9).

**Next:** [Chapter 10](10_quantum_computing.md) treats quantum computing with sectors: decoherence-free
subspaces and noiseless subsystems, operator-algebra error correction, and mixers that respect a conserved
Hamming weight — the algebra $`\bigoplus_wM_{\binom nw}`$ of this chapter in a new role.
