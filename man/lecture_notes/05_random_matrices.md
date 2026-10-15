**English** · [Español](es/05_random_matrices.md) · [Français](fr/05_random_matrices.md) · [Deutsch](de/05_random_matrices.md) · [中文](zh/05_random_matrices.md) · [日本語](ja/05_random_matrices.md) · [Русский](ru/05_random_matrices.md)

# 5. Random matrices in an algebra

Random matrix theory asks what a "typical" operator looks like. In an algebra
$`\bigoplus_cM_{k_c}`$ the answer has to respect the sectors: the natural random unitaries are
those of the algebra, not those of the big matrix algebra $`M_D`$, and the natural random
Hamiltonians have independent blocks. This chapter constructs the Haar measure on the unitary
group of the algebra, computes its first and second moments exactly, and then shows what goes
wrong in spectral statistics when the sectors are ignored.

**You will learn** why the QR decomposition of a Gaussian matrix gives a Haar unitary (and why
the phases must be fixed), that averaging over the unitary group of the algebra is the
conditional expectation onto the centre, the Weingarten formula of order two, the Wigner
surmise, the spacing ratio and its exact value for $`3\times3`$ Gaussian matrices, and why level statistics computed on a spectrum with mixed
sectors look like those of independent levels, the more so the more sectors are mixed.

*Computation:* the numerical checks of this chapter — [Tutorial 5](../tutorial/05_random_matrices.md).

Throughout, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ with $`D=\sum_ck_c`$, $`P_c`$ are the
minimal central projections, and $`E`$ is the conditional expectation onto the centre of
Proposition 2.8. $`U(k)`$ is the group of unitary $`k\times k`$ matrices and
$`U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}`$ the unitary group of $`\mathcal M`$. Random variables
live on a probability space with expectation $`\mathbb E`$; the expectation of a random matrix is
taken entrywise. A *standard complex Gaussian* variable has independent real and imaginary parts
$`\sim N(0,\tfrac12)`$, so that $`\mathbb E|z|^2=1`$; a *Ginibre matrix* has independent standard complex Gaussian entries.

<a id="s-5-1"></a>

## 5.1 Haar measure

Let $`G`$ be a compact group of matrices (a closed subgroup of some $`U(n)`$). A Borel probability
measure $`\mu`$ on $`G`$ is **left-invariant** if $`\mu(gB)=\mu(B)`$ for all $`g\in G`$ and Borel sets
$`B`$, equivalently $`\int f(gh)\,d\mu(h)=\int f\,d\mu`$ for every bounded Borel function $`f`$;
**right-invariant** is defined with $`Bg`$ and $`f(hg)`$.

<a id="r-5-1"></a>

**Lemma 5.1 (uniqueness).** If $`\mu`$ is a left-invariant and $`\nu`$ a right-invariant Borel
probability measure on $`G`$, then $`\mu=\nu`$. Hence $`G`$ carries at most one left-invariant Borel
probability measure, and if it exists it is also right-invariant and invariant under $`g\mapsto g^{-1}`$.

*Proof.* Let $`f:G\to\mathbb R`$ be continuous (hence bounded). By Fubini's theorem,
$`\int\!\big(\int f(gh)\,d\mu(h)\big)d\nu(g)=\int\!\big(\int f(gh)\,d\nu(g)\big)d\mu(h)`$. The inner
integral on the left is $`\int f\,d\mu`$ by left invariance, the one on the right is $`\int f\,d\nu`$
by right invariance; so $`\int f\,d\mu=\int f\,d\nu`$. Two Borel probability measures on a metric
space with the same integrals of continuous functions are equal (Corollary A.2: the continuous functions
$`\max(0,1-n\mathrm{dist}(\cdot,K))`$ decrease to the indicator of a closed set $`K`$, and the closed
sets determine the measure by Dynkin's lemma, Lemma A.1). If $`\mu_1,\mu_2`$ are left-invariant, then
$`\nu(B)=\mu_2(B^{-1})`$ is right-invariant, since $`(Bg)^{-1}=g^{-1}B^{-1}`$; hence $`\mu_1=\nu=\mu_2`$,
so $`\mu_1`$ is right-invariant and inversion-invariant. $`\square`$

**Definition.** The **Haar measure** of $`G`$ is its left-invariant Borel probability measure. A
random element with this law is **Haar distributed**.

Lemma 5.1 says nothing about existence. The course needs the Haar measures of $`U(k)`$, $`U(1)`$,
$`U(\mathcal M)`$, $`SU(k)`$, $`O(k)`$ and of the diagonal groups $`T`$ of Lemma C.2 only, and each of them is
constructed explicitly: $`U(k)`$ in Theorem 5.3, $`U(\mathcal M)`$ in Corollary 5.4, $`SU(k)`$ after it,
$`O(k)`$ in Theorem B.17, $`T`$ in Lemma C.2 (see Section B.7).

<a id="r-5-2"></a>

**Lemma 5.2 (QR decomposition).** Every invertible $`Z\in M_k(\mathbb C)`$ can be written $`Z=QR`$
with $`Q`$ unitary and $`R`$ upper triangular with positive diagonal. $`Q`$ and $`R`$ are unique and
depend continuously on $`Z`$.

*Proof.* Apply Gram–Schmidt to the columns $`z_1,\dots,z_k`$ of $`Z`$, which are linearly independent:
$`q_j=w_j/\|w_j\|`$ with $`w_j=z_j-\sum_{i<j}\langle q_i,z_j\rangle q_i`$; here $`w_j\neq0`$ because
$`z_j\notin\mathrm{span}\{z_1,\dots,z_{j-1}\}=\mathrm{span}\{q_1,\dots,q_{j-1}\}`$. Then
$`z_j=\sum_{i<j}\langle q_i,z_j\rangle q_i+\|w_j\|q_j`$, i.e. $`Z=QR`$ with $`R_{jj}=\|w_j\|>0`$, and
every step is continuous on the invertible matrices. If $`Q_1R_1=Q_2R_2`$, then
$`T=Q_2^*Q_1=R_2R_1^{-1}`$ is unitary and upper triangular with positive diagonal. Its inverse
$`T^*`$ is lower triangular, and also upper triangular (the inverse of an upper triangular matrix
is upper triangular), so $`T`$ is diagonal; a diagonal unitary with positive diagonal is $`1`$. $`\square`$

<a id="r-5-3"></a>

**Theorem 5.3 (Haar unitaries from Gaussian matrices).** Let $`Z`$ be a $`k\times k`$
Ginibre matrix. Then $`Z`$ is invertible with probability one, and the factor $`Q`$ of Lemma 5.2 is
Haar distributed on $`U(k)`$. Moreover, if $`Z=Q'R'`$ is any QR decomposition with $`Q'`$ unitary and $`R'`$
upper triangular (for instance the one returned by `torch.linalg.qr`), then
$`Q=Q'\Lambda`$ with $`\Lambda=\mathrm{diag}(R'_{ii}/|R'_{ii}|)`$ is Haar distributed.

*Proof.* *Invertibility.* $`|\det Z|^2`$ is a real polynomial in the $`2k^2`$ real coordinates of $`Z`$, not
identically zero (it equals $`1`$ at $`Z=1`$), and $`Z`$ has the density $`\pi^{-k^2}e^{-\|Z\|_2^2}`$
(Proposition A.16(b)); so $`P(\det Z=0)=0`$ by Corollary A.7(b).

*Invariance of the Gaussian law.* For $`V\in U(k)`$ the map $`Z\mapsto VZ`$ is real-linear and
preserves $`\|Z\|_2`$ (Section 1.4), so it is an orthogonal map of $`M_k(\mathbb C)\cong\mathbb R^{2k^2}`$, which
preserves Lebesgue measure (Lemma A.3(b)) and the density: $`VZ`$ has the same law as $`Z`$
(Proposition A.16(b)).

*Equivariance.* If $`Z=QR`$, then $`VZ=(VQ)R`$ is the decomposition of Lemma 5.2 of $`VZ`$, so
$`Q(VZ)=VQ(Z)`$. As $`Q`$ is continuous on the invertible matrices (Lemma 5.2), $`Q(Z)`$ is a random
element of $`U(k)`$, and for every Borel $`B\subseteq U(k)`$,
$`P(VQ(Z)\in B)=P(Q(VZ)\in B)=P(Q(Z)\in B)`$: the law of $`Q(Z)`$ is left-invariant, and by
Lemma 5.1 it is the Haar measure (this also proves that the Haar measure of $`U(k)`$ exists).
Finally $`Q'R'=(Q'\Lambda)(\Lambda^*R')`$, and $`\Lambda^*R'`$ is upper triangular with diagonal
$`|R'_{ii}|>0`$, so $`Q'\Lambda=Q(Z)`$ by uniqueness. $`\square`$

The same proof with real Gaussian matrices gives the Haar measure of the orthogonal group $`O(k)`$,
with the correction $`\mathrm{diag}(\mathrm{sgn}\,R'_{ii})`$; this is written out in Theorem B.17.
Without the correction $`\Lambda`$ the law of $`Q'`$ depends on the sign convention of the QR routine and
need not be the Haar measure (the recipe with $`\Lambda`$ is due to F. Mezzadri, *How to
generate random matrices from the classical compact groups*, Notices AMS 54 (2007), 592–604). For example, if the routine makes $`R'_{11}`$ real with the sign opposite to
$`\mathrm{Re}\,Z_{11}`$ (as the Householder QR of LAPACK does), the first columns $`z_1=R'_{11}q'_1`$ of $`Z=Q'R'`$ give
$`\mathrm{Re}\,Q'_{11}=-|\mathrm{Re}\,Z_{11}|/\|z_1\|<0`$ almost surely ($`\mathrm{Re}\,Z_{11}\sim N(0,\frac12)`$ vanishes with probability $`0`$);
so $`\mathbb E\,\mathrm{Re}\,Q'_{11}<0`$, whereas $`\mathbb EU_{11}=0`$ for a Haar unitary (Corollary 5.5).

*Computation:* the library's sampler `TypeIAlgebra.random_unitary` is the corrected QR of Theorem 5.3, and the uncorrected $`Q'`$ of `torch.linalg.qr` is not Haar distributed — [tutorial, Section 5.1](../tutorial/05_random_matrices.md#s-5-1).

<a id="r-5-4"></a>

**Corollary 5.4 (Haar measure of the algebra).** $`U(\mathcal M)=\{\bigoplus_cU_c:U_c\in U(k_c)\}\cong\prod_cU(k_c)`$,
and $`U=\bigoplus_cU_c`$ is Haar distributed on $`U(\mathcal M)`$ iff the blocks $`U_c`$ are independent
and each $`U_c`$ is Haar distributed on $`U(k_c)`$.

*Proof.* $`U^*U=\bigoplus_cU_c^*U_c`$, so $`U`$ is unitary iff every block is. Let the $`U_c`$ be
independent and Haar distributed. For $`V=\bigoplus_cV_c\in U(\mathcal M)`$, the blocks $`V_cU_c`$ of
$`VU`$ are independent (functions of independent variables) and $`V_cU_c`$ has the law of $`U_c`$, so
$`VU`$ has the law of $`U`$. The law of $`U`$ is therefore left-invariant, and by Lemma 5.1 it is the Haar
measure of the compact group $`U(\mathcal M)`$. Conversely, by uniqueness every Haar distributed $`U`$ has
this law. $`\square`$

*In the library:* Haar unitaries of $`U(\mathcal M)`$, with independent blocks, are drawn by `TypeIAlgebra.random_unitary_operator` — [tutorial, Section 5.1](../tutorial/05_random_matrices.md#s-5-1).

The group $`U(\mathcal M)`$ is much smaller than $`U(D)`$: for $`C\ge2`$ a Haar unitary $`U`$ of $`U(D)`$ lies
in $`\mathcal M`$ with probability zero. Indeed, by Theorem 5.3 and Lemma 5.1, $`U`$ has the law of the
factor $`Q`$ of a $`D\times D`$ Ginibre matrix $`Z`$, and by the Gram–Schmidt construction in Lemma 5.2 the
first column of $`Q`$ is $`z/\|z\|`$ for the first column $`z`$ of $`Z`$, a standard complex Gaussian vector of
$`\mathbb C^D`$ with a density (Proposition A.16(a)). For each $`i`$, $`\{z_i=0\}`$ is a proper subspace of
$`\mathbb C^D`$, so $`P(z_i=0)=0`$ (Corollary A.7(c)), and all entries of the first column are nonzero almost
surely. An element of $`\mathcal M`$ has zero entries in the first column outside the rows of the
first block, and these rows exist when $`C\ge2`$.

The Haar measure of $`SU(k)`$ is obtained in the same spirit: if $`U`$ is Haar on $`U(k)`$ and $`\theta=\arg\det U\in(-\pi,\pi]`$ (a Borel function of $`U`$), then $`V=e^{-i\theta/k}U\in SU(k)`$, and for $`W\in SU(k)`$,
$`\det(WU)=\det U`$ gives $`V(WU)=WV(U)`$; since $`WU`$ has the law of $`U`$, the law of $`V`$ is
left-invariant on $`SU(k)`$, hence the Haar measure of $`SU(k)`$ by Lemma 5.1.

<a id="r-5-5"></a>

**Corollary 5.5 (first moments of the entries).** If $`U`$ is Haar distributed on $`U(k)`$, then
$`\mathbb EU_{ij}=0`$ and $`\mathbb E|U_{ij}|^2=1/k`$.

*Proof.* $`-1\in U(k)`$, so $`-U`$ has the law of $`U`$ and $`\mathbb EU=-\mathbb EU`$. For a permutation
matrix $`P`$, $`PU`$ has the law of $`U`$, so $`\mathbb E|U_{ij}|^2`$ does not depend on $`i`$; as
$`\sum_i|U_{ij}|^2=1`$, it equals $`1/k`$. $`\square`$

*Computation:* Corollaries 5.4 and 5.5 by Monte Carlo, and the Haar measure of $`SU(k)`$ (`measure='haar_su'`) — [tutorial, Section 5.1](../tutorial/05_random_matrices.md#s-5-1).

<a id="s-5-2"></a>

## 5.2 Moments and the Weingarten calculus

<a id="r-5-6"></a>

**Lemma 5.6.** If $`T\in M_k(\mathbb C)`$ commutes with all diagonal unitaries and all permutation
matrices, then $`T\in\mathbb C1`$. In particular the commutant of $`U(k)`$ is $`\mathbb C1`$.

*Proof.* For $`\Delta=\mathrm{diag}(e^{i\theta_1},\dots,e^{i\theta_k})`$ with distinct
$`e^{i\theta_j}`$, $`\Delta T=T\Delta`$ gives $`(e^{i\theta_i}-e^{i\theta_j})T_{ij}=0`$, so $`T`$ is
diagonal. For the permutation matrix $`P_\sigma e_j=e_{\sigma(j)}`$, $`P_\sigma TP_\sigma^{-1}=T`$ gives
$`T_{\sigma(j)\sigma(j)}=T_{jj}`$, so all diagonal entries are equal. $`\square`$

<a id="r-5-7"></a>

**Proposition 5.7 (the Haar twirl is the conditional expectation onto the centre).** Let $`U`$ be
Haar distributed on $`U(\mathcal M)`$. For every $`X\in\mathcal M`$,
```math
\mathbb E\,[UXU^*]=E(X)=\bigoplus_c\frac{\mathrm{Tr}\,X_c}{k_c}\,1_{k_c}.
```

*Proof.* $`(UXU^*)_c=U_cX_cU_c^*`$, and $`U_c`$ is Haar on $`U(k_c)`$ (Corollary 5.4). Let
$`T_c=\mathbb E[U_cX_cU_c^*]`$ (the entries are bounded). For $`V\in U(k_c)`$, $`VU_c`$ has the law of
$`U_c`$, so $`VT_cV^*=\mathbb E[(VU_c)X_c(VU_c)^*]=T_c`$. By Lemma 5.6, $`T_c=\lambda1`$, and
$`k_c\lambda=\mathrm{Tr}\,T_c=\mathbb E\mathrm{Tr}(U_cX_cU_c^*)=\mathrm{Tr}\,X_c`$. $`\square`$

The normalisation is that of Proposition 2.8 and of `channels.center_expectation`. If
$`X\in\mathcal M`$ commutes with $`U(\mathcal M)`$, then $`UXU^*=X`$, so $`X=E(X)`$ is central; the Haar
average is the continuous analogue of the group average of Chapter 2, Exercise 3. If instead $`U`$ is
Haar on $`U(D)`$, the same proof (with $`C=1`$) gives $`\mathbb E[UXU^*]=\frac{\mathrm{Tr}\,X}{D}1`$: the
sector weights are lost, and $`UXU^*`$ need not even lie in $`\mathcal M`$ (for $`C\ge2`$, $`X=P_1`$ and the rotation $`U`$
by $`\pi/4`$ in the plane of $`e_1`$ and $`e_{k_1+1}`$, $`UP_1U^*`$ has the entry $`\frac12`$ in position $`(1,k_1+1)`$).

*Computation:* Proposition 5.7 by Monte Carlo, against `channels.center_expectation` — [tutorial, Section 5.2](../tutorial/05_random_matrices.md#s-5-2).

> **Physical meaning** (physical interpretation, not a theorem of this course). A Haar-random unitary of the algebra is a random evolution that respects
> the superselection rule. Averaging over it scrambles everything inside each sector and keeps
> exactly the classical part: the weight of each charge. This is the "twirl" of quantum
> information, and it is the conditional expectation onto the centre of Chapter 2.

For second moments we need the commutant of $`\{V\otimes V\}`$ on $`\mathbb C^k\otimes\mathbb C^k`$.
Let $`F`$ be the **swap**, $`F(x\otimes y)=y\otimes x`$; it is self-adjoint and unitary, and
$`\mathrm{Tr}\,1=k^2`$, $`\mathrm{Tr}\,F=\sum_{i,j}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k`$, $`F^2=1`$.

<a id="r-5-8"></a>

**Lemma 5.8.** Let $`k\ge2`$. A matrix $`X`$ on $`\mathbb C^k\otimes\mathbb C^k`$ commutes with $`V\otimes V`$
for all $`V\in U(k)`$ iff $`X=a1+bF`$ with $`a,b\in\mathbb C`$; and $`1,F`$ are linearly independent.

*Proof.* $`F(V\otimes V)(x\otimes y)=Vy\otimes Vx=(V\otimes V)F(x\otimes y)`$, so $`1`$ and $`F`$ commute
with all $`V\otimes V`$. Conversely let $`X`$ commute with them, and write
$`X_{ij,pq}=\langle e_i\otimes e_j,X(e_p\otimes e_q)\rangle`$.
*Diagonal unitaries:* $`\Delta\otimes\Delta`$ multiplies $`e_i\otimes e_j`$ by $`e^{i(\theta_i+\theta_j)}`$,
so $`(e^{i(\theta_i+\theta_j)}-e^{i(\theta_p+\theta_q)})X_{ij,pq}=0`$. If the multisets $`\{i,j\}`$ and
$`\{p,q\}`$ differ, some index $`m`$ occurs $`n\in\{\pm1,\pm2\}`$ times more often in the first; with
$`\theta_m=\pi/2`$ and all other $`\theta=0`$ the bracket is $`e^{i(\theta_p+\theta_q)}(e^{in\pi/2}-1)\neq0`$,
so $`X_{ij,pq}=0`$. Hence $`X(e_i\otimes e_i)=a_ie_i\otimes e_i`$, and for $`i\neq j`$,
$`X(e_i\otimes e_j)=b_{ij}\,e_i\otimes e_j+c_{ij}\,e_j\otimes e_i`$.
*Permutations:* applying $`X(P_\sigma\otimes P_\sigma)=(P_\sigma\otimes P_\sigma)X`$ to $`e_i\otimes e_j`$ gives
$`a_{\sigma(i)}=a_i`$, $`b_{\sigma(i)\sigma(j)}=b_{ij}`$, $`c_{\sigma(i)\sigma(j)}=c_{ij}`$. The permutations act
transitively on ordered pairs of distinct indices, so $`a_i=a'`$, $`b_{ij}=b'`$, $`c_{ij}=c'`$; thus
$`X=b'1+c'F`$ on the span of the $`e_i\otimes e_j`$ with $`i\ne j`$.
*A rotation:* let $`R`$ be the unitary with $`Re_1=(e_1+e_2)/\sqrt2`$, $`Re_2=(-e_1+e_2)/\sqrt2`$, $`Re_m=e_m`$ for $`m\ge3`$.
Then $`(R\otimes R)(e_1\otimes e_1)=\frac12(e_1\otimes e_1+e_1\otimes e_2+e_2\otimes e_1+e_2\otimes e_2)`$, and
$`X(R\otimes R)(e_1\otimes e_1)=\frac12\big(a'e_1\otimes e_1+(b'+c')(e_1\otimes e_2+e_2\otimes e_1)+a'e_2\otimes e_2\big)`$
must equal $`(R\otimes R)X(e_1\otimes e_1)=a'(R\otimes R)(e_1\otimes e_1)`$. Comparing the coefficients of
$`e_1\otimes e_2`$ gives $`a'=b'+c'`$, which is the value of $`b'1+c'F`$ on $`e_i\otimes e_i`$. So $`X=b'1+c'F`$.
Finally $`1`$ and $`F`$ map $`e_1\otimes e_2`$ to the independent vectors $`e_1\otimes e_2`$ and $`e_2\otimes e_1`$. $`\square`$

*Computation:* Lemma 5.8 by linear algebra, with the commutant of three sampled $`V\otimes V`$ — [tutorial, Section 5.2](../tutorial/05_random_matrices.md#s-5-2).

<a id="r-5-9"></a>

**Corollary 5.9 (Schur–Weyl duality of order two).** Let $`k\ge2`$ and $`P_\pm=\frac12(1\pm F)`$.
They are the projections onto the symmetric and antisymmetric subspaces
$`\mathrm{Sym}=\{x:Fx=x\}`$ and $`\mathrm{Alt}=\{x:Fx=-x\}`$, of dimensions $`k(k+1)/2`$ and $`k(k-1)/2`$.
The representation $`V\mapsto V\otimes V`$ of $`U(k)`$ is irreducible on $`\mathrm{Sym}`$ and on $`\mathrm{Alt}`$,
and the two representations are inequivalent.

*Proof.* $`F=F^*`$ and $`F^2=1`$ make $`P_\pm`$ projections with $`P_++P_-=1`$ and ranges $`\mathrm{Sym}`$,
$`\mathrm{Alt}`$; their dimensions are $`\mathrm{Tr}\,P_\pm=(k^2\pm k)/2`$. The span $`\mathcal A`$ of
$`\{V\otimes V\}`$ is a $`*`$-algebra ($`1\otimes1`$, products $`VW\otimes VW`$, adjoints $`V^*\otimes V^*`$) with
commutant $`\mathcal A'=\mathrm{span}\{1,F\}=\mathrm{span}\{P_+,P_-\}`$ by Lemma 5.8. Let
$`K\subseteq\mathrm{Sym}`$ be invariant under all $`V\otimes V`$. By Lemma 2.1 its projection is in $`\mathcal A'`$, so it is
$`xP_++yP_-`$ with $`x,y\in\{0,1\}`$ (it is a projection); it vanishes on $`\mathrm{Alt}`$, so $`y=0`$, and
$`K=0`$ or $`K=\mathrm{Sym}`$. The same argument applies to $`\mathrm{Alt}`$. If
$`T:\mathrm{Sym}\to\mathrm{Alt}`$ intertwines the two representations, then $`TP_+`$ (viewed as an
operator on $`\mathbb C^k\otimes\mathbb C^k`$) commutes with $`\mathcal A`$, so $`TP_+=xP_++yP_-`$; it
vanishes on $`\mathrm{Alt}`$, so $`y=0`$, and its range lies in $`\mathrm{Alt}\cap\mathrm{Sym}=0`$, so $`T=0`$. $`\square`$

<a id="r-5-10"></a>

**Theorem 5.10 (Weingarten formula of order two).** Let $`k\ge2`$, $`U`$ Haar distributed on $`U(k)`$,
and $`X`$ a matrix on $`\mathbb C^k\otimes\mathbb C^k`$. Then
```math
\mathbb E\,[(U\otimes U)X(U\otimes U)^*]=\alpha\,1+\beta F,\qquad
\alpha=\frac{\mathrm{Tr}\,X-\mathrm{Tr}(FX)/k}{k^2-1},\quad
\beta=\frac{\mathrm{Tr}(FX)-\mathrm{Tr}\,X/k}{k^2-1}.
```
Equivalently, with $`\mathrm{Wg}(e)=\frac1{k^2-1}`$, $`\mathrm{Wg}(\tau)=-\frac1{k(k^2-1)}`$ on the
group $`S_2=\{e,\tau\}`$,
```math
\mathbb E\big[U_{i_1j_1}U_{i_2j_2}\overline{U_{i'_1j'_1}U_{i'_2j'_2}}\big]
=\sum_{\sigma,\pi\in S_2}\delta_{i_1i'_{\sigma(1)}}\delta_{i_2i'_{\sigma(2)}}\,
\delta_{j_1j'_{\pi(1)}}\delta_{j_2j'_{\pi(2)}}\,\mathrm{Wg}(\sigma\pi^{-1}).
```

*Proof.* Let $`T=\mathbb E[(U\otimes U)X(U\otimes U)^*]`$. As in Proposition 5.7, $`T`$ commutes with
all $`V\otimes V`$, so $`T=\alpha1+\beta F`$ by Lemma 5.8. Taking traces,
$`\mathrm{Tr}\,T=\mathrm{Tr}\,X`$, and, since $`F`$ commutes with $`U\otimes U`$,
$`\mathrm{Tr}(FT)=\mathbb E\mathrm{Tr}((U\otimes U)^*F(U\otimes U)X)=\mathrm{Tr}(FX)`$.
This gives $`k^2\alpha+k\beta=\mathrm{Tr}\,X`$ and $`k\alpha+k^2\beta=\mathrm{Tr}(FX)`$, a
linear system with determinant $`k^2(k^2-1)\neq0`$, whose solution is the one stated. For the
entries take $`X=|e_{j_1}\otimes e_{j_2}\rangle\langle e_{j'_1}\otimes e_{j'_2}|`$: then
$`\langle e_{i_1}\otimes e_{i_2},T\,e_{i'_1}\otimes e_{i'_2}\rangle`$ is the left side, $`\mathrm{Tr}\,X=\delta_{j_1j'_1}\delta_{j_2j'_2}`$
(the term $`\pi=e`$), $`\mathrm{Tr}(FX)=\delta_{j_1j'_2}\delta_{j_2j'_1}`$ (the term $`\pi=\tau`$), and
$`1`$, $`F`$ contribute $`\delta_{i_1i'_1}\delta_{i_2i'_2}`$ (the term $`\sigma=e`$) and
$`\delta_{i_1i'_2}\delta_{i_2i'_1}`$ ($`\sigma=\tau`$). Collecting terms,
$`\alpha=\sum_\pi\mathrm{Wg}(e\pi^{-1})[\pi]`$ and $`\beta=\sum_\pi\mathrm{Wg}(\tau\pi^{-1})[\pi]`$,
where $`[\pi]`$ denotes the $`j`$-deltas of $`\pi`$. $`\square`$

The matrix $`\big(\mathrm{Wg}(\sigma^{-1}\pi)\big)_{\sigma,\pi}`$ is the inverse of the Gram matrix
$`\begin{pmatrix}k^2&k\\k&k^2\end{pmatrix}`$ of $`1`$ and $`F`$ with respect to
$`\langle A,B\rangle=\mathrm{Tr}(A^*B)`$; this is the form that generalises.

*Computation:* Theorem 5.10 by Monte Carlo for $`k=3`$ — [tutorial, Section 5.2](../tutorial/05_random_matrices.md#s-5-2).

<a id="r-5-11"></a>

**Corollary 5.11.** For $`k\ge2`$:
```math
\mathbb E|U_{11}|^4=\frac2{k(k+1)},\qquad
\mathbb E|U_{11}|^2|U_{12}|^2=\frac1{k(k+1)},\qquad
\mathbb E\big[U_{11}U_{22}\overline{U_{12}U_{21}}\big]=-\frac1{(k-1)k(k+1)}.
```

*Proof.* In Theorem 5.10, (i) all indices equal $`1`$: every product of deltas is $`1`$, so the
moment is $`2(\mathrm{Wg}(e)+\mathrm{Wg}(\tau))=2\frac{k-1}{k(k^2-1)}`$. (ii) $`i=(1,1)`$, $`j=(1,2)`$,
$`i'=(1,1)`$, $`j'=(1,2)`$: only $`\pi=e`$ survives, giving $`\mathrm{Wg}(e)+\mathrm{Wg}(\tau)`$.
(iii) $`i=(1,2)`$, $`j=(1,2)`$, $`i'=(1,2)`$, $`j'=(2,1)`$: only $`\sigma=e`$, $`\pi=\tau`$ survive,
giving $`\mathrm{Wg}(\tau)`$. $`\square`$

*Computation:* the three moments of Corollary 5.11 by Monte Carlo for $`k=3`$ — [tutorial, Section 5.2](../tutorial/05_random_matrices.md#s-5-2).

The last moment is negative: the entries of a Haar unitary are not independent, and the
correlations are of relative order $`1/k`$. For the algebra, entries of different blocks are
independent (Corollary 5.4), so mixed moments across sectors factorise into moments of single blocks.

For moments of any order the same reasoning applies once the commutant of
$`\{V^{\otimes n}\}`$ is known; it is spanned by the permutations of the $`n`$ tensor factors
(Schur–Weyl duality, Theorem B.6(c)). Appendix B carries this out.

<a id="r-5-12"></a>

**Theorem 5.12 (Weingarten formula of every order).** Let $`n\ge1`$ and $`U`$ Haar on $`U(k)`$. Then
```math
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\pi\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\pi(m)}}\;\mathrm{Wg}(\sigma\pi^{-1},k),
```
where $`\mathrm{Wg}(\cdot,k)`$ is a function on $`S_n`$ that depends only on the cycle type. For $`n\le k`$
the matrix $`\big(\mathrm{Wg}(\sigma^{-1}\pi,k)\big)_{\sigma,\pi}`$ is the inverse of the Gram
matrix $`\big(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)}\big)_{\sigma,\pi}`$; for $`n>k`$ this matrix is singular
and its pseudo-inverse is taken. Moments with different numbers of $`U`$ and $`\overline U`$ vanish.

*Proof.* This is Theorem B.13 with $`d=k`$: part (a) gives $`W_{\sigma\pi}=\mathrm{Wg}(\sigma^{-1}\pi)`$ with $`\mathrm{Wg}`$ a
function of the cycle type, part (b) the invertibility of the Gram matrix iff $`n\le k`$, part (d) the
formula for the moments and part (e) the vanishing of unbalanced moments. $`\square`$

(The formula is due to D. Weingarten, J. Math. Phys. 19 (1978), 999–1001, and B. Collins, Int. Math. Res.
Not. 2003, 953–982; see also B. Collins, P. Śniady, Comm. Math. Phys. 264 (2006), 773–795.)

*Computation:* the Weingarten matrix of order $`3`$ for $`k=3`$ and $`\mathbb E|U_{11}|^6=\frac1{10}`$ (Exercise 2) — [tutorial, Exercises](../tutorial/05_random_matrices.md#x-5).

<a id="s-5-3"></a>

## 5.3 Level statistics

Let $`H`$ be a random self-adjoint $`n\times n`$ matrix with eigenvalues
$`\lambda_1\le\dots\le\lambda_n`$. Its **spacings** are $`s_i=\lambda_{i+1}-\lambda_i`$. Two ensembles
are standard: the **GOE**, $`H=(G+G^T)/2`$ with $`G`$ real with independent $`N(0,1)`$ entries, and the
**GUE**, $`H=(G+G^*)/\sqrt2`$ with $`G`$ Ginibre (Section 1.7). In both, the diagonal entries are
$`N(0,1)`$, the entries above the diagonal are $`N(0,\tfrac12)`$ in the GOE and have independent real and
imaginary parts $`N(0,\tfrac12)`$ in the GUE, the entries $`h_{ij}`$, $`i\le j`$, are independent, and the
density of $`H`$ is proportional to $`e^{-\mathrm{Tr}\,H^2/2}`$ (Proposition A.17; the proportionality
comes from $`\mathrm{Tr}\,H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2`$).

<a id="r-5-13"></a>

**Proposition 5.13 (Wigner surmise).** Let $`n=2`$, $`s=\lambda_2-\lambda_1`$ and $`S=s/\mathbb Es`$.
For the GOE ($`\beta=1`$) and the GUE ($`\beta=2`$), $`S`$ has the densities
```math
p_1(S)=\frac\pi2\,S\,e^{-\pi S^2/4},\qquad p_2(S)=\frac{32}{\pi^2}\,S^2\,e^{-4S^2/\pi},
```
with distribution functions $`P(S\le t)=1-e^{-\pi t^2/4}`$ and
$`P(S\le t)=\mathrm{erf}(2t/\sqrt\pi)-\frac{4t}\pi e^{-4t^2/\pi}`$. In particular
$`P(S\le\varepsilon)`$ is of order $`\varepsilon^{\beta+1}`$ as $`\varepsilon\to0`$.

*Proof.* For $`H=\begin{pmatrix}a&b\\\bar b&c\end{pmatrix}`$ the eigenvalues are
$`\frac{a+c}2\pm\frac12\sqrt{(a-c)^2+4|b|^2}`$, so $`s=\sqrt{(a-c)^2+4|b|^2}`$. By Proposition A.17, $`a`$, $`c`$ and
the real coordinates of $`b`$ are independent, $`a,c\sim N(0,1)`$ and each real coordinate of $`b`$ is $`N(0,\frac12)`$; so
$`a-c\sim N(0,2)`$ (Proposition A.14(c)) and $`2\mathrm{Re}\,b`$, $`2\mathrm{Im}\,b`$ are $`N(0,2)`$, all independent.
*GOE:* $`x=a-c`$ and $`y=2b`$ are independent $`N(0,2)`$ variables, and in polar coordinates (Corollary A.19 with
$`n=2`$, $`s_2=2\pi`$)
$`P(s\le t)=\int_{x^2+y^2\le t^2}\frac{e^{-(x^2+y^2)/4}}{4\pi}\,dx\,dy=\int_0^t\frac r2e^{-r^2/4}dr=1-e^{-t^2/4}`$.
Then $`\mathbb Es=\int_0^\infty P(s>t)\,dt=\int_0^\infty e^{-t^2/4}dt=\sqrt\pi`$, and
$`P(S\le t)=P(s\le\sqrt\pi t)=1-e^{-\pi t^2/4}`$; differentiate. *GUE:* $`x=a-c`$, $`y=2\mathrm{Re}\,b`$,
$`y'=2\mathrm{Im}\,b`$ are independent $`N(0,2)`$, and by Theorem A.20 ($`n=3`$, $`s_3=4\pi`$) the density of
$`s=\sqrt{x^2+y^2+y'^2}`$ is $`4\pi r^2(4\pi)^{-3/2}e^{-r^2/4}=\frac{r^2e^{-r^2/4}}{2\sqrt\pi}`$. Its mean is
$`\frac1{2\sqrt\pi}\int_0^\infty r^3e^{-r^2/4}dr=\frac{8}{2\sqrt\pi}=\frac4{\sqrt\pi}`$. Substituting
$`r=4S/\sqrt\pi`$ gives $`p_2`$, and integrating $`p_2`$ by parts
($`\int x^2e^{-ax^2}dx=-\frac{x}{2a}e^{-ax^2}+\frac1{2a}\int e^{-ax^2}dx`$ with $`a=4/\pi`$) gives the
distribution function. The small-$`\varepsilon`$ behaviour follows from $`p_\beta(S)\sim cS^\beta`$. $`\square`$

*Computation:* Proposition 5.13 by Monte Carlo for the $`2\times2`$ GOE and GUE — [tutorial, Section 5.3](../tutorial/05_random_matrices.md#s-5-3).

The exponent counts dimensions: a real symmetric $`2\times2`$ matrix is degenerate iff two real
conditions hold ($`a=c`$, $`b=0`$), a complex Hermitian one iff three hold. For every $`n`$ the ordered
eigenvalues $`\lambda_1<\dots<\lambda_n`$ of the GOE ($`\beta=1`$) and the GUE ($`\beta=2`$) have the joint density
```math
Z_{n,\beta}^{-1}\prod_{i<j}(\lambda_j-\lambda_i)^\beta\,e^{-\sum_i\lambda_i^2/2}
```
(Theorem C.7, proved in Appendix C with the Weyl integration formula), and the smallest spacing $`\delta`$
satisfies $`c\,\varepsilon^{\beta+1}\le P(\delta<\varepsilon)\le C\varepsilon^{\beta+1}`$ for $`0<\varepsilon\le1`$, with constants
$`0<c\le C`$ depending on $`n`$ and $`\beta`$ (Proposition C.9). The limit $`n\to\infty`$ is not treated in this course
(see M. L. Mehta, *Random Matrices*, 3rd ed., Elsevier, 2004).

> **Physical meaning** (physical interpretation, not a theorem of this course). Wigner proposed that
> the energy levels of complex quantum systems (heavy nuclei) are statistically like the eigenvalues of a random Hamiltonian with the right symmetry:
> real symmetric if there is a time-reversal symmetry ($`\beta=1`$), complex Hermitian otherwise
> ($`\beta=2`$). Levels then *repel*: two levels of the same symmetry almost never come close
> (von Neumann and Wigner, 1929). Repulsion is the typical signature of quantum chaos, independent
> levels (no repulsion) that of integrable or localised systems; these are the Bohigas–Giannoni–Schmit
> and Berry–Tabor conjectures of physics, which are neither proved nor used in this course.

Spacings depend on the local density of levels, which varies along the spectrum; comparing them
with $`p_\beta`$ requires **unfolding**, the change of variables $`\lambda\mapsto\bar N(\lambda)`$ by a
smooth approximation $`\bar N`$ of the counting function, which makes the mean spacing $`1`$. The result
depends on the choice of $`\bar N`$. The ratio of consecutive spacings avoids this choice.

**Definition.** The **spacing ratio** is $`r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})\in[0,1]`$
(V. Oganesyan, D. A. Huse, Phys. Rev. B 75 (2007), 155111). It is unchanged by
$`\lambda\mapsto a\lambda+b`$, $`a>0`$, so it does not depend on the local mean spacing; this is why it is
used without unfolding.

*In the library:* the mean spacing ratio of each row of levels is `SpinChain.level_spacing_ratio` — [tutorial, Section 5.3](../tutorial/05_random_matrices.md#s-5-3).

A ratio needs three levels, so the smallest Gaussian matrices with a spacing ratio are $`3\times3`$; for them
the mean ratio can be computed exactly, as the surmise was for $`2\times2`$.

<a id="r-5-14"></a>

**Proposition 5.14 (mean spacing ratios).** (a) *Independent levels.* Let $`s,s'`$ be independent with the
same exponential distribution (as two consecutive spacings of the levels $`e_0+\dots+e_{k-1}`$ built from
independent exponential $`e_j`$, Section A.7). Then
$`r=\min(s,s')/\max(s,s')`$ has $`P(r\le t)=\frac{2t}{1+t}`$ on $`[0,1]`$, density $`\frac2{(1+t)^2}`$, and
```math
\mathbb Er=2\ln2-1\approx0.3863 .
```
(b) *$`3\times3`$ GOE and GUE.* Let $`H`$ be a $`3\times3`$ GOE ($`\beta=1`$) or GUE ($`\beta=2`$) matrix with eigenvalues
$`\lambda_1<\lambda_2<\lambda_3`$, spacings $`s_1=\lambda_2-\lambda_1`$, $`s_2=\lambda_3-\lambda_2`$, and $`r=\min(s_1,s_2)/\max(s_1,s_2)`$.
Then $`x=s_2/s_1`$ has a density proportional to
```math
f_\beta(x)=\frac{(x+x^2)^\beta}{(1+x+x^2)^{1+3\beta/2}}\qquad\text{on }(0,\infty),
```
$`r`$ has the density $`f_\beta/I_0`$ on $`(0,1)`$, where $`I_0=\int_0^1f_\beta`$, and
```math
\mathbb Er=4-2\sqrt3\approx0.5359\quad(\text{GOE}),\qquad
\mathbb Er=\frac{2\sqrt3}\pi-\frac12\approx0.6027\quad(\text{GUE}).
```

*Proof.* (a) By scaling, the rate is $`1`$. For $`x=s/s'`$, by independence and Tonelli's theorem,
$`P(x\le t)=\int_0^\infty e^{-y}(1-e^{-ty})\,dy=1-\frac1{1+t}=\frac t{1+t}`$.
For $`t<1`$ the events $`\{x\le t\}`$ and $`\{1/x\le t\}`$ are disjoint, $`x`$ and $`1/x`$ have the same law
($`s,s'`$ are exchangeable), and $`\{r\le t\}`$ is their union, so $`P(r\le t)=\frac{2t}{1+t}`$. Then
$`\mathbb Er=\int_0^1\frac{2t}{(1+t)^2}dt=2\Big[\ln(1+t)+\frac1{1+t}\Big]_0^1=2\ln2-1`$.

(b) *The law of $`x`$.* By Theorem C.7 the vector $`\lambda=(\lambda_1,\lambda_2,\lambda_3)`$ has the density
$`Z^{-1}\Delta(\lambda)^\beta e^{-|\lambda|^2/2}`$ on $`C=\{\lambda_1<\lambda_2<\lambda_3\}`$, $`\Delta(\lambda)=\prod_{i<j}(\lambda_j-\lambda_i)`$. The linear map
```math
(m,s_1,s_2)\mapsto\Big(m-\tfrac{2s_1+s_2}3,\ m+\tfrac{s_1-s_2}3,\ m+\tfrac{s_1+2s_2}3\Big)
```
has determinant $`1`$ and maps $`\mathbb R\times(0,\infty)^2`$ onto $`C`$; its inverse is $`m=\frac13\sum_i\lambda_i`$,
$`s_1=\lambda_2-\lambda_1`$, $`s_2=\lambda_3-\lambda_2`$. In these variables $`\Delta=s_1s_2(s_1+s_2)`$ and
$`|\lambda|^2=3m^2+\frac23(s_1^2+s_1s_2+s_2^2)`$ (expand). By Lemma A.3(c) and Fubini's theorem (integrate out $`m`$),
$`(s_1,s_2)`$ has the density $`c\,\big(s_1s_2(s_1+s_2)\big)^\beta e^{-(s_1^2+s_1s_2+s_2^2)/3}`$ on $`(0,\infty)^2`$ for a constant
$`c>0`$. Let $`q(x)=1+x+x^2`$. For Borel $`A\subseteq(0,\infty)`$, substituting $`s_2=xs_1`$ in the inner integral and
using Tonelli's theorem,
```math
P(x\in A)=c\int_0^\infty\!\!\int_A s_1^{3\beta+1}(x+x^2)^\beta e^{-s_1^2q(x)/3}\,dx\,ds_1
=c\,\gamma_\beta\int_A(x+x^2)^\beta\Big(\frac3{q(x)}\Big)^{1+3\beta/2}dx ,
```
where the substitution $`s_1=(3/q(x))^{1/2}u`$ gives the factor, with $`\gamma_\beta=\int_0^\infty u^{3\beta+1}e^{-u^2}du\in(0,\infty)`$.
So $`x`$ has the density $`f_\beta/\int_0^\infty f_\beta`$.

*Reduction to $`[0,1]`$.* $`f_\beta(1/x)=x^2f_\beta(x)`$ (multiply numerator and denominator by powers of $`x`$). With
$`x=1/y`$, $`\int_{1/A}f_\beta(x)\,dx=\int_Af_\beta(1/y)y^{-2}dy=\int_Af_\beta`$ for Borel $`A\subseteq(0,1)`$, where
$`1/A=\{1/y:y\in A\}`$; in particular $`\int_0^\infty f_\beta=2I_0`$ with $`I_j=\int_0^1x^jf_\beta(x)\,dx`$. As $`r=\min(x,1/x)`$ and
$`P(x=1)=0`$, $`\{r\in A\}=\{x\in A\}\cup\{x\in1/A\}`$ (disjoint), so
```math
P(r\in A)=\frac{2\int_Af_\beta}{2I_0}=\int_A\frac{f_\beta}{I_0},\qquad \mathbb Er=\frac{I_1}{I_0} .
```
Since $`x+x^2=q-1`$, $`f_\beta=g_\beta(q)`$ with $`g_1(u)=u^{-3/2}-u^{-5/2}`$ and $`g_2(u)=(u-1)^2u^{-4}`$. Write
$`x=\frac12(q'(x)-1)`$ with $`q'=2x+1`$, and substitute $`u=q(x)`$ (from $`q(0)=1`$ to $`q(1)=3`$):
```math
I_1=\frac12\int_1^3g_\beta(u)\,du-\frac12I_0,\qquad \mathbb Er=\frac{\int_1^3g_\beta}{2I_0}-\frac12 .
```
Here $`\int_1^3g_1=\big[-2u^{-1/2}+\frac23u^{-3/2}\big]_1^3=\frac43-\frac{16}{9\sqrt3}`$ and
$`\int_1^3g_2=\big[-u^{-1}+u^{-2}-\frac13u^{-3}\big]_1^3=\frac8{81}`$.

*The integrals $`I_0`$.* Let $`L_a=\int_0^1q^{-a}dx`$. From $`(2x+1)^2=4q-3`$,
$`\frac d{dx}\big[(2x+1)q^{-a}\big]=2q^{-a}-a(2x+1)^2q^{-a-1}=(2-4a)q^{-a}+3aq^{-a-1}`$, and integrating over $`[0,1]`$,
```math
3a\,L_{a+1}=3^{1-a}-1+(4a-2)L_a .
```
GOE: $`f_1=q^{-3/2}-q^{-5/2}`$. With $`a=\frac12`$, $`L_{3/2}=\frac23(\sqrt3-1)`$; with $`a=\frac32`$,
$`L_{5/2}=\frac29\big(\frac1{\sqrt3}-1+4L_{3/2}\big)=\frac{18\sqrt3-22}{27}`$. So $`I_0=L_{3/2}-L_{5/2}=\frac4{27}`$ and
$`\mathbb Er=\frac{27}8\big(\frac43-\frac{16}{9\sqrt3}\big)-\frac12=\frac92-2\sqrt3-\frac12=4-2\sqrt3`$.
GUE: $`f_2=q^{-2}-2q^{-3}+q^{-4}`$. Since $`q=(x+\frac12)^2+\frac34`$,
$`L_1=\frac2{\sqrt3}\big[\arctan\frac{2x+1}{\sqrt3}\big]_0^1=\frac2{\sqrt3}\big(\frac\pi3-\frac\pi6\big)=\frac\pi{3\sqrt3}`$, and the
recursion with $`a=1,2,3`$ gives $`L_2=\frac23L_1=\frac{2\pi}{9\sqrt3}`$, $`L_3=L_2-\frac19`$, $`L_4=\frac19\big(10L_3-\frac89\big)`$. So
$`I_0=L_2-2L_3+L_4=\frac{2\pi}{81\sqrt3}`$ and $`\mathbb Er=\frac{81\sqrt3}{4\pi}\cdot\frac8{81}-\frac12=\frac{2\sqrt3}\pi-\frac12`$. $`\square`$

The values in (b) were found by Y. Y. Atas, E. Bogomolny, O. Giraud and G. Roux (Phys. Rev. Lett. 110
(2013), 084101). For larger $`n`$ the mean ratio is not computed in this course.

Part (a) holds exactly for $`n`$ independent uniform levels on an interval: their $`n+1`$ spacings (including the
two end gaps) have the law of $`(e_0,\dots,e_n)/\sum_je_j`$ for independent exponential $`e_j`$, so each ratio
$`s_i/s_{i+1}`$ of consecutive spacings has the law of $`e_i/e_{i+1}`$, and each spacing ratio $`r_i`$ has the law of $`r`$ in (a)
(Theorem A.25 and Corollary A.26(b),(d)).

*Computation:* Proposition 5.14 by Monte Carlo and by quadrature of $`f_\beta`$, and the mean ratio of the GOE with $`n=64`$, observed numerically — [tutorial, Section 5.3](../tutorial/05_random_matrices.md#s-5-3).

<a id="s-5-4"></a>

## 5.4 Why sectors must not be mixed

Let $`H=\bigoplus_cH_c\in\mathcal M`$ be a random self-adjoint element with **independent blocks**,
for example GOE blocks. As an operator on $`\mathbb C^D`$ its spectrum is the union of the spectra
of the blocks. Inside a block the levels repel; we show that levels of different blocks do not, and that with many
independent blocks the merged spectrum has the spacing ratios of independent levels.

Assume that each block has a bounded continuous **level density** $`\rho_c`$:
$`\mathbb E\sum_i\varphi(\lambda_i(H_c))=\int\varphi\,\rho_c`$ for bounded Borel $`\varphi`$ (so
$`\int\rho_c=k_c`$). For GOE and GUE blocks this holds by Corollary C.8, which derives $`\rho_c`$ from the joint
density of Theorem C.7.

<a id="r-5-15"></a>

**Theorem 5.15 (no repulsion between sectors).** For $`c\neq d`$ and every bounded Borel
$`\psi:\mathbb R^2\to\mathbb R`$,
```math
\mathbb E\sum_{i,j}\psi\big(\lambda_i(H_c),\lambda_j(H_d)\big)=\iint\psi(x,y)\,\rho_c(x)\rho_d(y)\,dx\,dy .
```
Consequently, if $`N_\varepsilon`$ is the number of pairs of levels from different blocks at distance
$`<\varepsilon`$, then $`\mathbb EN_\varepsilon/(2\varepsilon)\to\sum_{c<d}\int\rho_c\rho_d`$ as
$`\varepsilon\to0`$, and the smallest spacing $`\delta`$ of the merged spectrum satisfies
```math
P(\delta<\varepsilon)\ \ge\ \frac{\mathbb EN_\varepsilon}{\sum_{c<d}k_ck_d}
=\frac{2\varepsilon\sum_{c<d}\int\rho_c\rho_d}{\sum_{c<d}k_ck_d}+o(\varepsilon).
```

*Proof.* The ordered eigenvalues are continuous functions of the matrix (Corollary A.11), so
$`\Lambda_c=(\lambda_i(H_c))_i`$ and $`\Lambda_d`$ are independent random vectors. By Fubini's theorem (their joint law is
the product of their laws), $`\mathbb E\sum_{i,j}\psi(\Lambda_{c,i},\Lambda_{d,j})=\mathbb E\sum_j g(\Lambda_{d,j})`$ with
$`g(y)=\mathbb E\sum_i\psi(\Lambda_{c,i},y)=\int\psi(x,y)\rho_c(x)\,dx`$, which is bounded and Borel; apply the
definition of $`\rho_d`$ to $`g`$. With $`\psi(x,y)=1\{|x-y|<\varepsilon\}`$,
```math
\frac{\mathbb EN_\varepsilon}{2\varepsilon}=\sum_{c<d}\int\rho_c(x)\Big(\frac1{2\varepsilon}\int_{x-\varepsilon}^{x+\varepsilon}\rho_d(y)\,dy\Big)dx ;
```
the bracket tends to $`\rho_d(x)`$ by continuity and is bounded by $`\sup\rho_d`$, and $`\rho_c`$ is integrable,
so dominated convergence gives the limit. Finally $`N_\varepsilon\le\sum_{c<d}k_ck_d`$, so
$`\mathbb EN_\varepsilon\le\big(\sum_{c<d}k_ck_d\big)P(N_\varepsilon\ge1)`$, and $`N_\varepsilon\ge1`$ implies $`\delta<\varepsilon`$. $`\square`$

So for two GOE or GUE blocks $`P(\delta<\varepsilon)`$ is at least of order $`\varepsilon`$: their level densities are
positive everywhere (in Corollary C.8, $`q(x,\cdot)>0`$ outside the null set where two coordinates of $`(x,y)`$
coincide), so $`\int\rho_c\rho_d>0`$. Inside a $`2\times2`$ GOE block it is $`1-e^{-\varepsilon^2/4}\le\varepsilon^2/4`$
(Proposition 5.13, unnormalised spacing), and inside an $`n\times n`$ GOE block it is at most a constant times
$`\varepsilon^2`$ (Proposition C.9). The level-repulsion test fails for the merged spectrum although
every block is a perfect GOE.

With many sectors the merged spectrum loses every trace of the blocks: its spacing ratios become those of
independent levels. We prove this, with a rate, for $`C`$ independent copies $`H_1,\dots,H_C`$ of one random
self-adjoint $`k\times k`$ matrix $`H`$, with $`k`$ fixed and $`C\to\infty`$, for the levels in a window where the level
density is positive. Two hypotheses on $`H`$ are used:

- **(H1)** its level density $`\rho`$ is bounded and Lipschitz: $`\rho\le\rho_{\max}`$ and $`|\rho(x)-\rho(y)|\le L|x-y|`$;
- **(H2)** its smallest spacing $`\delta=\min_i(\lambda_{i+1}-\lambda_i)`$ ($`\delta=+\infty`$ if $`k=1`$) satisfies
  $`P(\delta\le\varepsilon)\le K\varepsilon^2`$ for all $`\varepsilon\ge0`$.

(H2) is a weak form of level repulsion inside the block; with $`\varepsilon=0`$ it says that the spectrum of $`H`$ is
simple almost surely.

<a id="r-5-16"></a>

**Lemma 5.16 (GOE and GUE blocks).** GOE and GUE matrices of every size $`n\ge1`$ satisfy (H1) and (H2), and their
level density is positive everywhere.

*Proof.* Positivity was shown after Theorem 5.15. *(H2).* For $`n=1`$ there is nothing to prove. For $`n\ge2`$ and
$`0\le\varepsilon<1`$, Proposition C.9 gives a constant $`K_0`$ with $`P(\delta\le\varepsilon)\le P(\delta<\varepsilon')\le K_0\varepsilon'^{\beta+1}\le K_0\varepsilon'^2`$
for every $`\varepsilon'\in(\varepsilon,1]`$; let $`\varepsilon'\downarrow\varepsilon`$. With $`K=\max(K_0,1)`$ the bound also holds for $`\varepsilon\ge1`$.
*(H1).* $`\rho`$ is bounded by Corollary C.8. In its notation, with $`x\in\mathbb R`$ and $`y\in\mathbb R^{n-1}`$,
$`|\Delta(x,y)|=|\Delta(y)|\prod_i|x-y_i|`$, so
```math
\rho(x)=\frac{n}{n!\,Z_{n,\beta}}\int_{\mathbb R^{n-1}}|\Delta(y)|^\beta e^{-|y|^2/2}\,g_y(x)\,dy,\qquad
g_y(x)=f_0(x)\prod_{i=1}^{n-1}f_i(x),
```
with $`f_0(x)=e^{-x^2/(2n)}`$ and $`f_i(x)=|x-y_i|^\beta e^{-x^2/(2n)}`$ (for $`n=1`$ there is no $`y`$, and $`\rho=f_0/Z_{1,\beta}`$). Let
$`a=(\beta+1)\sup_x(1+|x|)^{\beta+1}e^{-x^2/(2n)}<\infty`$. From $`|x-y_i|\le(1+|x|)(1+|y_i|)`$ we get $`|f_i|\le a(1+|y_i|)^\beta`$, and $`f_i`$ is
continuous, continuously differentiable off $`y_i`$, with $`|f_i'(x)|\le\big(\beta|x-y_i|^{\beta-1}+|x||x-y_i|^\beta\big)e^{-x^2/(2n)}\le a(1+|y_i|)^\beta`$.
By the mean value theorem on each side of $`y_i`$, $`f_i`$ is Lipschitz with constant $`a(1+|y_i|)^\beta`$; likewise $`f_0`$ with
constant $`a`$, and $`|f_0|\le1\le a`$. Writing the difference of two products as a telescoping sum,
$`|g_y(x)-g_y(x')|\le\sum_i|f_i(x)-f_i(x')|\prod_{j\neq i}\sup|f_j|\le na^n\prod_i(1+|y_i|)^\beta\,|x-x'|`$. Hence $`\rho`$ is
Lipschitz with $`L=\frac{n^2a^n}{n!\,Z_{n,\beta}}\int|\Delta(y)|^\beta\prod_i(1+|y_i|)^\beta e^{-|y|^2/2}dy`$, which is finite by the bound
$`|\Delta(y)|\le\prod_i(1+|y_i|)^{n-2}`$ of the proof of Theorem C.7. $`\square`$

Let $`\Lambda_c`$ be the set of levels of $`H_c`$; the **merged spectrum of $`m`$ blocks** is $`\Lambda_1\cup\dots\cup\Lambda_m`$. For a point
$`x`$ let $`L_m(x)`$ and $`R_m(x)`$ be the distances from $`x`$ to the nearest level of the merged spectrum of $`m`$ blocks below and above
$`x`$ ($`+\infty`$ if there is none), and, when both are finite, $`r_m(x)=\min(L_m(x),R_m(x))/\max(L_m(x),R_m(x))`$. For an interval
$`J`$ of length $`|J|`$ let $`N(J)`$ be the number of levels of $`H`$ in $`J`$, $`\rho(J)=\int_J\rho=\mathbb EN(J)`$ and $`\sigma(J)=P(N(J)=0)`$.

<a id="r-5-17"></a>

**Lemma 5.17 (the two gaps around a fixed point).** Assume (H1) and (H2), and let $`\ell,\rho_*>0`$ and $`x\in\mathbb R`$ with
$`\rho\ge\rho_*`$ on $`[x-\ell,x+\ell]`$. There is a constant $`A`$, depending only on $`k,K,L,\rho_{\max},\rho_*,\ell`$, such that for
all $`m\ge1`$ and $`0\le t\le1`$
```math
\Big|P\big(L_m(x),R_m(x)<\infty,\ r_m(x)\le t\big)-\frac{2t}{1+t}\Big|\le\frac Am .
```

*Proof.* For $`m\le2`$ any $`A\ge2`$ will do; let $`m\ge3`$. *One block.* Since $`1\{N\ge1\}\le N\le k\,1\{N\ge1\}`$, taking
expectations gives (1): $`1-\rho(J)\le\sigma(J)\le1-\rho(J)/k\le e^{-\rho(J)/k}`$. If $`N(J)\ge2`$, two levels lie in $`J`$ and
$`\delta\le|J|`$; as $`N-1\{N\ge1\}\le k\,1\{N\ge2\}`$, (H2) gives (2): $`\sigma(J)\le1-\rho(J)+kK|J|^2`$. For $`c>0`$, by Tonelli's
theorem and (H2),
```math
\mathbb Ee^{-c\delta}=\int_0^1P\big(e^{-c\delta}>u\big)\,du\le\int_0^1\frac{K\ln^2u}{c^2}\,du=\frac{2K}{c^2}.\qquad(3)
```
*Simplicity.* By (H2) with $`\varepsilon=0`$ each block has a simple spectrum almost surely; by Theorem 5.15 with
$`\psi(u,v)=1\{u=v\}`$ (the diagonal is a null set) two blocks have no common level almost surely; and $`x`$ is almost surely not
a level, since $`\mathbb EN(\{x\})=0`$.

*A formula.* Write $`\Lambda`$ for the merged spectrum of $`m`$ blocks, $`L,R`$ for $`L_m(x),R_m(x)`$, fix $`t`$, and let
$`J_b=[x-tb,x+b)`$ and $`\sigma_b=\sigma(J_b)`$ for $`b>0`$; $`J_b`$ grows with $`b`$, so $`\sigma_b`$ is nonincreasing. Almost surely
```math
1\{R<\infty,\ L>tR\}=\sum_{y\in\Lambda,\ y>x}1\{\Lambda\cap J_{y-x}=\{y\}\},
```
since for $`y=x+R`$ the other levels in $`J_R`$ are those in $`[x-tR,x)`$, which exist iff $`L\le tR`$, and for $`y>x+R`$ the level
$`x+R`$ lies in $`J_{y-x}`$. Split the sum according to the block of $`y`$. The blocks are independent and have the same law, and
given the levels of the block of $`y`$, the other $`m-1`$ blocks miss $`J_{y-x}`$ with probability $`\sigma_{y-x}^{m-1}`$; by Fubini's
theorem as in the proof of Theorem 5.15,
```math
P(R<\infty,\ L>tR)=m\,\mathbb E\sum_{y\in\Lambda_1,\ y>x}1\{\Lambda_1\cap J_{y-x}=\{y\}\}\,\sigma_{y-x}^{m-1}=M-E_1,
```
where $`M=m\,\mathbb E\sum_{y\in\Lambda_1,y>x}\sigma_{y-x}^{m-1}=m\int_0^\infty\rho(x+b)\sigma_b^{m-1}db`$ (definition of $`\rho`$) and $`E_1\ge0`$
is the same expectation with $`1\{\Lambda_1\cap J_{y-x}\neq\{y\}\}`$.

*$`E_1`$ is small.* On that event $`H_1`$ has a second level within distance $`|J_{y-x}|\le2(y-x)`$ of $`y`$, so $`y-x\ge\delta_1/2`$,
where $`\delta_1`$ is the smallest spacing of $`H_1`$. By (1), $`\sigma_b\le e^{-\rho_*\min(b,\ell)/k}`$, as $`J_b\supseteq[x,x+\min(b,\ell))`$. With at
most $`k`$ terms, $`\mu=(m-1)\rho_*/k`$ and (3),
$`E_1\le mk\,\mathbb Ee^{-\mu\min(\delta_1/2,\ell)}\le mk\big(8Kk^2/((m-1)\rho_*)^2+e^{-(m-1)\rho_*\ell/k}\big)`$.

*$`M`$ is close to $`\frac1{1+t}`$.* Let $`\rho_0=\rho(x)`$, $`a=(1+t)\rho_0`$ and $`b_1=\min(\ell,1/a)\ge b_0=\min(\ell,1/(2\rho_{\max}))`$.
Then $`m\int_0^{1/a}\rho_0(1-ab)^{m-1}db=\rho_0/a=\frac1{1+t}`$. For $`0<b\le b_1`$: $`|\rho(J_b)-ab|\le\int_{J_b}|\rho-\rho_0|\le L\frac{(1+t^2)b^2}2\le Lb^2`$,
so by (1) and (2) $`|\sigma_b-(1-ab)|\le K_1b^2`$ with $`K_1=L+4kK`$; both numbers lie in $`[0,e^{-\rho_*b/k}]`$; and
$`|\rho(x+b)-\rho_0|\le Lb`$. By $`|u^n-v^n|\le n|u-v|\max(u,v)^{n-1}`$ for $`u,v\ge0`$,
```math
\big|\rho(x+b)\sigma_b^{m-1}-\rho_0(1-ab)^{m-1}\big|\le\rho_{\max}(m-1)K_1b^2e^{-(m-2)\rho_*b/k}+Lb\,e^{-(m-1)\rho_*b/k}.
```
With $`\int_0^\infty b^je^{-\lambda b}db=j!/\lambda^{j+1}`$, $`m`$ times the integral of this over $`(0,b_1]`$ is at most
$`2m(m-1)\rho_{\max}K_1k^3/((m-2)\rho_*)^3+mLk^2/((m-1)\rho_*)^2`$. The remaining parts are
$`m\int_{b_1}^\infty\rho(x+b)\sigma_b^{m-1}db\le mk\,e^{-(m-1)\rho_*b_0/k}`$ (as $`\int\rho=k`$ and $`\sigma_b\le\sigma_{b_1}`$) and
$`m\int_{b_1}^{1/a}\rho_0(1-ab)^{m-1}db\le(1-ab_1)^m\le e^{-m\rho_*b_0}`$.

*Both sides.* The same argument for the matrices $`-H_c`$ (level density $`\rho(-\,\cdot)`$, same hypotheses) at the point $`-x`$
gives $`P(L<\infty,\ R>tL)=M'-E_1'`$ with the same bounds. For $`t<1`$ the events $`\{R<\infty,\ L\le tR\}`$ and
$`\{L<\infty,\ R\le tL\}`$ are disjoint (together they give $`L\le t^2L`$, while $`L>0`$), and their union is
$`\{L,R<\infty,\ r_m(x)\le t\}`$, so
```math
P\big(L,R<\infty,\ r_m(x)\le t\big)-\frac{2t}{1+t}=-P(R=\infty)-P(L=\infty)-\Big(M-\frac1{1+t}\Big)-\Big(M'-\frac1{1+t}\Big)+E_1+E_1'.
```
By (1), $`P(R=\infty)\le\sigma((x,x+\ell))^m\le e^{-m\rho_*\ell/k}`$, and likewise for $`L`$; for $`t=1`$ the probability is
$`1-P(L=\infty\text{ or }R=\infty)`$. Every term is at most a constant over $`m`$, because $`m(m-1)/(m-2)^3\le18/m`$,
$`m/(m-1)^2\le4/m`$ and $`me^{-(m-1)c}\le4/(mc^2)`$ for $`c>0`$ (from $`u^2e^{-u}\le1`$). $`\square`$

<a id="r-5-18"></a>

**Theorem 5.18 (many independent sectors look like independent levels).** Assume (H1) and (H2). Let $`I=[\alpha,\gamma]`$
and $`\ell,\rho_*>0`$ with $`\rho\ge\rho_*`$ on $`[\alpha-\ell,\gamma+\ell]`$, let $`C\ge2`$, and let $`\Lambda`$ be the merged spectrum of
$`C`$ blocks (simple almost surely). For a level $`y`$ with neighbours $`y_-<y<y_+`$ in $`\Lambda`$ let
$`r(y)=\min(y-y_-,y_+-y)/\max(y-y_-,y_+-y)`$. Let $`N(I)`$ be the number of levels in $`I`$, $`N_t(I)`$ the number of those
with two neighbours and $`r(y)\le t`$, and $`S(I)`$ the sum of $`r(y)`$ over the levels in $`I`$ with two neighbours. Then
$`\mathbb EN(I)=C\rho(I)`$, and there is a constant $`A`$, depending only on $`k,K,L,\rho_{\max},\rho_*,\ell`$ and $`\rho(I)`$, such
that for $`0\le t\le1`$
```math
\Big|\frac{\mathbb EN_t(I)}{\mathbb EN(I)}-\frac{2t}{1+t}\Big|\le\frac AC,\qquad
\Big|\frac{\mathbb ES(I)}{\mathbb EN(I)}-(2\ln2-1)\Big|\le\frac AC .
```

The ratios of expectations are what a computation pooled over many samples estimates: the fraction of the levels in
$`I`$ with $`r\le t`$, and their mean spacing ratio. By Lemma 5.16 the theorem applies to independent GOE or GUE blocks of
any fixed size and to every compact interval $`I`$.

*Proof.* $`\mathbb EN(I)=\sum_c\rho(I)`$. The blocks have the same law, so
$`\mathbb EN_t(I)=C\,\mathbb E\sum_{y\in\Lambda_1\cap I}1\{y\text{ has two neighbours},\ r(y)\le t\}`$. Let $`m=C-1`$, and for
$`y\in\Lambda_1`$ let $`L'(y),R'(y)`$ be its gaps in the merged spectrum $`\Lambda'`$ of $`H_2,\dots,H_C`$ (as in Lemma 5.17) and $`r'(y)`$
their ratio. The other levels of $`H_1`$ are at distance $`\ge\delta_1`$ from $`y`$; so if $`\max(L'(y),R'(y))<\delta_1`$, the neighbours of
$`y`$ in $`\Lambda`$ are those in $`\Lambda'`$, and $`r(y)=r'(y)`$. Now $`\Lambda'`$ is independent of $`H_1`$, and $`\max(L'(y),R'(y))\ge s`$
requires $`\Lambda'`$ to miss $`(y-s,y)`$ or $`(y,y+s)`$, which by (1) has probability at most $`2e^{-m\rho_*\min(s,\ell)/k}`$ for
$`y\in I`$. By Fubini's theorem (as in the proof of Theorem 5.15) and (3),
```math
\Big|\mathbb EN_t(I)-C\,\mathbb E\sum_{y\in\Lambda_1\cap I}p_m(y,t)\Big|\le2Ck\big(\mathbb Ee^{-m\rho_*\delta_1/k}+e^{-m\rho_*\ell/k}\big)
\le2Ck\Big(\frac{2Kk^2}{(m\rho_*)^2}+e^{-m\rho_*\ell/k}\Big),
```
where $`p_m(x,t)=P\big(L_m(x),R_m(x)<\infty,\ r_m(x)\le t\big)`$, a Borel function of $`x`$ by Fubini's theorem. By the definition of
$`\rho`$, $`\mathbb E\sum_{y\in\Lambda_1\cap I}p_m(y,t)=\int_I\rho(x)p_m(x,t)\,dx`$, and Lemma 5.17 applies at every $`x\in I`$ with one
constant $`A_0`$. Dividing by $`C\rho(I)`$,
$`\big|\mathbb EN_t(I)/\mathbb EN(I)-\frac{2t}{1+t}\big|\le A_0/m+2k\big(2Kk^2/(m\rho_*)^2+e^{-m\rho_*\ell/k}\big)/\rho(I)`$, which is at
most a constant over $`C`$ since $`m\ge C/2`$. For the mean, $`r(y)=\int_0^11\{r(y)>t\}\,dt`$, so by Tonelli's theorem
$`\mathbb ES(I)=\int_0^1\big(\mathbb EN_1(I)-\mathbb EN_t(I)\big)dt`$, while $`\int_0^1\big(1-\frac{2t}{1+t}\big)dt=\int_0^1\frac{1-t}{1+t}dt=2\ln2-1`$
(Proposition 5.14(a)); so the second deviation is at most twice the bound of the first. $`\square`$

The proof shows the mechanism. A level $`y`$ of one block sees the other $`C-1`$ blocks as an independent configuration,
whose levels next to $`y`$ lie at distances of order $`1/(C\rho)`$; at such distances the other levels of the block of $`y`$ are
absent with high probability (this is where $`C\to\infty`$ and (H2) enter), so $`r(y)`$ is the ratio in which $`y`$ divides the
gap of the other blocks that contains it. For a constant density and at most one level per block near $`y`$ the formula of
the proof gives exactly the value of independent levels, $`m\int_0^{1/a}\rho_0(1-ab)^{m-1}db=\frac1{1+t}`$, for every $`m`$; the
errors come from the variation of $`\rho`$ over a gap and from two levels of one block in a gap. The rate $`1/C`$ is what the proof
gives; for GOE blocks the deviation found numerically decreases like $`1/C^2`$ (a numerical observation, in the tutorial).
That superpositions of many independent sparse point processes are close to Poisson processes is a classical theme of
probability (B. Grigelionis, Theory Probab. Appl. 8 (1963), 177–182); Theorem 5.18 is the case needed here, with a proof
within this course. Independence of the blocks is essential to the proof: the theorem does not apply to the sectors of the
Heisenberg chain with random fields (Chapter 9), which all depend on the same fields.

*Computation:* Theorem 5.18 for $`C`$ independent GOE blocks of size $`64`$, with a check of the rate $`1/C`$; and the Heisenberg chain with random fields, whose blocks are not independent, where the same effect is a numerical observation — [tutorial, Section 5.4](../tutorial/05_random_matrices.md#s-5-4).

> **Physical meaning** (physical interpretation, not a theorem of this course). Levels with different
> conserved quantum numbers do not interact, so they cross freely; only levels of the same symmetry repel. A spectrum that mixes symmetry sectors
> therefore shows Poisson-like statistics even when the dynamics in each sector is fully chaotic,
> and an unresolved symmetry can be mistaken for integrability or localisation. The rule of the
> trade — resolve every symmetry before computing level statistics — is the statement that one
> should compute in the algebra $`\bigoplus_cM_{k_c}`$, not in $`M_D`$.

<a id="s-5-5"></a>

## 5.5 What to compare with

*In the library:* comparison with other software — [tutorial, Section 5.5](../tutorial/05_random_matrices.md#s-5-5).

<a id="x-5"></a>

## Exercises

1. Let $`U`$ be Haar on $`U(\mathcal M)`$ and $`\mathrm{Tr}\,U=\sum_c\mathrm{Tr}\,U_c`$. Show that
   $`\mathbb E|\mathrm{Tr}\,U|^2=C`$ and $`\mathbb E|\mathrm{Tr}\,U|^4=\sum_cm(k_c)+2C(C-1)`$ with
   $`m(1)=1`$, $`m(k)=2`$ for $`k\ge2`$; compare with $`1`$ and $`2`$ for a Haar unitary of $`U(D)`$, $`D\ge2`$, and
   evaluate for $`\mathcal M=M_1\oplus M_2\oplus M_3`$.
2. Using Theorem 5.12, show that $`\mathbb E|U_{11}|^{2n}=\frac{n!}{k(k+1)\cdots(k+n-1)}`$ for $`n\le k`$
   (Corollary B.14 extends this to $`n>k`$).
3. The **circular orthogonal ensemble** is $`S=W^TW`$ with $`W`$ Haar on $`U(k)`$. Show
   that $`S`$ is unitary and symmetric, and that its law is invariant under $`S\mapsto V^TSV`$ for every
   $`V\in U(k)`$. Let $`\lambda_1,\dots,\lambda_k`$ be the eigenvalues of $`S`$ (with multiplicity) and
   $`N_B(S)=\#\{j:\lambda_j\in B\}`$ for Borel sets $`B`$ of the unit circle $`\mathbb T`$. Show that $`N_B(S)`$ is a
   random variable, that $`\nu(B)=\frac1k\mathbb EN_B(S)`$ is the normalised arc length (the mean density of the
   eigenphases is uniform), and that $`\mathbb E\mathrm{Tr}\,S^m=0`$ for every integer $`m\neq0`$.

<details>
<summary>Solutions</summary>

1. Write $`t_c=\mathrm{Tr}\,U_c`$; the $`t_c`$ are independent (Corollary 5.4). Multiplying $`U_c`$ by
   $`e^{i\theta}`$ does not change its law and multiplies $`t_c`$ by $`e^{i\theta}`$; so a moment
   $`\mathbb E[t_at_b\bar t_c\bar t_d]`$ vanishes unless every sector occurs as often in $`(a,b)`$ as in
   $`(c,d)`$, i.e. $`\{a,b\}=\{c,d\}`$ as multisets, and likewise $`\mathbb E[t_a\bar t_b]=0`$ for $`a\neq b`$. By
   Proposition 5.7, $`\mathbb E[(U_c)_{ii}\overline{(U_c)_{jj}}]=\langle e_i,\mathbb E[U_ce_ie_j^*U_c^*]e_j\rangle=\delta_{ij}/k_c`$,
   so $`\mathbb E|t_c|^2=1`$ and $`\mathbb E|\mathrm{Tr}\,U|^2=C`$. In $`\mathbb E|\sum_at_a|^4`$ the terms with
   $`a=b`$ give $`\sum_a\mathbb E|t_a|^4`$, those with $`a\neq b`$ give
   $`2\sum_{a\neq b}\mathbb E|t_a|^2\,\mathbb E|t_b|^2=2C(C-1)`$ (two choices of $`(c,d)`$). For $`k=1`$, $`|t|=1`$.
   For $`k\ge2`$, $`\mathbb E|t|^4=\sum_{i,j,i',j'}\mathbb E[U_{ii}U_{jj}\overline{U_{i'i'}U_{j'j'}}]`$; in
   Theorem 5.10 the two terms with $`\sigma=\pi`$ force $`\{i',j'\}=\{i,j\}`$ in one way and give
   $`k^2\mathrm{Wg}(e)`$ each after summation, the two with $`\sigma\neq\pi`$ force $`i=j=i'=j'`$ and give
   $`k\,\mathrm{Wg}(\tau)`$ each: in total $`2\big(\frac{k^2}{k^2-1}-\frac1{k^2-1}\big)=2`$. The case $`C=1`$,
   $`k=D`$ gives $`1`$ and $`2`$ for $`U(D)`$. For $`M_1\oplus M_2\oplus M_3`$:
   $`\mathbb E|\mathrm{Tr}\,U|^4=1+2+2+12=17`$, against $`2`$ for a Haar unitary of $`U(6)`$.

2. With all indices equal to $`1`$, every delta is $`1`$, so
   $`\mathbb E|U_{11}|^{2n}=\sum_{\sigma,\pi}\mathrm{Wg}(\sigma\pi^{-1})=\mathbf 1^TW\mathbf 1`$, where
   $`W=\Gamma^{-1}=(\mathrm{Wg}(\sigma^{-1}\pi))`$ and $`\Gamma=(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)})`$. Here we used that
   $`\mathrm{Wg}`$ depends only on the cycle type, and that $`\sigma\pi^{-1}=\sigma(\pi^{-1}\sigma)\sigma^{-1}`$ is
   conjugate to $`\pi^{-1}\sigma=(\sigma^{-1}\pi)^{-1}`$, which has the cycle type of $`\sigma^{-1}\pi`$. Every
   row sum of $`\Gamma`$ is $`c_n(k)=\sum_{\rho\in S_n}k^{\#\mathrm{cycles}(\rho)}`$, so $`\Gamma\mathbf 1=c_n(k)\mathbf 1`$,
   $`W\mathbf 1=\mathbf 1/c_n(k)`$ and $`\mathbf 1^TW\mathbf 1=n!/c_n(k)`$. Finally
   $`c_n(x)=x(x+1)\cdots(x+n-1)`$ by induction: a permutation of $`\{1,\dots,n\}`$ arises from one of
   $`\{1,\dots,n-1\}`$ either by adding $`n`$ as a fixed point (one more cycle, factor $`x`$) or by inserting $`n`$
   after one of the $`n-1`$ elements in its cycle ($`n-1`$ ways, same number of cycles); this is a bijection, so
   $`c_n(x)=(x+n-1)c_{n-1}(x)`$.

3. $`S^*S=W^*\bar WW^TW=1`$ since $`W^T`$ is unitary, and $`S^T=W^TW=S`$. For $`V\in U(k)`$, $`WV`$ has the law
   of $`W`$ (right invariance, Lemma 5.1), and $`(WV)^T(WV)=V^TSV`$. With $`V=e^{i\alpha/2}1`$: $`e^{i\alpha}S`$
   has the law of $`S`$ for every real $`\alpha`$.

   *Measurability.* By Theorem 1.3, $`S=Y\mathrm{diag}(\lambda_1,\dots,\lambda_k)Y^*`$ with $`Y`$ unitary. If $`v`$ is a
   unit vector with $`Sv=\mu v`$, then writing $`v=\sum_jc_jYe_j`$ gives $`c_j(\lambda_j-\mu)=0`$, so $`v`$ lies in
   $`E_K=\mathrm{span}\{Ye_j:\lambda_j\in K\}`$ whenever $`\mu\in K`$; and $`\dim E_K=N_K(S)`$. Hence, for $`K\subseteq\mathbb T`$
   and $`m\ge1`$: $`N_K(S)\ge m`$ iff $`S`$ has $`m`$ orthonormal eigenvectors with eigenvalues in $`K`$. If $`K`$ is
   closed, the set $`F_{K,m}`$ of unitaries with this property is closed in $`U(k)`$: if $`S_l\to S`$ and $`S_l`$ has
   orthonormal eigenvectors $`v_i^{(l)}`$ with eigenvalues $`\mu_i^{(l)}\in K`$ ($`i\le m`$), a subsequence has
   $`v_i^{(l)}\to v_i`$ and $`\mu_i^{(l)}\to\mu_i\in K`$ for all $`i\le m`$ (Bolzano–Weierstrass in the compact set
   $`(\text{unit sphere}\times K)^m`$), and the limits are orthonormal with $`Sv_i=\mu_iv_i`$.
   Let $`0\le\theta_1(S)\le\dots\le\theta_k(S)<2\pi`$ be the eigenphases, so that the eigenvalues are the
   $`e^{i\theta_j(S)}`$ (Corollary 1.4: $`|\lambda_j|=1`$). For $`0\le t<2\pi`$
   the arc $`K_t=\{e^{i\varphi}:0\le\varphi\le t\}`$ is closed and $`\{\theta_j\le t\}=\{N_{K_t}\ge j\}=F_{K_t,j}`$, so each
   $`\theta_j`$ is a Borel function on $`U(k)`$. Hence $`N_B(S)=\sum_j1_B(e^{i\theta_j(S)})`$ is a random variable, and
   $`\nu(B)=\frac1k\sum_jP(e^{i\theta_j(S)}\in B)`$ is a Borel probability measure on $`\mathbb T`$, the average of the laws
   of the $`k`$ eigenvalues.

   *Uniformity.* The eigenvalues of $`e^{i\alpha}S`$ are the $`e^{i\alpha}\lambda_j`$, so $`N_B(e^{i\alpha}S)=N_{e^{-i\alpha}B}(S)`$,
   and since $`e^{i\alpha}S`$ has the law of $`S`$, $`\nu(e^{-i\alpha}B)=\frac1k\mathbb EN_B(e^{i\alpha}S)=\nu(B)`$. So $`\nu`$ is
   invariant under rotations, and by Proposition A.22 (with $`N=1`$) it is the normalised arc length:
   $`\mathbb EN_B(S)=k\,|B|/2\pi`$ for every arc of length $`|B|`$. Likewise $`\mathrm{Tr}\,S^m=\sum_j\lambda_j^m`$ (also
   for $`m<0`$, as $`S^{-1}=S^*`$) is a bounded continuous function of $`S`$, and
   $`\mathbb E\mathrm{Tr}\,S^m=\mathbb E\mathrm{Tr}(e^{i\alpha}S)^m=e^{im\alpha}\,\mathbb E\mathrm{Tr}\,S^m`$ for all $`\alpha`$, so it
   vanishes for $`m\neq0`$: these are $`k`$ times the Fourier coefficients of $`\nu`$. In particular the mean
   eigenphase density is constant, and the spacings around the circle need no unfolding.

</details>

*Computation:* the numerical parts of Exercises 1–3 — [tutorial, Exercises](../tutorial/05_random_matrices.md#x-5).

**Next:** [Chapter 6](06_quantum_channels.md) treats quantum channels: Kraus and Stinespring
representations, duals with respect to a trace, and maps between sectors. The conditional
expectation $`E`$, met here as a Haar average, is the first example.
