**English** · [Español](es/05_random_matrices.md) · [Français](fr/05_random_matrices.md) · [Deutsch](de/05_random_matrices.md) · [中文](zh/05_random_matrices.md) · [日本語](ja/05_random_matrices.md) · [Русский](ru/05_random_matrices.md)

# 5. Random matrices in an algebra

Random matrix theory asks what a "typical" operator looks like. In an algebra
$\bigoplus_cM_{k_c}$ the answer has to respect the sectors: the natural random unitaries are
those of the algebra, not those of the big matrix algebra $M_D$, and the natural random
Hamiltonians have independent blocks. This chapter constructs the Haar measure on the unitary
group of the algebra, computes its first and second moments exactly, and then shows what goes
wrong in spectral statistics when the sectors are ignored.

**You will learn** why the QR decomposition of a Gaussian matrix gives a Haar unitary (and why
the phases must be fixed), that averaging over the unitary group of the algebra is the
conditional expectation onto the centre, the Weingarten formula of order two, the Wigner
surmise and the spacing ratio, and why level statistics computed on a spectrum with mixed
sectors look like those of independent levels.

Throughout, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ with $D=\sum_ck_c$, $P_c$ are the
minimal central projections, and $E$ is the conditional expectation onto the centre of
Proposition 2.8. $U(k)$ is the group of unitary $k\times k$ matrices and
$U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}$ the unitary group of $\mathcal M$. Random variables
live on a probability space with expectation $\mathbb E$; the expectation of a random matrix is
taken entrywise. A *standard complex Gaussian* variable has independent real and imaginary parts
$\sim N(0,\tfrac12)$, so that $\mathbb E|z|^2=1$ (this is what `torch.randn` draws for a complex
dtype); a *Ginibre matrix* has independent standard complex Gaussian entries.

## 5.1 Haar measure

Let $G$ be a compact group of matrices (a closed subgroup of some $U(n)$). A Borel probability
measure $\mu$ on $G$ is **left-invariant** if $\mu(gB)=\mu(B)$ for all $g\in G$ and Borel sets
$B$, equivalently $\int f(gh)\,d\mu(h)=\int f\,d\mu$ for every bounded Borel function $f$;
**right-invariant** is defined with $Bg$ and $f(hg)$.

**Lemma 5.1 (uniqueness).** If $\mu$ is a left-invariant and $\nu$ a right-invariant Borel
probability measure on $G$, then $\mu=\nu$. Hence $G$ carries at most one left-invariant Borel
probability measure, and if it exists it is also right-invariant and invariant under $g\mapsto g^{-1}$.

*Proof.* Let $f:G\to\mathbb R$ be continuous (hence bounded). By Fubini's theorem,
$\int\!\big(\int f(gh)\,d\mu(h)\big)d\nu(g)=\int\!\big(\int f(gh)\,d\nu(g)\big)d\mu(h)$. The inner
integral on the left is $\int f\,d\mu$ by left invariance, the one on the right is $\int f\,d\nu$
by right invariance; so $\int f\,d\mu=\int f\,d\nu$. Two Borel probability measures on a metric
space with the same integrals of continuous functions are equal: the continuous functions
$\max(0,1-n\operatorname{dist}(\cdot,K))$ decrease to the indicator of a closed set $K$, so the measures
agree on closed sets (dominated convergence), which generate the Borel sets and are closed under
intersections (Dynkin's $\pi$–$\lambda$ theorem). If $\mu_1,\mu_2$ are left-invariant, then
$\nu(B)=\mu_2(B^{-1})$ is right-invariant, since $(Bg)^{-1}=g^{-1}B^{-1}$; hence $\mu_1=\nu=\mu_2$,
so $\mu_1$ is right-invariant and inversion-invariant. $\square$

**Definition.** The **Haar measure** of $G$ is its left-invariant Borel probability measure. A
random element with this law is **Haar distributed**.

Every compact group has a Haar measure (see G. B. Folland, *A Course in Abstract Harmonic
Analysis*, 2nd ed., CRC Press, 2016, Section 2.2; stated here without proof). For $U(k)$ we do
not need this general theorem: Theorem 5.3 constructs the measure.

**Lemma 5.2 (QR decomposition).** Every invertible $Z\in M_k(\mathbb C)$ can be written $Z=QR$
with $Q$ unitary and $R$ upper triangular with positive diagonal. $Q$ and $R$ are unique and
depend continuously on $Z$.

*Proof.* Apply Gram–Schmidt to the columns $z_1,\dots,z_k$ of $Z$, which are linearly independent:
$q_j=w_j/\|w_j\|$ with $w_j=z_j-\sum_{i<j}\langle q_i,z_j\rangle q_i$; here $w_j\neq0$ because
$z_j\notin\operatorname{span}\{z_1,\dots,z_{j-1}\}=\operatorname{span}\{q_1,\dots,q_{j-1}\}$. Then
$z_j=\sum_{i<j}\langle q_i,z_j\rangle q_i+\|w_j\|q_j$, i.e. $Z=QR$ with $R_{jj}=\|w_j\|>0$, and
every step is continuous on the invertible matrices. If $Q_1R_1=Q_2R_2$, then
$T=Q_2^*Q_1=R_2R_1^{-1}$ is unitary and upper triangular with positive diagonal. Its inverse
$T^*$ is lower triangular, and also upper triangular (the inverse of an upper triangular matrix
is upper triangular), so $T$ is diagonal; a diagonal unitary with positive diagonal is $1$. $\square$

**Theorem 5.3 (Haar unitaries from Gaussian matrices).** Let $Z$ be a $k\times k$
Ginibre matrix. Then $Z$ is invertible with probability one, and the factor $Q$ of Lemma 5.2 is
Haar distributed on $U(k)$. Moreover, if $Z=Q'R'$ is any QR decomposition with $Q'$ unitary and $R'$
upper triangular (for instance the one returned by LAPACK, whose diagonal is real of either sign), then
$Q=Q'\Lambda$ with $\Lambda=\operatorname{diag}(R'_{ii}/|R'_{ii}|)$ is Haar distributed.

*Proof.* *Invertibility.* $\det Z$ is a nonzero polynomial in the $2k^2$ real coordinates of $Z$,
and the zero set of a nonzero real polynomial has Lebesgue measure zero (induction on the number
of variables: write $p=\sum_jp_j(x')x_n^j$ with some $p_j\neq0$; for $x'$ outside the null zero set of
$p_j$ there are finitely many roots $x_n$; conclude by Fubini). $Z$ has a density, so
$P(\det Z=0)=0$.

*Invariance of the Gaussian law.* The density of $Z$ with respect to Lebesgue measure on
$M_k(\mathbb C)\cong\mathbb R^{2k^2}$ is $\pi^{-k^2}e^{-\|Z\|_2^2}$. For $V\in U(k)$ the map
$Z\mapsto VZ$ is real-linear and preserves $\|Z\|_2$ (Section 1.4), so it is an orthogonal map of
$\mathbb R^{2k^2}$ with Jacobian $1$, and it preserves the density: $VZ$ has the same law as $Z$.

*Equivariance.* If $Z=QR$, then $VZ=(VQ)R$ is the decomposition of Lemma 5.2 of $VZ$, so
$Q(VZ)=VQ(Z)$. Hence for every Borel $B\subseteq U(k)$,
$P(VQ(Z)\in B)=P(Q(VZ)\in B)=P(Q(Z)\in B)$: the law of $Q(Z)$ is left-invariant, and by
Lemma 5.1 it is the Haar measure (this also proves that the Haar measure of $U(k)$ exists).
Finally $Q'R'=(Q'\Lambda)(\Lambda^*R')$, and $\Lambda^*R'$ is upper triangular with diagonal
$|R'_{ii}|>0$, so $Q'\Lambda=Q(Z)$ by uniqueness. $\square$

The same proof with real Gaussian matrices gives the Haar measure of the orthogonal group $O(k)$.
Without the correction $\Lambda$ the result is *not* Haar distributed: LAPACK's convention for
the signs of $R'_{ii}$ introduces a bias, which the code below exhibits
(F. Mezzadri, *How to generate random matrices from the classical compact groups*, Notices AMS 54
(2007), 592–604).

**Corollary 5.4 (Haar measure of the algebra).** $U(\mathcal M)=\{\bigoplus_cU_c:U_c\in U(k_c)\}\cong\prod_cU(k_c)$,
and $U=\bigoplus_cU_c$ is Haar distributed on $U(\mathcal M)$ iff the blocks $U_c$ are independent
and each $U_c$ is Haar distributed on $U(k_c)$.

*Proof.* $U^*U=\bigoplus_cU_c^*U_c$, so $U$ is unitary iff every block is. Let the $U_c$ be
independent and Haar distributed. For $V=\bigoplus_cV_c\in U(\mathcal M)$, the blocks $V_cU_c$ of
$VU$ are independent (functions of independent variables) and $V_cU_c$ has the law of $U_c$, so
$VU$ has the law of $U$. The law of $U$ is therefore left-invariant, and by Lemma 5.1 it is the Haar
measure of the compact group $U(\mathcal M)$. Conversely, by uniqueness every Haar distributed $U$ has
this law. $\square$

The group $U(\mathcal M)$ is much smaller than $U(D)$: for $C\ge2$ a Haar unitary of $U(D)$ lies
in $\mathcal M$ with probability zero (by Theorem 5.3 and the Gram–Schmidt construction in
Lemma 5.2, its first column is $z/\|z\|$ for a Gaussian vector $z$, so all its entries are nonzero
almost surely). The Haar measure of $SU(k)$ is obtained in the same spirit: if $U$ is Haar on
$U(k)$ and $\theta=\arg\det U$, then $V=e^{-i\theta/k}U\in SU(k)$, and for $W\in SU(k)$,
$\det(WU)=\det U$ gives $V(WU)=WV(U)$; since $WU$ has the law of $U$, the law of $V$ is
left-invariant on $SU(k)$, hence the Haar measure of $SU(k)$ by Lemma 5.1.

**Corollary 5.5 (first moments of the entries).** If $U$ is Haar distributed on $U(k)$, then
$\mathbb EU_{ij}=0$ and $\mathbb E|U_{ij}|^2=1/k$.

*Proof.* $-1\in U(k)$, so $-U$ has the law of $U$ and $\mathbb EU=-\mathbb EU$. For a permutation
matrix $P$, $PU$ has the law of $U$, so $\mathbb E|U_{ij}|^2$ does not depend on $i$; as
$\sum_i|U_{ij}|^2=1$, it equals $1/k$. $\square$

The library draws $U(\mathcal M)$-Haar unitaries with `random_unitary_operator`, one batched QR
decomposition per block size for all sectors and samples at once (on a GPU, Gram–Schmidt, which
gives the same $Q$ by the uniqueness in Lemma 5.2). The measure `'haar'` is Theorem 5.3 in every
block; for a complex algebra `'haar_su'` multiplies by $e^{-i\theta/k}$ as above; `'coe'` and
`'cse'` are Dyson's circular ensembles (the first is Exercise 3). We check Corollaries 5.4 and 5.5,
and that the library's sampler is exactly the corrected QR of Theorem 5.3:

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinChain, channels

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def mc(samples, exact, name):
    """Monte Carlo mean of real samples with its standard error; must agree within 5 standard errors."""
    mean, err = samples.mean().item(), samples.std().item() / math.sqrt(samples.numel())
    print(f"{name}: {mean:.4f} +- {err:.4f}   (theory {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12
    return mean

alg = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')        # M_1 (+) M_2 (+) M_3, D = 6
B = 100_000
U = alg.random_unitary_operator(B)                                   # Haar on U(1) x U(2) x U(3)
assert close(U.adjoint() @ U, alg.identity())                        # unitary elements of the algebra
u = [U.matrix[:, c, :k, :k] for c, k in enumerate(alg.k_factors)]   # the blocks U_c
for c, k in enumerate(alg.k_factors):                                # Corollary 5.5
    mc(u[c][:, 0, 0].real, 0.0, f"E Re (U_{c})_11")
    mc(u[c][:, 0, 0].abs() ** 2, 1 / k, f"E |(U_{c})_11|^2")
mc(u[1][:, 0, 0].abs() ** 2 * u[2][:, 0, 0].abs() ** 2, 1 / 6, "E |(U_1)_11|^2 |(U_2)_11|^2")   # independent blocks

SU = alg.random_unitary_operator(4, measure='haar_su').matrix
assert torch.allclose(torch.linalg.det(SU[:, 2, :3, :3]), torch.ones(4, dtype=dt))

torch.manual_seed(1)
Z = torch.randn(B, 2, 2, dtype=dt)                                   # Ginibre matrices
torch.manual_seed(1)
V = alg.random_unitary(2, batch_size=B)                              # the library, same random numbers
Q, R = torch.linalg.qr(Z)                                            # LAPACK: diag(R) real, either sign
d = R.diagonal(dim1=-2, dim2=-1)
assert close(V, Q * (d / d.abs()).unsqueeze(-2))                     # Theorem 5.3: V = Q Lambda
raw = Q.diagonal(dim1=-2, dim2=-1).sum(-1).real                      # Re Tr Q without the correction
print(f"E Re Tr Q without the phase correction: {raw.mean():.4f} +- {raw.std() / math.sqrt(B):.4f}")
assert abs(raw.mean()) > 20 * raw.std() / math.sqrt(B)               # far from 0: not Haar
```

Without the correction, $\mathbb E\operatorname{Re}\operatorname{Tr}Q\approx-0.84$ for $k=2$, whereas a
Haar unitary has $\mathbb EU=0$.

## 5.2 Moments and the Weingarten calculus

**Lemma 5.6.** If $T\in M_k(\mathbb C)$ commutes with all diagonal unitaries and all permutation
matrices, then $T\in\mathbb C1$. In particular the commutant of $U(k)$ is $\mathbb C1$.

*Proof.* For $\Delta=\operatorname{diag}(e^{i\theta_1},\dots,e^{i\theta_k})$ with distinct
$e^{i\theta_j}$, $\Delta T=T\Delta$ gives $(e^{i\theta_i}-e^{i\theta_j})T_{ij}=0$, so $T$ is
diagonal. For the permutation matrix $P_\sigma e_j=e_{\sigma(j)}$, $P_\sigma TP_\sigma^{-1}=T$ gives
$T_{\sigma(j)\sigma(j)}=T_{jj}$, so all diagonal entries are equal. $\square$

**Proposition 5.7 (the Haar twirl is the conditional expectation onto the centre).** Let $U$ be
Haar distributed on $U(\mathcal M)$. For every $X\in\mathcal M$,
$$
\mathbb E\,[UXU^*]=E(X)=\bigoplus_c\frac{\operatorname{Tr}X_c}{k_c}\,1_{k_c}.
$$

*Proof.* $(UXU^*)_c=U_cX_cU_c^*$, and $U_c$ is Haar on $U(k_c)$ (Corollary 5.4). Let
$T_c=\mathbb E[U_cX_cU_c^*]$ (the entries are bounded). For $V\in U(k_c)$, $VU_c$ has the law of
$U_c$, so $VT_cV^*=\mathbb E[(VU_c)X_c(VU_c)^*]=T_c$. By Lemma 5.6, $T_c=\lambda1$, and
$k_c\lambda=\operatorname{Tr}T_c=\mathbb E\operatorname{Tr}(U_cX_cU_c^*)=\operatorname{Tr}X_c$. $\square$

The normalisation is that of Proposition 2.8 and of `channels.center_expectation`. If
$X\in\mathcal M$ commutes with $U(\mathcal M)$, then $UXU^*=X$, so $X=E(X)$ is central; the Haar
average is the continuous analogue of the group average of Chapter 2, Exercise 3. If instead $U$ is
Haar on $U(D)$, the same proof (with $C=1$) gives $\mathbb E[UXU^*]=\frac{\operatorname{Tr}X}{D}1$: the
sector weights are lost, and $UXU^*$ is not even in $\mathcal M$.

> **Physical meaning.** A Haar-random unitary of the algebra is a random evolution that respects
> the superselection rule. Averaging over it scrambles everything inside each sector and keeps
> exactly the classical part: the weight of each charge. This is the "twirl" of quantum
> information, and it is the conditional expectation onto the centre of Chapter 2.

For second moments we need the commutant of $\{V\otimes V\}$ on $\mathbb C^k\otimes\mathbb C^k$.
Let $F$ be the **swap**, $F(x\otimes y)=y\otimes x$; it is self-adjoint and unitary, and
$\operatorname{Tr}1=k^2$, $\operatorname{Tr}F=\sum_{i,j}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k$, $F^2=1$.

**Lemma 5.8.** Let $k\ge2$. A matrix $X$ on $\mathbb C^k\otimes\mathbb C^k$ commutes with $V\otimes V$
for all $V\in U(k)$ iff $X=a1+bF$ with $a,b\in\mathbb C$; and $1,F$ are linearly independent.

*Proof.* $F(V\otimes V)(x\otimes y)=Vy\otimes Vx=(V\otimes V)F(x\otimes y)$, so $1$ and $F$ commute
with all $V\otimes V$. Conversely let $X$ commute with them, and write
$X_{ij,pq}=\langle e_i\otimes e_j,X(e_p\otimes e_q)\rangle$.
*Diagonal unitaries:* $\Delta\otimes\Delta$ multiplies $e_i\otimes e_j$ by $e^{i(\theta_i+\theta_j)}$,
so $(e^{i(\theta_i+\theta_j)}-e^{i(\theta_p+\theta_q)})X_{ij,pq}=0$. If the multisets $\{i,j\}$ and
$\{p,q\}$ differ, some index $m$ occurs $n\in\{\pm1,\pm2\}$ times more often in the first; with
$\theta_m=\pi/2$ and all other $\theta=0$ the bracket is $e^{i(\theta_p+\theta_q)}(e^{in\pi/2}-1)\neq0$,
so $X_{ij,pq}=0$. Hence $X(e_i\otimes e_i)=a_ie_i\otimes e_i$, and for $i\neq j$,
$X(e_i\otimes e_j)=b_{ij}\,e_i\otimes e_j+c_{ij}\,e_j\otimes e_i$.
*Permutations:* applying $X(P_\sigma\otimes P_\sigma)=(P_\sigma\otimes P_\sigma)X$ to $e_i\otimes e_j$ gives
$a_{\sigma(i)}=a_i$, $b_{\sigma(i)\sigma(j)}=b_{ij}$, $c_{\sigma(i)\sigma(j)}=c_{ij}$. The permutations act
transitively on ordered pairs of distinct indices, so $a_i=a'$, $b_{ij}=b'$, $c_{ij}=c'$; thus
$X=b'1+c'F$ on the span of the $e_i\otimes e_j$ with $i\ne j$.
*A rotation:* let $R$ be the unitary with $Re_1=(e_1+e_2)/\sqrt2$, $Re_2=(-e_1+e_2)/\sqrt2$, $Re_m=e_m$ for $m\ge3$.
Then $(R\otimes R)(e_1\otimes e_1)=\frac12(e_1\otimes e_1+e_1\otimes e_2+e_2\otimes e_1+e_2\otimes e_2)$, and
$X(R\otimes R)(e_1\otimes e_1)=\frac12\big(a'e_1\otimes e_1+(b'+c')(e_1\otimes e_2+e_2\otimes e_1)+a'e_2\otimes e_2\big)$
must equal $(R\otimes R)X(e_1\otimes e_1)=a'(R\otimes R)(e_1\otimes e_1)$. Comparing the coefficients of
$e_1\otimes e_2$ gives $a'=b'+c'$, which is the value of $b'1+c'F$ on $e_i\otimes e_i$. So $X=b'1+c'F$.
Finally $1$ and $F$ map $e_1\otimes e_2$ to the independent vectors $e_1\otimes e_2$ and $e_2\otimes e_1$. $\square$

**Corollary 5.9 (Schur–Weyl duality of order two).** Let $k\ge2$ and $P_\pm=\frac12(1\pm F)$.
They are the projections onto the symmetric and antisymmetric subspaces
$\mathrm{Sym}=\{x:Fx=x\}$ and $\mathrm{Alt}=\{x:Fx=-x\}$, of dimensions $k(k+1)/2$ and $k(k-1)/2$.
The representation $V\mapsto V\otimes V$ of $U(k)$ is irreducible on $\mathrm{Sym}$ and on $\mathrm{Alt}$,
and the two representations are inequivalent.

*Proof.* $F=F^*$ and $F^2=1$ make $P_\pm$ projections with $P_++P_-=1$ and ranges $\mathrm{Sym}$,
$\mathrm{Alt}$; their dimensions are $\operatorname{Tr}P_\pm=(k^2\pm k)/2$. The span $\mathcal A$ of
$\{V\otimes V\}$ is a $*$-algebra ($1\otimes1$, products $VW\otimes VW$, adjoints $V^*\otimes V^*$) with
commutant $\mathcal A'=\operatorname{span}\{1,F\}=\operatorname{span}\{P_+,P_-\}$ by Lemma 5.8. Let
$K\subseteq\mathrm{Sym}$ be invariant under all $V\otimes V$. By Lemma 2.1 its projection is in $\mathcal A'$, so it is
$xP_++yP_-$ with $x,y\in\{0,1\}$ (it is a projection); it vanishes on $\mathrm{Alt}$, so $y=0$, and
$K=0$ or $K=\mathrm{Sym}$. The same argument applies to $\mathrm{Alt}$. If
$T:\mathrm{Sym}\to\mathrm{Alt}$ intertwines the two representations, then $TP_+$ (viewed as an
operator on $\mathbb C^k\otimes\mathbb C^k$) commutes with $\mathcal A$, so $TP_+=xP_++yP_-$; it
vanishes on $\mathrm{Alt}$, so $y=0$, and its range lies in $\mathrm{Alt}\cap\mathrm{Sym}=0$, so $T=0$. $\square$

**Theorem 5.10 (Weingarten formula of order two).** Let $k\ge2$, $U$ Haar distributed on $U(k)$,
and $X$ a matrix on $\mathbb C^k\otimes\mathbb C^k$. Then
$$
\mathbb E\,[(U\otimes U)X(U\otimes U)^*]=\alpha\,1+\beta F,\qquad
\alpha=\frac{\operatorname{Tr}X-\operatorname{Tr}(FX)/k}{k^2-1},\quad
\beta=\frac{\operatorname{Tr}(FX)-\operatorname{Tr}X/k}{k^2-1}.
$$
Equivalently, with $\mathrm{Wg}(e)=\frac1{k^2-1}$, $\mathrm{Wg}(\tau)=-\frac1{k(k^2-1)}$ on the
group $S_2=\{e,\tau\}$,
$$
\mathbb E\big[U_{i_1j_1}U_{i_2j_2}\overline{U_{i'_1j'_1}U_{i'_2j'_2}}\big]
=\sum_{\sigma,\pi\in S_2}\delta_{i_1i'_{\sigma(1)}}\delta_{i_2i'_{\sigma(2)}}\,
\delta_{j_1j'_{\pi(1)}}\delta_{j_2j'_{\pi(2)}}\,\mathrm{Wg}(\sigma\pi^{-1}).
$$

*Proof.* Let $T=\mathbb E[(U\otimes U)X(U\otimes U)^*]$. As in Proposition 5.7, $T$ commutes with
all $V\otimes V$, so $T=\alpha1+\beta F$ by Lemma 5.8. Taking traces,
$\operatorname{Tr}T=\operatorname{Tr}X$, and, since $F$ commutes with $U\otimes U$,
$\operatorname{Tr}(FT)=\mathbb E\operatorname{Tr}((U\otimes U)^*F(U\otimes U)X)=\operatorname{Tr}(FX)$.
This gives $k^2\alpha+k\beta=\operatorname{Tr}X$ and $k\alpha+k^2\beta=\operatorname{Tr}(FX)$, a
linear system with determinant $k^2(k^2-1)\neq0$, whose solution is the one stated. For the
entries take $X=|e_{j_1}\otimes e_{j_2}\rangle\langle e_{j'_1}\otimes e_{j'_2}|$: then
$\langle e_{i_1}\otimes e_{i_2},T\,e_{i'_1}\otimes e_{i'_2}\rangle$ is the left side, $\operatorname{Tr}X=\delta_{j_1j'_1}\delta_{j_2j'_2}$
(the term $\pi=e$), $\operatorname{Tr}(FX)=\delta_{j_1j'_2}\delta_{j_2j'_1}$ (the term $\pi=\tau$), and
$1$, $F$ contribute $\delta_{i_1i'_1}\delta_{i_2i'_2}$ (the term $\sigma=e$) and
$\delta_{i_1i'_2}\delta_{i_2i'_1}$ ($\sigma=\tau$). Collecting terms,
$\alpha=\sum_\pi\mathrm{Wg}(e\pi^{-1})[\pi]$ and $\beta=\sum_\pi\mathrm{Wg}(\tau\pi^{-1})[\pi]$,
where $[\pi]$ denotes the $j$-deltas of $\pi$. $\square$

The matrix $\big(\mathrm{Wg}(\sigma^{-1}\pi)\big)_{\sigma,\pi}$ is the inverse of the Gram matrix
$\begin{pmatrix}k^2&k\\k&k^2\end{pmatrix}$ of $1$ and $F$ with respect to
$\langle A,B\rangle=\operatorname{Tr}(A^*B)$; this is the form that generalises.

**Corollary 5.11.** For $k\ge2$:
$$
\mathbb E|U_{11}|^4=\frac2{k(k+1)},\qquad
\mathbb E|U_{11}|^2|U_{12}|^2=\frac1{k(k+1)},\qquad
\mathbb E\big[U_{11}U_{22}\overline{U_{12}U_{21}}\big]=-\frac1{(k-1)k(k+1)}.
$$

*Proof.* In Theorem 5.10, (i) all indices equal $1$: every product of deltas is $1$, so the
moment is $2(\mathrm{Wg}(e)+\mathrm{Wg}(\tau))=2\frac{k-1}{k(k^2-1)}$. (ii) $i=(1,1)$, $j=(1,2)$,
$i'=(1,1)$, $j'=(1,2)$: only $\pi=e$ survives, giving $\mathrm{Wg}(e)+\mathrm{Wg}(\tau)$.
(iii) $i=(1,2)$, $j=(1,2)$, $i'=(1,2)$, $j'=(2,1)$: only $\sigma=e$, $\pi=\tau$ survive,
giving $\mathrm{Wg}(\tau)$. $\square$

The last moment is negative: the entries of a Haar unitary are not independent, and the
correlations are of relative order $1/k$. For the algebra, entries of different blocks are
independent (Corollary 5.4), so mixed moments across sectors factorise into moments of single blocks.

For moments of any order the same reasoning applies once the commutant of
$\{V^{\otimes n}\}$ is known; it is spanned by the permutations of the $n$ tensor factors
(Schur–Weyl duality, see R. Goodman, N. R. Wallach, *Symmetry, Representations, and Invariants*,
Springer, 2009, Chapter 9). We state the result without proof.

**Theorem 5.12 (Weingarten formula; stated without proof).** Let $n\ge1$ and $U$ Haar on $U(k)$. Then
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\pi\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\pi(m)}}\;\mathrm{Wg}(\sigma\pi^{-1},k),
$$
where $\mathrm{Wg}(\cdot,k)$ is a function on $S_n$ that depends only on the cycle type. For $n\le k$
the matrix $\big(\mathrm{Wg}(\sigma^{-1}\pi,k)\big)_{\sigma,\pi}$ is the inverse of the Gram
matrix $\big(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)}\big)_{\sigma,\pi}$; for $n>k$ this matrix is singular
and its pseudo-inverse is taken. Moments with different numbers of $U$ and $\overline U$ vanish. (D. Weingarten,
J. Math. Phys. 19 (1978), 999–1001; B. Collins, Int. Math. Res. Not. 2003, 953–982; B. Collins,
P. Śniady, Comm. Math. Phys. 264 (2006), 773–795.)

We check Proposition 5.7 against `channels.center_expectation`, the commutant of Lemma 5.8 by
exact linear algebra, and Theorem 5.10 and Corollary 5.11 by Monte Carlo with the $U(3)$ block
of the samples above:

```python
X = alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])      # a fixed element
T = (U @ X @ U.adjoint()).matrix                                    # B samples of U X U*
mean, err = T.mean(0), T.std(0) / math.sqrt(B)
EX = channels.center_expectation(alg)(X).matrix[0]
random = err > 1e-12                                                 # entries that are not constant
z = ((mean - EX).abs()[random] / err[random]).max().item()
print(f"twirl over U(M) vs E(X): largest deviation {z:.2f} standard errors")
assert ((mean - EX).abs() <= 5 * err + 1e-12).all()                 # Proposition 5.7

def commutant_dim(ops, tol=1e-9):
    """Dimension of {X : XA = AX for all A in ops}, as the null space of X -> (AX - XA)_A."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    ops = [A.contiguous() for A in ops]                              # torch.kron needs contiguous input
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return int((torch.linalg.svdvals(L) < tol).sum())

k = 3
F = torch.zeros(k * k, k * k, dtype=dt)
for i in range(k):
    for j in range(k):
        F[j * k + i, i * k + j] = 1                                  # F (e_i (x) e_j) = e_j (x) e_i
VV = [torch.kron(W, W) for W in alg.random_unitary(k, batch_size=3).contiguous()]
assert all(close(F @ W, W @ F) for W in VV)
assert commutant_dim(VV) == 2                                        # Lemma 5.8: span{1, F}

u3 = u[2]                                                            # Haar on U(3)
mc(u3[:, 0, 0].abs() ** 4, 2 / (k * (k + 1)), "E|U_11|^4")
mc(u3[:, 0, 0].abs() ** 2 * u3[:, 0, 1].abs() ** 2, 1 / (k * (k + 1)), "E|U_11|^2 |U_12|^2")
mc((u3[:, 0, 0] * u3[:, 1, 1] * (u3[:, 0, 1] * u3[:, 1, 0]).conj()).real,
   -1 / ((k - 1) * k * (k + 1)), "E U_11 U_22 conj(U_12 U_21)")

n_s = 20_000                                                         # Theorem 5.10 for a random X
X2 = torch.randn(k * k, k * k, dtype=dt)
trX, trFX = X2.diagonal().sum(), (F @ X2).diagonal().sum()
alpha, beta = (trX - trFX / k) / (k ** 2 - 1), (trFX - trX / k) / (k ** 2 - 1)
UU = torch.einsum('bij,bkl->bikjl', u3[:n_s], u3[:n_s]).reshape(n_s, k * k, k * k)   # U (x) U
T2 = UU @ X2 @ UU.mH
mean, err = T2.mean(0), T2.std(0) / math.sqrt(n_s)
assert ((mean - alpha * torch.eye(k * k, dtype=dt) - beta * F).abs() <= 5 * err).all()
```

## 5.3 Level statistics

Let $H$ be a random self-adjoint $n\times n$ matrix with eigenvalues
$\lambda_1\le\dots\le\lambda_n$. Its **spacings** are $s_i=\lambda_{i+1}-\lambda_i$. Two ensembles
are standard: the **GOE**, $H=(G+G^T)/2$ with $G$ real with independent $N(0,1)$ entries, and the
**GUE**, $H=(G+G^*)/\sqrt2$ with $G$ Ginibre (Section 1.7). In both, the diagonal entries are
$N(0,1)$, the real and imaginary parts of the entries above the diagonal are $N(0,\tfrac12)$, and the
density of $H$ is proportional to $e^{-\operatorname{Tr}H^2/2}$.

**Proposition 5.13 (Wigner surmise).** Let $n=2$, $s=\lambda_2-\lambda_1$ and $S=s/\mathbb Es$.
For the GOE ($\beta=1$) and the GUE ($\beta=2$), $S$ has the densities
$$
p_1(S)=\frac\pi2\,S\,e^{-\pi S^2/4},\qquad p_2(S)=\frac{32}{\pi^2}\,S^2\,e^{-4S^2/\pi},
$$
with distribution functions $P(S\le t)=1-e^{-\pi t^2/4}$ and
$P(S\le t)=\operatorname{erf}(2t/\sqrt\pi)-\frac{4t}\pi e^{-4t^2/\pi}$. In particular
$P(S\le\varepsilon)$ is of order $\varepsilon^{\beta+1}$ as $\varepsilon\to0$.

*Proof.* For $H=\begin{pmatrix}a&b\\\bar b&c\end{pmatrix}$ the eigenvalues are
$\frac{a+c}2\pm\frac12\sqrt{(a-c)^2+4|b|^2}$, so $s=\sqrt{(a-c)^2+4|b|^2}$. *GOE:* $x=a-c$ and
$y=2b$ are independent $N(0,2)$ variables, and in polar coordinates
$P(s\le t)=\int_{x^2+y^2\le t^2}\frac{e^{-(x^2+y^2)/4}}{4\pi}\,dx\,dy=\int_0^t\frac r2e^{-r^2/4}dr=1-e^{-t^2/4}$.
Then $\mathbb Es=\int_0^\infty P(s>t)\,dt=\int_0^\infty e^{-t^2/4}dt=\sqrt\pi$, and
$P(S\le t)=P(s\le\sqrt\pi t)=1-e^{-\pi t^2/4}$; differentiate. *GUE:* $x=a-c$, $y=2\operatorname{Re}b$,
$y'=2\operatorname{Im}b$ are independent $N(0,2)$, and in spherical coordinates the density of
$s=\sqrt{x^2+y^2+y'^2}$ is $4\pi r^2(4\pi)^{-3/2}e^{-r^2/4}=\frac{r^2e^{-r^2/4}}{2\sqrt\pi}$. Its mean is
$\frac1{2\sqrt\pi}\int_0^\infty r^3e^{-r^2/4}dr=\frac{8}{2\sqrt\pi}=\frac4{\sqrt\pi}$. Substituting
$r=4S/\sqrt\pi$ gives $p_2$, and integrating $p_2$ by parts
($\int x^2e^{-ax^2}dx=-\frac{x}{2a}e^{-ax^2}+\frac1{2a}\int e^{-ax^2}dx$ with $a=4/\pi$) gives the
distribution function. The small-$\varepsilon$ behaviour follows from $p_\beta(S)\sim cS^\beta$. $\square$

The exponent counts dimensions: a real symmetric $2\times2$ matrix is degenerate iff two real
conditions hold ($a=c$, $b=0$), a complex Hermitian one iff three hold. For every $n$ the joint
density of the eigenvalues is proportional to $\prod_{i<j}|\lambda_i-\lambda_j|^\beta e^{-\sum_i\lambda_i^2/2}$, and
as $n\to\infty$ the spacing distribution in the bulk converges to a limit (Gaudin–Mehta) that is close
to, but not equal to, the surmise (M. L. Mehta, *Random Matrices*, 3rd ed., Elsevier, 2004; stated without
proof).

> **Physical meaning.** Wigner proposed that the energy levels of complex quantum systems (heavy
> nuclei) are statistically like the eigenvalues of a random Hamiltonian with the right symmetry:
> real symmetric if there is a time-reversal symmetry ($\beta=1$), complex Hermitian otherwise
> ($\beta=2$). Levels then *repel*: two levels of the same symmetry almost never come close
> (von Neumann and Wigner, 1929). Repulsion is the typical signature of quantum chaos, independent
> levels (no repulsion) that of integrable or localised systems (the Bohigas–Giannoni–Schmit and
> Berry–Tabor conjectures).

Spacings depend on the local density of levels, which varies along the spectrum; comparing them
with $p_\beta$ requires **unfolding**, the change of variables $\lambda\mapsto\bar N(\lambda)$ by a
smooth approximation $\bar N$ of the counting function, which makes the mean spacing $1$. The choice
of $\bar N$ is a source of bias. The ratio of consecutive spacings avoids it.

**Definition.** The **spacing ratio** is $r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})\in[0,1]$
(V. Oganesyan, D. A. Huse, Phys. Rev. B 75 (2007), 155111). It is unchanged by
$\lambda\mapsto a\lambda+b$, $a>0$, so it needs no unfolding where the density is roughly constant
over two spacings.

**Proposition 5.14 (independent levels).** Let $s,s'$ be independent with the same exponential
distribution (as two consecutive spacings of a Poisson process). Then $r=\min(s,s')/\max(s,s')$ has
$P(r\le t)=\frac{2t}{1+t}$ on $[0,1]$, density $\frac2{(1+t)^2}$, and
$$
\mathbb Er=2\ln2-1\approx0.3863 .
$$

*Proof.* By scaling, the rate is $1$. For $x=s/s'$, $P(x\le t)=\mathbb E\,P(s\le ts'\mid s')=\mathbb E(1-e^{-ts'})=1-\frac1{1+t}=\frac t{1+t}$.
For $t<1$ the events $\{x\le t\}$ and $\{1/x\le t\}$ are disjoint, $x$ and $1/x$ have the same law
($s,s'$ are exchangeable), and $\{r\le t\}$ is their union, so $P(r\le t)=\frac{2t}{1+t}$. Then
$\mathbb Er=\int_0^1\frac{2t}{(1+t)^2}dt=2\Big[\ln(1+t)+\frac1{1+t}\Big]_0^1=2\ln2-1$. $\square$

The same holds exactly for $n$ independent uniform levels on an interval, because their spacings
have the law of $(e_1,\dots,e_{n+1})/\sum_je_j$ for independent exponential $e_j$, so that each
ratio $s_i/s_{i+1}=e_i/e_{i+1}$ (W. Feller, *An Introduction to Probability Theory and Its
Applications*, Vol. II, Wiley, 1971, Chapter I; stated without proof). For the GOE and the GUE in the
limit $n\to\infty$, $\langle r\rangle=0.5307$ and $0.5996$ (numerical values; the $3\times3$ analogues
of the surmise give $4-2\sqrt3\approx0.5359$ and $2\sqrt3/\pi-\frac12\approx0.6027$; Y. Y. Atas,
E. Bogomolny, O. Giraud, G. Roux, Phys. Rev. Lett. 110 (2013), 084101; stated without proof).

The function `SpinChain.level_spacing_ratio` (in `chains`) returns the mean ratio of each row of
levels. We check Propositions 5.13 and 5.14 and the GOE value:

```python
ratio = SpinChain.level_spacing_ratio
n2 = 200_000
alg_r = TypeIAlgebra([2], [2], complex_valued=False, precision='double')     # real symmetric 2 x 2
alg_c = TypeIAlgebra([2], [2], precision='double')
G = torch.randn(n2, 1, 2, 2, dtype=torch.float64)
Gc = torch.randn(n2, 1, 2, 2, dtype=dt)
goe2 = alg_r.operator((G + G.mT) / 2, is_self_adjoint=True)
gue2 = alg_c.operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
cdf = {1: lambda t: 1 - math.exp(-math.pi * t ** 2 / 4),
       2: lambda t: math.erf(2 * t / math.sqrt(math.pi)) - 4 * t / math.pi * math.exp(-4 * t ** 2 / math.pi)}
for beta, H, mean_s in ((1, goe2, math.sqrt(math.pi)), (2, gue2, 4 / math.sqrt(math.pi))):
    lam = H.eigenvalues()[0]                                         # (n2, 2), ascending
    S = (lam[:, 1] - lam[:, 0]) / mean_s
    mc(S, 1.0, f"beta={beta}: E S")
    for t in (0.25, 1.0, 2.0):
        mc((S <= t).double(), cdf[beta](t), f"beta={beta}: P(S <= {t})")

levels = torch.rand(4000, 50, dtype=torch.float64)                   # independent levels
mc(ratio(levels), 2 * math.log(2) - 1, "<r>, independent uniform levels")

n = 64
alg64 = TypeIAlgebra([n], [n], complex_valued=False, precision='double')
G = torch.randn(800, 1, n, n, dtype=torch.float64)
lam = alg64.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()[0]
mc(ratio(lam[:, n // 4: 3 * n // 4]), 0.5307, "<r>, GOE(64), middle half of the spectrum")
```

## 5.4 Why sectors must not be mixed

Let $H=\bigoplus_cH_c\in\mathcal M$ be a random self-adjoint element with **independent blocks**,
for example GOE blocks. As an operator on $\mathbb C^D$ its spectrum is the union of the spectra
of the blocks. Inside a block the levels repel; we show that levels of different blocks do not.

Assume that each block has a bounded continuous **level density** $\rho_c$:
$\mathbb E\sum_i\varphi(\lambda_i(H_c))=\int\varphi\,\rho_c$ for bounded Borel $\varphi$ (so
$\int\rho_c=k_c$). For the GOE and the GUE this follows from the joint density quoted above.

**Theorem 5.15 (no repulsion between sectors).** For $c\neq d$ and every bounded Borel
$\psi:\mathbb R^2\to\mathbb R$,
$$
\mathbb E\sum_{i,j}\psi\big(\lambda_i(H_c),\lambda_j(H_d)\big)=\iint\psi(x,y)\,\rho_c(x)\rho_d(y)\,dx\,dy .
$$
Consequently, if $N_\varepsilon$ is the number of pairs of levels from different blocks at distance
$<\varepsilon$, then $\mathbb EN_\varepsilon/(2\varepsilon)\to\sum_{c<d}\int\rho_c\rho_d$ as
$\varepsilon\to0$, and the smallest spacing $\delta$ of the merged spectrum satisfies
$$
P(\delta<\varepsilon)\ \ge\ \frac{\mathbb EN_\varepsilon}{\sum_{c<d}k_ck_d}
=\frac{2\varepsilon\sum_{c<d}\int\rho_c\rho_d}{\sum_{c<d}k_ck_d}+o(\varepsilon).
$$

*Proof.* The ordered eigenvalues are continuous functions of the matrix, so
$\Lambda_c=(\lambda_i(H_c))_i$ and $\Lambda_d$ are independent random vectors. By Fubini (conditioning on
$\Lambda_d$), $\mathbb E\sum_{i,j}\psi(\Lambda_{c,i},\Lambda_{d,j})=\mathbb E\sum_j g(\Lambda_{d,j})$ with
$g(y)=\mathbb E\sum_i\psi(\Lambda_{c,i},y)=\int\psi(x,y)\rho_c(x)\,dx$, which is bounded and Borel; apply the
definition of $\rho_d$ to $g$. With $\psi(x,y)=1\{|x-y|<\varepsilon\}$,
$$
\frac{\mathbb EN_\varepsilon}{2\varepsilon}=\sum_{c<d}\int\rho_c(x)\Big(\frac1{2\varepsilon}\int_{x-\varepsilon}^{x+\varepsilon}\rho_d(y)\,dy\Big)dx ;
$$
the bracket tends to $\rho_d(x)$ by continuity and is bounded by $\sup\rho_d$, and $\rho_c$ is integrable,
so dominated convergence gives the limit. Finally $N_\varepsilon\le\sum_{c<d}k_ck_d$, so
$\mathbb EN_\varepsilon\le\big(\sum_{c<d}k_ck_d\big)P(N_\varepsilon\ge1)$, and $N_\varepsilon\ge1$ implies $\delta<\varepsilon$. $\square$

So for two blocks with overlapping spectra, $P(\delta<\varepsilon)$ is at least of order $\varepsilon$,
while inside a $2\times2$ GOE block it is $1-e^{-\varepsilon^2/4}\le\varepsilon^2/4$
(Proposition 5.13, unnormalised spacing), and for an $n\times n$ GOE block it is of order
$\varepsilon^2$ (from the joint density; not proved here). The level-repulsion test fails for the
merged spectrum although every block is a perfect GOE.

*Many sectors (a derivation, not a proof).* For stationary spectra, the probability $\mathcal G(s)$ of
finding no level in an interval of length $s$ determines the spacing density, $p(s)=\mathcal G''(s)$ for
unit mean spacing (Mehta, *op. cit.*). Independent spectra superpose by multiplying these gap
probabilities. If $C$ sectors each contribute a fraction $1/C$ of the levels and each has the gap
probability $\mathcal G_1$ in its own unfolded units, then $\mathcal G(s)=\mathcal G_1(s/C)^C$, with
$\mathcal G_1(x)=1-x+\int_0^x(x-t)p_1(t)\,dt$. Differentiating twice at $0$,
$$
p(0)=1-\frac1C+\frac{p_1(0)}C=1-\frac1C ,
$$
the analogue of Theorem 5.15 for unfolded spectra, and $\mathcal G_1(s/C)^C\to e^{-s}$, the gap probability
of a Poisson process: as $C\to\infty$ the merged spectrum looks like independent levels
(N. Rosenzweig, C. E. Porter, Phys. Rev. 120 (1960), 1698–1714; M. V. Berry, M. Robnik, J. Phys. A 17
(1984), 2413–2421; rigorous versions of the Poisson limit are limit theorems for superpositions of
sparse point processes, see B. Grigelionis, Theory Probab. Appl. 8 (1963)).

In the library the blocks are separate from the start: `eigenvalues()` returns one spectrum per
sector, and mixing them is an explicit `torch.cat`. We put $C$ independent GOE blocks of size $64$
in one algebra and compare the mean ratio per sector with that of the merged spectrum:

```python
k, merged_r = 64, {}
for C in (1, 2, 4, 8, 16):
    alg_C = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double')
    G = torch.randn(800 // C, C, k, k, dtype=torch.float64)
    lam = alg_C.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()    # one spectrum per sector
    per_sector = torch.stack([ratio(l[:, k // 4: 3 * k // 4]) for l in lam])
    merged = torch.cat(lam, dim=-1).sort(dim=-1).values                       # the spectrum on C^D
    r_m = ratio(merged[:, k * C // 4: 3 * k * C // 4])
    merged_r[C] = r_m.mean().item()
    print(f"C = {C:2d}: <r> per sector {per_sector.mean():.4f}, merged {merged_r[C]:.4f}"
          f" +- {r_m.std() / math.sqrt(len(r_m)):.4f}")
    assert abs(per_sector.mean() - 0.5307) < 0.01                    # every sector is GOE
assert merged_r[1] > merged_r[2] > merged_r[4] > merged_r[16]
assert abs(merged_r[16] - (2 * math.log(2) - 1)) < 0.01              # close to independent levels
```

Every sector has $\langle r\rangle\approx0.53$; the merged spectrum drops to $0.42$ for two sectors and
to $0.39$, the value for independent levels, for $8$ and $16$.

The same happens in a physical model. The Heisenberg chain with random fields (Chapter 9) conserves
the magnetisation, so its Hamiltonian lies in $\bigoplus_{N=0}^LM_{\binom LN}$; for weak disorder the
levels inside each sector have GOE statistics (numerically; A. Pal, D. A. Huse, Phys. Rev. B 82 (2010),
174411). Its blocks are *not* independent — all of them depend on the same random fields — so
Theorem 5.15 does not apply as stated; the experiment shows the same effect.

```python
chain = SpinChain(10, complex_valued=False)                         # sectors N = 0, ..., 10
lam = chain.random_field_heisenberg(W=0.5, batch_size=20).eigenvalues()
k5 = chain.sector_dim(5)                                            # 252
r_sector = ratio(lam[5][:, k5 // 4: 3 * k5 // 4]).mean().item()
spectrum = torch.cat(lam, dim=-1).sort(dim=-1).values               # all 1024 levels mixed
r_all = ratio(spectrum[:, 256:768]).mean().item()
print(f"Heisenberg chain, L = 10: <r> in sector N = 5: {r_sector:.4f}, all sectors mixed: {r_all:.4f}")
assert r_sector > 0.49 and r_all < 0.42
```

Inside the sector $N=5$ the ratio is $0.52$ (a finite-size value slightly below $0.5307$); the
spectrum with all sectors mixed gives $0.39$.

> **Physical meaning.** Levels with different conserved quantum numbers do not interact, so they
> cross freely; only levels of the same symmetry repel. A spectrum that mixes symmetry sectors
> therefore shows Poisson-like statistics even when the dynamics in each sector is fully chaotic,
> and an unresolved symmetry can be mistaken for integrability or localisation. The rule of the
> trade — resolve every symmetry before computing level statistics — is the statement that one
> should compute in the algebra $\bigoplus_cM_{k_c}$, not in $M_D$.

## 5.5 What to compare with

Haar unitaries: `scipy.stats.unitary_group` (the corrected QR of Theorem 5.3) and `ortho_group`,
QuTiP's Haar-random unitaries, and any `numpy`/`torch` QR followed by the phase correction.
Weingarten calculus: the symbolic packages IntU (Mathematica; Z. Puchała, J. A. Miszczak), RTNI
(Mathematica and Python; M. Fukuda, R. König, I. Nechita, J. Phys. A 52 (2019), 425303) and the
Python package Haarpy compute moments of every order symbolically; the formulas proved in this
chapter are their order-two case. Level statistics are usually computed with a few lines of `numpy`,
and exact-diagonalisation codes such as QuSpin supply spectra sector by sector. What the library
adds is the unitary group of an algebra as one object: independent Haar blocks of all sizes drawn in
batched calls, the twirl identified with the conditional expectation
`channels.center_expectation`, and spectra returned per sector, so that mixing sectors is never
the default.

## Exercises

1. Let $U$ be Haar on $U(\mathcal M)$ and $\operatorname{Tr}U=\sum_c\operatorname{Tr}U_c$. Show that
   $\mathbb E|\operatorname{Tr}U|^2=C$ and $\mathbb E|\operatorname{Tr}U|^4=\sum_cm(k_c)+2C(C-1)$ with
   $m(1)=1$, $m(k)=2$ for $k\ge2$; compare with $1$ and $2$ for a Haar unitary of $U(D)$, $D\ge2$. Check
   for $\mathcal M=M_1\oplus M_2\oplus M_3$.
2. Using Theorem 5.12, show that $\mathbb E|U_{11}|^{2n}=\frac{n!}{k(k+1)\cdots(k+n-1)}$ for $n\le k$.
   Compute the Weingarten matrix of order $3$ for $k=3$ numerically and check $\mathbb E|U_{11}|^6=\frac1{10}$.
3. The **circular orthogonal ensemble** (`measure='coe'`) is $S=W^TW$ with $W$ Haar on $U(k)$. Show
   that $S$ is unitary and symmetric, that its law is invariant under $S\mapsto V^TSV$ for every
   $V\in U(k)$, and that the mean density of its eigenphases is uniform. Check numerically that its
   eigenphases have $\langle r\rangle\approx0.53$, like the GOE.

<details>
<summary>Solutions</summary>

1. Write $t_c=\operatorname{Tr}U_c$; the $t_c$ are independent (Corollary 5.4). Multiplying $U_c$ by
   $e^{i\theta}$ does not change its law and multiplies $t_c$ by $e^{i\theta}$; so a moment
   $\mathbb E[t_at_b\bar t_c\bar t_d]$ vanishes unless every sector occurs as often in $(a,b)$ as in
   $(c,d)$, i.e. $\{a,b\}=\{c,d\}$ as multisets, and likewise $\mathbb E[t_a\bar t_b]=0$ for $a\neq b$. By
   Proposition 5.7, $\mathbb E[(U_c)_{ii}\overline{(U_c)_{jj}}]=\langle e_i,\mathbb E[U_ce_ie_j^*U_c^*]e_j\rangle=\delta_{ij}/k_c$,
   so $\mathbb E|t_c|^2=1$ and $\mathbb E|\operatorname{Tr}U|^2=C$. In $\mathbb E|\sum_at_a|^4$ the terms with
   $a=b$ give $\sum_a\mathbb E|t_a|^4$, those with $a\neq b$ give
   $2\sum_{a\neq b}\mathbb E|t_a|^2\,\mathbb E|t_b|^2=2C(C-1)$ (two choices of $(c,d)$). For $k=1$, $|t|=1$.
   For $k\ge2$, $\mathbb E|t|^4=\sum_{i,j,i',j'}\mathbb E[U_{ii}U_{jj}\overline{U_{i'i'}U_{j'j'}}]$; in
   Theorem 5.10 the two terms with $\sigma=\pi$ force $\{i',j'\}=\{i,j\}$ in one way and give
   $k^2\mathrm{Wg}(e)$ each after summation, the two with $\sigma\neq\pi$ force $i=j=i'=j'$ and give
   $k\,\mathrm{Wg}(\tau)$ each: in total $2\big(\frac{k^2}{k^2-1}-\frac1{k^2-1}\big)=2$. The case $C=1$,
   $k=D$ gives $1$ and $2$ for $U(D)$. For $M_1\oplus M_2\oplus M_3$:
   $\mathbb E|\operatorname{Tr}U|^4=1+2+2+12=17$, against $2$ for a Haar unitary of $U(6)$.

```python
tr = U.trace                                                        # Tr U = sum_c Tr U_c, shape (B,)
mc(tr.abs() ** 2, 3.0, "E|Tr U|^2")
mc(tr.abs() ** 4, 17.0, "E|Tr U|^4")
```

2. With all indices equal to $1$, every delta is $1$, so
   $\mathbb E|U_{11}|^{2n}=\sum_{\sigma,\pi}\mathrm{Wg}(\sigma\pi^{-1})=\mathbf 1^TW\mathbf 1$, where
   $W=\Gamma^{-1}=(\mathrm{Wg}(\sigma^{-1}\pi))$ and $\Gamma=(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)})$. Here we used that
   $\mathrm{Wg}$ depends only on the cycle type, and that $\sigma\pi^{-1}=\sigma(\pi^{-1}\sigma)\sigma^{-1}$ is
   conjugate to $\pi^{-1}\sigma=(\sigma^{-1}\pi)^{-1}$, which has the cycle type of $\sigma^{-1}\pi$. Every
   row sum of $\Gamma$ is $c_n(k)=\sum_{\rho\in S_n}k^{\#\mathrm{cycles}(\rho)}$, so $\Gamma\mathbf 1=c_n(k)\mathbf 1$,
   $W\mathbf 1=\mathbf 1/c_n(k)$ and $\mathbf 1^TW\mathbf 1=n!/c_n(k)$. Finally
   $c_n(x)=x(x+1)\cdots(x+n-1)$ by induction: a permutation of $\{1,\dots,n\}$ arises from one of
   $\{1,\dots,n-1\}$ either by adding $n$ as a fixed point (one more cycle, factor $x$) or by inserting $n$
   after one of the $n-1$ elements in its cycle ($n-1$ ways, same number of cycles); this is a bijection, so
   $c_n(x)=(x+n-1)c_{n-1}(x)$.

```python
import itertools

def n_cycles(p):                                                    # number of cycles of a permutation
    seen, count = set(), 0
    for i in range(len(p)):
        count += i not in seen
        while i not in seen:
            seen.add(i)
            i = p[i]
    return count

perms = list(itertools.permutations(range(3)))                      # S_3; s.index(j) = s^{-1}(j)
Gram = torch.tensor([[3.0 ** n_cycles([s.index(p[i]) for i in range(3)]) for p in perms] for s in perms],
                    dtype=torch.float64)                            # k^{#cycles(s^{-1} p)}, k = 3
Wg3 = torch.linalg.inv(Gram)                                        # the Weingarten matrix, n = k = 3
assert abs(Wg3.sum().item() - 6 / (3 * 4 * 5)) < 1e-12
mc(u3[:, 0, 0].abs() ** 6, 0.1, "E|U_11|^6")
```

3. $S^*S=W^*\bar WW^TW=1$ since $W^T$ is unitary, and $S^T=W^TW=S$. For $V\in U(k)$, $WV$ has the law
   of $W$ (right invariance, Lemma 5.1), and $(WV)^T(WV)=V^TSV$. With $V=e^{i\theta/2}1$, $e^{i\theta}S$
   has the law of $S$, so the law of the set of eigenphases is invariant under rotations, and so is the
   mean density, which is therefore uniform. In particular the spacings need no unfolding. The local
   statistics of the COE coincide with those of the GOE for large $k$ (F. J. Dyson, J. Math. Phys. 3
   (1962), 140–156; Mehta, *op. cit.*; stated without proof).

```python
alg_s = TypeIAlgebra([32], [32], precision='double')
S = alg_s.random_unitary_operator(400, measure='coe').matrix[:, 0]
assert close(S, S.mT) and close(S @ S.mH, torch.eye(32, dtype=dt))
phases = torch.linalg.eigvals(S).angle().sort(dim=-1).values
s = torch.diff(torch.cat([phases, phases[:, :1] + 2 * math.pi], dim=-1), dim=-1)   # around the circle
r = torch.minimum(s, s.roll(-1, -1)) / torch.maximum(s, s.roll(-1, -1))
mc(r.mean(-1), 0.5307, "<r>, COE(32) eigenphases")
```

</details>

**Next:** [Chapter 6](06_quantum_channels.md) treats quantum channels: Kraus and Stinespring
representations, duals with respect to a trace, and maps between sectors. The conditional
expectation $E$, met here as a Haar average, is the first example.
