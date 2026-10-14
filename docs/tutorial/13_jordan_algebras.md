**English** · [Español](es/13_jordan_algebras.md) · [Français](fr/13_jordan_algebras.md) · [Deutsch](de/13_jordan_algebras.md) · [中文](zh/13_jordan_algebras.md) · [日本語](ja/13_jordan_algebras.md) · [Русский](ru/13_jordan_algebras.md)

# 13. Jordan algebras: why quantum mechanics uses *-algebras

The product $xy$ of two observables is not an observable unless $x$ and $y$ commute, but the symmetrised product
$x\circ y=\frac12(xy+yx)$ always is. In 1933 P. Jordan proposed to take the observables with this product as the
basic structure of quantum mechanics, and Jordan, von Neumann and Wigner classified the finite-dimensional algebras
that arise. This chapter develops the finite-dimensional theory, with complete proofs here or, for the four long
ones (the classification, the Koecher–Vinberg theorem, Kadison's theorem and the Jordan automorphisms of
$\bigoplus_cM_{k_c}$), in Appendices E and F, and asks what the Jordan structure misses. The answer is the sign of the
commutator: the Jordan algebra of one system determines the commutator only up to sign (Section 13.5), but composite
systems need the sign: it decides which joint operators are positive, and the PPT test built from it detects
entanglement (Proposition 13.12).

**You will learn** the Jordan product and the quadratic representation $P(y)x=yxy$, which, unlike $\circ$,
preserves positivity; spin factors, their spectral theory, and the Lorentz cone; Jordan subalgebras of Hermitian
matrices that are not self-adjoint parts of $*$-algebras; the Jordan–von Neumann–Wigner classification; why the
transpose is a Jordan automorphism and not a $*$-automorphism, and what it does on a composite system; and a trace
inequality of Chapter 11 written in Jordan language and proved in a spin factor.

Throughout, $\operatorname{Herm}_n(\mathbb K)$ is the real vector space of Hermitian $n\times n$ matrices over
$\mathbb K=\mathbb R$ or $\mathbb C$ (Section 13.4 also uses the quaternions $\mathbb H$ and the octonions
$\mathbb O$, constructed in Section E.6, and $\operatorname{Herm}_n$ over them, Section E.7), $M_{sa}$ the self-adjoint part of $\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$,
$X,Y,Z$ the Pauli matrices and $\sigma=(X,Y,Z)$, so $u\cdot\sigma=u_1X+u_2Y+u_3Z$ for $u\in\mathbb R^3$.
$[a,b]=ab-ba$, and $|u|$ is the Euclidean norm of $u\in\mathbb R^m$.

## 13.1 Observables and the Jordan product

**Definition.** A **Jordan algebra** is a real vector space $J$ with a bilinear, commutative product $x\circ y$
that satisfies the **Jordan identity** $(x^2\circ y)\circ x=x^2\circ(y\circ x)$, where $x^2=x\circ x$. It is
**special** if it is isomorphic to a subspace of an associative algebra closed under $\frac12(xy+yx)$. On $J$,
$L(y)x=y\circ x$ and the **quadratic representation** is $P(y)=2L(y)^2-L(y^2)$.

**Lemma 13.1.** Let $\mathcal A$ be an associative real or complex algebra and $x\circ y=\frac12(xy+yx)$. Then
$\circ$ satisfies the Jordan identity, and the Jordan powers $x^{\circ1}=x$, $x^{\circ(k+1)}=x\circ x^{\circ k}$ equal
the associative powers $x^k$. In particular $M_{sa}$, $\operatorname{Herm}_n(\mathbb C)$ and
$\operatorname{Herm}_n(\mathbb R)$ are special Jordan algebras.

*Proof.* Expanding, $4(x^2\circ y)\circ x=x^2yx+yx^3+x^3y+xyx^2$ and
$4x^2\circ(y\circ x)=x^2yx+x^3y+yx^3+xyx^2$, the same four terms. By induction,
$x\circ x^k=\frac12(x^{k+1}+x^{k+1})=x^{k+1}$. If $x,y$ are Hermitian, $(xy+yx)^*=yx+xy$, so the three spaces are closed
under $\circ$; for real symmetric $x,y$ the product is real. $\square$

**Proposition 13.2.** Let $x,y,z\in M_{sa}$.

1. $P(y)x=yxy$, and the **fundamental formula** $P(P(y)x)=P(y)P(x)P(y)$ holds.
2. $\operatorname{Tr}P(y)x=\operatorname{Tr}y^2x$.
3. If $x\ge0$ then $P(y)x\ge0$. If $y$ is invertible, $P(y)$ is a bijection of the positive cone onto itself, with
   inverse $P(y^{-1})$.
4. The Jordan product of positive elements need not be positive: for $p=\begin{pmatrix}1&0\\0&0\end{pmatrix}$ and
   $q=\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}$, $\det(p\circ q)=-\frac1{16}$.

*Proof.* (1) $2y\circ(y\circ x)=\frac12(y^2x+2yxy+xy^2)$ and $y^2\circ x=\frac12(y^2x+xy^2)$; subtract. Then
$P(P(y)x)z=yxy\,z\,yxy=P(y)P(x)P(y)z$ by associativity. (2) Cyclicity of the trace. (3)
$yxy=(x^{1/2}y)^*(x^{1/2}y)\ge0$ (Corollary 1.7), and $P(y)P(y^{-1})=P(y^{-1})P(y)=\mathrm{id}$. (4)
$p\circ q=\frac14\begin{pmatrix}2&1\\1&0\end{pmatrix}$. $\square$

With $P(y)$ written as $2L(y)^2-L(y^2)$, every identity in (1) is expressed with $\circ$ alone. The proof of (1) uses
only the associativity of the matrix product, so in every special Jordan algebra $J\subseteq\mathcal A$ ($\mathcal A$
associative) it gives $P(y)x=yxy$ and the fundamental formula $P(P(y)x)=P(y)P(x)P(y)$ on $J$. (This covers every
algebra in which the course uses it; the fundamental formula is not used for non-special algebras.)

Among the products $\alpha xy+\beta yx$ ($\alpha,\beta\in\mathbb R$) on $\operatorname{Herm}_n(\mathbb C)$, $n\ge2$,
$\frac12(xy+yx)$ is the only one that is symmetric in $x,y$ and equals $x^2$ for $y=x$: symmetry means
$(\alpha-\beta)(xy-yx)=0$ for all $x,y$, and $XZ\ne ZX$ (in the upper left $2\times2$ block) forces $\alpha=\beta$;
then $(\alpha+\beta)1=1$ for $x=y=1$ gives $\alpha=\beta=\frac12$. The module `torch_vn_algebra.jordan` implements these maps on operators of a `TypeIAlgebra`, lazily
and batched over all sectors:

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinFactor
from torch_vn_algebra import jordan as jd

torch.manual_seed(0)
dt = torch.complex128
jp, P = jd.jordan_product, jd.quadratic_representation

def close(A, B, tol=1e-10):
    A = A.matrix if hasattr(A, 'matrix') else A
    B = B.matrix if hasattr(B, 'matrix') else B
    return (A - B).abs().max().item() < tol

alg = TypeIAlgebra([2, 3], [2, 3], precision='double')           # M_2 (+) M_3
def sa(B=4):
    return alg.operator_from_eigenvalues(lambda d: torch.randn(B, d), batch_size=B, force_self_adjoint=True)
x, y, z = sa(), sa(), sa()
x2 = jp(x, x)
assert close(jp(jp(x2, y), x), jp(x2, jp(y, x)))                  # Jordan identity (Lemma 13.1)
assert close(jd.jordan_power(x, 5), x @ x @ x @ x @ x)            # power associativity
assert close(P(y, x), 2 * jp(y, jp(y, x)) - jp(jp(y, y), x))      # P(y) = 2L(y)^2 - L(y^2)
assert close(P(P(y, x), z), P(y, P(x, P(y, z))))                  # fundamental formula
assert close(P(y, x).trace, (jp(y, y) @ x).trace)                 # Tr P(y)x = Tr y^2 x
pos = alg.operator_from_eigenvalues(lambda d: torch.rand(4, d), batch_size=4, force_positive=True)
assert all((e > -1e-12).all() for e in P(y, pos).eigenvalues())  # Proposition 13.2(3)

M2 = TypeIAlgebra([2], [2], precision='double')
p = M2.operator(torch.tensor([[[[1.0, 0], [0, 0]]]]), is_self_adjoint=True, is_positive=True)
q = M2.operator(torch.tensor([[[[0.5, 0.5], [0.5, 0.5]]]]), is_self_adjoint=True, is_positive=True)
print(jp(p, q).eigenvalues()[0])                                  # (1 -+ sqrt 2)/4: not positive
assert abs(torch.linalg.det(jp(p, q).matrix[0, 0]).real.item() + 1 / 16) < 1e-14
```

> **Physical meaning** (physical interpretation, not a theorem of this course). $x\circ y$ is the symmetrised product
> of observables (the unique such combination, as shown above), and it is again an observable. The quadratic
> representation is the **Lüders rule**: measuring a yes/no question $p$ and keeping the outcome "yes" changes a
> density $\rho$ into $P(p)\rho/\operatorname{Tr}P(p)\rho=p\rho p/\operatorname{Tr}(p\rho)$. Proposition 13.2(4)
> says that $\circ$ itself has no such operational meaning for positive elements.

## 13.2 Spin factors and the Lorentz cone

**Definition.** For $m\ge1$, the **spin factor** $V_m=\mathbb R\oplus\mathbb R^m$ has elements $x=(a,u)$, the
product $(a,u)\circ(b,v)=(ab+u\cdot v,\;av+bu)$ and the unit $1=(1,0)$. Its **trace** is $\operatorname{tr}x=2a$.

The product is bilinear and commutative; the Jordan identity follows from a matrix representation.

**Lemma 13.3 (anticommuting symmetries).** Let $p=\lfloor m/2\rfloor$ and $D=2^p$ ($D=2$ if $m=1$). On
$(\mathbb C^2)^{\otimes p}$ put $s_{2j-1}=Z^{\otimes(j-1)}\otimes X\otimes1$, $s_{2j}=Z^{\otimes(j-1)}\otimes Y\otimes1$
for $j\le p$, and $s_m=Z^{\otimes p}$ if $m$ is odd ($s_1=Z$ if $m=1$). Then $s_i=s_i^*$ and
$s_is_j+s_js_i=2\delta_{ij}1$. The map $\phi(a,u)=a1+\sum_iu_is_i$ is injective and satisfies
$\phi(x\circ y)=\phi(x)\circ\phi(y)$; hence $V_m$ is a special Jordan algebra. For $m\ge2$ the $s_i$ generate
$M_D(\mathbb C)$ as a $*$-algebra.

*Proof.* Pauli matrices square to $1$ and pairwise anticommute; so does $Z$ with $X$ and $Y$. Two different strings
$s_i,s_j$ differ at exactly one site in an anticommuting way: for $s_{2j-1},s_{2j}$ at site $j$ ($X$ against
$Y$); for strings with indices $j<j'$, at site $j$ ($X$ or $Y$ against $Z$), while at the other sites they carry equal
factors or a factor $1$; the same holds for $s_m=Z^{\otimes p}$. Hence $s_is_j=-s_js_i$ for $i\ne j$, and
$s_i^2=1$. So $(u\cdot s)\circ(v\cdot s)=\sum_{ij}u_iv_j\,s_i\circ s_j=(u\cdot v)1$, which gives
$\phi(x)\circ\phi(y)=\phi(x\circ y)$. Each $s_i$ is traceless (a tensor product with a traceless factor) and
$\operatorname{Tr}(s_is_j)=\operatorname{Tr}(s_i\circ s_j)=D\delta_{ij}$, so $1,s_1,\dots,s_m$ are orthogonal for
$\operatorname{Tr}(ab)$ and $\phi$ is injective. The Jordan identity holds in $\operatorname{Herm}_D$ (Lemma 13.1) and is
transported back by $\phi$. Finally $s_{2j-1}s_{2j}=1\otimes\cdots\otimes XY\otimes\cdots\otimes1=iZ_j$ ($Z$ at site
$j$), so the $*$-algebra generated contains $Z_j$, then $X_j=Z_1\cdots Z_{j-1}s_{2j-1}$ and $Y_j$, hence all
tensor products of Pauli matrices, a basis of $M_D(\mathbb C)$. $\square$

Up to signs, $s_{2j-1}$ and $s_{2j}$ are the Majorana operators $c_j+c_j^*$ and $i(c_j-c_j^*)$ of the Jordan–Wigner
fermions $c_j$ of Section 9.4, with site $j-1$ of Chapter 9 as the $j$-th tensor factor: there $s^a=\frac12\sigma_a$,
$n=s^z+\frac12$ and $s^\pm=s^x\pm is^y$, so $(-1)^n=1-2n=-Z$ and $s^-=\frac12(X-iY)$, hence
$c_j=(-1)^{j-1}Z^{\otimes(j-1)}\otimes s^-\otimes1$. As $s^-+s^+=X$ and $i(s^--s^+)=i(-iY)=Y$,
$$
s_{2j-1}=(-1)^{j-1}(c_j+c_j^*),\qquad s_{2j}=(-1)^{j-1}\,i(c_j-c_j^*).
$$

**Theorem 13.4 (spectral theory of $V_m$).** Let $x=(a,u)\in V_m$ and $\lambda_\pm=a\pm|u|$.

1. The idempotents of $V_m$ are $0$, $1$ and $c_w=\frac12(1,w)$ with $|w|=1$; and $c_w\circ c_{w'}=0$ iff $w'=-w$.
2. If $u\ne0$, $x=\lambda_+e_++\lambda_-e_-$ with $e_\pm=c_{\pm u/|u|}$; if $u=0$, $x=a1$. If $x=\mu_1c_w+\mu_2c_{-w}$,
   then $\{\mu_1,\mu_2\}=\{\lambda_+,\lambda_-\}$. We call $\lambda_\pm$ the **eigenvalues** of $x$.
3. For $f:\{\lambda_+,\lambda_-\}\to\mathbb R$ put $f(x)=f(\lambda_+)e_++f(\lambda_-)e_-$ ($f(a)1$ if $u=0$). Then
   $x^{\circ k}=\lambda_+^ke_++\lambda_-^ke_-$, $\phi(f(x))=f(\phi(x))$, the eigenvalues of $\phi(x)$ are $\lambda_\pm$,
   each with multiplicity $D/2$, and $\operatorname{Tr}\phi(x)=\frac D2\operatorname{tr}x$. Also
   $\operatorname{tr}x=\lambda_++\lambda_-$ and $\det x:=\lambda_+\lambda_-=a^2-|u|^2$.
4. The cone of squares $Q=\{x^2:x\in V_m\}$ is the **Lorentz cone** $\{(a,u):a\ge|u|\}=\{x:\lambda_-\ge0\}=\{x:\phi(x)\ge0\}$.
   It is **self-dual**: $y\in Q$ iff $\operatorname{tr}(x\circ y)\ge0$ for all $x\in Q$.
5. $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$; hence $P(y)Q\subseteq Q$ for every $y\in V_m$.

*Proof.* (1) $(a,u)^2=(a^2+|u|^2,2au)$ equals $(a,u)$ iff $u=0$, $a\in\{0,1\}$, or $a=\frac12$, $|u|=\frac12$.
$c_w\circ c_{w'}=\frac14(1+w\cdot w',\,w+w')$ vanishes iff $w'=-w$. (2)
$\lambda_+e_++\lambda_-e_-=(a,0)+|u|(0,u/|u|)=x$, and $\mu_1c_w+\mu_2c_{-w}=(\frac{\mu_1+\mu_2}2,\frac{\mu_1-\mu_2}2w)$, so
$\mu_1+\mu_2=2a$ and $|\mu_1-\mu_2|=2|u|$. (3) By (1), $(\alpha e_++\beta e_-)\circ(\gamma e_++\delta e_-)=\alpha\gamma
e_++\beta\delta e_-$, which gives the powers by induction. With $w=u/|u|$, $(w\cdot s)^2=1$, so
$\phi(e_\pm)=\frac12(1\pm w\cdot s)$ are orthogonal projections with sum $1$ and product $0$, each of trace $D/2$; thus
$\phi(x)=\lambda_+\phi(e_+)+\lambda_-\phi(e_-)$ is a spectral decomposition, and
$f(\phi(x))=f(\lambda_+)\phi(e_+)+f(\lambda_-)\phi(e_-)=\phi(f(x))$ (Theorem 1.5). (4) $x^2$ has eigenvalues
$\lambda_\pm^2\ge0$; conversely, if $\lambda_-\ge0$ then $x=(\sqrt x)^2$. $\lambda_-\ge0$ iff $a\ge|u|$, iff
$\phi(x)\ge0$ by (3). For self-duality, $\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$. If $a\ge|u|$ and $b\ge|v|$, then
$ab\ge|u||v|\ge-u\cdot v$. If $b<|v|$, then $x=(|v|,-v)\in Q$ gives $2(|v|b-|v|^2)<0$ (for $v=0$, $b<0$, take $x=1$).
(5) $P(y)x=2y\circ(y\circ x)-y^2\circ x$ and $\phi$ preserves $\circ$, so $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ by
Proposition 13.2(1); apply (4) and Proposition 13.2(3). $\square$

**Proposition 13.5.** (a) For $m=3$, $\phi(a,u)=a1+u\cdot\sigma$ is a Jordan isomorphism
$V_3\cong\operatorname{Herm}_2(\mathbb C)$. (b) The Jordan automorphisms of $V_m$ are the maps $(a,u)\mapsto(a,Ou)$,
$O\in O(m)$. (c) For $m\ge2$, $\phi(V_m)$ has dimension $m+1$, and the self-adjoint part of the $*$-algebra it
generates has dimension $D^2=4^{\lfloor m/2\rfloor}$; they are equal only for $m=3$.

*Proof.* (a) For $m=3$, $(s_1,s_2,s_3)=(X,Y,Z)$, and $\phi$ is injective between spaces of dimension $4$. (b) The maps
$(a,u)\mapsto(a,Ou)$ preserve $u\cdot v$, hence $\circ$. Conversely let $\alpha$ be a Jordan automorphism. Then
$\alpha(1)$ is a unit of $V_m$, so $\alpha(1)=1$. By the formula for $x^2$, $x^2=1$ iff $x=\pm1$ or $x=(0,w)$ with
$|w|=1$. For $|w|=1$, $\alpha(0,w)^2=1$ and $\alpha(0,w)\ne\pm1=\alpha(\pm1)$, so $\alpha(0,w)=(0,w')$ with $|w'|=1$.
Writing $\alpha(0,u)=(\ell(u),Ou)$ with linear $\ell,O$, this says $\ell=0$ and $|Ou|=|u|$, i.e. $O\in O(m)$.
(c) Lemma 13.3; $4^{\lfloor m/2\rfloor}>m+1$ for $m=2$ and $m\ge4$. $\square$

```python
V = SpinFactor(5)
xs, ys = V.random(1000), V.random(1000)                          # (1000, 6): batches of (a, u)
x2 = V.square(xs)
assert (V.product(V.product(x2, ys), xs) - V.product(x2, V.product(ys, xs))).abs().max() < 1e-12
lam, e = V.spectral_decomposition(xs)                            # (lambda-, lambda+), (e-, e+)
assert (lam[:, 1, None] * e[:, 1] + lam[:, 0, None] * e[:, 0] - xs).abs().max() < 1e-12
Xm, Ym, Dm = V.to_matrix(xs), V.to_matrix(ys), V.representation_dim  # phi into Herm_4(C)
assert (V.to_matrix(V.product(xs, ys)) - 0.5 * (Xm @ Ym + Ym @ Xm)).abs().max() < 1e-12
assert (torch.linalg.eigvalsh(Xm) - lam.repeat_interleave(Dm // 2, dim=-1)).abs().max() < 1e-12
assert (V.to_matrix(V.exp(xs)) - torch.linalg.matrix_exp(Xm)).abs().max() < 1e-9       # phi(f(x)) = f(phi(x))
cone = V.random(1000, cone=True)                                 # squares
assert V.in_cone(cone).all() and (torch.linalg.eigvalsh(V.to_matrix(cone)) > -1e-12).all()
assert (V.square(V.sqrt(cone)) - cone).abs().max() < 1e-10
O, _ = torch.linalg.qr(torch.randn(5, 5, dtype=torch.float64))  # Proposition 13.5(b)
rot = lambda w: torch.cat([w[..., :1], w[..., 1:] @ O.T], dim=-1)
assert (rot(V.product(xs, ys)) - V.product(rot(xs), rot(ys))).abs().max() < 1e-12

for m in range(1, 7):                                            # Proposition 13.5(c)
    r = jd.compare_with_star_algebra(jd.anticommuting_symmetries(m))
    print(f"m = {m}: dim V_m = {r['jordan_dim']}, dim of the self-adjoint part of the *-algebra = {r['star_dim']}")
    assert r['jordan_dim'] == m + 1 and r['equal'] == (m in (1, 3))

pq = 3                                                           # Majorana operators of Section 9.4, p = 3 sites
Xp, Yp, Zp = jd.anticommuting_symmetries(3)
sminus, I2 = 0.5 * (Xp - 1j * Yp), torch.eye(2, dtype=dt)
parity = I2 - 2 * (0.5 * Zp + 0.5 * I2)                          # (-1)^n with n = s^z + 1/2
def kron_all(fs):
    out = torch.ones(1, 1, dtype=dt)
    for f_ in fs:
        out = torch.kron(out, f_)
    return out
s6 = jd.anticommuting_symmetries(2 * pq)
for j in range(1, pq + 1):
    c = kron_all([parity] * (j - 1) + [sminus] + [I2] * (pq - j))
    sign = (-1) ** (j - 1)
    assert close(s6[2 * j - 2], sign * (c + c.mH)) and close(s6[2 * j - 1], sign * 1j * (c - c.mH))
```

For $m=1$, $s_1=Z$ generates the diagonal matrices, of dimension $2=m+1$, and $V_1\cong\mathbb R\oplus\mathbb R$ is not
simple. The printed dimensions are $3,4$ for $m=2$; $5,16$ for $m=4$; $7,64$ for $m=6$.

**States of $V_m$.** A **state** of $V_m$ is a linear functional $\omega$ with $\omega(Q)\ge0$ and $\omega(1)=1$. The
states are exactly the maps $x\mapsto\operatorname{tr}(\rho\circ x)$ with $\rho=\frac12(1,r)$, $|r|\le1$; so the state
space is a ball of dimension $m$. *Proof.* $\langle x,y\rangle=\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$ is an inner
product on $V_m$, so every linear $\omega$ is $\operatorname{tr}(\rho\circ\cdot)$ for exactly one $\rho=(b,v)$ (finite-dimensional
Riesz representation: in an orthonormal basis, $\rho=\sum_k\omega(e_k)e_k$). By the self-duality in Theorem 13.4(4),
$\omega(Q)\ge0$ iff $\rho\in Q$, i.e. $b\ge|v|$; and $\omega(1)=\operatorname{tr}\rho=2b=1$. So $\rho=\frac12(1,r)$ with
$r=2v$, $|r|\le1$; conversely every such $\rho$ gives a state. $\square$

> **Physical meaning** (physical interpretation, not a theorem of this course). For $m=3$ the state space of $V_m$
> is the Bloch ball of a qubit (Proposition 13.5(a)). The Lorentz cone is the forward
> light cone of special relativity in $1+m$ dimensions, and Proposition 13.5(b) says that the Jordan automorphisms
> are the orthogonal maps of space that fix the time axis.

The cones $Q$ are the constraint sets of **second-order cone programming** (SOCP): linear optimisation over
products of Lorentz cones. The quantities of Theorem 13.4 (eigenvalues, $\det$, and $P(y)$) are those in which
interior-point methods for SOCP are usually formulated; these algorithms are not treated in this course (further
reading: F. Alizadeh, D. Goldfarb, *Second-order cone programming*, Math. Program. 95 (2003) 3–51).

## 13.3 Jordan subalgebras that are not *-subalgebras

A Jordan subalgebra $J\subseteq\operatorname{Herm}_n(\mathbb C)$ generates a complex $*$-algebra $\mathcal A$ with
$J\subseteq\mathcal A_{sa}$. Proposition 13.5(c) gave examples with $J\ne\mathcal A_{sa}$; here are two more.

**Proposition 13.6.** (a) $\operatorname{Herm}_n(\mathbb R)$ is a Jordan subalgebra of $\operatorname{Herm}_n(\mathbb C)$ of
dimension $n(n+1)/2$; the $*$-algebra it generates is $M_n(\mathbb C)$, with self-adjoint part of dimension $n^2$.
(b) For $x,y\in\operatorname{Herm}_2(\mathbb C)$ there is a unitary $U$ such that $UxU^*$ and $UyU^*$ are real
symmetric. The Jordan algebra generated by $1,x,y$ is $\operatorname{span}\{1,x,y\}$; if $x,y,1$ are linearly
independent, it is $U^*\operatorname{Herm}_2(\mathbb R)U\cong V_2$, of dimension $3$, while the $*$-algebra generated by
$1,x,y$ is $M_2(\mathbb C)$, of dimension $4$.

*Proof.* (a) Lemma 13.1. For $i\ne j$, $E_{ij}=E_{ii}(E_{ij}+E_{ji})$ is a product of real symmetric matrices, so the
generated algebra contains all matrix units. (b) Write $x=a1+u\cdot\sigma$, $y=b1+v\cdot\sigma$ and pick a unit vector
$n\perp u,v$. The matrices $n\cdot\sigma$ and $Y$ are Hermitian, traceless, and square to $1$, so both have the
eigenvalues $1,-1$, and $U(n\cdot\sigma)U^*=Y$ for some unitary $U$ (Theorem 1.3). Conjugation by $U$ preserves
$\operatorname{Tr}$ and $\circ$, so it maps $u\cdot\sigma$ to $(Ru)\cdot\sigma$ with $R$ linear and
$Ru\cdot Rv=\frac12\operatorname{Tr}\big((Ru\cdot\sigma)(Rv\cdot\sigma)\big)=u\cdot v$: $R\in O(3)$, and $Rn=e_2$.
Hence $Ru,Rv\perp e_2$, and $UxU^*=a1+(Ru)_1X+(Ru)_3Z$ is real, since $X,Z$ are; likewise $UyU^*$. Next,
$x\circ y=(ab+u\cdot v)1+av\cdot\sigma+bu\cdot\sigma\in\operatorname{span}\{1,x,y\}$, and similarly $x^2,y^2$, so the span is
closed. If it has dimension $3$, its conjugate by $U$ is all of $\operatorname{Herm}_2(\mathbb R)=\operatorname{span}\{1,X,Z\}$,
which is $\cong V_2$ via $(a,u)\mapsto a1+u_1X+u_2Z$ (the proof of Lemma 13.3 uses only that the $s_i$ are
anticommuting Hermitian symmetries, and $X,Z$ are such), and generates $M_2(\mathbb C)$ by (a). $\square$

```python
for n in range(2, 6):                                            # Proposition 13.6(a), a random real pair
    G = torch.randn(2, n, n, dtype=torch.float64)
    r = jd.compare_with_star_algebra(G + G.mT)
    assert (r['jordan_dim'], r['star_dim']) == (n * (n + 1) // 2, n * n)

H = torch.randn(3, 2, 2, dtype=dt)
H = H + H.mH                                                     # three random Hermitian 2 x 2 matrices
r2, r3 = jd.compare_with_star_algebra(H[:2]), jd.compare_with_star_algebra(H)
print("two generators:", r2['jordan_dim'], r2['star_dim'], "  three:", r3['jordan_dim'], r3['star_dim'])
assert (r2['jordan_dim'], r2['star_dim'], r3['jordan_dim']) == (3, 4, 4)

V3 = SpinFactor(3)
u, v = V3.from_matrix(H[0])[1:], V3.from_matrix(H[1])[1:]        # Bloch vectors
nrm = torch.linalg.cross(u, v)
nrm = nrm / nrm.norm()
_, W1 = torch.linalg.eigh(V3.to_matrix(torch.cat([torch.zeros(1, dtype=torch.float64), nrm])))
_, W2 = torch.linalg.eigh(jd.anticommuting_symmetries(3)[1])     # Y
Um = W2 @ W1.mH                                                  # U (n . sigma) U* = Y
assert all((Um @ A @ Um.mH).imag.abs().max() < 1e-12 for A in H[:2])
```

Three generic Hermitian $2\times2$ matrices generate all of $\operatorname{Herm}_2(\mathbb C)$: dimension $4$ in both
columns. In a block $M_2$ of $\mathcal M$, therefore, two observables never reveal the complex structure; a third one
can.

## 13.4 The classification

A finite-dimensional Jordan algebra is **formally real** if $\sum_ix_i^2=0$ implies that all $x_i=0$, and
**Euclidean** if it has a unit $1$ and an inner product with $\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$. The two
notions coincide (Theorem E.5(4); the direction "Euclidean $\Rightarrow$ formally real" is one line:
$0=\sum_i\langle x_i^2,1\rangle=\sum_i\langle x_i,x_i\rangle$). $M_{sa}$ with $\langle x,y\rangle=\operatorname{Tr}(xy)$ is
Euclidean: $\operatorname{Tr}(x^2)=\sum_{ij}|x_{ij}|^2$, and by cyclicity
$\operatorname{Tr}((x\circ y)z)=\frac12\operatorname{Tr}(xyz+yxz)=\operatorname{Tr}(y(x\circ z))$. So is $V_m$ with
$\langle x,y\rangle=\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$: for $z=(c,w)$,
$\operatorname{tr}((x\circ y)\circ z)=2(abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w)$ is symmetric in $x,y,z$.

**Theorem 13.7 (Jordan–von Neumann–Wigner).** Every finite-dimensional formally real Jordan
algebra is a direct sum of simple ideals, and every simple one is isomorphic to exactly one of: $\mathbb R$;
$\operatorname{Herm}_r(\mathbb K)$ for $\mathbb K=\mathbb R,\mathbb C,\mathbb H$ and $r\ge3$; the spin factors
$V_m$, $m\ge2$; the **Albert algebra** $\operatorname{Herm}_3(\mathbb O)$ of Hermitian $3\times3$ matrices over the
octonions (dimension $27$). Moreover $\operatorname{Herm}_1(\mathbb K)\cong\mathbb R$ and
$\operatorname{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}$ for $\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O$.

*Proof.* Appendix E proves it completely: the decomposition is Theorem E.9, the list and "exactly one" are Theorem
E.27(1),(2), and the isomorphisms in rank one and two are Theorem E.24(3). $\square$

(P. Jordan, J. von Neumann, E. Wigner, *On an algebraic generalization of the quantum mechanical formalism*, Ann. of
Math. 35 (1934) 29–64; the proof in Appendix E follows [FK94, Chapters II–V].)

All of these except the Albert algebra are special: $\mathbb R=\operatorname{Herm}_1(\mathbb R)$ and
$\operatorname{Herm}_r(\mathbb K)$, $\mathbb K=\mathbb R,\mathbb C$, by Lemma 13.1, $\operatorname{Herm}_r(\mathbb H)$
because writing quaternions as $2\times2$ complex matrices embeds it into $\operatorname{Herm}_{2r}(\mathbb C)$
(Theorem E.24(1)), and $V_m$ by Lemma 13.3. The Albert algebra is a Jordan algebra (Theorem E.24(2)) and is
**exceptional**: it is not isomorphic to a Jordan subalgebra of any associative real or complex algebra (Theorem E.26,
an argument of the kind used by A. A. Albert, *On a certain algebra of quantum mechanics*, Ann. of Math. 35 (1934)
65–73: a Jordan embedding into an associative algebra would make the octonion product associative, which it is not,
Lemma E.15(2)). For $r\ge4$, $\operatorname{Herm}_r(\mathbb O)$ does not satisfy the Jordan identity (Corollary E.25),
and a simple algebra with a Jordan frame of $r\ge4$ elements is $\operatorname{Herm}_r(D)$ with an associative
composition algebra $D\ne\mathbb O$ (Theorem E.22), which is why the octonions occur only for $r\le3$. So $\operatorname{Herm}_3(\mathbb O)$
cannot be represented faithfully by real or complex matrices (or inside any associative algebra) with the product
$\frac12(xy+yx)$; a computation in it has to use the octonion product (or its structure constants) instead of a
matrix product.

Let $\Omega$ be a nonempty open convex cone in a Euclidean space $V$ (a finite-dimensional real inner product space
$\ne\{0\}$). Its **dual cone** is $\Omega^*=\{y:\langle x,y\rangle>0\text{ for all }x\in\bar\Omega\setminus\{0\}\}$; $\Omega$ is
**self-dual** if $\Omega^*=\Omega$, **homogeneous** if the linear maps $g$ with $g\Omega=\Omega$ act transitively on
$\Omega$, and **symmetric** if it is both.

**Theorem 13.8 (Koecher–Vinberg).** A nonempty open convex cone $\Omega$ in a Euclidean space $V$ is symmetric if and
only if there is a product $\circ$ on $V$ such that $(V,\circ,\langle\cdot,\cdot\rangle)$ is a Euclidean Jordan algebra
and $\Omega$ is the interior of its cone of squares.

*Proof.* Theorem F.12 ("if": Theorems F.5(2) and F.6(2)). $\square$

(M. Koecher, 1957; E. B. Vinberg, 1960; the proof in Appendix F follows [FK94, Chapters I–III].) Theorem 13.4(4)
proves the closed form of self-duality for the Lorentz cone ($y\in Q$ iff $\operatorname{tr}(x\circ y)\ge0$ for all
$x\in Q$), and Exercise 1 its homogeneity. The cone of squares of a direct sum of Euclidean Jordan algebras (with the
orthogonal sum of the inner products) is the product of their cones of squares, since $(x_1,\dots,x_k)^2=(x_1^2,\dots,x_k^2)$.
So, by the "if" direction, the interiors of the cones of semidefinite programming (products of
$\operatorname{Herm}_n(\mathbb R)^+$ or $\operatorname{Herm}_n(\mathbb C)^+$, the cones of squares by Theorem 1.6; for real $x\ge0$ the square root is real, since its complex conjugate is
also a positive square root of $x$), of SOCP
(products of Lorentz cones) and of linear programming ($\mathbb R_+^n$, the cone of squares of
$\mathbb R^n=\mathbb R\oplus\dots\oplus\mathbb R$) are symmetric cones.

[FK94]: J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, Oxford University Press, 1994.

## 13.5 Jordan maps versus *-maps

A **Jordan $*$-homomorphism** of von Neumann algebras is a complex-linear map $\Phi$ with $\Phi(x^*)=\Phi(x)^*$ and
$\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$; on self-adjoint parts it is a Jordan homomorphism, and every Jordan
homomorphism of self-adjoint parts extends to one by $\Phi(x+iy)=\Phi(x)+i\Phi(y)$ ($x,y$ self-adjoint). (This
$\Phi$ is complex linear and satisfies $\Phi(x-iy)=\Phi(x+iy)^*$; and since $\circ$ is complex bilinear and maps pairs
of self-adjoint elements to self-adjoint elements, $(x+iy)\circ(x'+iy')=(x\circ x'-y\circ y')+i(x\circ y'+y\circ x')$ is
mapped to $\Phi(x+iy)\circ\Phi(x'+iy')$.)

**Proposition 13.9 (the transpose).** Let $T(A)=\bigoplus_cA_c^T$ on $\mathcal M$. (a) $T$ is a unital Jordan
$*$-automorphism with $T(AB)=T(B)T(A)$; it preserves spectra, so $T$ and $T^{-1}=T$ are positive. (b) On $M_k$,
$k\ge2$, the Choi matrix $\sum_{ij}E_{ij}\otimes T(E_{ij})$ is the swap $F(e_a\otimes e_b)=e_b\otimes e_a$, with
eigenvalues $+1$ (multiplicity $k(k+1)/2$) and $-1$ (multiplicity $k(k-1)/2$). So $T$ is not completely positive.

*Proof.* (a) $(A^T)^*=(A^*)^T$, $(AB)^T=B^TA^T$, hence $(A\circ B)^T=B^T\circ A^T=A^T\circ B^T$;
$\det(A^T-\lambda)=\det(A-\lambda)$. (b) $\sum_{ij}E_{ij}e_a\otimes E_{ji}e_b=e_b\otimes e_a$. $F^2=1$, and $F$ is $+1$
on the symmetric and $-1$ on the antisymmetric tensors, which span. By Theorem 6.2 (and Proposition 6.1) $T$ is not
CP. $\square$

So a Jordan automorphism need not be a $*$-automorphism. Up to this phenomenon, positivity determines the Jordan
structure:

**Theorem 13.10 (Kadison, finite dimensions).** Let $\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$ and
$\mathcal N=\bigoplus_dM_{k'_d}(\mathbb C)$, and let $\Phi:\mathcal M\to\mathcal N$ be a linear bijection with
$\Phi(1)=1$ such that $\Phi$ and $\Phi^{-1}$ are positive. Then $\Phi$ is a Jordan $*$-isomorphism. Conversely, every
Jordan $*$-isomorphism $\mathcal M\to\mathcal N$ is unital, and it and its inverse are positive.

*Proof.* Theorem F.15. (The idea: $\Phi$ maps the convex set $\{0\le a\le1\}$ onto the corresponding set of
$\mathcal N$, hence extreme points, which are the projections, onto projections; it preserves $p+q\le1$, which for
projections means $pq=0$; so it maps spectral decompositions to spectral decompositions and $\Phi(x^2)=\Phi(x)^2$.)
$\square$

R. V. Kadison proved the theorem for all unital C\*-algebras (*A generalized Schwarz inequality and algebraic
invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503); only the finite-dimensional case is proved and
used in this course. For the factors $M_n$, the Jordan automorphisms are known exactly:
every Jordan $*$-automorphism of $M_n(\mathbb C)$ is $A\mapsto UAU^*$ or $A\mapsto UA^TU^*$ with $U$ unitary
(Theorem F.17, for all $n$; the result is due to N. Jacobson, C. E. Rickart, Trans. Amer. Math. Soc. 69 (1950)
479–502, and I. N. Herstein, Trans. Amer. Math. Soc. 81 (1956) 331–341). The first maps are multiplicative and the
second antimultiplicative, so a Jordan $*$-automorphism of $M_n$ either preserves all commutators or reverses all of
them: $\Phi([x,y])=\pm[\Phi(x),\Phi(y)]$. On $\bigoplus_cM_{k_c}$ a Jordan $*$-isomorphism permutes blocks of equal
size and acts on each block in one of the two ways (Theorem F.18); it is unital and positive in both directions
(Theorem 13.10), and by Corollary F.19 it is completely positive iff it
transposes no block of size $\ge2$, iff it is a $*$-isomorphism. Here is a direct proof of the case $n=2$ in the
language of spin factors.

**Proposition 13.11.** Every Jordan automorphism of $\operatorname{Herm}_2(\mathbb C)$ is $A\mapsto UAU^*$ or
$A\mapsto UA^TU^*$ with $U\in U(2)$.

*Proof.* By Proposition 13.5(a,b), a Jordan automorphism acts as $a1+u\cdot\sigma\mapsto a1+(Ou)\cdot\sigma$ with
$O\in O(3)$. Write $R_U\in O(3)$ for the matrix of $A\mapsto UAU^*$ in this sense; $R_{UV}=R_UR_V$. (i) $\det R_U=1$:
$U=W\operatorname{diag}(e^{i\theta_1},e^{i\theta_2})W^*$ is joined to $1$ by the path
$U_t=W\operatorname{diag}(e^{it\theta_1},e^{it\theta_2})W^*$, and $\det R_{U_t}\in\{\pm1\}$ depends continuously on $t$.
(ii) Every $R\in SO(3)$ is some $R_U$: as in Proposition 13.6, there is $U_1$ with
$U_1ZU_1^*=(Re_3)\cdot\sigma$, i.e. $R_{U_1}e_3=Re_3$. Then $R_{U_1}^{-1}R\in SO(3)$ fixes $e_3$, so it is a rotation
by some angle $\beta$ about $e_3$. A direct computation gives, for $U_2=\operatorname{diag}(e^{-i\beta/2},e^{i\beta/2})$,
$U_2XU_2^*=\cos\beta\,X+\sin\beta\,Y$ and $U_2YU_2^*=-\sin\beta\,X+\cos\beta\,Y$, so $R=R_{U_1U_2}$. (iii) The
transpose fixes $X,Z$ and maps $Y\mapsto-Y$: $O_T=\operatorname{diag}(1,-1,1)$, $\det O_T=-1$. Hence $O=R_U$ if
$\det O=1$, and $OO_T=R_U$, i.e. the automorphism is $A\mapsto UA^TU^*$, if $\det O=-1$. $\square$

**Composite systems.** The Jordan structure of a factor cannot tell $T$ apart from the identity. For two systems this
ambiguity becomes visible:

**Proposition 13.12.** Let $x,x'\in\operatorname{Herm}_n(\mathbb C)$, $y,y'\in\operatorname{Herm}_m(\mathbb C)$. In
$\operatorname{Herm}_{nm}(\mathbb C)$, acting on $\mathbb C^n\otimes\mathbb C^m$,
$$
(x\otimes y)\circ(x'\otimes y')=(x\circ x')\otimes(y\circ y')+\tfrac14[x,x']\otimes[y,y'].
$$
Let $\theta=\mathrm{id}\otimes T$ (the **partial transpose**) and $A\circ_\theta B=\theta(\theta(A)\circ\theta(B))$. Then
$\circ_\theta$ is a Jordan product on $\operatorname{Herm}_{nm}(\mathbb C)$, isomorphic to $\circ$ via $\theta$, that
agrees with $\circ$ on all pairs of local observables $x\otimes1$, $1\otimes y$, and
$(x\otimes y)\circ_\theta(x'\otimes y')=(x\circ x')\otimes(y\circ y')-\frac14[x,x']\otimes[y,y']$. Its cone of squares is
$\theta(\operatorname{Herm}_{nm}^+)$. For $n=m=2$ and $\Omega=e_1\otimes e_1+e_2\otimes e_2$, the state
$\omega=\frac12|\Omega\rangle\langle\Omega|\ge0$ is not in this cone: $\theta(\omega)=\frac12F$ has the eigenvalue
$-\frac12$.

*Proof.* With $j=x\circ x'$, $c=\frac12[x,x']$, $k=y\circ y'$, $d=\frac12[y,y']$: $xx'=j+c$, $x'x=j-c$, and likewise for
$y$, so $\frac12(xx'\otimes yy'+x'x\otimes y'y)=\frac12\big((j+c)\otimes(k+d)+(j-c)\otimes(k-d)\big)=j\otimes k+c\otimes d$.
$\theta$ is a real-linear involution of $\operatorname{Herm}_{nm}$: it maps $x\otimes y$ to $x\otimes y^T$, and the
$x\otimes y$ span $\operatorname{Herm}_{nm}$ over $\mathbb R$ (products $h_\alpha\otimes g_\beta$ of real bases of
$\operatorname{Herm}_n,\operatorname{Herm}_m$ form a complex basis of $M_{nm}$, and a Hermitian matrix has real
coefficients in it). So $\circ_\theta$ is a Jordan product and $\theta$ an isomorphism. Since $(y^T\circ y'^T)^T=y\circ y'$ and
$[y^T,y'^T]^T=-[y,y']$, applying $\theta$ to the formula for $(x\otimes y^T)\circ(x'\otimes y'^T)$ gives the formula for
$\circ_\theta$. On local observables the commutator term vanishes ($[1,1]=0$, $[x,1]=0$). The squares are
$A\circ_\theta A=\theta(\theta(A)^2)$, and $\omega$ lies in $\theta(\operatorname{Herm}_{nm}^+)$ iff $\theta(\omega)\ge0$,
since $\theta^2=\mathrm{id}$. Finally $\theta(\omega)=\frac12\sum_{ij}E_{ij}\otimes E_{ji}=\frac12F$
(Proposition 13.9(b); this is the computation in the proof of Proposition 6.1). $\square$

The two products are equally good Jordan algebras, they agree on everything that can be measured on one subsystem
at a time, and they disagree on which joint operators are positive, that is, on what the states of the
composite system are. What distinguishes them is the sign of the commutator term, and this sign is fixed by the
$*$-algebra $M_n\otimes M_m$ and not by the Jordan algebras of the factors. ($\circ_\theta$ is the Jordan product of the
associative product $A\star B=\theta(\theta(A)\theta(B))$, for which $(x\otimes y)\star(x'\otimes y')=xx'\otimes y'y$:
the algebra $M_n\otimes M_m^{\rm op}$.) The cone $\theta(\operatorname{Herm}_{nm}^+)$ is
what gives the **PPT criterion** (A. Peres, Phys. Rev. Lett. 77 (1996) 1413). Call a state $\rho$ of
$M_n\otimes M_m$ **separable** if $\rho=\sum_kp_ka_k\otimes b_k$ (a finite sum) with $p_k\ge0$, $a_k,b_k\ge0$, and
**entangled** otherwise. A separable state lies in both cones, since
$\theta(\rho)=\sum_kp_ka_k\otimes b_k^T$ with $b_k^T\ge0$ (Proposition 13.9(a)) and
$a\otimes b=(a^{1/2}\otimes b^{1/2})^2\ge0$. So $\omega$, which is not in $\theta(\operatorname{Herm}_4^+)$, is
entangled. (Further reading on composites of Jordan algebras, not used here: H. Hanche-Olsen, *On the structure and
tensor products of JC-algebras*, Canad. J. Math. 35 (1983) 1059–1074; H. Barnum, M. Graydon, A. Wilce,
*Composites and categories of Euclidean Jordan algebras*, Quantum 4 (2020) 359.)

**Dynamics.** One might hope that time evolution forces the associative product. It does not:

**Proposition 13.13.** Let $H\in\operatorname{Herm}_n(\mathbb C)$ and $\delta(z)=i[H,z]$. (a) $\delta$ maps
$\operatorname{Herm}_n$ to itself and is a derivation of $\circ$; $e^{t\delta}(z)=e^{itH}ze^{-itH}$ is a Jordan
automorphism. (b) $[L(x),L(y)]z=\frac14[[x,y],z]$ for $x,y,z\in\operatorname{Herm}_n$. (c) $\delta=4\sum_k[L(x_k),L(y_k)]$
for some $x_k,y_k\in\operatorname{Herm}_n$: Hamiltonian dynamics is generated by the Jordan product alone.

*Proof.* (a) $(i[H,z])^*=-i[z,H]=i[H,z]$; $[H,\cdot]$ is a derivation of the associative product, hence of $\circ$.
$z(t)=e^{itH}ze^{-itH}$ satisfies $z'=i[H,z]$ and $z(0)=z$ (Proposition 8.2(3) with $\mathcal M=M_n$), so
$z(t)=e^{t\delta}(z)$ by the uniqueness in Lemma 8.1; conjugation by the unitary $e^{itH}$ preserves products and
adjoints, hence $\circ$ and $\operatorname{Herm}_n$.
(b) $4x\circ(y\circ z)=xyz+xzy+yzx+zyx$; subtract the same with $x,y$ exchanged. (c) By (b) it suffices to write
$iH-\frac{\operatorname{Tr}H}n i1$ (which has the same commutators as $iH$) as a sum of $[x_k,y_k]$. For $i\ne j$ and
$a=E_{ij}+E_{ji}$, $b=i(E_{ij}-E_{ji})$: $[a,E_{ii}]=E_{ji}-E_{ij}$, $[b,E_{ii}]=-i(E_{ij}+E_{ji})$ and
$[a,b]=2i(E_{jj}-E_{ii})$; these span the traceless anti-Hermitian matrices over $\mathbb R$. $\square$

What the Jordan algebra does not contain is the **correspondence** $H\mapsto i[H,\cdot]$ between an observable and the
dynamics it generates. In $\operatorname{Herm}_n(\mathbb R)$, for instance, the derivations
$[L(x),L(y)]=\frac14[[x,y],\cdot]$ come from real antisymmetric $[x,y]$, which is not $i$ times a real symmetric
matrix (Exercise 4). The correspondence is not determined by $\circ$: on $\operatorname{Herm}_n(\mathbb C)$ both
$H\mapsto i[H,\cdot]$ and $H\mapsto-i[H,\cdot]$ send every observable to a derivation of $\circ$ (Proposition 13.13(a)),
and the transpose, a Jordan automorphism, exchanges them, since $T(i[H,T(z)])=-i[H^T,z]$. The second is $i$ times the
commutator of the opposite algebra ($a\cdot_{\rm op}b=ba$): $-i[H,z]=i(H\cdot_{\rm op}z-z\cdot_{\rm op}H)$ — the two signs
of Proposition 13.12. (Further reading on such "dynamical correspondences", not used here: E. M. Alfsen, F. W. Shultz, *Geometry of
State Spaces of Operator Algebras*, Birkhäuser, 2003.)

```python
x, y = sa(), sa()
T = jd.transpose_map
assert close(T(jp(x, y)), jp(T(x), T(y))) and close(T(x @ y), T(y) @ T(x))   # Proposition 13.9(a)
assert not close(T(x @ y), T(x) @ T(y))

k = 3
Mk = TypeIAlgebra([k], [k], precision='double')
units = torch.eye(k * k, dtype=dt).reshape(k * k, 1, k, k)                     # E_ij, batch index i*k + j
TE = T(Mk.operator(units)).matrix[:, 0]
choi = TE.reshape(k, k, k, k).permute(0, 2, 1, 3).reshape(k * k, k * k)       # block (i, j) = T(E_ij)
swap = torch.eye(k * k, dtype=dt).reshape(k, k, k, k).permute(0, 1, 3, 2).reshape(k * k, k * k)
w = torch.linalg.eigvalsh(choi)
assert torch.equal(choi, swap) and (w < 0).sum() == k * (k - 1) // 2          # Proposition 13.9(b)

sig = jd.anticommuting_symmetries(3)                                           # X, Y, Z
O_T = torch.stack([V3.from_matrix(s.T)[1:] for s in sig], dim=1)               # transpose on Bloch vectors
assert torch.allclose(O_T, torch.diag(torch.tensor([1.0, -1.0, 1.0], dtype=torch.float64)))

def ptranspose(R, n, m):                                                       # theta = id (x) T
    return R.reshape(n, m, n, m).permute(0, 3, 2, 1).reshape(n * m, n * m)
def herm(n):
    A = torch.randn(n, n, dtype=dt)
    return A + A.mH
jm = lambda A, B: 0.5 * (A @ B + B @ A)
cm = lambda A, B: A @ B - B @ A
a, a2, b, b2 = herm(2), herm(2), herm(3), herm(3)
local = torch.kron(jm(a, a2), jm(b, b2))
assert close(jm(torch.kron(a, b), torch.kron(a2, b2)), local + 0.25 * torch.kron(cm(a, a2), cm(b, b2)))
theta = lambda R: ptranspose(R, 2, 3)
assert close(theta(jm(theta(torch.kron(a, b)), theta(torch.kron(a2, b2)))),
             local - 0.25 * torch.kron(cm(a, a2), cm(b, b2)))                  # Proposition 13.12
Om = torch.tensor([1.0, 0, 0, 1.0], dtype=dt) / math.sqrt(2)
omega = torch.outer(Om, Om.conj())                                             # Bell state
print("spectrum of the partial transpose of the Bell state:", torch.linalg.eigvalsh(ptranspose(omega, 2, 2)))
assert abs(torch.linalg.eigvalsh(ptranspose(omega, 2, 2))[0] + 0.5) < 1e-12

Hm, z1, z2 = herm(3), herm(3), herm(3)
delta = lambda A: 1j * cm(Hm, A)
assert close(delta(jm(z1, z2)), jm(delta(z1), z2) + jm(z1, delta(z2)))         # Proposition 13.13(a)
assert close(jm(z1, jm(z2, Hm)) - jm(z2, jm(z1, Hm)), 0.25 * cm(cm(z1, z2), Hm))   # (b)
```

> **Physical meaning** (physical interpretation, not a theorem of this course). For a spinless particle the transpose (complex conjugation in a real basis) is **time
> reversal**. It is a symmetry of the observables of a single system, but it is antiunitary, and applied to one
> subsystem only it is not a physical operation: Proposition 13.12 shows that it turns the Bell state into an
> operator with a negative eigenvalue. Quantum mechanics uses complex $*$-algebras because composite systems
> need the commutator with its sign, not only the symmetrised product.

## 13.6 Back to trace inequalities

Proposition 11.9 (Exp. 3) states $\operatorname{Tr}|YXY|\le\operatorname{Tr}(Y|X|Y)$ for $X=X^*$ and $Y\ge0$ (here
$X,Y$ are operators of $\mathcal M$, not Pauli matrices). Its proof never uses $Y\ge0$: $YX_\pm Y=T^*T$ with
$T=X_\pm^{1/2}Y$ needs only $Y=Y^*$, and Lemma 11.8 is applied to $YX_+Y$ and $YX_-Y$. So it holds for all $Y=Y^*$,
and in Jordan language it reads
$$
\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|,
$$
which makes sense in every Euclidean Jordan algebra once a trace is fixed. We prove it in $V_m$ by an argument that uses three properties
only.

**Proposition 13.14 (Exp. 3 in a spin factor).** For $x,y\in V_m$,
$\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|$, with equality if $x\in Q$ or $-x\in Q$.

*Proof.* Let $x_\pm=f_\pm(x)$ with $f_+(t)=\max(t,0)$, $f_-(t)=\max(-t,0)$ (Theorem 13.4(3)): $x_\pm\in Q$, $x=x_+-x_-$ and
$|x|=x_++x_-$. Put $A=P(y)x_+$ and $B=P(y)x_-$; they lie in $Q$ by Theorem 13.4(5), and $P(y)x=A-B$. Let
$P(y)x=\mu_1c_1+\mu_2c_2$ be a spectral decomposition (Theorem 13.4(2); $c_1=c_w$, $c_2=c_{-w}$, so $c_1+c_2=1$ and
$c_1\circ c_2=0$), $c$ the sum of the $c_i$ with $\mu_i>0$, and
$s=2c-1=\sum_i\operatorname{sgn}'(\mu_i)c_i$ with $\operatorname{sgn}'(t)=1$ for $t>0$ and $-1$ otherwise. Then
$s\circ P(y)x=\sum_i|\mu_i|c_i=|P(y)x|$ (by Theorem 13.4(2,3) this does not depend on the choice of $w$ when
$\mu_1=\mu_2$), and $c,1-c\in Q$ (idempotents are squares). By the easy half of self-duality (Theorem 13.4(4):
$\operatorname{tr}(a\circ b)\ge0$ for $a,b\in Q$), applied to $\operatorname{tr}(A\circ(1-c))$ and $\operatorname{tr}(B\circ c)$,
$$
\operatorname{tr}|P(y)x|=\operatorname{tr}(A\circ c)-\operatorname{tr}(A\circ(1-c))-\operatorname{tr}(B\circ c)+\operatorname{tr}(B\circ(1-c))
\le\operatorname{tr}A+\operatorname{tr}B=\operatorname{tr}P(y)|x|.
$$
If $x\in Q$, then $x_-=0$, $B=0$ and $|P(y)x|=A$, so both sides equal $\operatorname{tr}A$; likewise for $-x\in Q$.
$\square$

The proof uses only: (i) $P(y)Q\subseteq Q$ for every $y$; (ii) the spectral theorem: every $z$ is
$\sum_i\mu_ic_i$ with orthogonal idempotents $c_i$ summing to $1$, and $f(z)=\sum_if(\mu_i)c_i$ (so $x_\pm$ and the
idempotents are squares); (iii) $\operatorname{tr}$ is linear, $\operatorname{tr}(a\circ1)=\operatorname{tr}a$, and
$\operatorname{tr}(a\circ b)\ge0$ for $a,b\in Q$. All three hold in every Euclidean Jordan algebra $V$ with unit $e$
for the linear form $\tau(x)=\langle x,e\rangle$, for which $\tau(a\circ b)=\langle b,a\circ e\rangle=\langle a,b\rangle$:
(ii) is Theorem F.3, (i) is Theorem F.6(4), and $\langle a,b\rangle\ge0$ for $a,b\in Q$ is Theorem F.5(1). Corollary F.7
carries out the argument: $\tau(|P(y)x|)\le\tau(P(y)|x|)$ in every Euclidean Jordan algebra, with equality if $x\in Q$
or $-x\in Q$. On $M_{sa}$ with $\langle x,y\rangle=\operatorname{Tr}(xy)$, $\tau=\operatorname{Tr}$ and $P(y)x=yxy$;
the idempotents of $(M_{sa},\circ)$ are the projections ($c\circ c=c^2$), and for projections $c\circ c'=0$ iff
$cc'=0$ (Lemma F.14(2)), so by the uniqueness in Theorem F.3 its decomposition of $x$ is the spectral decomposition into
eigenprojections and $|x|$ is the matrix absolute value. So this reproves Proposition 11.9 for all $Y=Y^*$. On $V_m$
with $\langle x,y\rangle=\operatorname{tr}(x\circ y)$, $\tau=\operatorname{tr}$, the decomposition of Theorem F.3 is
that of Theorem 13.4(2) (with $x=a1$ if $u=0$), and this is Proposition 13.14. In a simple Euclidean Jordan algebra $J$ of rank $r$, the trace
$\operatorname{tr}x=\frac r{\dim J}\operatorname{Tr}L(x)$ of Corollary E.28 satisfies
$\operatorname{tr}(\sum_i\mu_ic_i)=\sum_i\mu_i$ for every Jordan frame $(c_i)$ (Theorem E.7), and
$\langle x,y\rangle=\operatorname{tr}(x\circ y)$ is an associative inner product; for it $\tau=\operatorname{tr}$, so the
inequality holds for this trace. This applies to $\operatorname{Herm}_r(\mathbb H)$ and to the Albert algebra, which are
simple Euclidean Jordan algebras (Theorem E.24). For $r\ge3$ they are not Jordan isomorphic to $M_{sa}$ for any
$\mathcal M=\bigoplus_cM_{k_c}$: $M_{sa}=\bigoplus_c\operatorname{Herm}_{k_c}(\mathbb C)$ is a direct sum of ideals, so it
is simple only with a single block, and $\operatorname{Herm}_k(\mathbb C)$ ($\cong\mathbb R$ for $k=1$, $\cong V_3$ for
$k=2$) is not isomorphic to $\operatorname{Herm}_r(\mathbb H)$, $r\ge3$, or to $\operatorname{Herm}_3(\mathbb O)$, by the
uniqueness in Theorem 13.7. In $V_m$ the inequality also follows from Proposition 11.9 (extended above to $Y=Y^*$)
through $\phi$, since $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ and $\phi(|x|)=|\phi(x)|$; the code checks both:

```python
for m in (2, 4, 7):
    Vm = SpinFactor(m)
    xs, ys = Vm.random(20_000), Vm.random(20_000)
    lhs = Vm.trace(Vm.apply_function(Vm.quadratic_representation(ys, xs), torch.abs))
    rhs = Vm.trace(Vm.quadratic_representation(ys, Vm.apply_function(xs, torch.abs)))
    definite = Vm.determinant(xs) >= 0                                 # x or -x in the cone
    print(f"m = {m}: min gap {(rhs - lhs).min().item():.1e}, max gap {(rhs - lhs).max().item():.1f}, "
          f"definite x: {definite.double().mean().item():.3f}")
    assert (rhs - lhs >= -1e-12 * rhs).all()                            # Proposition 13.14
    assert ((rhs - lhs)[definite].abs() <= 1e-12 * rhs[definite]).all()
    Xs, Ys = Vm.to_matrix(xs[:200]), Vm.to_matrix(ys[:200])
    tn = torch.linalg.svdvals(Ys @ Xs @ Ys).sum(-1)                     # Tr |phi(y) phi(x) phi(y)|
    assert (tn - Vm.representation_dim / 2 * lhs[:200]).abs().max() < 1e-10 * tn.max()
    wx, Ux = torch.linalg.eigh(Xs)
    absX = Ux @ torch.diag_embed(wx.abs().to(dt)) @ Ux.mH                  # |phi(x)|
    tr_rhs = (Ys @ absX @ Ys).diagonal(dim1=-2, dim2=-1).sum(-1).real    # Tr phi(y) |phi(x)| phi(y)
    assert (torch.linalg.eigvalsh(Ys)[:, 0] < 0).any()                  # some phi(y) are indefinite
    assert (tr_rhs - tn >= -1e-10 * tr_rhs).all()                       # Proposition 11.9 for Y = Y*
```

The gap $\operatorname{tr}P(y)|x|-\operatorname{tr}|P(y)x|$ vanishes for definite $x$ (Proposition 13.14). Numerical
observation (not a theorem): for these Gaussian samples the printed maximal gaps are of order $10^1$–$10^2$, so the
inequality is far from an equality for indefinite $x$.

## 13.7 What to compare with

To our knowledge there is no mainstream numerical library for Jordan algebras; the operations of this chapter are
usually written by hand in terms of matrices. Computer algebra covers the algebraic side: SageMath's `JordanAlgebra` builds
special Jordan algebras, Jordan algebras of symmetric bilinear forms (which include the spin factors) and, from an
octonion algebra, the 27-dimensional exceptional Jordan algebra; Magma also has constructors for Jordan algebras,
including the exceptional one. These are exact algebraic tools rather than batched numerical ones. The Lorentz cones are standard in
optimisation: SOCP solvers (ECOS, SCS, Clarabel, MOSEK, accessible through cvxpy's `SOC` constraint) use the spectral
decomposition of $V_m$ (to project onto the cone) or its quadratic representation (for the Nesterov–Todd scaling of
interior-point methods), and semidefinite solvers do the same for $\operatorname{Herm}_n$. The library adds Jordan products, the quadratic
representation and the transpose as batched, differentiable operators on $\bigoplus_cM_{k_c}$; spin factors with their
functional calculus and matrix representation; and a randomised computation of the Jordan algebra generated by given
matrices, compared with the $*$-algebra they generate.

## Exercises

1. (**Homogeneity of the Lorentz cone.**) Let $y\in V_m$ with $\det y\ne0$. Show that $y^{-1}=(a,-u)/\det y$ satisfies
   $y\circ y^{-1}=1$, that $P(y)P(y^{-1})=\mathrm{id}$, that $P(y)$ maps the interior of $Q$ onto itself, and that
   $P(x^{1/2})1=x$ for $x$ in the interior. Deduce that the interior of $Q$ is homogeneous.
2. (**Isotropic states.**) For $p\in[0,1]$ let $\rho_p=p\,\omega+(1-p)\frac14$ on $\mathbb C^2\otimes\mathbb C^2$
   (Proposition 13.12). Compute the spectrum of $\theta(\rho_p)$ and show that $\rho_p$ is entangled for $p>\frac13$.
3. (**Products of boundary points.**) Let $|w|=|w'|=1$ in $\mathbb R^m$. Show that $(1,w)\circ(1,w')\in Q$ iff
   $w'=w$ or $w'=-w$: the Jordan product of two boundary rays of the Lorentz cone is almost never in the cone.
4. (**Real quantum mechanics.**) Show that for $x,y\in\operatorname{Herm}_n(\mathbb R)$ the derivation
   $[L(x),L(y)]$ of $\operatorname{Herm}_n(\mathbb R)$ is $z\mapsto[K,z]$ with $K$ real antisymmetric, and that for
   $n\ge2$ such a map is never $z\mapsto i[H,z]$ with $H\in\operatorname{Herm}_n(\mathbb R)$ unless it is $0$.

<details>
<summary>Solutions</summary>

1. $(a,u)\circ(a,-u)=(a^2-|u|^2,0)$, which gives $y\circ y^{-1}=1$. $\phi(y)$ and $\phi(y^{-1})$ commute (both lie in
   $\operatorname{span}\{1,u\cdot s\}$), so $\phi(y)\phi(y^{-1})=\phi(y)\circ\phi(y^{-1})=1$, and by Theorem 13.4(5)
   $\phi(P(y)P(y^{-1})x)=\phi(x)$; $\phi$ is injective. $P(y)$ maps $Q$ into $Q$, and so does its inverse $P(y^{-1})$, so
   $P(y)Q=Q$; being a linear homeomorphism, it maps the interior onto the interior. Finally
   $P(y)1=2y\circ y-y^2=y^2$, so $P(x^{1/2})1=x$, and $P(x^{1/2})$ is invertible when $x$ is interior
   ($\det x^{1/2}=\sqrt{\lambda_+\lambda_-}>0$; the interior of $Q$ is $\{a>|u|\}=\{\lambda_->0\}$, since
   $a-|u|$ is continuous and $(a-\varepsilon,u)\notin Q$ when $a=|u|$). Any two interior points are joined by
   $P(x'^{1/2})P(x^{1/2})^{-1}$.

```python
V4 = SpinFactor(4)
xi = V4.random(100, cone=True)                                   # interior points (almost surely)
r = V4.sqrt(xi)
assert (V4.quadratic_representation(r, V4.unit((100,))) - xi).abs().max() < 1e-10
yv = V4.random(100)
back = V4.quadratic_representation(V4.inverse(yv), V4.quadratic_representation(yv, xi))
assert (back - xi).abs().max() < 1e-8 * xi.abs().max()
assert V4.in_cone(V4.quadratic_representation(yv, xi)).all()
```

2. $\theta(\rho_p)=\frac p2F+\frac{1-p}4$ by Proposition 13.12. $F$ has eigenvalues $1$ (three times) and $-1$ (once), so
   the spectrum is $\frac{1+p}4$ (three times) and $\frac{1-3p}4$, which is negative iff $p>\frac13$; then $\rho_p$ is not
   separable by the PPT criterion.

```python
for pp in (0.2, 1 / 3, 0.5):
    rho_p = pp * omega + (1 - pp) * torch.eye(4, dtype=dt) / 4
    lam_min = torch.linalg.eigvalsh(ptranspose(rho_p, 2, 2))[0].item()
    assert abs(lam_min - (1 - 3 * pp) / 4) < 1e-12
```

3. $(1,w)\circ(1,w')=(1+t,\,w+w')$ with $t=w\cdot w'\in[-1,1]$ and $|w+w'|=\sqrt{2+2t}$. It lies in $Q$ iff
   $1+t\ge\sqrt{2(1+t)}$, i.e. iff $t=-1$ (the product is $0$) or $1+t\ge2$, i.e. $t=1$.

```python
V3c = SpinFactor(3)
wv = torch.nn.functional.normalize(torch.randn(2, 1000, 3, dtype=torch.float64), dim=-1)
prod = V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[1]))
assert not V3c.in_cone(prod).any()
assert V3c.in_cone(V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[0]))).all()
```

4. By Proposition 13.13(b), $[L(x),L(y)]z=[K,z]$ with $K=\frac14[x,y]$, and $K^T=\frac14[y^T,x^T]=-K$ is real. For real $H$
   and $z$, $i[H,z]$ is purely imaginary, while $[K,z]$ is real, so they agree only if both vanish. Thus in real
   quantum mechanics the generators of the dynamics, $K$, are not observables.

</details>

**Where to go from here.** This is the last chapter of the course; Appendices E and F complete its proofs. For
further reading, the standard reference for the finite-dimensional theory is [FK94]; K. McCrimmon, *A Taste of Jordan Algebras* (Springer 2004), covers the algebra and the history, and
H. Hanche-Olsen and E. Størmer, *Jordan Operator Algebras* (Pitman 1984), and Alfsen–Shultz (2003) the
infinite-dimensional theory and its relation to C\*-algebras. For the role of composite systems in reconstructions of
quantum theory from operational axioms, see the papers of Barnum and Wilce cited above. For trace inequalities beyond
Chapter 11, E. Carlen, *Trace inequalities and quantum entropy* (Contemp. Math. 529, 2010). The quadratic
representation and the Lorentz cone reappear in conic optimisation: Alizadeh–Goldfarb (2003) and the Jordan-algebraic
interior-point methods surveyed there.
