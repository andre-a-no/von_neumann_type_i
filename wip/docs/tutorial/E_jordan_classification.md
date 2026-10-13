**English** · [Español](es/E_jordan_classification.md) · [Français](fr/E_jordan_classification.md) · [Deutsch](de/E_jordan_classification.md) · [中文](zh/E_jordan_classification.md) · [日本語](ja/E_jordan_classification.md) · [Русский](ru/E_jordan_classification.md)

# Appendix E. Classification of Euclidean Jordan algebras

Chapter 13 used the Jordan–von Neumann–Wigner classification of finite-dimensional formally real Jordan algebras.
This appendix proves it completely, from the Jordan identity and finite-dimensional linear algebra only. The route is
the classical one (J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, Chapters II–V; K. McCrimmon, *A Taste of
Jordan Algebras*): power associativity and a spectral lemma give idempotents and an associative positive trace form;
the Peirce decomposition with respect to a Jordan frame splits the algebra into pieces $J_{ij}$ whose products obey
strict rules; in a simple algebra of rank $r\ge3$ these rules turn the Peirce space $J_{13}$ into a composition algebra, which by
Hurwitz's theorem is $\mathbb R,\mathbb C,\mathbb H$ or $\mathbb O$, and the whole algebra is then
$\operatorname{Herm}_r$ of it. Section E.6 also proves Frobenius' theorem on real associative division algebras, and
Section E.8 proves that the Albert algebra $\operatorname{Herm}_3(\mathbb O)$ is a Jordan algebra that does not embed
into any associative algebra.

**You will learn** why power associativity, positivity of the trace form and the existence of a unit follow from
formal reality alone; the Peirce rules and how they force the coordinates of a simple algebra to multiply like
numbers; the proof of Hurwitz's theorem by Cayley–Dickson doubling; and why the octonions give a new simple algebra only in rank $3$.

**Conventions.** All algebras are finite-dimensional real vector spaces with a bilinear product. A **Jordan algebra**
is as in Section 13.1: $x\circ y=y\circ x$ and $(x^2\circ y)\circ x=x^2\circ(y\circ x)$, $x^2=x\circ x$; a unit is not
assumed. $L(x)y=x\circ y$, $[A,B]=AB-BA$, and $\operatorname{Tr}$ is the trace of a linear map. A **subalgebra** is a
subspace closed under $\circ$, an **ideal** a subspace $I$ with $J\circ I\subseteq I$, and $J$ is **simple** if
$J\ne0$ and its only ideals are $0$ and $J$. An **idempotent** is $c$ with $c^2=c$; idempotents $c,c'$ are
**orthogonal** if $c\circ c'=0$. $J$ is **formally real** if $\sum_ix_i^2=0$ implies that all $x_i=0$, and
**Euclidean** if it has a unit $1$ and an inner product with $\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$.
Polynomial identities are **linearised** as follows: if $p(t)=\sum_kp_kt^k$ has vector coefficients and $p(t)=0$ for
all real $t$, then all $p_k=0$ (evaluate at $\deg p+1$ distinct points and invert the Vandermonde matrix).

## E.1 Linearised Jordan identities and power associativity

**Lemma E.1.** In every Jordan algebra, for all $x,y,z,w$:

1. $[L(x),L(x^2)]=0$ and $[L(x^2),L(z)]+2[L(x\circ z),L(x)]=0$;
2. $L(x^2\circ y)=2L(x\circ y)L(x)+L(x^2)L(y)-2L(x)L(y)L(x)$;
3. $L((x\circ w)\circ y)=L(w\circ y)L(x)+L(x\circ y)L(w)+L(x\circ w)L(y)-L(w)L(y)L(x)-L(x)L(y)L(w)$;
4. $L((x\circ w)\circ y)-L(x\circ(w\circ y))=[L(w),[L(x),L(y)]]$; in particular
   $\operatorname{Tr}L((x\circ w)\circ y)=\operatorname{Tr}L(x\circ(w\circ y))$.

*Proof.* (1) The Jordan identity says $L(x)L(x^2)y=L(x^2)L(x)y$. Replace $x$ by $x+tz$; the coefficient of $t$ is
$[L(z),L(x^2)]+2[L(x),L(x\circ z)]=0$. (2) Replace $x$ by $x+tz$ in the Jordan identity itself; the coefficient of $t$ is
$$
2((x\circ z)\circ y)\circ x+(x^2\circ y)\circ z=2(x\circ z)\circ(y\circ x)+x^2\circ(y\circ z).
$$
Read as linear maps of $z$, the four terms are $2L(x)L(y)L(x)z$, $L(x^2\circ y)z$, $2L(x\circ y)L(x)z$ and
$L(x^2)L(y)z$. (3) Replace $x$ by $x+sw$ in (2) and take the coefficient of $s$, divided by $2$. (4) Write (3) once
as it stands and once with $x$ and $y$ exchanged. The first three terms of the right-hand sides agree, and
$(y\circ w)\circ x=x\circ(w\circ y)$, so the difference is
$-L(w)L(y)L(x)-L(x)L(y)L(w)+L(w)L(x)L(y)+L(y)L(x)L(w)=[L(w),[L(x),L(y)]]$, a commutator, whose trace is $0$.
$\square$

**Lemma E.2 (power associativity).** Put $x^1=x$ and $x^{k+1}=x\circ x^k$. Then $x^p\circ x^q=x^{p+q}$ for all
$p,q\ge1$, and the operators $L(x^k)$ commute pairwise. Hence, for real polynomials $f,g$ without constant term and
$f(x):=\sum_kf_kx^k$, $(fg)(x)=f(x)\circ g(x)$, and $\mathbb R[x]=\{f(x)\}$ is an associative commutative subalgebra.

*Proof.* $L(x)$ and $L(x^2)$ commute by Lemma E.1(1); let $\mathcal A$ be the commutative algebra of operators they
generate. We show by induction on $n\ge2$: $x^p\circ x^q=x^{p+q}$ whenever $p+q\le n$, and $L(x^k)\in\mathcal A$ for
$k\le n$. For $n=2$ this is the definition. Assume it for $n$, and let $p+q=n+1$. If $q=1$, $x^p\circ x=x^{n+1}$ by
commutativity. If $q\ge2$, then
$x^p\circ x^q=L(x^p)L(x)x^{q-1}=L(x)L(x^p)x^{q-1}=L(x)x^{n}=x^{n+1}$. Next, Lemma E.1(2) with $y=x^{n-1}$ and
$x^2\circ x^{n-1}=x^{n+1}$ gives $L(x^{n+1})=2L(x^n)L(x)+L(x^2)L(x^{n-1})-2L(x)^2L(x^{n-1})\in\mathcal A$. The rest
follows by bilinearity. $\square$

## E.2 Idempotents and the Peirce decomposition

**Lemma E.3 (spectral lemma).** Let $J$ be formally real and $x\in J$. There are $k\ge0$, distinct non-zero reals
$\lambda_1,\dots,\lambda_k$ and non-zero, pairwise orthogonal idempotents $c_1,\dots,c_k\in\mathbb R[x]$ with
$$
x=\sum_{i=1}^k\lambda_ic_i,\qquad x^n=\sum_{i=1}^k\lambda_i^nc_i\quad(n\ge1),\qquad\mathbb R[x]=\operatorname{span}\{c_1,\dots,c_k\}.
$$
In particular every non-zero formally real Jordan algebra contains a non-zero idempotent.

*Proof.* Formal reality with one term says that $z^2=0$ implies $z=0$. Let $K=\{q\in\mathbb R[t]:(tq)(x)=0\}$. By
Lemma E.2, $K$ is an ideal of $\mathbb R[t]$: if $q\in K$ and $s=s_0+s_1$ with $s_1(0)=0$, then
$(tqs)(x)=s_0(tq)(x)+(tq)(x)\circ s_1(x)=0$. Since $x,x^2,\dots,x^{\dim J+1}$ are linearly dependent, $K\ne0$, so
$K=\mu\,\mathbb R[t]$ for a monic $\mu$. (i) $\mu(0)\neq0$: if $\mu=th$, then $w=(th)(x)$ has
$w^2=(t\cdot\mu h)(x)=0$, as $\mu h\in K$; so $w=0$, $h\in K$, and $\mu\mid h$, impossible by degree. (ii) $\mu$ has no
repeated factor: if $\mu=g^2h$ with $\deg g\ge1$, then $w=(tgh)(x)$ has $w^2=(t\cdot\mu\cdot th)(x)=0$, so $gh\in K$,
impossible by degree. (iii) $\mu$ has no non-real root: if $\mu=((t-\alpha)^2+\beta^2)h$ with $\beta\ne0$, then
$a=(t(t-\alpha)h)(x)$ and $b=(\beta th)(x)$ have $a^2+b^2=(t\cdot\mu\cdot th)(x)=0$, so $b=0$ and $h\in K$, impossible.
Hence $\mu=\prod_{i=1}^k(t-\lambda_i)$ with distinct non-zero real $\lambda_i$. Put
$$
p_i(t)=\frac t{\lambda_i}\prod_{j\ne i}\frac{t-\lambda_j}{\lambda_i-\lambda_j},\qquad c_i=p_i(x).
$$
The polynomial $p_ip_l-\delta_{il}p_i$ is $t$ times a polynomial vanishing at every $\lambda_j$, hence $t$ times an
element of $K$; so $c_i\circ c_l=\delta_{il}c_i$. Likewise $t-\sum_i\lambda_ip_i=t\big(1-\sum_i\lambda_ip_i/t\big)$ and the
bracket vanishes at every $\lambda_j$, so $x=\sum_i\lambda_ic_i$. If $c_i=0$ then $p_i/t\in K$, of degree $k-1<\deg\mu$;
so $c_i\ne0$. The formula for $x^n$ follows from orthogonality, and it shows $\mathbb R[x]\subseteq\operatorname{span}\{c_i\}$.
$\square$

**Proposition E.4 (Peirce decomposition for an idempotent).** Let $c$ be an idempotent of a Jordan algebra $J$.

1. $L(c)(2L(c)-1)(L(c)-1)=0$. Hence $J=J_1(c)\oplus J_{1/2}(c)\oplus J_0(c)$ with
   $J_\alpha(c)=\{x:c\circ x=\alpha x\}$, and $c\in J_1(c)$.
2. $[L(c),L(z)]=0$ for $z\in J_0(c)\cup J_1(c)$, and $[L(c),L(y^2)]=0$ for $y\in J_{1/2}(c)$.
3. $J_1\circ J_1\subseteq J_1$, $J_0\circ J_0\subseteq J_0$, $J_1\circ J_0=0$, $(J_1+J_0)\circ J_{1/2}\subseteq J_{1/2}$
   and $J_{1/2}\circ J_{1/2}\subseteq J_1+J_0$ (all with respect to $c$).
4. If $c'$ is an idempotent orthogonal to $c$, then $c+c'$ is an idempotent, $[L(c),L(c')]=0$ and
   $J_1(c')\subseteq J_0(c)$.

*Proof.* (1) Lemma E.1(2) with $x=y=c$ gives $L(c)=3L(c)^2-2L(c)^3$. A linear map annihilated by a polynomial with
simple roots is diagonalisable with eigenvalues among the roots. (2) Lemma E.1(1) with $x=c$ gives
$[L(c),L(z)]=2[L(c),L(c\circ z)]$, i.e. $(1-2\alpha)[L(c),L(z)]=0$ for $z\in J_\alpha(c)$. Lemma E.1(1) with $x=y$
and $z=c$ gives $[L(y^2),L(c)]+2[L(y\circ c),L(y)]=[L(y^2),L(c)]+[L(y),L(y)]=[L(y^2),L(c)]=0$. (3) By (2), $L(z)$ maps
each eigenspace of $L(c)$ into itself if $z\in J_0\cup J_1$. This gives the first, second and fourth rule, and
$J_1\circ J_0\subseteq J_1\cap J_0=0$. For $y\in J_{1/2}$, $L(y^2)$ commutes with $L(c)$, so $c\circ y^2=y^2\circ c$ is
an eigenvector of $L(c)$ for $1$, as $c$ is; writing $y^2=a_1+a_{1/2}+a_0$, $c\circ y^2=a_1+\frac12a_{1/2}$ lies in $J_1$,
so $a_{1/2}=0$. Polarising ($y\to y+y'$) gives the last rule. (4) $(c+c')^2=c+2c\circ c'+c'=c+c'$, and
$c'\in J_0(c)$, so $L(c)$ and $L(c')$ commute by (2) and are simultaneously diagonalisable; in particular $J_1(c')$ is
spanned by common eigenvectors. On a common eigenvector
with eigenvalues $\alpha$ for $L(c)$ and $1$ for $L(c')$, $L(c+c')$ has the eigenvalue $\alpha+1\in\{0,\frac12,1\}$, so
$\alpha=0$. $\square$

## E.3 The trace form; Euclidean = formally real

**Theorem E.5.** Let $J$ be formally real and $\tau(x,y)=\operatorname{Tr}L(x\circ y)$.

1. $\tau$ is symmetric, bilinear and **associative**: $\tau(x\circ y,z)=\tau(x,y\circ z)$.
2. $\tau$ is positive definite; for every non-zero idempotent $c$, $\tau(c,c)=\operatorname{Tr}L(c)\ge1$.
3. $J$ has a unit element.
4. $J$ is Euclidean, with the inner product $\tau$. Conversely, every Euclidean Jordan algebra is formally real.

*Proof.* (1) Lemma E.1(4). (2) By Proposition E.4(1), $\operatorname{Tr}L(c)=\dim J_1(c)+\frac12\dim J_{1/2}(c)\ge1$, as
$0\ne c\in J_1(c)$. For $x=\sum_i\lambda_ic_i$ as in Lemma E.3, $x^2=\sum_i\lambda_i^2c_i$ and
$\tau(x,x)=\operatorname{Tr}L(x^2)=\sum_i\lambda_i^2\operatorname{Tr}L(c_i)$, which is $>0$ unless $k=0$, i.e. $x=0$.
(3) Let $J\neq0$ and let $e$ be a non-zero idempotent (Lemma E.3) with $\dim J_1(e)$ maximal. If $J_0(e)\ne0$, it is a
formally real subalgebra (Proposition E.4(3)) and contains a non-zero idempotent $f$; $f$ is orthogonal to $e$, so
$e+f$ is an idempotent, $J_1(e)\subseteq J_1(e+f)$ (as $f\circ J_1(e)\subseteq J_0\circ J_1=0$), and
$f\in J_1(e+f)\setminus J_1(e)$, contradicting maximality. So $J_0(e)=0$. For $y\in J_{1/2}(e)$, then
$y^2\in J_1(e)$, and associativity gives
$\frac12\tau(y,y)=\tau(e\circ y,y)=\tau(e,y^2)=\operatorname{Tr}L(e\circ y^2)=\operatorname{Tr}L(y^2)=\tau(y,y)$, so $y=0$.
Hence $J=J_1(e)$ and $e$ is a unit. (4) By (1)–(3), $\tau(x\circ y,z)=\tau(y,x\circ z)$. Conversely, if $J$ is Euclidean and
$\sum_ix_i^2=0$, then $0=\sum_i\langle x_i\circ x_i,1\rangle=\sum_i\langle x_i,x_i\rangle$. $\square$

So "formally real" and "Euclidean" are the same class of finite-dimensional algebras. From now on $J$ is formally real,
$1$ is its unit, and $\tau$ its trace form. The code checks Lemma E.1(4), Lemma E.3 (the idempotents $p_i(x)$ computed with
Jordan powers only), Proposition E.4 and Theorem E.5 in $\operatorname{Herm}_3(\mathbb C)$, given by structure constants;
the Cayley–Dickson product `cd` used to build it is the one of Section E.6.

```python
import math
import torch

torch.manual_seed(0)
dt = torch.float64

def conj(x):                                      # Cayley-Dickson conjugation (a, b)* = (a*, -b) on R^(2^k)
    if x.shape[-1] == 1:
        return x
    h = x.shape[-1] // 2
    return torch.cat([conj(x[..., :h]), -x[..., h:]], -1)

def cd(x, y):                                     # (a, b)(c, d) = (ac - d*b, da + bc*): R, C, H, O for 1, 2, 4, 8
    if x.shape[-1] == 1:
        return x * y
    h = x.shape[-1] // 2
    a, b, c, d = x[..., :h], x[..., h:], y[..., :h], y[..., h:]
    return torch.cat([cd(a, c) - cd(conj(d), b), cd(d, a) + cd(b, conj(c))], -1)

def herm(n, d):
    """Herm_n(D), D = R, C, H, O for d = 1, 2, 4, 8; coordinates: diagonal, then the entries a_ij, i < j."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    def to_mat(v):
        X = torch.zeros(v.shape[:-1] + (n, n, d), dtype=dt)
        X[..., range(n), range(n), 0] = v[..., :n]
        for p, (i, j) in enumerate(pairs):
            X[..., i, j, :] = v[..., n + d * p:n + d * (p + 1)]
            X[..., j, i, :] = conj(X[..., i, j, :])
        return X
    def from_mat(X):
        return torch.cat([X[..., range(n), range(n), 0]] + [X[..., i, j, :] for i, j in pairs], -1)
    def matmul(X, Y):
        return torch.stack([torch.stack([sum(cd(X[..., i, k, :], Y[..., k, j, :]) for k in range(n))
                                         for j in range(n)], -2) for i in range(n)], -3)
    def prod(x, y):
        X, Y = to_mat(x), to_mat(y)
        return from_mat(0.5 * (matmul(X, Y) + matmul(Y, X)))
    return n + d * len(pairs), prod

def spin(m):                                      # V_m: (a, u) o (b, v) = (ab + u.v, av + bu)
    def prod(x, y):
        return torch.cat([x[..., :1] * y[..., :1] + (x[..., 1:] * y[..., 1:]).sum(-1, keepdim=True),
                          x[..., :1] * y[..., 1:] + y[..., :1] * x[..., 1:]], -1)
    return m + 1, prod

def structure(dim, prod):                         # C[a, b] = e_a o e_b
    E = torch.eye(dim, dtype=dt)
    return prod(E[:, None, :].expand(dim, dim, dim), E[None, :, :].expand(dim, dim, dim))

mul = lambda C, x, y: torch.einsum('...a,...b,abc->...c', x, y, C)
Lop = lambda C, x: torch.einsum('...a,abc->...cb', x, C)          # matrix of L(x)
comm = lambda A, B: A @ B - B @ A

C3 = structure(*herm(3, 2))                       # Herm_3(C), dimension 9
x, y, w = torch.randn(3, 9, dtype=dt)
lhs = Lop(C3, mul(C3, mul(C3, x, w), y)) - Lop(C3, mul(C3, x, mul(C3, w, y)))
assert torch.allclose(lhs, comm(Lop(C3, w), comm(Lop(C3, x), Lop(C3, y))))     # Lemma E.1(4)

def jpow(C, x, k):                                # Jordan powers x^1 = x, x^(k+1) = x o x^k
    p = x
    for _ in range(k - 1):
        p = mul(C, x, p)
    return p

def coords(Z):                                    # complex Hermitian 3 x 3 -> coordinates of Herm_3(C)
    return torch.cat([Z.diagonal().real] + [torch.stack([Z[i, j].real, Z[i, j].imag]) for i, j in ((0, 1), (0, 2), (1, 2))])

G = torch.randn(3, 3, dtype=torch.complex128)
Z = G + G.mH
lam = torch.linalg.eigvalsh(Z)                    # distinct and non-zero (almost surely)
xz = coords(Z)
powers = torch.stack([jpow(C3, xz, k) for k in range(1, 4)])  # x, x^2, x^3
cs = []
for i in range(3):                                # Lagrange idempotents p_i(x) of Lemma E.3
    others = [lam[j] for j in range(3) if j != i]
    coef = torch.tensor([others[0] * others[1], -(others[0] + others[1]), 1.0], dtype=dt)
    den = lam[i] * (lam[i] - others[0]) * (lam[i] - others[1])
    cs.append((coef[:, None] * powers).sum(0) / den)
for i in range(3):
    for j in range(3):
        assert torch.allclose(mul(C3, cs[i], cs[j]), cs[i] if i == j else torch.zeros(9, dtype=dt), atol=1e-10)
assert torch.allclose(sum(l * c for l, c in zip(lam, cs)), xz)

c = cs[0] + cs[1]                                 # the idempotent c_1 + c_2: Proposition E.4
Lc = Lop(C3, c)
assert torch.allclose(Lc @ (2 * Lc - torch.eye(9, dtype=dt)) @ (Lc - torch.eye(9, dtype=dt)), torch.zeros(9, 9, dtype=dt), atol=1e-10)
ev = torch.linalg.eigvals(Lc).real.sort().values
print("eigenvalues of L(c):", ev.numpy().round(6))
tau = lambda C: torch.einsum('abc,c->ab', C, torch.einsum('cbb->c', Lop(C, torch.eye(C.shape[0], dtype=dt))))  # Tr L(e_a o e_b)
T3 = tau(C3)
assert torch.linalg.eigvalsh(T3).min() > 0                                     # Theorem E.5(2)
assert torch.allclose(mul(C3, x, y) @ T3 @ w, x @ T3 @ mul(C3, y, w))         # Theorem E.5(1)
```

## E.4 Jordan frames and the spectral theorem

A non-zero idempotent is **primitive** if it is not the sum of two non-zero orthogonal idempotents. A **Jordan frame**
is a family $(c_1,\dots,c_r)$ of pairwise orthogonal primitive idempotents with $\sum_ic_i=1$.

**Lemma E.6.** (1) If $c$ is primitive, then $J_1(c)=\mathbb Rc$. (2) Pairwise orthogonal non-zero idempotents are
linearly independent; every non-zero idempotent is a sum of pairwise orthogonal primitive idempotents; and $J$ has a
Jordan frame.

*Proof.* (1) $J_1(c)$ is a formally real subalgebra with unit $c$ (Proposition E.4(3)). For $x\in J_1(c)$ write
$x=\sum_i\lambda_ic_i$ by Lemma E.3, with $c_i\in\mathbb R[x]\subseteq J_1(c)$. Then $c-c_i$ is an idempotent orthogonal
to $c_i$ ($(c-c_i)^2=c-2c_i+c_i$ and $c_i\circ(c-c_i)=0$), so $c-c_i=0$ by primitivity. Distinct $c_i$ are orthogonal,
so there is at most one, and $x\in\mathbb Rc$. (2) If $\sum_ia_if_i=0$, the product with $f_j$ gives $a_jf_j=0$. So an
idempotent $e$ is a sum of a largest possible number of pairwise orthogonal non-zero idempotents. If one of them, $f$,
were $g+g'$ with $g,g'$ orthogonal non-zero idempotents, then $g=f\circ g\in J_1(f)\subseteq J_0(f')$ for every other
member $f'$ (Proposition E.4(4)), and likewise $g'$, so $f$ could be replaced by $g,g'$. Take $e=1$. $\square$

**Theorem E.7 (spectral theorem).** For every $x\in J$ there are a Jordan frame $(c_1,\dots,c_r)$ and real numbers
$\lambda_1,\dots,\lambda_r$ (not necessarily distinct) with $x=\sum_i\lambda_ic_i$.

*Proof.* Write $x=\sum_{i=1}^k\lambda_ie_i$ by Lemma E.3 and add $e_0=1-\sum_ie_i$, $\lambda_0=0$; $e_0$ is an
idempotent orthogonal to all $e_i$. Split every non-zero $e_i$ into pairwise orthogonal primitive idempotents
(Lemma E.6(2)). A piece $f$ of $e_i$ lies in $J_1(e_i)\subseteq J_0(e_l)$ for $l\ne i$, and a piece of $e_l$ lies in
$J_1(e_l)$; as $J_0(e_l)\circ J_1(e_l)=0$, all pieces are pairwise orthogonal. Their sum is $\sum_ie_i=1$. $\square$

That all Jordan frames of $J$ have the same number of elements is proved in Corollary E.28.

**Theorem E.8 (Peirce decomposition for a frame).** Let $(c_1,\dots,c_r)$ be a Jordan frame, $t_i=\tau(c_i,c_i)$ and,
for $i\ne j$, $J_{ij}=J_{ji}=J_{1/2}(c_i)\cap J_{1/2}(c_j)$.

1. The $L(c_k)$ commute, $J=\bigoplus_i\mathbb Rc_i\oplus\bigoplus_{i<j}J_{ij}$, and
   $c_k\circ x=\frac12(\delta_{ki}+\delta_{kj})x$ for $x\in J_{ij}$.
2. For distinct $i,j,k$: $J_{ij}\circ J_{jk}\subseteq J_{ik}$; and $J_{ij}\circ J_{kl}=0$ if $\{i,j\}\cap\{k,l\}=\emptyset$.
3. There is a positive definite quadratic form $N$ on $J_{ij}$ with $x^2=N(x)(c_i+c_j)$; with its polar form
   $\langle\cdot,\cdot\rangle$, $x\circ y=\langle x,y\rangle(c_i+c_j)$ and $\tau(x,y)=(t_i+t_j)\langle x,y\rangle$ for
   $x,y\in J_{ij}$. If $J_{ij}\ne0$, then $t_i=t_j$.

*Proof.* (1) The $L(c_k)$ commute by Proposition E.4(4), are diagonalisable and add up to $L(1)=\mathrm{id}$. A common
eigenvector has eigenvalues $\alpha_k\in\{0,\frac12,1\}$ with $\sum_k\alpha_k=1$: either one $\alpha_i=1$ and the others
vanish, and then the vector lies in $J_1(c_i)=\mathbb Rc_i$ (Lemma E.6(1)); or $\alpha_i=\alpha_j=\frac12$ for one pair and
the others vanish. So $J_{ij}$ is exactly the common eigenspace of the second kind, and $J$ is the direct sum of the common
eigenspaces. (2) Let $x\in J_{ij}$ and $y\in J_{jk}$. With respect to $c_i$, $x\in J_{1/2}$ and $y\in J_0$, so
$x\circ y\in J_{1/2}(c_i)$; with respect to $c_k$, likewise $x\circ y\in J_{1/2}(c_k)$. Hence $x\circ y\in J_{ik}$. For
$y\in J_{kl}$ instead, $x\circ y$ lies in $J_{1/2}(c_p)$ for $p=i,j,k,l$; a common eigenvector with four eigenvalues
$\frac12$ does not exist, so $x\circ y=0$. (3) For $x,y\in J_{ij}$, $x\circ y\in J_0(c_p)+J_1(c_p)$ for $p=i,j$ and
$x\circ y\in J_0(c_k)$ for $k\notin\{i,j\}$ (Proposition E.4(3)); by (1), $x\circ y\in\mathbb Rc_i+\mathbb Rc_j$. Write
$x^2=\alpha c_i+\beta c_j$. The Jordan identity with $y=c_i$ gives $(x^2\circ c_i)\circ x=\frac12\alpha x$ and
$x^2\circ(c_i\circ x)=\frac12x^2\circ x=\frac14(\alpha+\beta)x$, so $\alpha=\beta$ if $x\ne0$; put $N(x)=\alpha$. Since
$\tau(c_i,1)=\operatorname{Tr}L(c_i)=t_i$, $\tau(x,x)=\tau(x^2,1)=N(x)(t_i+t_j)$; so $N$ is a quadratic form, positive
definite by Theorem E.5(2), and polarisation gives the formulas for $x\circ y$ and $\tau$. Finally
$N(x)t_i=\tau(x^2,c_i)=\tau(x,x\circ c_i)=\frac12\tau(x,x)$, and the same with $j$; so $t_i=t_j$ if some $x\ne0$. $\square$

## E.5 Simple algebras and rank two

**Theorem E.9.** $J$ is a direct sum $J^{(1)}\oplus\dots\oplus J^{(s)}$ of simple ideals with
$J^{(p)}\circ J^{(q)}=0$ for $p\ne q$. Each $J^{(p)}$ is a formally real Jordan algebra.

*Proof.* If $I$ is an ideal, so is its orthogonal complement $I^\perp$ for $\tau$:
$\tau(J\circ I^\perp,I)=\tau(I^\perp,J\circ I)=0$. Hence $J=I\oplus I^\perp$ and $I\circ I^\perp\subseteq I\cap I^\perp=0$, so
$J$ is the direct sum of the algebras $I$ and $I^\perp$, and an ideal of one summand is an ideal of $J$. Let $I$ be a
non-zero ideal of smallest dimension; then $I$ is simple, and we conclude by induction on the dimension, applied to the
formally real algebra $I^\perp$. $\square$

**Lemma E.10.** Let $(c_i)$ be a Jordan frame, $i,j,k,l$ distinct, $u\in J_{ij}$, $w\in J_{jk}$ and $z\in J_{kl}$. Then
(1) $u\circ(u\circ w)=\frac14N(u)w$, and (2) $(u\circ w)\circ z=u\circ(w\circ z)$.

*Proof.* (1) Apply Lemma E.1(1) with $x=u$ to $c_j$. With $L(u^2)=N(u)L(c_i+c_j)$, $c_j\circ w=\frac12w$,
$c_j\circ u=\frac12u$ and $c_j\circ(u\circ w)=0$ (as $u\circ w\in J_{ik}$):
$[L(u^2),L(w)]c_j=\frac14N(u)w-\frac12N(u)w$ and $2[L(u\circ w),L(u)]c_j=(u\circ w)\circ u$; their sum vanishes. (2) Apply
Lemma E.1(3), with $(x,w,y)$ replaced by $(u,z,w)$, to $c_j$. Here $u\circ z=0$ (Theorem E.8(2)), $c_j\circ u=\frac12u$
and $c_j\circ z=0$, so $0=\frac12(z\circ w)\circ u-\frac12z\circ(w\circ u)$. $\square$

**Proposition E.11.** Let $(c_1,\dots,c_r)$ be a Jordan frame of $J$, and join $i\neq j$ by an edge if $J_{ij}\ne0$.

1. $J$ is simple if and only if this graph is connected.
2. If $J$ is simple and $r\ge2$, all $J_{ij}$, $i\ne j$, have the same dimension $d\ge1$, all $t_i$ equal
   $t=1+\frac{(r-1)d}2$, and $\dim J=r+d\,\frac{r(r-1)}2=rt$.

*Proof.* (1) If the graph is not connected, let $S$ be a connected component and $I$ the span of the $c_i$ and $J_{ij}$
with $i,j\in S$, $I'$ the same for the complement. Then $J=I\oplus I'$ (all $J_{ij}$ with $i\in S$, $j\notin S$ vanish),
$I\circ I\subseteq I$ and $I\circ I'=0$ by Theorem E.8; so $I$ is an ideal, $\ne0,J$. Conversely let the graph be
connected and $I\ne0$ an ideal. $I$ is invariant under the commuting diagonalisable $L(c_k)$, so it is the sum of its
intersections with the common eigenspaces. If $0\ne x\in I\cap J_{ij}$, then $c_i\circ x^2=N(x)c_i\in I$; so some $c_i$
lies in $I$. If $c_i\in I$ and $0\ne x\in J_{ij}$, then $x=2c_i\circ x\in I$ and $c_j=c_j\circ x^2/N(x)\in I$. By
connectedness all $c_j\in I$, so $1\in I$ and $I=J$. (2) Let $u\in J_{ij}$ with $N(u)=1$ and $k\notin\{i,j\}$. By
Theorem E.8(2) $2L(u)$ maps $J_{jk}$ to $J_{ik}$ and $J_{ik}$ to $J_{jk}$, and by Lemma E.10(1) (also with $i,j$
exchanged) its square is the identity on both. So $\dim J_{ik}=\dim J_{jk}$ whenever $J_{ij}\ne0$. Thus edges $ij$, $jk$
force an edge $ik$; a connected graph is then complete, and changing one index at a time shows that all $\dim J_{ij}$
are equal. Finally $L(c_i)$ has the eigenvalue $1$ on $\mathbb Rc_i$, $\frac12$ on the $r-1$ spaces $J_{ij}$ and $0$
elsewhere, so $t_i=\operatorname{Tr}L(c_i)=1+\frac{(r-1)d}2$. $\square$

For the rest of the section $J$ is simple, $(c_1,\dots,c_r)$ is a Jordan frame with $r\ge2$, and for $x\in J_{ij}$,
$y\in J_{jk}$ ($i,j,k$ distinct) we write
$$
m(x,y)=2\,x\circ y\in J_{ik}.
$$
By Theorem E.8(3) and Proposition E.11(2), $\tau(x,y)=2t\langle x,y\rangle$ on every $J_{ij}$.

**Lemma E.12.** For distinct $i,j,k$ and $x\in J_{ij}$, $y\in J_{jk}$, $z\in J_{ik}$:

1. $N(m(x,y))=N(x)N(y)$ and $m(x,m(x,y))=N(x)y$;
2. $\langle m(x,y),z\rangle=\langle y,m(x,z)\rangle=\langle x,m(y,z)\rangle$;
3. $m(m(x,y),v)=m(x,m(y,v))$ for $l\notin\{i,j,k\}$ and $v\in J_{kl}$.

*Proof.* (2) The three numbers are $\tau(x\circ y,z)/t$, $\tau(y,x\circ z)/t$ and $\tau(x,y\circ z)/t$, equal by
associativity and commutativity. (1) The second formula is Lemma E.10(1). Then
$2tN(m(x,y))=4\tau(x\circ y,x\circ y)=4\tau(y,x\circ(x\circ y))=N(x)\tau(y,y)=2tN(x)N(y)$. (3) Lemma E.10(2). $\square$

**Theorem E.13 (rank one and two).** Let $J$ be simple with a Jordan frame of $r$ elements. If $r=1$, then
$J=\mathbb Rc_1\cong\mathbb R$. If $r=2$, then $J\cong V_{d+1}$ with $d=\dim J_{12}\ge1$ (spin factors as in Section 13.2).

*Proof.* $r=1$: Theorem E.8(1). $r=2$: $J=\mathbb R1\oplus W$ with $W=\mathbb Rf\oplus J_{12}$, $f=c_1-c_2$. By
Theorem E.8, $f\circ f=1$, $f\circ x=\frac12x-\frac12x=0$ and $x\circ y=\langle x,y\rangle1$ for $x,y\in J_{12}$. So
$w\circ w'=B(w,w')1$ for $w,w'\in W$, where $B(\alpha f+x,\beta f+y)=\alpha\beta+\langle x,y\rangle$ is an inner product.
If $U:W\to\mathbb R^{d+1}$ is an isometry, $a1+w\mapsto(a,Uw)$ maps
$(a1+w)\circ(b1+w')=(ab+B(w,w'))1+aw'+bw$ to the product of $V_{d+1}$. $\square$

## E.6 Composition algebras: the theorems of Hurwitz and Frobenius

A **composition algebra** is an algebra $A$ with unit $1$ and a positive definite quadratic form $N$ such that
$N(xy)=N(x)N(y)$. We write $\langle\cdot,\cdot\rangle$ for the polar form ($N(x)=\langle x,x\rangle$),
$\operatorname{Re}x=\langle x,1\rangle$ and $\bar x=2\operatorname{Re}(x)1-x$. Since $N(1)=1$ (Lemma E.14(1), whose
proof does not use the conjugation), this is the reflection in the line $\mathbb R1$, so $x\mapsto\bar x$ is an isometric
involution and $\operatorname{Re}\bar x=\operatorname{Re}x$. An isomorphism of composition algebras is an isometric algebra
isomorphism. The **associator** is $(x,y,z)=(xy)z-x(yz)$.

**Lemma E.14.** In a composition algebra, for all $x,y,z,w$:

1. $N(1)=1$, $\langle xy,xz\rangle=N(x)\langle y,z\rangle$, $\langle yx,zx\rangle=\langle y,z\rangle N(x)$ and
   $\langle xy,wz\rangle+\langle wy,xz\rangle=2\langle x,w\rangle\langle y,z\rangle$;
2. $\langle xy,z\rangle=\langle y,\bar xz\rangle=\langle x,z\bar y\rangle$;
3. $\bar x(xy)=N(x)y=(yx)\bar x$; hence $\bar xx=x\bar x=N(x)1$, $x^2-2\operatorname{Re}(x)x+N(x)1=0$,
   $xy+yx=2\operatorname{Re}(x)y+2\operatorname{Re}(y)x-2\langle x,y\rangle1$,
   $\bar x(wy)+\bar w(xy)=2\langle x,w\rangle y$ and $(yx)\bar w+(yw)\bar x=2\langle x,w\rangle y$;
4. $\overline{xy}=\bar y\bar x$, $\operatorname{Re}(xy)=\operatorname{Re}(yx)$ and
   $\operatorname{Re}((xy)z)=\operatorname{Re}(x(yz))$;
5. $x(xy)=x^2y$ and $(yx)x=yx^2$; the associator is alternating (it changes sign when two arguments are exchanged).

*Proof.* (1) $N(1)=N(1)^2>0$. Polarise $N(xy)=N(x)N(y)$ in $y$, respectively in $x$, and then $\langle xy,xz\rangle=N(x)\langle y,z\rangle$ in $x$.
(2) The last identity of (1) with $w=1$ reads $\langle xy,z\rangle+\langle y,xz\rangle=2\operatorname{Re}(x)\langle y,z\rangle$,
i.e. $\langle xy,z\rangle=\langle y,\bar xz\rangle$; the polarisation of $\langle yx,zx\rangle=\langle y,z\rangle N(x)$ in $x$,
with one of the two variables equal to $1$, gives $\langle yx,z\rangle=\langle y,z\bar x\rangle$. (3) By (2) and (1),
$\langle\bar x(xy),z\rangle=\langle xy,xz\rangle=N(x)\langle y,z\rangle$ for all $z$, and likewise on the right. Put $y=1$
and substitute $\bar x=2\operatorname{Re}(x)1-x$; then polarise. (4) By (2) and the isometry of the conjugation,
$\langle\bar y\bar x,z\rangle=\langle\bar x,yz\rangle=\langle\bar x\bar z,y\rangle=\langle\bar z,xy\rangle=\langle z,\overline{xy}\rangle$.
Next $\operatorname{Re}(xy)=\langle y,\bar x\rangle=\langle\bar y,x\rangle=\operatorname{Re}(yx)$, and
$\operatorname{Re}((xy)z)=\langle xy,\bar z\rangle=\langle y,\bar x\bar z\rangle=\langle yz,\bar x\rangle=\operatorname{Re}(x(yz))$.
(5) By (3), $x(xy)=(2\operatorname{Re}(x)-\bar x)(xy)=(2\operatorname{Re}(x)x-N(x))y=x^2y$, and likewise on the right. So
$(x,x,y)=(y,x,x)=0$; polarising, the associator changes sign when the first two or the last two arguments are
exchanged, and these transpositions generate all permutations. $\square$

**The Cayley–Dickson construction.** For an algebra $B$ with unit, involution $b\mapsto\bar b$ and quadratic form $N$, let
$\operatorname{CD}(B)=B\times B$ with
$$
(a,b)(c,d)=(ac-\bar db,\;da+b\bar c),\qquad N(a,b)=N(a)+N(b),\qquad\overline{(a,b)}=(\bar a,-b).
$$
Starting from $\mathbb R$ ($\bar x=x$, $N(x)=x^2$) we define $\mathbb C=\operatorname{CD}(\mathbb R)$,
$\mathbb H=\operatorname{CD}(\mathbb C)$ and $\mathbb O=\operatorname{CD}(\mathbb H)$, of dimensions $1,2,4,8$, each doubling
using the involution and the form just defined (by Lemma E.15(1) below, they are the conjugation and the norm of the
composition algebra). In
$\mathbb C$, $(a,b)(c,d)=(ac-bd,ad+bc)$: these are the complex numbers $a+bi$. $B$ is the subalgebra $B\times0$.

**Lemma E.15.** (1) If $B$ is an associative composition algebra, then $\operatorname{CD}(B)$ is a composition algebra
with unit $(1,0)$ and conjugation $\overline{(a,b)}=(\bar a,-b)$. (2) $\mathbb C$ is commutative and associative.
$M(a,b)=\begin{pmatrix}a&-b\\\bar b&\bar a\end{pmatrix}$ is an injective algebra homomorphism $\mathbb H\to M_2(\mathbb C)$
with $M(\bar q)=M(q)^*$; so $\mathbb H$ is associative; it is not commutative. $\mathbb O$ is not associative. (3)
$\mathbb R,\mathbb C,\mathbb H,\mathbb O$ are composition algebras, and $\mathbb R,\mathbb C,\mathbb H$ are division algebras:
every $x\ne0$ has the two-sided inverse $\bar x/N(x)$, so $xy=0$ implies $x=0$ or $y=0$.

*Proof.* (1) $N((a,b)(c,d))=N(a)N(c)+N(d)N(b)+N(d)N(a)+N(b)N(c)-2\langle ac,\bar db\rangle+2\langle da,b\bar c\rangle$, and by
Lemma E.14(2) $\langle ac,\bar db\rangle=\langle d(ac),b\rangle$ and $\langle da,b\bar c\rangle=\langle(da)c,b\rangle$, which
agree when $B$ is associative. $(1,0)$ is a unit, and $2\langle(a,b),(1,0)\rangle(1,0)-(a,b)=(\bar a,-b)$. (2) Multiplying
out and using that $\mathbb C$ is commutative, $M(a,b)M(c,d)=M(ac-\bar db,da+b\bar c)$; $M$ is linear and injective, and
$M(\bar a,-b)=M(a,b)^*$. In $\mathbb H$, $(i,0)(0,1)=(0,i)$ but $(0,1)(i,0)=(0,-i)$. In $\mathbb O=\mathbb H\times\mathbb H$ let
$p=((i,0),0)$, $q=((0,1),0)$, $s=(0,1)$; then $(pq)s=(0,(0,i))$ and $p(qs)=(0,(0,-i))$. (3) By (1), (2) and induction. If
$x\ne0$, $x(\bar x/N(x))=(\bar x/N(x))x=1$ by Lemma E.14(3); in an associative algebra, $xy=0$ then gives
$y=(\bar x/N(x))(xy)=0$. $\square$

**Lemma E.16 (doubling).** Let $A$ be a composition algebra, $B\subseteq A$ a subalgebra with $1\in B$, and $i\in A$ with
$i\perp B$ and $N(i)=1$. For $a,b,c,d\in B$: $Bi\perp B$, $N(a+bi)=N(a)+N(b)$ and
$$
(a+bi)(c+di)=(ac-\bar db)+(da+b\bar c)\,i .
$$
So $B\oplus Bi$ is a subalgebra of dimension $2\dim B$, and $a+bi\mapsto(h(a),h(b))$ is an isomorphism
$B\oplus Bi\cong\operatorname{CD}(D)$ whenever $h:B\to D$ is an isomorphism of composition algebras. Moreover $B$ is
associative; if $B\oplus Bi$ is associative, $B$ is commutative; if $B\oplus Bi$ is commutative, $B=\mathbb R1$.

*Proof.* $B$ contains $\bar b=2\operatorname{Re}(b)1-b$ with $b$. By Lemma E.14, $\langle bi,a\rangle=\langle i,\bar ba\rangle=0$,
$N(bi)=N(b)$, $\bar i=-i$, and $ci+ic=2\operatorname{Re}(c)i$, i.e. $ci=i\bar c$, for $c\in B$. The identities of
Lemma E.14(3) give, in turn: with $(x,w,y)=(\bar a,i,b)$, $a(ib)=i(\bar ab)$; since $ib=\bar bi$ and
$i(\bar ab)=\overline{\bar ab}\,i=(\bar ba)i$, replacing $b$ by $\bar b$ gives $a(bi)=(ba)i$. With $(y,x,w)=(a,i,\bar c)$,
$(ai)c=(a\bar c)i$. With $x=i$, $i(iy)=-y$. With $(x,w,y)=(ai,i,\bar d)$, where $\overline{ai}=-ai$ (as $ai\perp1$) and
$\langle ai,i\rangle=\operatorname{Re}a$ (Lemma E.14(1)),
$$
-(ai)(i\bar d)-i((ai)\bar d)=2\operatorname{Re}(a)\bar d .
$$
Here $i\bar d=di$ and $i((ai)\bar d)=i((ad)i)=i(i\,\overline{ad})=-\bar d\bar a$, so
$(ai)(di)=\bar d\bar a-2\operatorname{Re}(a)\bar d=-\bar da$. Expanding $(a+bi)(c+di)$ with these four rules gives the formula,
and $N$ and the products match those of $\operatorname{CD}(D)$ under $h$. Since $N$ is multiplicative on
$B\oplus Bi\subseteq A$, the computation in the proof of Lemma E.15(1) gives $\langle d(ac)-(da)c,b\rangle=0$ for all
$b\in B$, so $B$ is associative. If $B\oplus Bi$ is associative, $(ab)i=a(bi)=(ba)i$, so $ab=ba$. If it is commutative,
$ai=ia=\bar ai$, so $a=\bar a\in\mathbb R1$. $\square$

**Theorem E.17 (Hurwitz).** Every composition algebra is isomorphic to exactly one of $\mathbb R,\mathbb C,\mathbb H,\mathbb O$.

*Proof.* Put $B_0=\mathbb R1\cong\mathbb R$. While $B_k\ne A$, choose $i_{k+1}\perp B_k$ with $N(i_{k+1})=1$ and put
$B_{k+1}=B_k\oplus B_ki_{k+1}$; by Lemma E.16, $B_k\cong\operatorname{CD}^k(\mathbb R)$ has dimension $2^k$. If $B_4$
existed, Lemma E.16 would make $B_3$ associative, hence $B_2$ commutative, hence $B_1=\mathbb R1$, which has dimension
$1\ne2$. So $A=B_k$ for some $k\le3$. The four algebras have different dimensions. $\square$

**Lemma E.18.** Let $A$ be a composition algebra. (1) Any two elements of $A$ lie in an associative subalgebra.
(2) (**Moufang identity**) $(ux)(yu)=u((xy)u)=(u(xy))u$ for all $u,x,y$.

*Proof.* (1) Let $p,q\in A$; we may assume $p\notin\mathbb R1$ (else exchange them, and if both are real use $\mathbb R1$).
$B=\operatorname{span}\{1,p\}$ is a subalgebra ($p^2=2\operatorname{Re}(p)p-N(p)1$), a composition algebra of dimension $2$,
so $B\cong\mathbb C$ (Theorem E.17). If $q\in B$ we are done; otherwise let $i$ be the normalised component of $q$
orthogonal to $B$. By Lemma E.16, $B\oplus Bi\cong\operatorname{CD}(\mathbb C)=\mathbb H$ is associative and contains $p,q$.
(2) The associator is alternating (Lemma E.14(5)), so $(u,xy,u)=0$, which is the second equality. Let $X=(u,x,y)$.
$\operatorname{Re}X=0$ by Lemma E.14(4), and $\langle u,X\rangle=\langle u\bar y,ux\rangle-N(u)\langle1,xy\rangle
=N(u)(\langle\bar y,x\rangle-\langle\bar x,y\rangle)=0$; so the polarised quadratic identity of Lemma E.14(3) gives
$uX+Xu=2\operatorname{Re}(u)X$, i.e. $\bar uX=Xu$. Expanding the definitions shows that in every algebra
$$
(wx,y,z)-(w,xy,z)+(w,x,yz)=w(x,y,z)+(w,x,y)z .
$$
With $(w,x,y,z)=(u,u,x,y)$, and $(u^2,x,y)=2\operatorname{Re}(u)X$: $(u,ux,y)=2\operatorname{Re}(u)X-uX=\bar uX=Xu$. With
$(w,x,y,z)=(u,x,y,u)$, $(ux,y,u)=(u,ux,y)=Xu$ and $(x,y,u)=X$ (cyclic permutations): $(u,x,yu)=uX+Xu-Xu=uX$. Hence
$(ux)(yu)-u((xy)u)=(u,x,yu)-u(x,y,u)=uX-uX=0$. $\square$

**Theorem E.19 (Frobenius).** Every finite-dimensional associative real algebra $A$ with unit $1\ne0$ and without zero
divisors is isomorphic to $\mathbb R$, $\mathbb C$ or $\mathbb H$.

*Proof.* (i) For $x\in A$ the powers $1,x,x^2,\dots$ are linearly dependent; let $\mu$ be the monic real polynomial of least
degree with $\mu(x)=0$. If $\mu=fg$ with $\deg f,\deg g\ge1$, then $f(x)g(x)=0$ forces $f(x)=0$ or $g(x)=0$; so $\mu$ is
irreducible, hence (fundamental theorem of algebra) $\mu=t-\alpha$ or $\mu=(t-\alpha)^2+\beta^2$, $\beta>0$. Let
$V=\{v\in A:v^2\in(-\infty,0]\cdot1\}$. Then $A=\mathbb R1+V$ (in the second case $x-\alpha\in V$), and
$\mathbb R1\cap V=0$. (ii) $V$ is a subspace: it is closed under scalars, and for $u,v\in V$ we show $u+v\in V$. If
$u,v$ are linearly dependent, $u+v$ is a multiple of one of them; so let them be independent. Then $1,u,v$ are independent: if
$1=\alpha u+\gamma v$ with, say, $\alpha\ne0$, then $u^2=(1-2\gamma v+\gamma^2v^2)/\alpha^2\in\mathbb R1$ forces $\gamma=0$
and $u\in\mathbb R1\cap V=0$. By (i), $(u\pm v)^2=p_\pm(u\pm v)+q_\pm1$ with real $p_\pm,q_\pm$; adding,
$2u^2+2v^2=(p_++p_-)u+(p_+-p_-)v+(q_++q_-)1$, so $p_\pm=0$. If $q_+>0$, then
$(u+v-\sqrt{q_+})(u+v+\sqrt{q_+})=0$ and $u+v\in\mathbb R1$, which is impossible; so $u+v\in V$. (iii) For $u,v\in V$,
$uv+vu=(u+v)^2-u^2-v^2=-2\langle u,v\rangle1$ defines a symmetric bilinear form with $\langle u,u\rangle=-u^2\ge0$, and
$\langle u,u\rangle=0$ only if $u^2=0$, i.e. $u=0$: an inner product. Orthogonal $u,v\in V$ anticommute. (iv) If $V=0$,
$A=\mathbb R$. If $\dim V=1$, $V=\mathbb Ri$ with $i^2=-1$, and $A\cong\mathbb C$. If $\dim V\ge2$, take orthonormal
$i,j\in V$ and $k=ij$. Then $k^2=ijij=-i^2j^2=-1$, so $k\in V$, $ik+ki=-j+j=0$ and $jk+kj=i-i=0$: $i,j,k$ are
orthonormal, with $ij=-ji=k$, $jk=-kj=i$, $ki=-ik=j$. If there were a unit $l\in V$ orthogonal to $i,j,k$, then
$lk=(li)j=-i(lj)=(ij)l=kl$ and $lk=-kl$, so $kl=0$, impossible. Hence $A=\operatorname{span}\{1,i,j,k\}$. In $\mathbb H$, the
elements $I=(i,0)$, $J=(0,1)$, $K=(0,i)$ have the same multiplication table (Lemma E.15(2) gives $IJ=K=-JI$, and
$I^2=J^2=-1$, $K=IJ$ determine the rest by associativity), so $A\cong\mathbb H$. $\square$

The code builds $\mathbb O$ by Cayley–Dickson doubling and checks the identities of Lemma E.14, the Moufang identity
of Lemma E.18 and the non-associativity of Lemma E.15 on random octonions, and the homomorphism $M$:

```python
o = torch.randn(4, 1000, 8, dtype=dt)             # random octonions
a, b, cc, u = o
N = lambda z: (z * z).sum(-1)
ip = lambda p, q: (p * q).sum(-1)
assert torch.allclose(N(cd(a, b)), N(a) * N(b))                               # composition
assert torch.allclose(ip(cd(a, b), cc), ip(b, cd(conj(a), cc)))               # Lemma E.14(2)
assert torch.allclose(ip(cd(a, b), cc), ip(a, cd(cc, conj(b))))
assert torch.allclose(conj(cd(a, b)), cd(conj(b), conj(a)))                   # Lemma E.14(4)
assert torch.allclose(cd(a, cd(a, b)), cd(cd(a, a), b))                       # alternative
assert torch.allclose(cd(cd(u, a), cd(b, u)), cd(u, cd(cd(a, b), u)))         # Moufang
assoc = cd(cd(a, b), cc) - cd(a, cd(b, cc))
assert assoc[..., 0].abs().max() < 1e-12 and assoc.abs().max() > 1          # Re-associative, not associative
e = torch.eye(8, dtype=dt)
print("(e1 e2) e4 =", cd(cd(e[1], e[2]), e[4]).tolist(), "\ne1 (e2 e4) =", cd(e[1], cd(e[2], e[4])).tolist())

def M(q):                                         # H -> M_2(C): (a, b) -> [[a, -b], [b*, a*]]
    za, zb = torch.complex(q[..., 0], q[..., 1]), torch.complex(q[..., 2], q[..., 3])
    return torch.stack([torch.stack([za, -zb], -1), torch.stack([zb.conj(), za.conj()], -1)], -2)
p4, q4 = torch.randn(2, 1000, 4, dtype=dt)
assert torch.allclose(M(cd(p4, q4)), M(p4) @ M(q4)) and torch.allclose(M(conj(p4)), M(p4).mH)
```

## E.7 Coordinatisation

For $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$, $d=\dim D$ and $n\ge1$, let
$\operatorname{Herm}_n(D)=\{X\in M_n(D):X_{ji}=\overline{X_{ij}}\}$ with $X\circ Y=\frac12(XY+YX)$, where
$(XY)_{ik}=\sum_jX_{ij}Y_{jk}$ (no associativity of $D$ is needed to define this). Diagonal entries satisfy $\bar x=x$, so
they are real. Let $E_{ij}$ be the matrix units and, for $i\ne j$ and $a\in D$, $a[ij]=aE_{ij}+\bar aE_{ji}$; thus
$a[ji]=\bar a[ij]$.

**Lemma E.20.** $\operatorname{Herm}_n(D)$ is closed under $\circ$, has the unit $\sum_iE_{ii}$ and dimension
$n+d\frac{n(n-1)}2$, and for distinct $i,j,k,l$ and $a,b\in D$:
$$
E_{ii}\circ E_{jj}=\delta_{ij}E_{ii},\quad E_{kk}\circ a[ij]=\tfrac12(\delta_{ki}+\delta_{kj})a[ij],\quad
a[ij]\circ b[ij]=\langle a,b\rangle(E_{ii}+E_{jj}),\quad a[ij]\circ b[jk]=\tfrac12(ab)[ik],\quad a[ij]\circ b[kl]=0 .
$$
Moreover $(X|Y)=\operatorname{Re}\operatorname{tr}(XY)$ is an inner product with $(X\circ Y|Z)=(Y|X\circ Z)$.

*Proof.* Matrix multiplication: $a[ij]\,b[ij]=a\bar bE_{ii}+\bar abE_{jj}$, and $a\bar b+b\bar a=\bar ab+\bar ba=2\langle a,b\rangle$
(polarise $x\bar x=\bar xx=N(x)$); $a[ij]\,b[jk]=abE_{ik}$ and $b[jk]\,a[ij]=\bar b\bar aE_{ki}=\overline{ab}E_{ki}$; the
other products are immediate. All basis products lie in $\operatorname{Herm}_n(D)$, hence so do all products. Next,
$(X|X)=\sum_{i,j}\operatorname{Re}(X_{ij}\overline{X_{ij}})=\sum_{i,j}N(X_{ij})>0$ for $X\ne0$. Let
$T(X,Y,Z)=\operatorname{Re}\operatorname{tr}((XY)Z)=\sum_{i,j,k}\operatorname{Re}((X_{ij}Y_{jk})Z_{ki})$. By Lemma E.14(4),
$\operatorname{Re}\operatorname{tr}(AB)=\operatorname{Re}\operatorname{tr}(BA)$, $T(X,Y,Z)=T(Z,X,Y)$, and, taking conjugates
($\operatorname{Re}\bar p=\operatorname{Re}p$, $\overline{(pq)r}=\bar r(\bar q\bar p)$) and using that $X,Y,Z$ are Hermitian,
$T(X,Y,Z)=\sum\operatorname{Re}(Z_{ik}(Y_{kj}X_{ji}))=T(Z,Y,X)$. So $T$ is symmetric, and
$(X\circ Y|Z)=\frac12(T(X,Y,Z)+T(Y,X,Z))=T(X,Y,Z)$ is symmetric in $X,Y,Z$. $\square$

**Lemma E.21 (orientations).** Let $P_{12},P_{13},P_{23}$ be inner product spaces, $P_{ji}=P_{ij}$, and for every
ordered triple $(i,j,k)$ of distinct elements of $\{1,2,3\}$ let $m:P_{ij}\times P_{jk}\to P_{ik}$ be bilinear, such that
$m(x,y)=m(y,x)$ (the right side being the map for $(k,j,i)$) and $\langle m(x,y),z\rangle=\langle y,m(x,z)\rangle$ for
$x\in P_{ij}$, $y\in P_{jk}$, $z\in P_{ik}$. Let $D$ be a composition algebra and $\iota_{ij}:D\to P_{ij}$ ($i<j$) linear
isometric bijections with $m(\iota_{12}a,\iota_{23}b)=\iota_{13}(ab)$; put $\iota_{ji}(a)=\iota_{ij}(\bar a)$. Then
$m(\iota_{ij}a,\iota_{jk}b)=\iota_{ik}(ab)$ for all ordered triples and all $a,b\in D$.

The hypotheses hold (i) for three indices of a Jordan frame of a simple $J$, with the Peirce spaces and $m=2\circ$
(commutativity and Lemma E.12(2)); (ii) in $\operatorname{Herm}_3(D)$, with $P_{ij}=\{a[ij]\}$,
$\langle a[ij],b[ij]\rangle=\langle a,b\rangle$ (well defined, as conjugation is isometric) and $m=2\circ$, since
$m(a[ij],b[jk])=(ab)[ik]$ and, for $c[ik]$, $\langle ab,c\rangle=\langle b,\bar ac\rangle$ (Lemma E.14(2)).

*Proof.* If the claim holds for $(i,j,k)$, it holds for $(k,j,i)$:
$m(\iota_{kj}a,\iota_{ji}b)=m(\iota_{ij}\bar b,\iota_{jk}\bar a)=\iota_{ik}(\bar b\bar a)=\iota_{ik}(\overline{ab})=\iota_{ki}(ab)$.
So it suffices to treat $(2,1,3)$ and $(1,3,2)$, given $(1,2,3)$. For all $c\in D$, using the adjointness, $(1,2,3)$ and
Lemma E.14(2),
$$
\langle m(\iota_{21}a,\iota_{13}b),\iota_{23}c\rangle=\langle\iota_{13}b,m(\iota_{12}\bar a,\iota_{23}c)\rangle=\langle b,\bar ac\rangle=\langle ab,c\rangle=\langle\iota_{23}(ab),\iota_{23}c\rangle,
$$
$$
\langle m(\iota_{13}a,\iota_{32}b),\iota_{12}c\rangle=\langle\iota_{13}a,m(\iota_{12}c,\iota_{23}\bar b)\rangle=\langle a,c\bar b\rangle=\langle ab,c\rangle=\langle\iota_{12}(ab),\iota_{12}c\rangle,
$$
and $\iota_{23},\iota_{12}$ are onto. (In the second line the adjointness is used for $y=\iota_{32}b$ in the role of $x$,
after $m(x,y)=m(y,x)$.) $\square$

**Theorem E.22 (coordinatisation).** Let $J$ be simple with a Jordan frame $(c_1,\dots,c_r)$, $r\ge3$. Then there are a
composition algebra $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$ of dimension $d$ and an algebra isomorphism
$\Phi:\operatorname{Herm}_r(D)\to J$ with $\Phi(E_{ii})=c_i$. If $r\ge4$, $D$ is associative, i.e. $D\ne\mathbb O$.

*Proof.* *Step 1 (one triangle).* Choose $e_{12}\in J_{12}$, $e_{23}\in J_{23}$ with $N=1$ (Proposition E.11). On
$D'=J_{13}$ define $p\cdot q=m(m(p,e_{23}),m(e_{12},q))$, using the products $J_{13}\times J_{32}\to J_{12}$,
$J_{21}\times J_{13}\to J_{23}$ and $J_{12}\times J_{23}\to J_{13}$. By Lemma E.12(1), $1'=m(e_{12},e_{23})$ is a unit:
$p\cdot1'=m(m(p,e_{23}),e_{23})=p$ and $1'\cdot q=m(e_{12},m(e_{12},q))=q$; and $N(p\cdot q)=N(p)N(q)$. So $(D',N)$ is a
composition algebra, and by Theorem E.17 there is an isomorphism $h:D\to D'$ with $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$.
Put $\iota_{13}=h$, $\iota_{12}=m(h(\cdot),e_{23})$ and $\iota_{23}=m(e_{12},h(\cdot))$; these are isometric bijections
(Lemma E.12(1)) and $m(\iota_{12}a,\iota_{23}b)=h(a)\cdot h(b)=\iota_{13}(ab)$. By Lemma E.21, for $k=3$:

> $(\mathrm T_k)$ there are isometric bijections $\iota_{ij}:D\to J_{ij}$, $i<j\le k$, such that, with $\iota_{ji}(a)=\iota_{ij}(\bar a)$,
> $m(\iota_{ij}a,\iota_{jl}b)=\iota_{il}(ab)$ for all distinct $i,j,l\le k$ and all $a,b\in D$.

*Step 2 (adding an index).* Assume $(\mathrm T_{k-1})$ with $4\le k\le r$. Choose $v\in J_{1k}$ with $N(v)=1$ and put
$\iota_{ik}(a)=m(\iota_{i1}a,v)$ for $2\le i<k$ and $\iota_{1k}(a)=m(\iota_{12}a,\iota_{2k}1)$, isometric bijections by
Lemma E.12(1). Using Lemma E.12(3) (four distinct indices each time) and $(\mathrm T_{k-1})$: for distinct
$i,j\in\{2,\dots,k-1\}$,
$$
m(\iota_{ij}a,\iota_{jk}b)=m(\iota_{ij}a,m(\iota_{j1}b,v))=m(m(\iota_{ij}a,\iota_{j1}b),v)=m(\iota_{i1}(ab),v)=\iota_{ik}(ab).
$$
In particular $\iota_{2k}(b)=m(\iota_{23}b,\iota_{3k}1)$, and hence
$$
m(\iota_{12}a,\iota_{2k}b)=m(m(\iota_{12}a,\iota_{23}b),\iota_{3k}1)=m(\iota_{13}(ab),\iota_{3k}1)
=m(m(\iota_{12}(ab),\iota_{23}1),\iota_{3k}1)=m(\iota_{12}(ab),\iota_{2k}1)=\iota_{1k}(ab).
$$
For $3\le j<k$, $\iota_{1j}a=m(\iota_{12}a,\iota_{2j}1)$, so
$m(\iota_{1j}a,\iota_{jk}b)=m(\iota_{12}a,m(\iota_{2j}1,\iota_{jk}b))=m(\iota_{12}a,\iota_{2k}b)=\iota_{1k}(ab)$. So the relation
holds for all triples $i<j<k$, and Lemma E.21, applied to each triangle $\{i,j,k\}$, gives $(\mathrm T_k)$.

*Step 3 (associativity).* If $r\ge4$, by $(\mathrm T_4)$ and Lemma E.12(3),
$\iota_{14}((ab)c)=m(m(\iota_{12}a,\iota_{23}b),\iota_{34}c)=m(\iota_{12}a,m(\iota_{23}b,\iota_{34}c))=\iota_{14}(a(bc))$.

*Step 4 (the isomorphism).* By $(\mathrm T_r)$, define the linear map
$\Phi(\sum_i\alpha_iE_{ii}+\sum_{i<j}a_{ij}[ij])=\sum_i\alpha_ic_i+\sum_{i<j}\iota_{ij}(a_{ij})$; it is bijective by
Theorem E.8(1), and $\Phi(a[ij])=\iota_{ij}(a)$ for all $i\ne j$. By Theorem E.8 and $(\mathrm T_r)$,
$c_i\circ c_j=\delta_{ij}c_i$, $c_k\circ\iota_{ij}(a)=\frac12(\delta_{ki}+\delta_{kj})\iota_{ij}(a)$,
$\iota_{ij}(a)\circ\iota_{ij}(b)=\langle a,b\rangle(c_i+c_j)$, $\iota_{ij}(a)\circ\iota_{jk}(b)=\frac12\iota_{ik}(ab)$ and
$\iota_{ij}(a)\circ\iota_{kl}(b)=0$: the table of Lemma E.20. By bilinearity $\Phi(X\circ Y)=\Phi(X)\circ\Phi(Y)$. $\square$

## E.8 The Albert algebra

**Lemma E.23 (automorphisms of $\operatorname{Herm}_3(D)$).** Let $D$ be a composition algebra and $f_{12},f_{23},f_{13}$
isometric linear bijections of $D$ with $f_{13}(ab)=f_{12}(a)f_{23}(b)$ for all $a,b$. Then
$$
\Phi\Big(\sum_i\alpha_iE_{ii}+\sum_{i<j}a_{ij}[ij]\Big)=\sum_i\alpha_iE_{ii}+\sum_{i<j}f_{ij}(a_{ij})[ij]
$$
is an automorphism of $(\operatorname{Herm}_3(D),\circ)$. For every $u\in D$ with $N(u)=1$,
$(f_{12},f_{23},f_{13})=(L_u,R_u,B_u)$ with $L_ua=ua$, $R_ua=au$, $B_ua=uau$ is such a triple.

*Proof.* Apply Lemma E.21 in example (ii) to $\iota_{ij}(a)=f_{ij}(a)[ij]$: then $\Phi(a[ij])=\iota_{ij}(a)$ for all $i\ne j$,
and $2\iota_{ij}(a)\circ\iota_{jk}(b)=\iota_{ik}(ab)$ for all distinct $i,j,k$. As $\Phi$ fixes the $E_{ii}$ and the $f_{ij}$
are isometric, $\Phi$ preserves every product in the table of Lemma E.20, hence $\circ$. For the triple, $L_u,R_u,B_u$
are isometric ($N(ua)=N(a)$), bijective, and $u(ab)u=(ua)(bu)$ is the Moufang identity (Lemma E.18(2)). $\square$

**Theorem E.24.** (1) For $D\in\{\mathbb R,\mathbb C,\mathbb H\}$ and $n\ge1$, $\operatorname{Herm}_n(D)$ is a Euclidean Jordan
algebra, and $\operatorname{Herm}_n(\mathbb H)$ is isomorphic to a Jordan subalgebra of $\operatorname{Herm}_{2n}(\mathbb C)$.
(2) The **Albert algebra** $\operatorname{Herm}_3(\mathbb O)$ (dimension $27$) and $\operatorname{Herm}_2(\mathbb O)$ are
Euclidean Jordan algebras. (3) In each of these, $(E_{11},\dots,E_{nn})$ is a Jordan frame with Peirce spaces
$J_{ij}=\{a[ij]:a\in D\}$ of dimension $d=\dim D$; the algebra is simple; $\operatorname{Herm}_1(D)\cong\mathbb R$ and
$\operatorname{Herm}_2(D)\cong V_{d+1}$ for $D=\mathbb R,\mathbb C,\mathbb H,\mathbb O$.

*Proof.* (1) $M_n(D)$ is an associative real algebra (Lemma E.15(2)), so by Lemma E.20 and Lemma 13.1
$\operatorname{Herm}_n(D)$ is a Jordan algebra, Euclidean by Lemma E.20. Replacing each entry $X_{ij}\in\mathbb H$ by the
$2\times2$ block $M(X_{ij})$ is an injective algebra homomorphism $M_n(\mathbb H)\to M_{2n}(\mathbb C)$ (Lemma E.15(2)) that
maps Hermitian matrices to Hermitian matrices, as $M(\bar q)=M(q)^*$; so it preserves $\circ$.

(2) Let $J(x,y)=(x^2\circ y)\circ x-x^2\circ(y\circ x)$ on $\operatorname{Herm}_3(\mathbb O)$; an automorphism $\Phi$ satisfies
$J(\Phi x,\Phi y)=\Phi J(x,y)$. Let $x=\sum_i\alpha_iE_{ii}+p[12]+q[23]+s[13]$. We compose three automorphisms of the
form of Lemma E.23 so that the entries at $(1,3)$ and $(2,3)$ become real. (a) If $s\notin\mathbb R$, write
$s=|s|(\cos\theta+v\sin\theta)$ with $v\perp1$, $N(v)=1$, so $v^2=-1$ (Lemma E.14(3)); for
$u=\cos\frac\theta2-v\sin\frac\theta2$, computing in the commutative subalgebra $\operatorname{span}\{1,v\}\cong\mathbb C$,
$B_us=u^2s=|s|$. (b) Now $s$ is real. Choose $u\perp1$ with $N(u)=1$ and $u\perp q$ (possible as $\dim\mathbb O=8$);
then $B_us=su^2=-s$ is real and $\operatorname{Re}(qu)=\langle u,\bar q\rangle=0$, as $\bar q\in\operatorname{span}\{1,q\}$.
(c) Now $q\perp1$; if $q\ne0$, take $u=q/|q|$: $qu=q^2/|q|=-|q|$ and $B_us=-s$ are real. So there is an automorphism
$\Phi$ such that $\Phi x$ has real entries except possibly $p'$ at $(1,2)$. Since $J(\Phi x,\cdot)$ is linear, it suffices
to show $J(\Phi x,y)=0$ for $y=E_{ii}$ and $y=c[ij]$, $c\in\mathbb O$. By Lemma E.18(1), $p'$ and $c$ lie in an associative
subalgebra $Q$; then $\Phi x,y\in\operatorname{Herm}_3(Q)$, a subalgebra of $\operatorname{Herm}_3(\mathbb O)$ inside the
associative algebra $M_3(Q)$, and $J(\Phi x,y)=0$ by Lemma 13.1. Hence $J(x,\Phi^{-1}y)=\Phi^{-1}J(\Phi x,y)=0$ for all $y$:
the Jordan identity holds, and Lemma E.20 makes the algebra Euclidean. The matrices with zero third row and column form
a subalgebra (Lemma E.20) with unit $E_{11}+E_{22}$, which is $\operatorname{Herm}_2(\mathbb O)$.

(3) By Lemma E.20 the $E_{ii}$ are pairwise orthogonal idempotents with sum $1$; $L(E_{kk})$ acts on $E_{jj}$ by
$\delta_{kj}$ and on $a[ij]$ by $\frac12(\delta_{ki}+\delta_{kj})$, so $J_1(E_{ii})=\mathbb RE_{ii}$ and $J_{1/2}(E_{ii})\cap
J_{1/2}(E_{jj})=\{a[ij]\}$. $E_{ii}$ is primitive: if $E_{ii}=g+g'$ with orthogonal non-zero idempotents, then
$g=E_{ii}\circ g\in\mathbb RE_{ii}$, so $g=E_{ii}$ and $g'=0$. The graph of Proposition E.11 is complete, so the algebra is
simple. For $n=1,2$ apply Theorem E.13. $\square$

**Corollary E.25.** For $n\ge4$, $\operatorname{Herm}_n(\mathbb O)$ is not a Jordan algebra.

*Proof.* If it were, Lemma E.10(2) would apply: its proof uses only the Jordan identity and the products of Lemma E.20
(with $c_k=E_{kk}$). Then $\frac14((ab)c)[14]=(a[12]\circ b[23])\circ c[34]=a[12]\circ(b[23]\circ c[34])=\frac14(a(bc))[14]$
for all $a,b,c\in\mathbb O$, contradicting Lemma E.15(2). $\square$

**Theorem E.26 (Albert).** Let $D$ be a composition algebra. If there are an associative real algebra $A$ and an
injective linear map $\varphi:\operatorname{Herm}_3(D)\to A$ with $\varphi(x\circ y)=\frac12(\varphi(x)\varphi(y)+\varphi(y)\varphi(x))$,
then $D$ is associative. In particular the Albert algebra is **exceptional**: it is not isomorphic to a Jordan subalgebra of
any associative algebra, such as $M_N(\mathbb C)$.

*Proof.* Put $e_i=\varphi(E_{ii})$ and $e=e_1+e_2+e_3=\varphi(1)$. Then $e_i^2=e_i$, and $e_ie_j+e_je_i=0$ for $i\ne j$;
multiplying by $e_i$ on the left, respectively on the right, gives $e_ie_j+e_ie_je_i=0=e_ie_je_i+e_je_i$, so
$e_ie_j=e_je_i=0$. For $x=\varphi(X)$, $ex+xe=2x$ gives $exe=ex=xe=x$, so $x=\sum_{p,q}x_{pq}$ with
$x_{pq}=e_pxe_q$. For $X=a[ij]$ the relations $e_i\circ x=e_j\circ x=\frac12x$ and $e_k\circ x=0$ ($k\ne i,j$) say
$\frac12([p=i]+[q=i])x_{pq}=\frac12x_{pq}$ and so on, which leaves only $x_{ij}$ and $x_{ji}$. Define, for ordered pairs,
$\alpha_{ij}(a)=e_i\varphi(a[ij])e_j$. For distinct $i,j,k$, $\varphi(a[ij])\varphi(b[jk])=\alpha_{ij}(a)\alpha_{jk}(b)$ is the
only term in $e_iAe_k$ of $\varphi(a[ij])\varphi(b[jk])+\varphi(b[jk])\varphi(a[ij])=\varphi((ab)[ik])$, so
$$
\alpha_{ij}(a)\,\alpha_{jk}(b)=\alpha_{ik}(ab).
$$
From $\varphi(1[13])^2=\varphi(E_{11}+E_{33})$ we get $\alpha_{13}(1)\alpha_{31}(1)=e_1$, $\alpha_{31}(1)\alpha_{13}(1)=e_3$;
from $\varphi(a[13])^2=N(a)(e_1+e_3)$, $\alpha_{13}(a)\alpha_{31}(\bar a)=N(a)e_1$, so $\alpha_{13}$ is injective. Now
$$
\alpha_{13}(ay)=\alpha_{12}(a)\alpha_{23}(y)=\alpha_{13}(a)\alpha_{32}(1)\alpha_{21}(1)\alpha_{13}(y)=\alpha_{13}(a)\alpha_{31}(1)\alpha_{13}(y),
$$
and therefore, with $\pi=\alpha_{31}(1)$ and associativity in $A$,
$\alpha_{13}((ab)c)=\alpha_{13}(a)\pi\alpha_{13}(b)\pi\alpha_{13}(c)=\alpha_{13}(a(bc))$. So $(ab)c=a(bc)$. By Lemma E.15(2),
$\mathbb O$ is not associative. $\square$

Albert's original argument (A. A. Albert, Ann. of Math. 35 (1934) 65–73) is of this kind. The code checks the Jordan
identity of $\operatorname{Herm}_3(\mathbb O)$ and the associativity and positivity of its trace form on random elements,
shows numerically that the Jordan identity fails in $\operatorname{Herm}_4(\mathbb O)$ (Corollary E.25), checks the
automorphism of Lemma E.23, and runs the three steps (a)–(c) of the proof of Theorem E.24(2) on a random element:

```python
n3, prod3 = herm(3, 8)
X, Y, Z = torch.randn(3, 200, 27, dtype=dt)
X2 = prod3(X, X)
assert (prod3(prod3(X2, Y), X) - prod3(X2, prod3(Y, X))).abs().max() < 1e-10    # Jordan identity, Theorem E.24(2)
n4, prod4 = herm(4, 8)
X4, Y4 = torch.randn(2, 200, n4, dtype=dt)
X42 = prod4(X4, X4)
defect4 = (prod4(prod4(X42, Y4), X4) - prod4(X42, prod4(Y4, X4))).abs().max().item()
print("Herm_4(O): max |Jordan defect| =", defect4)
assert defect4 > 1                                                               # Corollary E.25
retr = lambda v: v[..., :3].sum(-1)               # Re tr
assert torch.allclose(retr(prod3(prod3(X, Y), Z)), retr(prod3(X, prod3(Y, Z))))
sq = (X[:, :3] ** 2).sum(-1) + 2 * (X[:, 3:] ** 2).sum(-1)
assert torch.allclose(retr(prod3(X, X)), sq)                                    # (X|X) = sum of N(X_ij)

def moufang_auto(v, uu):                          # (L_u, R_u, B_u) on the entries a_12, a_23, a_13
    a12, a13, a23 = v[..., 3:11], v[..., 11:19], v[..., 19:27]   # pairs order (0,1),(0,2),(1,2)
    uu = uu.expand_as(a12)
    return torch.cat([v[..., :3], cd(uu, a12), cd(cd(uu, a13), uu), cd(a23, uu)], -1)
uu = torch.randn(8, dtype=dt); uu = uu / uu.norm()
assert torch.allclose(moufang_auto(prod3(X, Y), uu), prod3(moufang_auto(X, uu), moufang_auto(Y, uu)))

x = torch.randn(27, dtype=dt)                     # steps (a)-(c) of the proof of Theorem E.24(2)
s = x[11:19]
im = torch.cat([torch.zeros(1, dtype=dt), s[1:]])
v, th = im / im.norm(), torch.atan2(im.norm(), s[0]).item()
x1 = moufang_auto(x, math.cos(th / 2) * e[0] - math.sin(th / 2) * v)            # (a): (1,3) entry -> |s|
imq = torch.cat([torch.zeros(1, dtype=dt), x1[20:27]])
g = torch.cat([torch.zeros(1, dtype=dt), torch.randn(7, dtype=dt)])
g = g - (g @ imq) / (imq @ imq) * imq                                             # u perpendicular to 1 and q
x2 = moufang_auto(x1, g / g.norm())                                               # (b): (2,3) entry imaginary
x3 = moufang_auto(x2, x2[19:27] / x2[19:27].norm())                              # (c): (2,3) entry real
print("entries (1,3), (2,3) after (a)-(c):", x3[11:19].numpy().round(3), x3[19:27].numpy().round(3))
assert x3[12:19].abs().max() < 1e-12 and x3[20:27].abs().max() < 1e-12
assert abs(abs(x3[11].item()) - s.norm().item()) < 1e-12
```

## E.9 The classification

**Theorem E.27 (Jordan–von Neumann–Wigner).**

1. Every finite-dimensional formally real (equivalently, Euclidean) Jordan algebra is a direct sum of simple ideals,
   and every simple one is isomorphic to $\operatorname{Herm}_n(\mathbb R)$, $\operatorname{Herm}_n(\mathbb C)$ or
   $\operatorname{Herm}_n(\mathbb H)$ for some $n\ge1$, to a spin factor $V_m$ with $m\ge2$, or to the Albert algebra
   $\operatorname{Herm}_3(\mathbb O)$.
2. The algebras $\mathbb R$; $V_m$ ($m\ge2$); $\operatorname{Herm}_n(D)$ ($D=\mathbb R,\mathbb C,\mathbb H$, $n\ge3$);
   $\operatorname{Herm}_3(\mathbb O)$ are simple, formally real and pairwise non-isomorphic, so every simple formally
   real Jordan algebra is isomorphic to exactly one of them. Moreover $\operatorname{Herm}_1(D)\cong\mathbb R$ and
   $\operatorname{Herm}_2(D)\cong V_{\dim D+1}$ for $D=\mathbb R,\mathbb C,\mathbb H,\mathbb O$.
3. In a simple formally real Jordan algebra all Jordan frames have the same number $r$ of elements and all their Peirce
   spaces $J_{ij}$ the same dimension $d$: $r=1$ for $\mathbb R$, $(r,d)=(2,m-1)$ for $V_m$ and $(n,\dim D)$ for
   $\operatorname{Herm}_n(D)$, $n\ge2$.

*Proof.* (1) By Theorem E.9, $J$ is a direct sum of simple formally real ideals. Such an ideal has a Jordan frame of some
size $r$ (Lemma E.6(2)); it is $\cong\mathbb R=\operatorname{Herm}_1(\mathbb R)$ if $r=1$, $\cong V_{d+1}$ if $r=2$
(Theorem E.13), and $\cong\operatorname{Herm}_r(D)$ with $D\ne\mathbb O$ if $r\ge4$ and $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$
if $r=3$ (Theorem E.22).

(2) For $\operatorname{Herm}_n(D)$ see Theorem E.24. $V_m$ is a Jordan algebra by Lemma 13.3 and Euclidean for
$\langle(a,u),(b,v)\rangle=ab+u\cdot v$, since $\langle x\circ y,z\rangle=abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w$ is symmetric
in $x=(a,u)$, $y=(b,v)$, $z=(c,w)$. For a unit vector $w$, $c_{\pm w}=\frac12(1,\pm w)$ is a Jordan frame
(Theorem 13.4(1); $c_{\pm w}$ is primitive because a sum $c_{w'}+c_{-w'}$ of two orthogonal non-zero idempotents is $1$)
with $J_{12}=\{(0,v):v\perp w\}$, of dimension $m-1\ge1$ (the dimensions add up to $m+1$); so $V_m$
is simple by Proposition E.11. To separate the algebras, let $s(J)$ be the largest number of pairwise orthogonal non-zero
idempotents (finite by Lemma E.6(2)); $s$ and $\dim$ are invariant under isomorphism. Clearly $s(\mathbb R)=1$. In
$V_m$ the non-zero idempotents are $1$ and the $c_w$, with $1\circ c_w=c_w$ and $c_w\circ c_{w'}=0$ only for $w'=-w$
(Theorem 13.4(1)); so $s(V_m)=2$. For $D\in\{\mathbb R,\mathbb C,\mathbb H\}$, $s(\operatorname{Herm}_n(D))=n$: the
$E_{ii}$ give $s\ge n$. Conversely, orthogonal idempotents $p_k\ne0$ satisfy $p_kp_l=0$ for $k\ne l$ (as in the proof of
Theorem E.26), so the right $D$-submodules $V_k=p_kD^n$ of $D^n$ are non-zero and independent: if $\sum_kv_k=0$ with
$v_k\in V_k$, then $v_l=p_l\sum_kv_k=0$. A non-zero $V_k$ contains some $v\ne0$ and then $vD$, of real dimension $d$
($D$ has no zero divisors); so $sd\le nd$. For $\operatorname{Herm}_3(\mathbb O)$, the $E_{ii}$ give $s\ge3$; and $s$
pairwise orthogonal non-zero idempotents, together with $1$ minus their sum, refine as in the proof of Theorem E.7 to a
Jordan frame of $r'\ge s$ elements. By Proposition E.11(2) and Theorem E.22, $27=r'+d'\frac{r'(r'-1)}2$ with $d'\ge1$, and
$d'\in\{1,2,4\}$ if $r'\ge4$. For $r'=4,5,6$ the right side lies in $\{10,16,28\}$, $\{15,25,45\}$, $\{21,36,66\}$, and
for $r'\ge7$ it exceeds $27$; so $s\le r'=3$. The pairs $(s,\dim)$ are $(1,1)$, $(2,m+1)$, $(n,\frac{n(n+1)}2)$,
$(n,n^2)$, $(n,n(2n-1))$ for $n\ge3$, and $(3,27)$, all different, as $\frac{n(n+1)}2<n^2<n(2n-1)$ and $6,9,15\ne27$.
The isomorphisms in rank $1$ and $2$ are Theorem E.24(3).

(3) If $F$ is a Jordan frame of $r_F$ elements of a simple $J$, the proof of (1), applied with $F$, gives an isomorphism of
$J$ with a model of rank $r_F$ in (2), and $s$ of that model is $r_F$. So $r_F=s(J)$ for every frame, and then
$d=(\dim J-r)/\binom r2$ for $r\ge2$ (Proposition E.11(2)). $\square$

**Corollary E.28 (rank and trace).** In every formally real Jordan algebra all Jordan frames have the same number of
elements, the **rank**. If $J$ is simple of rank $r$, then $\operatorname{tr}x=\frac r{\dim J}\operatorname{Tr}L(x)$ satisfies
$\operatorname{tr}c=1$ for every primitive idempotent $c$, hence $\operatorname{tr}(\sum_i\lambda_ic_i)=\sum_i\lambda_i$ for
every decomposition of Theorem E.7, and $\operatorname{tr}(x\circ y)$ is an associative inner product.

*Proof.* Let $J=\bigoplus_pJ^{(p)}$ as in Theorem E.9 and $c$ primitive. Its components $c^{(p)}\in J^{(p)}$ are
pairwise orthogonal idempotents ($c^2=\sum_p(c^{(p)})^2$), so exactly one is non-zero. Hence the members of a Jordan frame
of $J$ that lie in $J^{(p)}$ form a Jordan frame of $J^{(p)}$ (they add up to the component $1^{(p)}$ of $1$), and the rank
of $J$ is the sum of the ranks of the $J^{(p)}$ (Theorem E.27(3)). In a simple $J$, a primitive $c$ belongs to a Jordan frame
(refine $1-c$ as in the proof of Theorem E.7), so $\operatorname{Tr}L(c)=t=\dim J/r$ by Proposition E.11(2) (for $r=1$,
$J=\mathbb Rc$ and $\operatorname{Tr}L(c)=1$). Finally $\operatorname{tr}(x\circ y)=\frac r{\dim J}\tau(x,y)$ (Theorem E.5). $\square$

The code computes, for every type, the Jordan frame $(E_{ii})$ or $(c_w,c_{-w})$, the dimensions of the Peirce spaces as
common eigenspaces of the $L(c_i)$ (Theorem E.8), the dimension formula of Proposition E.11(2) and the positivity of the trace
form $\tau$ (Theorem E.5):

```python
def peirce(C, frame):
    """dimensions of the joint Peirce spaces J_ij = {x : c_i o x = c_j o x = x/2}"""
    I = torch.eye(C.shape[0], dtype=dt)
    Ls = [Lop(C, f) for f in frame]
    return [C.shape[0] - torch.linalg.matrix_rank(torch.cat([Ls[i] - I / 2, Ls[j] - I / 2])).item()
            for i in range(len(frame)) for j in range(i + 1, len(frame))]

cases = [(f"Herm_{n}({'RCHO'[k]})", herm(n, 2 ** k), n, 2 ** k) for n, k in ((2, 0), (3, 0), (4, 0), (3, 1), (4, 1), (3, 2), (2, 3), (3, 3))]
cases += [(f"V_{m}", spin(m), 2, m - 1) for m in (2, 3, 5)]
for name, (dim, prod), r, d in cases:
    C = structure(dim, prod)
    if name.startswith('V'):
        frame = [0.5 * torch.eye(dim, dtype=dt)[0] + s * 0.5 * torch.eye(dim, dtype=dt)[1] for s in (1, -1)]
    else:
        frame = list(torch.eye(dim, dtype=dt)[:r])
    for i, f in enumerate(frame):
        assert torch.allclose(mul(C, f, f), f) and all(torch.allclose(mul(C, f, g), torch.zeros(dim, dtype=dt)) for g in frame[i + 1:])
    dims = peirce(C, frame)
    t = torch.einsum('cbb->c', Lop(C, torch.eye(dim, dtype=dt)))
    gram = torch.einsum('abc,c->ab', C, t)
    print(f"{name:10s} dim {dim:2d}  rank {r}  Peirce dimensions {dims}  min eig of trace form {torch.linalg.eigvalsh(gram).min().item():.2f}")
    assert dims == [d] * (r * (r - 1) // 2) and dim == r + d * r * (r - 1) // 2
    assert torch.linalg.eigvalsh(gram).min() > 0
```

## E.10 What this appendix proves for Chapter 13

| statement in Chapter 13 | proved here |
|---|---|
| formally real $\Leftrightarrow$ Euclidean (Section 13.4) | Theorem E.5 |
| Theorem 13.7: decomposition into simple ideals, the list, "exactly one" | Theorems E.9, E.27 |
| $\operatorname{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}$ for $\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O$ | Theorem E.24(3) |
| $\operatorname{Herm}_r(\mathbb H)$ embeds into $\operatorname{Herm}_{2r}(\mathbb C)$ | Theorem E.24(1) |
| the Albert algebra is a Jordan algebra and is exceptional | Theorems E.24(2), E.26 |
| $\operatorname{Herm}_r(\mathbb O)$ is a Jordan algebra only for $r\le3$ | Theorem E.24(2), Corollary E.25 |
| spectral theorem with Jordan frames, rank, $\operatorname{tr}z=\sum_i\mu_i$, associativity of $\operatorname{tr}(x\circ y)$ (Section 13.6) | Theorem E.7, Corollary E.28 |

**Notes.** The classification is due to P. Jordan, J. von Neumann and E. Wigner, *On an algebraic generalization of the
quantum mechanical formalism*, Ann. of Math. 35 (1934) 29–64, and the exceptionality of $\operatorname{Herm}_3(\mathbb O)$
to A. A. Albert (same volume, 65–73). The proof above follows the structure of Faraut–Korányi, *Analysis on Symmetric
Cones* (Oxford 1994), Chapters II–V; the doubling proof of Hurwitz's theorem (Nachr. Ges. Wiss. Göttingen 1898, 309–316) is as in
J. H. Conway, D. A. Smith, *On Quaternions and Octonions* (2003) and T. A. Springer, F. D. Veldkamp, *Octonions, Jordan
Algebras and Exceptional Groups* (2000), and the proof of Frobenius' theorem (1878) is R. S. Palais', Amer. Math. Monthly
75 (1968) 366–368. Further reading (not used anywhere in the course): C. M. Glennie, Pacific J. Math. 16 (1966)
47–59, for a different approach to exceptionality through polynomial identities, and McCrimmon's *A Taste of Jordan
Algebras* (2004) for the history and the general theory.
