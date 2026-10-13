**English** · [Español](es/E_jordan_classification.md) · [Français](fr/E_jordan_classification.md) · [Deutsch](de/E_jordan_classification.md) · [中文](zh/E_jordan_classification.md) · [日本語](ja/E_jordan_classification.md) · [Русский](ru/E_jordan_classification.md)

# Appendix E. Classification of Euclidean Jordan algebras

Chapter 13 used the Jordan–von Neumann–Wigner classification of finite-dimensional formally real Jordan algebras.
This appendix proves it completely, from the Jordan identity and finite-dimensional linear algebra only. The route is
the classical one (J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, Chapters II–V; K. McCrimmon, *A Taste of
Jordan Algebras*): power associativity and a spectral lemma give idempotents and an associative positive trace form;
the Peirce decomposition with respect to a Jordan frame splits the algebra into pieces $J_{ij}$ whose products obey
strict rules; in a simple algebra of rank $r\ge3$ these rules turn $J_{12}$ into a composition algebra, which by
Hurwitz's theorem is $\mathbb R,\mathbb C,\mathbb H$ or $\mathbb O$, and the whole algebra is then
$\operatorname{Herm}_r$ of it. Section E.6 also proves Frobenius' theorem on real associative division algebras, and
Section E.8 proves that the Albert algebra $\operatorname{Herm}_3(\mathbb O)$ is a Jordan algebra that does not embed
into any associative algebra.

**You will learn** why power associativity, positivity of the trace form and the existence of a unit follow from
formal reality alone; the Peirce rules and how they force the coordinates of a simple algebra to multiply like
numbers; the proof of Hurwitz's theorem by Cayley–Dickson doubling; and why the octonions occur only in rank $3$.

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
$c'\in J_0(c)$, so $L(c)$ and $L(c')$ commute by (2) and are simultaneously diagonalisable. On a common eigenvector
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
assert torch.allclose(lhs, comm(Lop(C3, w), comm(Lop(C3, x), Lop(C3, y))))     # Lemma E.1(d)

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

c = cs[0] + cs[1]                                 # an idempotent of "rank 2": Proposition E.4
Lc = Lop(C3, c)
assert torch.allclose(Lc @ (2 * Lc - torch.eye(9, dtype=dt)) @ (Lc - torch.eye(9, dtype=dt)), torch.zeros(9, 9, dtype=dt), atol=1e-10)
ev = torch.linalg.eigvals(Lc).real.sort().values
print("eigenvalues of L(c):", ev.numpy().round(6))
tau = lambda C: torch.einsum('abc,c->ab', C, torch.einsum('cbb->c', Lop(C, torch.eye(C.shape[0], dtype=dt))))  # Tr L(e_a o e_b)
T3 = tau(C3)
assert torch.linalg.eigvalsh(T3).min() > 0                                     # Theorem E.5(b)
assert torch.allclose(mul(C3, x, y) @ T3 @ w, x @ T3 @ mul(C3, y, w))         # Theorem E.5(a)

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

