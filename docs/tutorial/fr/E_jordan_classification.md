[English](../E_jordan_classification.md) · [Español](../es/E_jordan_classification.md) · **Français** · [Deutsch](../de/E_jordan_classification.md) · [中文](../zh/E_jordan_classification.md) · [日本語](../ja/E_jordan_classification.md) · [Русский](../ru/E_jordan_classification.md)

# Annexe E. Classification des algèbres de Jordan euclidiennes

Le chapitre 13 a utilisé la classification de Jordan–von Neumann–Wigner des algèbres de Jordan formellement réelles de
dimension finie. Cette annexe la démontre entièrement, à partir de la seule identité de Jordan et de l'algèbre linéaire
en dimension finie. La voie suivie est la voie classique (J. Faraut, A. Korányi, *Analysis on Symmetric Cones*,
chapitres II–V ; K. McCrimmon, *A Taste of Jordan Algebras*) : l'associativité des puissances et un lemme spectral
fournissent des idempotents et une forme trace associative et positive ; la décomposition de Peirce relative à un repère
de Jordan découpe l'algèbre en morceaux $J_{ij}$ dont les produits obéissent à des règles strictes ; dans une algèbre
simple de rang $r\ge3$, ces règles font de l'espace de Peirce $J_{13}$ une algèbre de composition, qui d'après le
théorème de Hurwitz est $\mathbb R,\mathbb C,\mathbb H$ ou $\mathbb O$, et l'algèbre tout entière est alors le
$\operatorname{Herm}_r$ de celle-ci. La section E.6 démontre en outre le théorème de Frobenius sur les algèbres à
division réelles associatives, et la section E.8 démontre que l'algèbre d'Albert $\operatorname{Herm}_3(\mathbb O)$ est
une algèbre de Jordan qui ne se plonge dans aucune algèbre associative.

**Vous apprendrez** pourquoi l'associativité des puissances, la positivité de la forme trace et l'existence d'une unité
découlent de la seule réalité formelle ; les règles de Peirce et comment elles contraignent les coordonnées d'une algèbre
simple à se multiplier comme des nombres ; la démonstration du théorème de Hurwitz par le doublement de Cayley–Dickson ;
et pourquoi les octonions ne donnent une nouvelle algèbre simple qu'en rang $3$.

**Conventions.** Toutes les algèbres sont des espaces vectoriels réels de dimension finie munis d'un produit bilinéaire.
Une **algèbre de Jordan** est définie comme à la section 13.1 : $x\circ y=y\circ x$ et
$(x^2\circ y)\circ x=x^2\circ(y\circ x)$, $x^2=x\circ x$ ; on ne suppose pas l'existence d'une unité. $L(x)y=x\circ y$,
$[A,B]=AB-BA$, et $\operatorname{Tr}$ est la trace d'une application linéaire. Une **sous-algèbre** est un sous-espace
stable par $\circ$, un **idéal** un sous-espace $I$ tel que $J\circ I\subseteq I$, et $J$ est **simple** si $J\ne0$ et
si ses seuls idéaux sont $0$ et $J$. Un **idempotent** est un $c$ tel que $c^2=c$ ; des idempotents $c,c'$ sont
**orthogonaux** si $c\circ c'=0$. $J$ est **formellement réelle** si $\sum_ix_i^2=0$ entraîne que tous les $x_i=0$, et
**euclidienne** si elle possède une unité $1$ et un produit scalaire tel que
$\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$. Les identités polynomiales sont **linéarisées** ainsi : si
$p(t)=\sum_kp_kt^k$ a des coefficients vectoriels et si $p(t)=0$ pour tout $t$ réel, alors tous les $p_k=0$ (on évalue
en $\deg p+1$ points distincts et on inverse la matrice de Vandermonde).

## E.1 Identités de Jordan linéarisées et associativité des puissances

**Lemme E.1.** Dans toute algèbre de Jordan, pour tous $x,y,z,w$ :

1. $[L(x),L(x^2)]=0$ et $[L(x^2),L(z)]+2[L(x\circ z),L(x)]=0$ ;
2. $L(x^2\circ y)=2L(x\circ y)L(x)+L(x^2)L(y)-2L(x)L(y)L(x)$ ;
3. $L((x\circ w)\circ y)=L(w\circ y)L(x)+L(x\circ y)L(w)+L(x\circ w)L(y)-L(w)L(y)L(x)-L(x)L(y)L(w)$ ;
4. $L((x\circ w)\circ y)-L(x\circ(w\circ y))=[L(w),[L(x),L(y)]]$ ; en particulier
   $\operatorname{Tr}L((x\circ w)\circ y)=\operatorname{Tr}L(x\circ(w\circ y))$.

*Démonstration.* (1) L'identité de Jordan dit que $L(x)L(x^2)y=L(x^2)L(x)y$. Remplaçons $x$ par $x+tz$ ; le
coefficient de $t$ est $[L(z),L(x^2)]+2[L(x),L(x\circ z)]=0$. (2) Remplaçons $x$ par $x+tz$ dans l'identité de Jordan
elle-même ; le coefficient de $t$ est
$$
2((x\circ z)\circ y)\circ x+(x^2\circ y)\circ z=2(x\circ z)\circ(y\circ x)+x^2\circ(y\circ z).
$$
Lus comme applications linéaires de $z$, les quatre termes sont $2L(x)L(y)L(x)z$, $L(x^2\circ y)z$, $2L(x\circ y)L(x)z$
et $L(x^2)L(y)z$. (3) On remplace $x$ par $x+sw$ dans (2) et on prend le coefficient de $s$, divisé par $2$. (4) On
écrit (3) une fois tel quel et une fois en échangeant $x$ et $y$. Les trois premiers termes des membres de droite
coïncident, et $(y\circ w)\circ x=x\circ(w\circ y)$, donc la différence est
$-L(w)L(y)L(x)-L(x)L(y)L(w)+L(w)L(x)L(y)+L(y)L(x)L(w)=[L(w),[L(x),L(y)]]$, un commutateur, dont la trace est $0$.
$\square$

**Lemme E.2 (associativité des puissances).** Posons $x^1=x$ et $x^{k+1}=x\circ x^k$. Alors $x^p\circ x^q=x^{p+q}$ pour
tous $p,q\ge1$, et les opérateurs $L(x^k)$ commutent deux à deux. Par conséquent, pour des polynômes réels $f,g$ sans
terme constant et $f(x):=\sum_kf_kx^k$, on a $(fg)(x)=f(x)\circ g(x)$, et $\mathbb R[x]=\{f(x)\}$ est une sous-algèbre
associative et commutative.

*Démonstration.* $L(x)$ et $L(x^2)$ commutent d'après le lemme E.1(1) ; soit $\mathcal A$ l'algèbre commutative
d'opérateurs qu'ils engendrent. Montrons par récurrence sur $n\ge2$ que $x^p\circ x^q=x^{p+q}$ dès que $p+q\le n$, et
que $L(x^k)\in\mathcal A$ pour $k\le n$. Pour $n=2$, c'est la définition. Supposons-le pour $n$, et soit $p+q=n+1$. Si
$q=1$, $x^p\circ x=x^{n+1}$ par commutativité. Si $q\ge2$, alors
$x^p\circ x^q=L(x^p)L(x)x^{q-1}=L(x)L(x^p)x^{q-1}=L(x)x^{n}=x^{n+1}$. Ensuite, le lemme E.1(2) avec $y=x^{n-1}$ et
$x^2\circ x^{n-1}=x^{n+1}$ donne $L(x^{n+1})=2L(x^n)L(x)+L(x^2)L(x^{n-1})-2L(x)^2L(x^{n-1})\in\mathcal A$. Le reste
découle de la bilinéarité. $\square$

## E.2 Idempotents et décomposition de Peirce

**Lemme E.3 (lemme spectral).** Soient $J$ formellement réelle et $x\in J$. Il existe $k\ge0$, des réels non nuls
distincts $\lambda_1,\dots,\lambda_k$ et des idempotents non nuls deux à deux orthogonaux $c_1,\dots,c_k\in\mathbb R[x]$
tels que
$$
x=\sum_{i=1}^k\lambda_ic_i,\qquad x^n=\sum_{i=1}^k\lambda_i^nc_i\quad(n\ge1),\qquad\mathbb R[x]=\operatorname{span}\{c_1,\dots,c_k\}.
$$
En particulier, toute algèbre de Jordan formellement réelle non nulle contient un idempotent non nul.

*Démonstration.* La réalité formelle avec un seul terme dit que $z^2=0$ entraîne $z=0$. Soit
$K=\{q\in\mathbb R[t]:(tq)(x)=0\}$. D'après le lemme E.2, $K$ est un idéal de $\mathbb R[t]$ : si $q\in K$ et
$s=s_0+s_1$ avec $s_1(0)=0$, alors $(tqs)(x)=s_0(tq)(x)+(tq)(x)\circ s_1(x)=0$. Comme $x,x^2,\dots,x^{\dim J+1}$ sont
linéairement dépendants, $K\ne0$, donc $K=\mu\,\mathbb R[t]$ pour un $\mu$ unitaire. (i) $\mu(0)\neq0$ : si $\mu=th$,
alors $w=(th)(x)$ vérifie $w^2=(t\cdot\mu h)(x)=0$, car $\mu h\in K$ ; donc $w=0$, $h\in K$ et $\mu\mid h$, ce qui est
impossible pour des raisons de degré. (ii) $\mu$ n'a pas de facteur répété : si $\mu=g^2h$ avec $\deg g\ge1$, alors
$w=(tgh)(x)$ vérifie $w^2=(t\cdot\mu\cdot th)(x)=0$, donc $gh\in K$, impossible pour des raisons de degré. (iii) $\mu$
n'a pas de racine non réelle : si $\mu=((t-\alpha)^2+\beta^2)h$ avec $\beta\ne0$, alors $a=(t(t-\alpha)h)(x)$ et
$b=(\beta th)(x)$ vérifient $a^2+b^2=(t\cdot\mu\cdot th)(x)=0$, donc $b=0$ et $h\in K$, impossible. Ainsi
$\mu=\prod_{i=1}^k(t-\lambda_i)$ avec des $\lambda_i$ réels, non nuls et distincts. Posons
$$
p_i(t)=\frac t{\lambda_i}\prod_{j\ne i}\frac{t-\lambda_j}{\lambda_i-\lambda_j},\qquad c_i=p_i(x).
$$
Le polynôme $p_ip_l-\delta_{il}p_i$ est $t$ fois un polynôme qui s'annule en chaque $\lambda_j$, donc $t$ fois un élément
de $K$ ; ainsi $c_i\circ c_l=\delta_{il}c_i$. De même $t-\sum_i\lambda_ip_i=t\big(1-\sum_i\lambda_ip_i/t\big)$ et le
facteur entre parenthèses s'annule en chaque $\lambda_j$, donc $x=\sum_i\lambda_ic_i$. Si $c_i=0$, alors $p_i/t\in K$,
de degré $k-1<\deg\mu$ ; donc $c_i\ne0$. La formule pour $x^n$ découle de l'orthogonalité, et elle montre que
$\mathbb R[x]\subseteq\operatorname{span}\{c_i\}$.
$\square$

**Proposition E.4 (décomposition de Peirce relative à un idempotent).** Soit $c$ un idempotent d'une algèbre de Jordan
$J$.

1. $L(c)(2L(c)-1)(L(c)-1)=0$. Par conséquent $J=J_1(c)\oplus J_{1/2}(c)\oplus J_0(c)$ avec
   $J_\alpha(c)=\{x:c\circ x=\alpha x\}$, et $c\in J_1(c)$.
2. $[L(c),L(z)]=0$ pour $z\in J_0(c)\cup J_1(c)$, et $[L(c),L(y^2)]=0$ pour $y\in J_{1/2}(c)$.
3. $J_1\circ J_1\subseteq J_1$, $J_0\circ J_0\subseteq J_0$, $J_1\circ J_0=0$, $(J_1+J_0)\circ J_{1/2}\subseteq J_{1/2}$
   et $J_{1/2}\circ J_{1/2}\subseteq J_1+J_0$ (le tout relativement à $c$).
4. Si $c'$ est un idempotent orthogonal à $c$, alors $c+c'$ est un idempotent, $[L(c),L(c')]=0$ et
   $J_1(c')\subseteq J_0(c)$.

*Démonstration.* (1) Le lemme E.1(2) avec $x=y=c$ donne $L(c)=3L(c)^2-2L(c)^3$. Une application linéaire annulée par
un polynôme à racines simples est diagonalisable, de valeurs propres parmi ces racines. (2) Le lemme E.1(1) avec $x=c$
donne $[L(c),L(z)]=2[L(c),L(c\circ z)]$, c'est-à-dire $(1-2\alpha)[L(c),L(z)]=0$ pour $z\in J_\alpha(c)$. Le lemme
E.1(1) avec $x=y$ et $z=c$ donne $[L(y^2),L(c)]+2[L(y\circ c),L(y)]=[L(y^2),L(c)]+[L(y),L(y)]=[L(y^2),L(c)]=0$.
(3) D'après (2), $L(z)$ envoie chaque sous-espace propre de $L(c)$ dans lui-même si $z\in J_0\cup J_1$. Cela donne la
première, la deuxième et la quatrième règle, ainsi que $J_1\circ J_0\subseteq J_1\cap J_0=0$. Pour $y\in J_{1/2}$,
$L(y^2)$ commute avec $L(c)$, donc $c\circ y^2=y^2\circ c$ est un vecteur propre de $L(c)$ pour la valeur propre $1$,
comme $c$ ; en écrivant $y^2=a_1+a_{1/2}+a_0$, on voit que $c\circ y^2=a_1+\frac12a_{1/2}$ est dans $J_1$, donc
$a_{1/2}=0$. La polarisation ($y\to y+y'$) donne la dernière règle. (4) $(c+c')^2=c+2c\circ c'+c'=c+c'$, et
$c'\in J_0(c)$, donc $L(c)$ et $L(c')$ commutent d'après (2) et sont simultanément diagonalisables ; en particulier
$J_1(c')$ est engendré par des vecteurs propres communs. Sur un vecteur propre commun de valeurs propres $\alpha$ pour
$L(c)$ et $1$ pour $L(c')$, $L(c+c')$ a la valeur propre $\alpha+1\in\{0,\frac12,1\}$, donc $\alpha=0$. $\square$

## E.3 La forme trace ; euclidienne = formellement réelle

**Théorème E.5.** Soient $J$ formellement réelle et $\tau(x,y)=\operatorname{Tr}L(x\circ y)$.

1. $\tau$ est symétrique, bilinéaire et **associative** : $\tau(x\circ y,z)=\tau(x,y\circ z)$.
2. $\tau$ est définie positive ; pour tout idempotent non nul $c$, $\tau(c,c)=\operatorname{Tr}L(c)\ge1$.
3. $J$ possède un élément unité.
4. $J$ est euclidienne, pour le produit scalaire $\tau$. Réciproquement, toute algèbre de Jordan euclidienne est
   formellement réelle.

*Démonstration.* (1) Lemme E.1(4). (2) D'après la proposition E.4(1),
$\operatorname{Tr}L(c)=\dim J_1(c)+\frac12\dim J_{1/2}(c)\ge1$, car $0\ne c\in J_1(c)$. Pour $x=\sum_i\lambda_ic_i$ comme
dans le lemme E.3, $x^2=\sum_i\lambda_i^2c_i$ et
$\tau(x,x)=\operatorname{Tr}L(x^2)=\sum_i\lambda_i^2\operatorname{Tr}L(c_i)$, qui est $>0$ sauf si $k=0$, c'est-à-dire
$x=0$. (3) Soient $J\neq0$ et $e$ un idempotent non nul (lemme E.3) avec $\dim J_1(e)$ maximale. Si $J_0(e)\ne0$, c'est
une sous-algèbre formellement réelle (proposition E.4(3)), qui contient un idempotent non nul $f$ ; $f$ est orthogonal à
$e$, donc $e+f$ est un idempotent, $J_1(e)\subseteq J_1(e+f)$ (car $f\circ J_1(e)\subseteq J_0\circ J_1=0$) et
$f\in J_1(e+f)\setminus J_1(e)$, ce qui contredit la maximalité. Donc $J_0(e)=0$. Pour $y\in J_{1/2}(e)$, on a alors
$y^2\in J_1(e)$, et l'associativité donne
$\frac12\tau(y,y)=\tau(e\circ y,y)=\tau(e,y^2)=\operatorname{Tr}L(e\circ y^2)=\operatorname{Tr}L(y^2)=\tau(y,y)$, donc
$y=0$. Ainsi $J=J_1(e)$ et $e$ est une unité. (4) D'après (1)–(3), $\tau(x\circ y,z)=\tau(y,x\circ z)$. Réciproquement,
si $J$ est euclidienne et $\sum_ix_i^2=0$, alors $0=\sum_i\langle x_i\circ x_i,1\rangle=\sum_i\langle x_i,x_i\rangle$.
$\square$

Ainsi « formellement réelle » et « euclidienne » désignent la même classe d'algèbres de dimension finie. Désormais $J$
est formellement réelle, $1$ est son unité et $\tau$ sa forme trace. Le code vérifie le lemme E.1(4), le lemme E.3 (les
idempotents $p_i(x)$ calculés uniquement avec des puissances de Jordan), la proposition E.4 et le théorème E.5 dans
$\operatorname{Herm}_3(\mathbb C)$, donnée par ses constantes de structure ; le produit de Cayley–Dickson `cd` utilisé
pour la construire est celui de la section E.6.

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

## E.4 Repères de Jordan et théorème spectral

Un idempotent non nul est **primitif** s'il n'est pas la somme de deux idempotents orthogonaux non nuls. Un **repère de
Jordan** est une famille $(c_1,\dots,c_r)$ d'idempotents primitifs deux à deux orthogonaux telle que $\sum_ic_i=1$.

**Lemme E.6.** (1) Si $c$ est primitif, alors $J_1(c)=\mathbb Rc$. (2) Des idempotents non nuls deux à deux orthogonaux
sont linéairement indépendants ; tout idempotent non nul est une somme d'idempotents primitifs deux à deux
orthogonaux ; et $J$ possède un repère de Jordan.

*Démonstration.* (1) $J_1(c)$ est une sous-algèbre formellement réelle d'unité $c$ (proposition E.4(3)). Pour
$x\in J_1(c)$, écrivons $x=\sum_i\lambda_ic_i$ d'après le lemme E.3, avec $c_i\in\mathbb R[x]\subseteq J_1(c)$. Alors
$c-c_i$ est un idempotent orthogonal à $c_i$ ($(c-c_i)^2=c-2c_i+c_i$ et $c_i\circ(c-c_i)=0$), donc $c-c_i=0$ par
primitivité. Des $c_i$ distincts sont orthogonaux, il y en a donc au plus un, et $x\in\mathbb Rc$. (2) Si
$\sum_ia_if_i=0$, le produit par $f_j$ donne $a_jf_j=0$. Un idempotent $e$ est donc somme d'un nombre maximal
d'idempotents non nuls deux à deux orthogonaux. Si l'un d'eux, $f$, était égal à $g+g'$ avec $g,g'$ idempotents
orthogonaux non nuls, alors $g=f\circ g\in J_1(f)\subseteq J_0(f')$ pour tout autre membre $f'$ (proposition E.4(4)),
et de même pour $g'$, donc $f$ pourrait être remplacé par $g,g'$. On prend $e=1$. $\square$

**Théorème E.7 (théorème spectral).** Pour tout $x\in J$, il existe un repère de Jordan $(c_1,\dots,c_r)$ et des réels
$\lambda_1,\dots,\lambda_r$ (pas nécessairement distincts) tels que $x=\sum_i\lambda_ic_i$.

*Démonstration.* Écrivons $x=\sum_{i=1}^k\lambda_ie_i$ d'après le lemme E.3 et ajoutons $e_0=1-\sum_ie_i$,
$\lambda_0=0$ ; $e_0$ est un idempotent orthogonal à tous les $e_i$. Décomposons chaque $e_i$ non nul en idempotents
primitifs deux à deux orthogonaux (lemme E.6(2)). Un morceau $f$ de $e_i$ est dans $J_1(e_i)\subseteq J_0(e_l)$ pour
$l\ne i$, et un morceau de $e_l$ est dans $J_1(e_l)$ ; comme $J_0(e_l)\circ J_1(e_l)=0$, tous les morceaux sont deux à
deux orthogonaux. Leur somme est $\sum_ie_i=1$. $\square$

Le fait que tous les repères de Jordan de $J$ ont le même nombre d'éléments est démontré au corollaire E.28.

**Théorème E.8 (décomposition de Peirce relative à un repère).** Soient $(c_1,\dots,c_r)$ un repère de Jordan,
$t_i=\tau(c_i,c_i)$ et, pour $i\ne j$, $J_{ij}=J_{ji}=J_{1/2}(c_i)\cap J_{1/2}(c_j)$.

1. Les $L(c_k)$ commutent, $J=\bigoplus_i\mathbb Rc_i\oplus\bigoplus_{i<j}J_{ij}$, et
   $c_k\circ x=\frac12(\delta_{ki}+\delta_{kj})x$ pour $x\in J_{ij}$.
2. Pour $i,j,k$ distincts : $J_{ij}\circ J_{jk}\subseteq J_{ik}$ ; et $J_{ij}\circ J_{kl}=0$ si $\{i,j\}\cap\{k,l\}=\emptyset$.
3. Il existe une forme quadratique définie positive $N$ sur $J_{ij}$ telle que $x^2=N(x)(c_i+c_j)$ ; avec sa forme
   polaire $\langle\cdot,\cdot\rangle$, $x\circ y=\langle x,y\rangle(c_i+c_j)$ et $\tau(x,y)=(t_i+t_j)\langle x,y\rangle$
   pour $x,y\in J_{ij}$. Si $J_{ij}\ne0$, alors $t_i=t_j$.

*Démonstration.* (1) Les $L(c_k)$ commutent d'après la proposition E.4(4), sont diagonalisables et ont pour somme
$L(1)=\mathrm{id}$. Un vecteur propre commun a des valeurs propres $\alpha_k\in\{0,\frac12,1\}$ avec
$\sum_k\alpha_k=1$ : ou bien un $\alpha_i=1$ et les autres sont nuls, et le vecteur est alors dans
$J_1(c_i)=\mathbb Rc_i$ (lemme E.6(1)) ; ou bien $\alpha_i=\alpha_j=\frac12$ pour une paire et les autres sont nuls.
Ainsi $J_{ij}$ est exactement le sous-espace propre commun du second type, et $J$ est la somme directe des sous-espaces
propres communs. (2) Soient $x\in J_{ij}$ et $y\in J_{jk}$. Relativement à $c_i$, $x\in J_{1/2}$ et $y\in J_0$, donc
$x\circ y\in J_{1/2}(c_i)$ ; relativement à $c_k$, de même $x\circ y\in J_{1/2}(c_k)$. Par conséquent
$x\circ y\in J_{ik}$. Si au contraire $y\in J_{kl}$, $x\circ y$ est dans $J_{1/2}(c_p)$ pour $p=i,j,k,l$ ; il n'existe
pas de vecteur propre commun ayant quatre valeurs propres $\frac12$, donc $x\circ y=0$. (3) Pour $x,y\in J_{ij}$,
$x\circ y\in J_0(c_p)+J_1(c_p)$ pour $p=i,j$ et $x\circ y\in J_0(c_k)$ pour $k\notin\{i,j\}$ (proposition E.4(3)) ;
d'après (1), $x\circ y\in\mathbb Rc_i+\mathbb Rc_j$. Écrivons $x^2=\alpha c_i+\beta c_j$. L'identité de Jordan avec
$y=c_i$ donne $(x^2\circ c_i)\circ x=\frac12\alpha x$ et $x^2\circ(c_i\circ x)=\frac12x^2\circ x=\frac14(\alpha+\beta)x$,
donc $\alpha=\beta$ si $x\ne0$ ; posons $N(x)=\alpha$. Comme $\tau(c_i,1)=\operatorname{Tr}L(c_i)=t_i$, on a
$\tau(x,x)=\tau(x^2,1)=N(x)(t_i+t_j)$ ; donc $N$ est une forme quadratique, définie positive d'après le théorème E.5(2),
et la polarisation donne les formules pour $x\circ y$ et $\tau$. Enfin
$N(x)t_i=\tau(x^2,c_i)=\tau(x,x\circ c_i)=\frac12\tau(x,x)$, et de même avec $j$ ; donc $t_i=t_j$ s'il existe un
$x\ne0$. $\square$

## E.5 Algèbres simples et rang deux

**Théorème E.9.** $J$ est une somme directe $J^{(1)}\oplus\dots\oplus J^{(s)}$ d'idéaux simples avec
$J^{(p)}\circ J^{(q)}=0$ pour $p\ne q$. Chaque $J^{(p)}$ est une algèbre de Jordan formellement réelle.

*Démonstration.* Si $I$ est un idéal, son orthogonal $I^\perp$ pour $\tau$ en est un aussi :
$\tau(J\circ I^\perp,I)=\tau(I^\perp,J\circ I)=0$. Par conséquent $J=I\oplus I^\perp$ et
$I\circ I^\perp\subseteq I\cap I^\perp=0$, donc $J$ est la somme directe des algèbres $I$ et $I^\perp$, et un idéal de
l'un des facteurs est un idéal de $J$. Soit $I$ un idéal non nul de dimension minimale ; alors $I$ est simple, et l'on
conclut par récurrence sur la dimension, appliquée à l'algèbre formellement réelle $I^\perp$. $\square$

**Lemme E.10.** Soient $(c_i)$ un repère de Jordan, $i,j,k,l$ distincts, $u\in J_{ij}$, $w\in J_{jk}$ et
$z\in J_{kl}$. Alors (1) $u\circ(u\circ w)=\frac14N(u)w$, et (2) $(u\circ w)\circ z=u\circ(w\circ z)$.

*Démonstration.* (1) On applique le lemme E.1(1) avec $x=u$ à $c_j$. Avec $L(u^2)=N(u)L(c_i+c_j)$,
$c_j\circ w=\frac12w$, $c_j\circ u=\frac12u$ et $c_j\circ(u\circ w)=0$ (car $u\circ w\in J_{ik}$) :
$[L(u^2),L(w)]c_j=\frac14N(u)w-\frac12N(u)w$ et $2[L(u\circ w),L(u)]c_j=(u\circ w)\circ u$ ; leur somme est nulle.
(2) On applique le lemme E.1(3), avec $(x,w,y)$ remplacé par $(u,z,w)$, à $c_j$. Ici $u\circ z=0$ (théorème E.8(2)),
$c_j\circ u=\frac12u$ et $c_j\circ z=0$, donc $0=\frac12(z\circ w)\circ u-\frac12z\circ(w\circ u)$. $\square$

**Proposition E.11.** Soit $(c_1,\dots,c_r)$ un repère de Jordan de $J$, et relions $i\neq j$ par une arête si
$J_{ij}\ne0$.

1. $J$ est simple si et seulement si ce graphe est connexe.
2. Si $J$ est simple et $r\ge2$, tous les $J_{ij}$, $i\ne j$, ont la même dimension $d\ge1$, tous les $t_i$ valent
   $t=1+\frac{(r-1)d}2$, et $\dim J=r+d\,\frac{r(r-1)}2=rt$.

*Démonstration.* (1) Si le graphe n'est pas connexe, soient $S$ une composante connexe et $I$ le sous-espace engendré par
les $c_i$ et les $J_{ij}$ avec $i,j\in S$, et $I'$ l'analogue pour le complémentaire. Alors $J=I\oplus I'$ (tous les
$J_{ij}$ avec $i\in S$, $j\notin S$ sont nuls), $I\circ I\subseteq I$ et $I\circ I'=0$ d'après le théorème E.8 ; donc
$I$ est un idéal, $\ne0,J$. Réciproquement, supposons le graphe connexe et soit $I\ne0$ un idéal. $I$ est invariant par
les $L(c_k)$, qui commutent et sont diagonalisables, donc il est la somme de ses intersections avec les sous-espaces
propres communs. Si $0\ne x\in I\cap J_{ij}$, alors $c_i\circ x^2=N(x)c_i\in I$ ; donc un certain $c_i$ est dans $I$.
Si $c_i\in I$ et $0\ne x\in J_{ij}$, alors $x=2c_i\circ x\in I$ et $c_j=c_j\circ x^2/N(x)\in I$. Par connexité, tous
les $c_j\in I$, donc $1\in I$ et $I=J$. (2) Soient $u\in J_{ij}$ avec $N(u)=1$ et $k\notin\{i,j\}$. D'après le
théorème E.8(2), $2L(u)$ envoie $J_{jk}$ dans $J_{ik}$ et $J_{ik}$ dans $J_{jk}$, et d'après le lemme E.10(1) (appliqué
aussi avec $i,j$ échangés) son carré est l'identité sur les deux. Donc $\dim J_{ik}=\dim J_{jk}$ dès que $J_{ij}\ne0$.
Ainsi des arêtes $ij$, $jk$ imposent une arête $ik$ ; un graphe connexe est alors complet, et en changeant un indice à
la fois on voit que toutes les $\dim J_{ij}$ sont égales. Enfin $L(c_i)$ a la valeur propre $1$ sur $\mathbb Rc_i$,
$\frac12$ sur les $r-1$ espaces $J_{ij}$ et $0$ ailleurs, donc $t_i=\operatorname{Tr}L(c_i)=1+\frac{(r-1)d}2$. $\square$

Dans le reste de la section, $J$ est simple, $(c_1,\dots,c_r)$ est un repère de Jordan avec $r\ge2$, et pour
$x\in J_{ij}$, $y\in J_{jk}$ ($i,j,k$ distincts) on écrit
$$
m(x,y)=2\,x\circ y\in J_{ik}.
$$
D'après le théorème E.8(3) et la proposition E.11(2), $\tau(x,y)=2t\langle x,y\rangle$ sur chaque $J_{ij}$.

**Lemme E.12.** Pour $i,j,k$ distincts et $x\in J_{ij}$, $y\in J_{jk}$, $z\in J_{ik}$ :

1. $N(m(x,y))=N(x)N(y)$ et $m(x,m(x,y))=N(x)y$ ;
2. $\langle m(x,y),z\rangle=\langle y,m(x,z)\rangle=\langle x,m(y,z)\rangle$ ;
3. $m(m(x,y),v)=m(x,m(y,v))$ pour $l\notin\{i,j,k\}$ et $v\in J_{kl}$.

*Démonstration.* (2) Les trois nombres sont $\tau(x\circ y,z)/t$, $\tau(y,x\circ z)/t$ et $\tau(x,y\circ z)/t$, égaux par
associativité et commutativité. (1) La seconde formule est le lemme E.10(1). Ensuite
$2tN(m(x,y))=4\tau(x\circ y,x\circ y)=4\tau(y,x\circ(x\circ y))=N(x)\tau(y,y)=2tN(x)N(y)$. (3) Lemme E.10(2). $\square$

**Théorème E.13 (rang un et deux).** Soit $J$ simple, avec un repère de Jordan à $r$ éléments. Si $r=1$, alors
$J=\mathbb Rc_1\cong\mathbb R$. Si $r=2$, alors $J\cong V_{d+1}$ avec $d=\dim J_{12}\ge1$ (facteurs de spin comme à la
section 13.2).

*Démonstration.* $r=1$ : théorème E.8(1). $r=2$ : $J=\mathbb R1\oplus W$ avec $W=\mathbb Rf\oplus J_{12}$, $f=c_1-c_2$.
D'après le théorème E.8, $f\circ f=1$, $f\circ x=\frac12x-\frac12x=0$ et $x\circ y=\langle x,y\rangle1$ pour
$x,y\in J_{12}$. Ainsi $w\circ w'=B(w,w')1$ pour $w,w'\in W$, où $B(\alpha f+x,\beta f+y)=\alpha\beta+\langle x,y\rangle$
est un produit scalaire. Si $U:W\to\mathbb R^{d+1}$ est une isométrie, $a1+w\mapsto(a,Uw)$ envoie
$(a1+w)\circ(b1+w')=(ab+B(w,w'))1+aw'+bw$ sur le produit de $V_{d+1}$. $\square$

## E.6 Algèbres de composition : les théorèmes de Hurwitz et de Frobenius

Une **algèbre de composition** est une algèbre $A$ d'unité $1$ munie d'une forme quadratique définie positive $N$ telle
que $N(xy)=N(x)N(y)$. On note $\langle\cdot,\cdot\rangle$ la forme polaire ($N(x)=\langle x,x\rangle$),
$\operatorname{Re}x=\langle x,1\rangle$ et $\bar x=2\operatorname{Re}(x)1-x$. Comme $N(1)=1$ (lemme E.14(1), dont la
démonstration n'utilise pas la conjugaison), c'est la symétrie orthogonale par rapport à la droite $\mathbb R1$, donc
$x\mapsto\bar x$ est une involution isométrique et $\operatorname{Re}\bar x=\operatorname{Re}x$. Un isomorphisme
d'algèbres de composition est un isomorphisme d'algèbres isométrique. L'**associateur** est $(x,y,z)=(xy)z-x(yz)$.

**Lemme E.14.** Dans une algèbre de composition, pour tous $x,y,z,w$ :

1. $N(1)=1$, $\langle xy,xz\rangle=N(x)\langle y,z\rangle$, $\langle yx,zx\rangle=\langle y,z\rangle N(x)$ et
   $\langle xy,wz\rangle+\langle wy,xz\rangle=2\langle x,w\rangle\langle y,z\rangle$ ;
2. $\langle xy,z\rangle=\langle y,\bar xz\rangle=\langle x,z\bar y\rangle$ ;
3. $\bar x(xy)=N(x)y=(yx)\bar x$ ; par conséquent $\bar xx=x\bar x=N(x)1$, $x^2-2\operatorname{Re}(x)x+N(x)1=0$,
   $xy+yx=2\operatorname{Re}(x)y+2\operatorname{Re}(y)x-2\langle x,y\rangle1$,
   $\bar x(wy)+\bar w(xy)=2\langle x,w\rangle y$ et $(yx)\bar w+(yw)\bar x=2\langle x,w\rangle y$ ;
4. $\overline{xy}=\bar y\bar x$, $\operatorname{Re}(xy)=\operatorname{Re}(yx)$ et
   $\operatorname{Re}((xy)z)=\operatorname{Re}(x(yz))$ ;
5. $x(xy)=x^2y$ et $(yx)x=yx^2$ ; l'associateur est alterné (il change de signe quand on échange deux arguments).

*Démonstration.* (1) $N(1)=N(1)^2>0$. On polarise $N(xy)=N(x)N(y)$ en $y$, respectivement en $x$, puis
$\langle xy,xz\rangle=N(x)\langle y,z\rangle$ en $x$. (2) La dernière identité de (1) avec $w=1$ s'écrit
$\langle xy,z\rangle+\langle y,xz\rangle=2\operatorname{Re}(x)\langle y,z\rangle$, c'est-à-dire
$\langle xy,z\rangle=\langle y,\bar xz\rangle$ ; la polarisation de $\langle yx,zx\rangle=\langle y,z\rangle N(x)$ en
$x$, avec l'une des deux variables égale à $1$, donne $\langle yx,z\rangle=\langle y,z\bar x\rangle$. (3) D'après (2)
et (1), $\langle\bar x(xy),z\rangle=\langle xy,xz\rangle=N(x)\langle y,z\rangle$ pour tout $z$, et de même à droite. On
pose $y=1$ et on substitue $\bar x=2\operatorname{Re}(x)1-x$ ; puis on polarise. (4) D'après (2) et l'isométrie de la
conjugaison,
$\langle\bar y\bar x,z\rangle=\langle\bar x,yz\rangle=\langle\bar x\bar z,y\rangle=\langle\bar z,xy\rangle=\langle z,\overline{xy}\rangle$.
Ensuite $\operatorname{Re}(xy)=\langle y,\bar x\rangle=\langle\bar y,x\rangle=\operatorname{Re}(yx)$, et
$\operatorname{Re}((xy)z)=\langle xy,\bar z\rangle=\langle y,\bar x\bar z\rangle=\langle yz,\bar x\rangle=\operatorname{Re}(x(yz))$.
(5) D'après (3), $x(xy)=(2\operatorname{Re}(x)-\bar x)(xy)=(2\operatorname{Re}(x)x-N(x))y=x^2y$, et de même à droite.
Donc $(x,x,y)=(y,x,x)=0$ ; par polarisation, l'associateur change de signe quand on échange les deux premiers ou les
deux derniers arguments, et ces transpositions engendrent toutes les permutations. $\square$

**La construction de Cayley–Dickson.** Pour une algèbre $B$ munie d'une unité, d'une involution $b\mapsto\bar b$ et
d'une forme quadratique $N$, soit $\operatorname{CD}(B)=B\times B$ avec
$$
(a,b)(c,d)=(ac-\bar db,\;da+b\bar c),\qquad N(a,b)=N(a)+N(b),\qquad\overline{(a,b)}=(\bar a,-b).
$$
En partant de $\mathbb R$ ($\bar x=x$, $N(x)=x^2$), on définit $\mathbb C=\operatorname{CD}(\mathbb R)$,
$\mathbb H=\operatorname{CD}(\mathbb C)$ et $\mathbb O=\operatorname{CD}(\mathbb H)$, de dimensions $1,2,4,8$, chaque
doublement utilisant l'involution et la forme que l'on vient de définir (d'après le lemme E.15(1) ci-dessous, ce sont la
conjugaison et la norme de l'algèbre de composition). Dans $\mathbb C$, $(a,b)(c,d)=(ac-bd,ad+bc)$ : ce sont les
nombres complexes $a+bi$. $B$ est la sous-algèbre $B\times0$.

**Lemme E.15.** (1) Si $B$ est une algèbre de composition associative, alors $\operatorname{CD}(B)$ est une algèbre de
composition d'unité $(1,0)$ et de conjugaison $\overline{(a,b)}=(\bar a,-b)$. (2) $\mathbb C$ est commutative et
associative. $M(a,b)=\begin{pmatrix}a&-b\\\bar b&\bar a\end{pmatrix}$ est un homomorphisme d'algèbres injectif
$\mathbb H\to M_2(\mathbb C)$ avec $M(\bar q)=M(q)^*$ ; donc $\mathbb H$ est associative ; elle n'est pas commutative.
$\mathbb O$ n'est pas associative. (3) $\mathbb R,\mathbb C,\mathbb H,\mathbb O$ sont des algèbres de composition, et
$\mathbb R,\mathbb C,\mathbb H$ sont des algèbres à division : tout $x\ne0$ a pour inverse bilatère $\bar x/N(x)$, donc
$xy=0$ entraîne $x=0$ ou $y=0$.

*Démonstration.* (1) $N((a,b)(c,d))=N(a)N(c)+N(d)N(b)+N(d)N(a)+N(b)N(c)-2\langle ac,\bar db\rangle+2\langle da,b\bar c\rangle$,
et d'après le lemme E.14(2), $\langle ac,\bar db\rangle=\langle d(ac),b\rangle$ et
$\langle da,b\bar c\rangle=\langle(da)c,b\rangle$, qui coïncident quand $B$ est associative. $(1,0)$ est une unité, et
$2\langle(a,b),(1,0)\rangle(1,0)-(a,b)=(\bar a,-b)$. (2) En développant et en utilisant que $\mathbb C$ est
commutative, $M(a,b)M(c,d)=M(ac-\bar db,da+b\bar c)$ ; $M$ est linéaire et injective, et $M(\bar a,-b)=M(a,b)^*$. Dans
$\mathbb H$, $(i,0)(0,1)=(0,i)$ mais $(0,1)(i,0)=(0,-i)$. Dans $\mathbb O=\mathbb H\times\mathbb H$, soient
$p=((i,0),0)$, $q=((0,1),0)$, $s=(0,1)$ ; alors $(pq)s=(0,(0,i))$ et $p(qs)=(0,(0,-i))$. (3) D'après (1), (2) et par
récurrence. Si $x\ne0$, $x(\bar x/N(x))=(\bar x/N(x))x=1$ d'après le lemme E.14(3) ; dans une algèbre associative,
$xy=0$ donne alors $y=(\bar x/N(x))(xy)=0$. $\square$

**Lemme E.16 (doublement).** Soient $A$ une algèbre de composition, $B\subseteq A$ une sous-algèbre avec $1\in B$, et
$i\in A$ avec $i\perp B$ et $N(i)=1$. Pour $a,b,c,d\in B$ : $Bi\perp B$, $N(a+bi)=N(a)+N(b)$ et
$$
(a+bi)(c+di)=(ac-\bar db)+(da+b\bar c)\,i .
$$
Ainsi $B\oplus Bi$ est une sous-algèbre de dimension $2\dim B$, et $a+bi\mapsto(h(a),h(b))$ est un isomorphisme
$B\oplus Bi\cong\operatorname{CD}(D)$ dès que $h:B\to D$ est un isomorphisme d'algèbres de composition. De plus $B$ est
associative ; si $B\oplus Bi$ est associative, $B$ est commutative ; si $B\oplus Bi$ est commutative, $B=\mathbb R1$.

*Démonstration.* $B$ contient $\bar b=2\operatorname{Re}(b)1-b$ en même temps que $b$. D'après le lemme E.14,
$\langle bi,a\rangle=\langle i,\bar ba\rangle=0$, $N(bi)=N(b)$, $\bar i=-i$, et $ci+ic=2\operatorname{Re}(c)i$,
c'est-à-dire $ci=i\bar c$, pour $c\in B$. Les identités du lemme E.14(3) donnent successivement : avec
$(x,w,y)=(\bar a,i,b)$, $a(ib)=i(\bar ab)$ ; comme $ib=\bar bi$ et $i(\bar ab)=\overline{\bar ab}\,i=(\bar ba)i$,
remplacer $b$ par $\bar b$ donne $a(bi)=(ba)i$. Avec $(y,x,w)=(a,i,\bar c)$, $(ai)c=(a\bar c)i$. Avec $x=i$,
$i(iy)=-y$. Avec $(x,w,y)=(ai,i,\bar d)$, où $\overline{ai}=-ai$ (car $ai\perp1$) et $\langle ai,i\rangle=\operatorname{Re}a$
(lemme E.14(1)),
$$
-(ai)(i\bar d)-i((ai)\bar d)=2\operatorname{Re}(a)\bar d .
$$
Ici $i\bar d=di$ et $i((ai)\bar d)=i((ad)i)=i(i\,\overline{ad})=-\bar d\bar a$, donc
$(ai)(di)=\bar d\bar a-2\operatorname{Re}(a)\bar d=-\bar da$. En développant $(a+bi)(c+di)$ avec ces quatre règles, on
obtient la formule, et $N$ et les produits correspondent à ceux de $\operatorname{CD}(D)$ via $h$. Comme $N$ est
multiplicative sur $B\oplus Bi\subseteq A$, le calcul de la démonstration du lemme E.15(1) donne
$\langle d(ac)-(da)c,b\rangle=0$ pour tout $b\in B$, donc $B$ est associative. Si $B\oplus Bi$ est associative,
$(ab)i=a(bi)=(ba)i$, donc $ab=ba$. Si elle est commutative, $ai=ia=\bar ai$, donc $a=\bar a\in\mathbb R1$. $\square$

**Théorème E.17 (Hurwitz).** Toute algèbre de composition est isomorphe à exactement une des algèbres
$\mathbb R,\mathbb C,\mathbb H,\mathbb O$.

*Démonstration.* Posons $B_0=\mathbb R1\cong\mathbb R$. Tant que $B_k\ne A$, choisissons $i_{k+1}\perp B_k$ avec
$N(i_{k+1})=1$ et posons $B_{k+1}=B_k\oplus B_ki_{k+1}$ ; d'après le lemme E.16, $B_k\cong\operatorname{CD}^k(\mathbb R)$
est de dimension $2^k$. Si $B_4$ existait, le lemme E.16 rendrait $B_3$ associative, donc $B_2$ commutative, donc
$B_1=\mathbb R1$, qui est de dimension $1\ne2$. Ainsi $A=B_k$ pour un $k\le3$. Les quatre algèbres ont des dimensions
différentes. $\square$

**Lemme E.18.** Soit $A$ une algèbre de composition. (1) Deux éléments quelconques de $A$ sont contenus dans une
sous-algèbre associative. (2) (**Identité de Moufang**) $(ux)(yu)=u((xy)u)=(u(xy))u$ pour tous $u,x,y$.

*Démonstration.* (1) Soient $p,q\in A$ ; on peut supposer $p\notin\mathbb R1$ (sinon on les échange, et si tous deux
sont réels on prend $\mathbb R1$). $B=\operatorname{span}\{1,p\}$ est une sous-algèbre
($p^2=2\operatorname{Re}(p)p-N(p)1$), une algèbre de composition de dimension $2$, donc $B\cong\mathbb C$ (théorème
E.17). Si $q\in B$, c'est terminé ; sinon, soit $i$ la composante normalisée de $q$ orthogonale à $B$. D'après le lemme
E.16, $B\oplus Bi\cong\operatorname{CD}(\mathbb C)=\mathbb H$ est associative et contient $p,q$. (2) L'associateur est
alterné (lemme E.14(5)), donc $(u,xy,u)=0$, ce qui est la seconde égalité. Soit $X=(u,x,y)$. $\operatorname{Re}X=0$
d'après le lemme E.14(4), et $\langle u,X\rangle=\langle u\bar y,ux\rangle-N(u)\langle1,xy\rangle
=N(u)(\langle\bar y,x\rangle-\langle\bar x,y\rangle)=0$ ; donc l'identité quadratique polarisée du lemme E.14(3) donne
$uX+Xu=2\operatorname{Re}(u)X$, c'est-à-dire $\bar uX=Xu$. En développant les définitions, on voit que dans toute
algèbre
$$
(wx,y,z)-(w,xy,z)+(w,x,yz)=w(x,y,z)+(w,x,y)z .
$$
Avec $(w,x,y,z)=(u,u,x,y)$, et $(u^2,x,y)=2\operatorname{Re}(u)X$ : $(u,ux,y)=2\operatorname{Re}(u)X-uX=\bar uX=Xu$.
Avec $(w,x,y,z)=(u,x,y,u)$, $(ux,y,u)=(u,ux,y)=Xu$ et $(x,y,u)=X$ (permutations circulaires) :
$(u,x,yu)=uX+Xu-Xu=uX$. Par conséquent $(ux)(yu)-u((xy)u)=(u,x,yu)-u(x,y,u)=uX-uX=0$. $\square$

**Théorème E.19 (Frobenius).** Toute algèbre réelle associative $A$ de dimension finie, d'unité $1\ne0$ et sans
diviseurs de zéro, est isomorphe à $\mathbb R$, $\mathbb C$ ou $\mathbb H$.

*Démonstration.* (i) Pour $x\in A$, les puissances $1,x,x^2,\dots$ sont linéairement dépendantes ; soit $\mu$ le
polynôme réel unitaire de plus petit degré tel que $\mu(x)=0$. Si $\mu=fg$ avec $\deg f,\deg g\ge1$, alors
$f(x)g(x)=0$ impose $f(x)=0$ ou $g(x)=0$ ; donc $\mu$ est irréductible et par conséquent (théorème fondamental de
l'algèbre) $\mu=t-\alpha$ ou $\mu=(t-\alpha)^2+\beta^2$, $\beta>0$. Soit $V=\{v\in A:v^2\in(-\infty,0]\cdot1\}$. Alors
$A=\mathbb R1+V$ (dans le second cas $x-\alpha\in V$), et $\mathbb R1\cap V=0$. (ii) $V$ est un sous-espace : il est
stable par multiplication par les scalaires, et pour $u,v\in V$ nous montrons que $u+v\in V$. Si $u,v$ sont
linéairement dépendants, $u+v$ est un multiple de l'un d'eux ; supposons-les donc indépendants. Alors $1,u,v$ sont
indépendants : si $1=\alpha u+\gamma v$ avec, par exemple, $\alpha\ne0$, alors
$u^2=(1-2\gamma v+\gamma^2v^2)/\alpha^2\in\mathbb R1$ impose $\gamma=0$ et $u\in\mathbb R1\cap V=0$. D'après (i),
$(u\pm v)^2=p_\pm(u\pm v)+q_\pm1$ avec $p_\pm,q_\pm$ réels ; en additionnant,
$2u^2+2v^2=(p_++p_-)u+(p_+-p_-)v+(q_++q_-)1$, donc $p_\pm=0$. Si $q_+>0$, alors
$(u+v-\sqrt{q_+})(u+v+\sqrt{q_+})=0$ et $u+v\in\mathbb R1$, ce qui est impossible ; donc $u+v\in V$. (iii) Pour
$u,v\in V$, $uv+vu=(u+v)^2-u^2-v^2=-2\langle u,v\rangle1$ définit une forme bilinéaire symétrique avec
$\langle u,u\rangle=-u^2\ge0$, et $\langle u,u\rangle=0$ seulement si $u^2=0$, c'est-à-dire $u=0$ : un produit
scalaire. Des $u,v\in V$ orthogonaux anticommutent. (iv) Si $V=0$, $A=\mathbb R$. Si $\dim V=1$, $V=\mathbb Ri$ avec
$i^2=-1$, et $A\cong\mathbb C$. Si $\dim V\ge2$, prenons $i,j\in V$ orthonormés et $k=ij$. Alors
$k^2=ijij=-i^2j^2=-1$, donc $k\in V$, $ik+ki=-j+j=0$ et $jk+kj=i-i=0$ : $i,j,k$ sont orthonormés, avec $ij=-ji=k$,
$jk=-kj=i$, $ki=-ik=j$. S'il existait un vecteur unitaire $l\in V$ orthogonal à $i,j,k$, alors
$lk=(li)j=-i(lj)=(ij)l=kl$ et $lk=-kl$, donc $kl=0$, ce qui est impossible. Par conséquent
$A=\operatorname{span}\{1,i,j,k\}$. Dans $\mathbb H$, les éléments $I=(i,0)$, $J=(0,1)$, $K=(0,i)$ ont la même table de
multiplication (le lemme E.15(2) donne $IJ=K=-JI$, et $I^2=J^2=-1$, $K=IJ$ déterminent le reste par associativité),
donc $A\cong\mathbb H$. $\square$

Le code construit $\mathbb O$ par doublement de Cayley–Dickson et vérifie sur des octonions aléatoires les identités du
lemme E.14, l'identité de Moufang du lemme E.18 et la non-associativité du lemme E.15, ainsi que l'homomorphisme $M$ :

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

Pour $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$, $d=\dim D$ et $n\ge1$, soit
$\operatorname{Herm}_n(D)=\{X\in M_n(D):X_{ji}=\overline{X_{ij}}\}$ muni de $X\circ Y=\frac12(XY+YX)$, où
$(XY)_{ik}=\sum_jX_{ij}Y_{jk}$ (cette définition n'utilise pas l'associativité de $D$). Les coefficients diagonaux
vérifient $\bar x=x$, ils sont donc réels. Soient $E_{ij}$ les unités matricielles et, pour $i\ne j$ et $a\in D$,
$a[ij]=aE_{ij}+\bar aE_{ji}$ ; ainsi $a[ji]=\bar a[ij]$.

**Lemme E.20.** $\operatorname{Herm}_n(D)$ est stable par $\circ$, a pour unité $\sum_iE_{ii}$ et pour dimension
$n+d\frac{n(n-1)}2$, et, pour $a,b\in D$, avec $i,j,k$ quelconques dans les deux premières formules et $i,j,k,l$ distincts dans les autres ($a[ij]$ toujours avec $i\ne j$) :
$$
E_{ii}\circ E_{jj}=\delta_{ij}E_{ii},\quad E_{kk}\circ a[ij]=\tfrac12(\delta_{ki}+\delta_{kj})a[ij],\quad
a[ij]\circ b[ij]=\langle a,b\rangle(E_{ii}+E_{jj}),\quad a[ij]\circ b[jk]=\tfrac12(ab)[ik],\quad a[ij]\circ b[kl]=0 .
$$
De plus $(X|Y)=\operatorname{Re}\operatorname{tr}(XY)$ est un produit scalaire tel que $(X\circ Y|Z)=(Y|X\circ Z)$.

*Démonstration.* Produit matriciel : $a[ij]\,b[ij]=a\bar bE_{ii}+\bar abE_{jj}$, et
$a\bar b+b\bar a=\bar ab+\bar ba=2\langle a,b\rangle$ (polariser $x\bar x=\bar xx=N(x)$) ; $a[ij]\,b[jk]=abE_{ik}$ et
$b[jk]\,a[ij]=\bar b\bar aE_{ki}=\overline{ab}E_{ki}$ ; les autres produits sont immédiats. Tous les produits d'éléments
de la base sont dans $\operatorname{Herm}_n(D)$, donc tous les produits aussi. Ensuite,
$(X|X)=\sum_{i,j}\operatorname{Re}(X_{ij}\overline{X_{ij}})=\sum_{i,j}N(X_{ij})>0$ pour $X\ne0$. Soit
$T(X,Y,Z)=\operatorname{Re}\operatorname{tr}((XY)Z)=\sum_{i,j,k}\operatorname{Re}((X_{ij}Y_{jk})Z_{ki})$. D'après le
lemme E.14(4), $\operatorname{Re}\operatorname{tr}(AB)=\operatorname{Re}\operatorname{tr}(BA)$, $T(X,Y,Z)=T(Z,X,Y)$ et,
en passant aux conjugués ($\operatorname{Re}\bar p=\operatorname{Re}p$, $\overline{(pq)r}=\bar r(\bar q\bar p)$) et en
utilisant que $X,Y,Z$ sont hermitiennes, $T(X,Y,Z)=\sum\operatorname{Re}(Z_{ik}(Y_{kj}X_{ji}))=T(Z,Y,X)$. Donc $T$ est
symétrique, et $(X\circ Y|Z)=\frac12(T(X,Y,Z)+T(Y,X,Z))=T(X,Y,Z)$ est symétrique en $X,Y,Z$. $\square$

**Lemme E.21 (orientations).** Soient $P_{12},P_{13},P_{23}$ des espaces munis d'un produit scalaire, $P_{ji}=P_{ij}$,
et pour tout triplet ordonné $(i,j,k)$ d'éléments distincts de $\{1,2,3\}$ soit $m:P_{ij}\times P_{jk}\to P_{ik}$
bilinéaire, de sorte que $m(x,y)=m(y,x)$ (le membre de droite étant l'application pour $(k,j,i)$) et
$\langle m(x,y),z\rangle=\langle y,m(x,z)\rangle$ pour $x\in P_{ij}$, $y\in P_{jk}$, $z\in P_{ik}$. Soient $D$ une
algèbre de composition et $\iota_{ij}:D\to P_{ij}$ ($i<j$) des bijections linéaires isométriques telles que
$m(\iota_{12}a,\iota_{23}b)=\iota_{13}(ab)$ ; posons $\iota_{ji}(a)=\iota_{ij}(\bar a)$. Alors
$m(\iota_{ij}a,\iota_{jk}b)=\iota_{ik}(ab)$ pour tous les triplets ordonnés et tous $a,b\in D$.

Les hypothèses sont satisfaites (i) pour trois indices d'un repère de Jordan d'une algèbre $J$ simple, avec les espaces
de Peirce et $m=2\circ$ (commutativité et lemme E.12(2)) ; (ii) dans $\operatorname{Herm}_3(D)$, avec
$P_{ij}=\{a[ij]\}$, $\langle a[ij],b[ij]\rangle=\langle a,b\rangle$ (bien défini, car la conjugaison est isométrique)
et $m=2\circ$, puisque $m(a[ij],b[jk])=(ab)[ik]$ et, pour $c[ik]$, $\langle ab,c\rangle=\langle b,\bar ac\rangle$
(lemme E.14(2)).

*Démonstration.* Si l'assertion vaut pour $(i,j,k)$, elle vaut pour $(k,j,i)$ :
$m(\iota_{kj}a,\iota_{ji}b)=m(\iota_{ij}\bar b,\iota_{jk}\bar a)=\iota_{ik}(\bar b\bar a)=\iota_{ik}(\overline{ab})=\iota_{ki}(ab)$.
Il suffit donc de traiter $(2,1,3)$ et $(1,3,2)$, étant donné $(1,2,3)$. Pour tout $c\in D$, en utilisant l'adjonction,
$(1,2,3)$ et le lemme E.14(2),
$$
\langle m(\iota_{21}a,\iota_{13}b),\iota_{23}c\rangle=\langle\iota_{13}b,m(\iota_{12}\bar a,\iota_{23}c)\rangle=\langle b,\bar ac\rangle=\langle ab,c\rangle=\langle\iota_{23}(ab),\iota_{23}c\rangle,
$$
$$
\langle m(\iota_{13}a,\iota_{32}b),\iota_{12}c\rangle=\langle\iota_{13}a,m(\iota_{12}c,\iota_{23}\bar b)\rangle=\langle a,c\bar b\rangle=\langle ab,c\rangle=\langle\iota_{12}(ab),\iota_{12}c\rangle,
$$
et $\iota_{23},\iota_{12}$ sont surjectives. (Dans la seconde ligne, l'adjonction est utilisée pour $y=\iota_{32}b$
dans le rôle de $x$, après $m(x,y)=m(y,x)$.) $\square$

**Théorème E.22 (coordinatisation).** Soit $J$ simple, avec un repère de Jordan $(c_1,\dots,c_r)$, $r\ge3$. Il existe
alors une algèbre de composition $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$ de dimension $d$ et un isomorphisme
d'algèbres $\Phi:\operatorname{Herm}_r(D)\to J$ tel que $\Phi(E_{ii})=c_i$. Si $r\ge4$, $D$ est associative,
c'est-à-dire $D\ne\mathbb O$.

*Démonstration.* *Étape 1 (un triangle).* Choisissons $e_{12}\in J_{12}$, $e_{23}\in J_{23}$ avec $N=1$ (proposition
E.11). Sur $D'=J_{13}$, définissons $p\cdot q=m(m(p,e_{23}),m(e_{12},q))$, en utilisant les produits
$J_{13}\times J_{32}\to J_{12}$, $J_{21}\times J_{13}\to J_{23}$ et $J_{12}\times J_{23}\to J_{13}$. D'après le lemme
E.12(1), $1'=m(e_{12},e_{23})$ est une unité : $p\cdot1'=m(m(p,e_{23}),e_{23})=p$ et
$1'\cdot q=m(e_{12},m(e_{12},q))=q$ ; et $N(p\cdot q)=N(p)N(q)$. Ainsi $(D',N)$ est une algèbre de composition, et
d'après le théorème E.17 il existe un isomorphisme $h:D\to D'$ avec $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$.
Posons $\iota_{13}=h$, $\iota_{12}=m(h(\cdot),e_{23})$ et $\iota_{23}=m(e_{12},h(\cdot))$ ; ce sont des bijections
isométriques (lemme E.12(1)) et $m(\iota_{12}a,\iota_{23}b)=h(a)\cdot h(b)=\iota_{13}(ab)$. D'après le lemme E.21, pour
$k=3$ :

> $(\mathrm T_k)$ il existe des bijections isométriques $\iota_{ij}:D\to J_{ij}$, $i<j\le k$, telles que, avec $\iota_{ji}(a)=\iota_{ij}(\bar a)$,
> $m(\iota_{ij}a,\iota_{jl}b)=\iota_{il}(ab)$ pour tous $i,j,l\le k$ distincts et tous $a,b\in D$.

*Étape 2 (ajout d'un indice).* Supposons $(\mathrm T_{k-1})$ avec $4\le k\le r$. Choisissons $v\in J_{1k}$ avec
$N(v)=1$ et posons $\iota_{ik}(a)=m(\iota_{i1}a,v)$ pour $2\le i<k$ et $\iota_{1k}(a)=m(\iota_{12}a,\iota_{2k}1)$,
bijections isométriques d'après le lemme E.12(1). En utilisant le lemme E.12(3) (avec chaque fois quatre indices
distincts) et $(\mathrm T_{k-1})$ : pour $i,j\in\{2,\dots,k-1\}$ distincts,
$$
m(\iota_{ij}a,\iota_{jk}b)=m(\iota_{ij}a,m(\iota_{j1}b,v))=m(m(\iota_{ij}a,\iota_{j1}b),v)=m(\iota_{i1}(ab),v)=\iota_{ik}(ab).
$$
En particulier $\iota_{2k}(b)=m(\iota_{23}b,\iota_{3k}1)$, et par conséquent
$$
m(\iota_{12}a,\iota_{2k}b)=m(m(\iota_{12}a,\iota_{23}b),\iota_{3k}1)=m(\iota_{13}(ab),\iota_{3k}1)
=m(m(\iota_{12}(ab),\iota_{23}1),\iota_{3k}1)=m(\iota_{12}(ab),\iota_{2k}1)=\iota_{1k}(ab).
$$
Pour $3\le j<k$, $\iota_{1j}a=m(\iota_{12}a,\iota_{2j}1)$, donc
$m(\iota_{1j}a,\iota_{jk}b)=m(\iota_{12}a,m(\iota_{2j}1,\iota_{jk}b))=m(\iota_{12}a,\iota_{2k}b)=\iota_{1k}(ab)$. La
relation vaut donc pour tous les triplets $i<j<k$, et le lemme E.21, appliqué à chaque triangle $\{i,j,k\}$, donne
$(\mathrm T_k)$.

*Étape 3 (associativité).* Si $r\ge4$, d'après $(\mathrm T_4)$ et le lemme E.12(3),
$\iota_{14}((ab)c)=m(m(\iota_{12}a,\iota_{23}b),\iota_{34}c)=m(\iota_{12}a,m(\iota_{23}b,\iota_{34}c))=\iota_{14}(a(bc))$.

*Étape 4 (l'isomorphisme).* Grâce à $(\mathrm T_r)$, définissons l'application linéaire
$\Phi(\sum_i\alpha_iE_{ii}+\sum_{i<j}a_{ij}[ij])=\sum_i\alpha_ic_i+\sum_{i<j}\iota_{ij}(a_{ij})$ ; elle est bijective
d'après le théorème E.8(1), et $\Phi(a[ij])=\iota_{ij}(a)$ pour tous $i\ne j$. D'après le théorème E.8 et
$(\mathrm T_r)$, $c_i\circ c_j=\delta_{ij}c_i$, $c_k\circ\iota_{ij}(a)=\frac12(\delta_{ki}+\delta_{kj})\iota_{ij}(a)$,
$\iota_{ij}(a)\circ\iota_{ij}(b)=\langle a,b\rangle(c_i+c_j)$, $\iota_{ij}(a)\circ\iota_{jk}(b)=\frac12\iota_{ik}(ab)$ et
$\iota_{ij}(a)\circ\iota_{kl}(b)=0$ : c'est la table du lemme E.20. Par bilinéarité, $\Phi(X\circ Y)=\Phi(X)\circ\Phi(Y)$.
$\square$

## E.8 L'algèbre d'Albert

**Lemme E.23 (automorphismes de $\operatorname{Herm}_3(D)$).** Soient $D$ une algèbre de composition et
$f_{12},f_{23},f_{13}$ des bijections linéaires isométriques de $D$ telles que $f_{13}(ab)=f_{12}(a)f_{23}(b)$ pour tous
$a,b$. Alors
$$
\Phi\Big(\sum_i\alpha_iE_{ii}+\sum_{i<j}a_{ij}[ij]\Big)=\sum_i\alpha_iE_{ii}+\sum_{i<j}f_{ij}(a_{ij})[ij]
$$
est un automorphisme de $(\operatorname{Herm}_3(D),\circ)$. Pour tout $u\in D$ tel que $N(u)=1$,
$(f_{12},f_{23},f_{13})=(L_u,R_u,B_u)$ avec $L_ua=ua$, $R_ua=au$, $B_ua=uau$ est un tel triplet.

*Démonstration.* On applique le lemme E.21 dans l'exemple (ii) à $\iota_{ij}(a)=f_{ij}(a)[ij]$ : alors
$\Phi(a[ij])=\iota_{ij}(a)$ pour tous $i\ne j$, et $2\iota_{ij}(a)\circ\iota_{jk}(b)=\iota_{ik}(ab)$ pour tous $i,j,k$
distincts. Comme $\Phi$ fixe les $E_{ii}$ et que les $f_{ij}$ sont isométriques, $\Phi$ préserve chaque produit de la
table du lemme E.20, donc $\circ$. Quant au triplet, $L_u,R_u,B_u$ sont isométriques ($N(ua)=N(a)$) et bijectives, et
$u(ab)u=(ua)(bu)$ est l'identité de Moufang (lemme E.18(2)). $\square$

**Théorème E.24.** (1) Pour $D\in\{\mathbb R,\mathbb C,\mathbb H\}$ et $n\ge1$, $\operatorname{Herm}_n(D)$ est une
algèbre de Jordan euclidienne, et $\operatorname{Herm}_n(\mathbb H)$ est isomorphe à une sous-algèbre de Jordan de
$\operatorname{Herm}_{2n}(\mathbb C)$. (2) L'**algèbre d'Albert** $\operatorname{Herm}_3(\mathbb O)$ (de dimension
$27$) et $\operatorname{Herm}_2(\mathbb O)$ sont des algèbres de Jordan euclidiennes. (3) Dans chacune d'elles,
$(E_{11},\dots,E_{nn})$ est un repère de Jordan d'espaces de Peirce $J_{ij}=\{a[ij]:a\in D\}$ de dimension
$d=\dim D$ ; l'algèbre est simple ; $\operatorname{Herm}_1(D)\cong\mathbb R$ et $\operatorname{Herm}_2(D)\cong V_{d+1}$
pour $D=\mathbb R,\mathbb C,\mathbb H,\mathbb O$.

*Démonstration.* (1) $M_n(D)$ est une algèbre réelle associative (lemme E.15(2)), donc, d'après le lemme E.20 et le
lemme 13.1, $\operatorname{Herm}_n(D)$ est une algèbre de Jordan, euclidienne d'après le lemme E.20. Remplacer chaque
coefficient $X_{ij}\in\mathbb H$ par le bloc $2\times2$ $M(X_{ij})$ est un homomorphisme d'algèbres injectif
$M_n(\mathbb H)\to M_{2n}(\mathbb C)$ (lemme E.15(2)) qui envoie les matrices hermitiennes sur des matrices
hermitiennes, puisque $M(\bar q)=M(q)^*$ ; il préserve donc $\circ$.

(2) Soit $J(x,y)=(x^2\circ y)\circ x-x^2\circ(y\circ x)$ sur $\operatorname{Herm}_3(\mathbb O)$ ; un automorphisme
$\Phi$ vérifie $J(\Phi x,\Phi y)=\Phi J(x,y)$. Soit $x=\sum_i\alpha_iE_{ii}+p[12]+q[23]+s[13]$. Nous composons trois
automorphismes de la forme du lemme E.23 de sorte que les coefficients en $(1,3)$ et $(2,3)$ deviennent réels. (a) Si
$s\notin\mathbb R$, écrivons $s=|s|(\cos\theta+v\sin\theta)$ avec $v\perp1$, $N(v)=1$, de sorte que $v^2=-1$ (lemme
E.14(3)) ; pour $u=\cos\frac\theta2-v\sin\frac\theta2$, en calculant dans la sous-algèbre commutative
$\operatorname{span}\{1,v\}\cong\mathbb C$, on a $B_us=u^2s=|s|$. (b) Maintenant $s$ est réel. Choisissons $u\perp1$
avec $N(u)=1$ et $u\perp q$ (possible car $\dim\mathbb O=8$) ; alors $B_us=su^2=-s$ est réel et
$\operatorname{Re}(qu)=\langle u,\bar q\rangle=0$, car $\bar q\in\operatorname{span}\{1,q\}$. (c) Maintenant
$q\perp1$ ; si $q\ne0$, prenons $u=q/|q|$ : $qu=q^2/|q|=-|q|$ et $B_us=-s$ sont réels. Il existe donc un automorphisme
$\Phi$ tel que $\Phi x$ ait des coefficients réels, sauf peut-être $p'$ en $(1,2)$. Comme $J(\Phi x,\cdot)$ est
linéaire, il suffit de montrer que $J(\Phi x,y)=0$ pour $y=E_{ii}$ et $y=c[ij]$, $c\in\mathbb O$. D'après le lemme
E.18(1), $p'$ et $c$ sont dans une sous-algèbre associative $Q$ ; alors $\Phi x,y\in\operatorname{Herm}_3(Q)$, une
sous-algèbre de $\operatorname{Herm}_3(\mathbb O)$ contenue dans l'algèbre associative $M_3(Q)$, et $J(\Phi x,y)=0$
d'après le lemme 13.1. Par conséquent $J(x,\Phi^{-1}y)=\Phi^{-1}J(\Phi x,y)=0$ pour tout $y$ : l'identité de Jordan est
vérifiée, et le lemme E.20 rend l'algèbre euclidienne. Les matrices dont la troisième ligne et la troisième colonne sont
nulles forment une sous-algèbre (lemme E.20) d'unité $E_{11}+E_{22}$, qui est $\operatorname{Herm}_2(\mathbb O)$.

(3) D'après le lemme E.20, les $E_{ii}$ sont des idempotents deux à deux orthogonaux de somme $1$ ; $L(E_{kk})$ agit
sur $E_{jj}$ par $\delta_{kj}$ et sur $a[ij]$ par $\frac12(\delta_{ki}+\delta_{kj})$, donc $J_1(E_{ii})=\mathbb RE_{ii}$
et $J_{1/2}(E_{ii})\cap
J_{1/2}(E_{jj})=\{a[ij]\}$. $E_{ii}$ est primitif : si $E_{ii}=g+g'$ avec des idempotents orthogonaux non nuls, alors
$g=E_{ii}\circ g\in\mathbb RE_{ii}$, donc $g=E_{ii}$ et $g'=0$. Le graphe de la proposition E.11 est complet, donc
l'algèbre est simple. Pour $n=1,2$, on applique le théorème E.13. $\square$

**Corollaire E.25.** Pour $n\ge4$, $\operatorname{Herm}_n(\mathbb O)$ n'est pas une algèbre de Jordan.

*Démonstration.* Si c'en était une, le lemme E.10(2) s'appliquerait : sa démonstration n'utilise que l'identité de
Jordan et les produits du lemme E.20 (avec $c_k=E_{kk}$). On aurait alors
$\frac14((ab)c)[14]=(a[12]\circ b[23])\circ c[34]=a[12]\circ(b[23]\circ c[34])=\frac14(a(bc))[14]$ pour tous
$a,b,c\in\mathbb O$, ce qui contredit le lemme E.15(2). $\square$

**Théorème E.26 (Albert).** Soit $D$ une algèbre de composition. S'il existe une algèbre réelle associative $A$ et une
application linéaire injective $\varphi:\operatorname{Herm}_3(D)\to A$ telle que
$\varphi(x\circ y)=\frac12(\varphi(x)\varphi(y)+\varphi(y)\varphi(x))$, alors $D$ est associative. En particulier,
l'algèbre d'Albert est **exceptionnelle** : elle n'est isomorphe à aucune sous-algèbre de Jordan d'une algèbre
associative, telle que $M_N(\mathbb C)$.

*Démonstration.* Posons $e_i=\varphi(E_{ii})$ et $e=e_1+e_2+e_3=\varphi(1)$. Alors $e_i^2=e_i$, et $e_ie_j+e_je_i=0$
pour $i\ne j$ ; en multipliant par $e_i$ à gauche, respectivement à droite, on obtient
$e_ie_j+e_ie_je_i=0=e_ie_je_i+e_je_i$, donc $e_ie_j=e_je_i=0$. Pour $x=\varphi(X)$, $ex+xe=2x$ donne
$exe=ex=xe=x$, donc $x=\sum_{p,q}x_{pq}$ avec $x_{pq}=e_pxe_q$. Pour $X=a[ij]$, les relations
$e_i\circ x=e_j\circ x=\frac12x$ et $e_k\circ x=0$ ($k\ne i,j$) s'écrivent $\frac12([p=i]+[q=i])x_{pq}=\frac12x_{pq}$,
etc., ce qui ne laisse que $x_{ij}$ et $x_{ji}$. Définissons, pour des couples ordonnés,
$\alpha_{ij}(a)=e_i\varphi(a[ij])e_j$. Pour $i,j,k$ distincts, $\varphi(a[ij])\varphi(b[jk])=\alpha_{ij}(a)\alpha_{jk}(b)$
est le seul terme dans $e_iAe_k$ de
$\varphi(a[ij])\varphi(b[jk])+\varphi(b[jk])\varphi(a[ij])=\varphi((ab)[ik])$, donc
$$
\alpha_{ij}(a)\,\alpha_{jk}(b)=\alpha_{ik}(ab).
$$
De $\varphi(1[13])^2=\varphi(E_{11}+E_{33})$, on tire $\alpha_{13}(1)\alpha_{31}(1)=e_1$,
$\alpha_{31}(1)\alpha_{13}(1)=e_3$ ; de $\varphi(a[13])^2=N(a)(e_1+e_3)$, $\alpha_{13}(a)\alpha_{31}(\bar a)=N(a)e_1$,
donc $\alpha_{13}$ est injective. Or
$$
\alpha_{13}(ay)=\alpha_{12}(a)\alpha_{23}(y)=\alpha_{13}(a)\alpha_{32}(1)\alpha_{21}(1)\alpha_{13}(y)=\alpha_{13}(a)\alpha_{31}(1)\alpha_{13}(y),
$$
et par conséquent, avec $\pi=\alpha_{31}(1)$ et l'associativité dans $A$,
$\alpha_{13}((ab)c)=\alpha_{13}(a)\pi\alpha_{13}(b)\pi\alpha_{13}(c)=\alpha_{13}(a(bc))$. Donc $(ab)c=a(bc)$. D'après le
lemme E.15(2), $\mathbb O$ n'est pas associative. $\square$

L'argument original d'Albert (A. A. Albert, Ann. of Math. 35 (1934) 65–73) est de ce type. Le code vérifie sur des
éléments aléatoires l'identité de Jordan de $\operatorname{Herm}_3(\mathbb O)$ ainsi que l'associativité et la
positivité de sa forme trace, montre numériquement que l'identité de Jordan est fausse dans
$\operatorname{Herm}_4(\mathbb O)$ (corollaire E.25), vérifie l'automorphisme du lemme E.23, et exécute les trois
étapes (a)–(c) de la démonstration du théorème E.24(2) sur un élément aléatoire :

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

## E.9 La classification

**Théorème E.27 (Jordan–von Neumann–Wigner).**

1. Toute algèbre de Jordan formellement réelle (de manière équivalente, euclidienne) de dimension finie est somme
   directe d'idéaux simples, et toute algèbre simple de ce type est isomorphe à $\operatorname{Herm}_n(\mathbb R)$,
   $\operatorname{Herm}_n(\mathbb C)$ ou
   $\operatorname{Herm}_n(\mathbb H)$ pour un $n\ge1$, à un facteur de spin $V_m$ avec $m\ge2$, ou à l'algèbre d'Albert
   $\operatorname{Herm}_3(\mathbb O)$.
2. Les algèbres $\mathbb R$ ; $V_m$ ($m\ge2$) ; $\operatorname{Herm}_n(D)$ ($D=\mathbb R,\mathbb C,\mathbb H$,
   $n\ge3$) ; $\operatorname{Herm}_3(\mathbb O)$ sont simples, formellement réelles et deux à deux non isomorphes, de
   sorte que toute algèbre de Jordan simple formellement réelle est isomorphe à exactement l'une d'elles. De plus
   $\operatorname{Herm}_1(D)\cong\mathbb R$ et
   $\operatorname{Herm}_2(D)\cong V_{\dim D+1}$ pour $D=\mathbb R,\mathbb C,\mathbb H,\mathbb O$.
3. Dans une algèbre de Jordan simple formellement réelle, tous les repères de Jordan ont le même nombre $r$ d'éléments
   et tous leurs espaces de Peirce $J_{ij}$ la même dimension $d$ : $r=1$ pour $\mathbb R$, $(r,d)=(2,m-1)$ pour $V_m$
   et $(n,\dim D)$ pour $\operatorname{Herm}_n(D)$, $n\ge2$.

*Démonstration.* (1) D'après le théorème E.9, $J$ est somme directe d'idéaux simples formellement réels. Un tel idéal
possède un repère de Jordan d'un certain cardinal $r$ (lemme E.6(2)) ; il est
$\cong\mathbb R=\operatorname{Herm}_1(\mathbb R)$ si $r=1$, $\cong V_{d+1}$ si $r=2$ (théorème E.13), et
$\cong\operatorname{Herm}_r(D)$ avec $D\ne\mathbb O$ si $r\ge4$ et $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$ si
$r=3$ (théorème E.22).

(2) Pour $\operatorname{Herm}_n(D)$, voir le théorème E.24. $V_m$ est une algèbre de Jordan d'après le lemme 13.3, et
euclidienne pour $\langle(a,u),(b,v)\rangle=ab+u\cdot v$, puisque
$\langle x\circ y,z\rangle=abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w$ est symétrique en $x=(a,u)$, $y=(b,v)$, $z=(c,w)$.
Pour un vecteur unitaire $w$, $c_{\pm w}=\frac12(1,\pm w)$ est un repère de Jordan (théorème 13.4(1) ; $c_{\pm w}$ est
primitif car une somme $c_{w'}+c_{-w'}$ de deux idempotents orthogonaux non nuls vaut $1$) avec
$J_{12}=\{(0,v):v\perp w\}$, de dimension $m-1\ge1$ (la somme des dimensions est $m+1$) ; donc $V_m$ est simple d'après
la proposition E.11. Pour distinguer les algèbres, soit $s(J)$ le plus grand nombre d'idempotents non nuls deux à deux
orthogonaux (fini d'après le lemme E.6(2)) ; $s$ et $\dim$ sont invariants par isomorphisme. Clairement
$s(\mathbb R)=1$. Dans $V_m$, les idempotents non nuls sont $1$ et les $c_w$, avec $1\circ c_w=c_w$ et
$c_w\circ c_{w'}=0$ seulement pour $w'=-w$ (théorème 13.4(1)) ; donc $s(V_m)=2$. Pour
$D\in\{\mathbb R,\mathbb C,\mathbb H\}$, $s(\operatorname{Herm}_n(D))=n$ : les $E_{ii}$ donnent $s\ge n$.
Réciproquement, des idempotents orthogonaux $p_k\ne0$ vérifient $p_kp_l=0$ pour $k\ne l$ (comme dans la démonstration
du théorème E.26), donc les sous-$D$-modules à droite $V_k=p_kD^n$ de $D^n$ sont non nuls et indépendants : si
$\sum_kv_k=0$ avec $v_k\in V_k$, alors $v_l=p_l\sum_kv_k=0$. Un $V_k$ non nul contient un $v\ne0$ et donc $vD$, de
dimension réelle $d$ ($D$ n'a pas de diviseurs de zéro) ; donc $sd\le nd$. Pour $\operatorname{Herm}_3(\mathbb O)$, les
$E_{ii}$ donnent $s\ge3$ ; et $s$ idempotents non nuls deux à deux orthogonaux, avec $1$ moins leur somme, se raffinent
comme dans la démonstration du théorème E.7 en un repère de Jordan à $r'\ge s$ éléments. D'après la proposition
E.11(2) et le théorème E.22, $27=r'+d'\frac{r'(r'-1)}2$ avec $d'\ge1$, et $d'\in\{1,2,4\}$ si $r'\ge4$. Pour
$r'=4,5,6$, le membre de droite appartient à $\{10,16,28\}$, $\{15,25,45\}$, $\{21,36,66\}$, et pour $r'\ge7$ il
dépasse $27$ ; donc $s\le r'=3$. Les couples $(s,\dim)$ sont $(1,1)$, $(2,m+1)$, $(n,\frac{n(n+1)}2)$, $(n,n^2)$,
$(n,n(2n-1))$ pour $n\ge3$, et $(3,27)$, tous différents, car $\frac{n(n+1)}2<n^2<n(2n-1)$ et $6,9,15\ne27$. Les
isomorphismes en rang $1$ et $2$ sont le théorème E.24(3).

(3) Si $F$ est un repère de Jordan à $r_F$ éléments d'une algèbre $J$ simple, la démonstration de (1), appliquée avec
$F$, donne un isomorphisme de $J$ sur un modèle de rang $r_F$ de (2), et le $s$ de ce modèle vaut $r_F$. Ainsi
$r_F=s(J)$ pour tout repère, et alors $d=(\dim J-r)/\binom r2$ pour $r\ge2$ (proposition E.11(2)). $\square$

**Corollaire E.28 (rang et trace).** Dans toute algèbre de Jordan formellement réelle, tous les repères de Jordan ont le
même nombre d'éléments, le **rang**. Si $J$ est simple de rang $r$, alors
$\operatorname{tr}x=\frac r{\dim J}\operatorname{Tr}L(x)$ vérifie $\operatorname{tr}c=1$ pour tout idempotent primitif
$c$, donc $\operatorname{tr}(\sum_i\lambda_ic_i)=\sum_i\lambda_i$ pour toute décomposition du théorème E.7, et
$\operatorname{tr}(x\circ y)$ est un produit scalaire associatif.

*Démonstration.* Soient $J=\bigoplus_pJ^{(p)}$ comme dans le théorème E.9 et $c$ primitif. Ses composantes
$c^{(p)}\in J^{(p)}$ sont des idempotents deux à deux orthogonaux ($c^2=\sum_p(c^{(p)})^2$), donc exactement l'une
d'elles est non nulle. Par conséquent, les membres d'un repère de Jordan de $J$ qui sont dans $J^{(p)}$ forment un
repère de Jordan de $J^{(p)}$ (leur somme est la composante $1^{(p)}$ de $1$), et le rang de $J$ est la somme des rangs
des $J^{(p)}$ (théorème E.27(3)). Dans une algèbre $J$ simple, un $c$ primitif appartient à un repère de Jordan (on
raffine $1-c$ comme dans la démonstration du théorème E.7), donc $\operatorname{Tr}L(c)=t=\dim J/r$ d'après la
proposition E.11(2) (pour $r=1$, $J=\mathbb Rc$ et $\operatorname{Tr}L(c)=1$). Enfin
$\operatorname{tr}(x\circ y)=\frac r{\dim J}\tau(x,y)$ (théorème E.5). $\square$

Le code calcule, pour chaque type, le repère de Jordan $(E_{ii})$ ou $(c_w,c_{-w})$, les dimensions des espaces de
Peirce comme sous-espaces propres communs des $L(c_i)$ (théorème E.8), la formule de dimension de la proposition
E.11(2) et la positivité de la forme trace $\tau$ (théorème E.5) :

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

## E.10 Ce que cette annexe démontre pour le chapitre 13

| énoncé du chapitre 13 | démontré ici |
|---|---|
| formellement réelle $\Leftrightarrow$ euclidienne (section 13.4) | théorème E.5 |
| théorème 13.7 : décomposition en idéaux simples, la liste, « exactement une » | théorèmes E.9, E.27 |
| $\operatorname{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}$ pour $\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O$ | théorème E.24(3) |
| $\operatorname{Herm}_r(\mathbb H)$ se plonge dans $\operatorname{Herm}_{2r}(\mathbb C)$ | théorème E.24(1) |
| l'algèbre d'Albert est une algèbre de Jordan et elle est exceptionnelle | théorèmes E.24(2), E.26 |
| $\operatorname{Herm}_r(\mathbb O)$ n'est une algèbre de Jordan que pour $r\le3$ | théorème E.24(2), corollaire E.25 |
| théorème spectral avec repères de Jordan, rang, $\operatorname{tr}z=\sum_i\mu_i$, associativité de $\operatorname{tr}(x\circ y)$ (section 13.6) | théorème E.7, corollaire E.28 |

**Notes.** La classification est due à P. Jordan, J. von Neumann et E. Wigner, *On an algebraic generalization of the
quantum mechanical formalism*, Ann. of Math. 35 (1934) 29–64, et l'exceptionnalité de
$\operatorname{Herm}_3(\mathbb O)$ à A. A. Albert (même volume, 65–73). La démonstration ci-dessus suit la structure de
Faraut–Korányi, *Analysis on Symmetric Cones* (Oxford 1994), chapitres II–V ; la démonstration du théorème de Hurwitz
par doublement (Nachr. Ges. Wiss. Göttingen 1898, 309–316) est celle de J. H. Conway, D. A. Smith, *On Quaternions and
Octonions* (2003) et de T. A. Springer, F. D. Veldkamp, *Octonions, Jordan Algebras and Exceptional Groups* (2000), et
la démonstration du théorème de Frobenius (1878) est celle de R. S. Palais, Amer. Math. Monthly 75 (1968) 366–368.
Pour aller plus loin (sans que cela soit utilisé nulle part dans le cours) : C. M. Glennie, Pacific J. Math. 16 (1966)
47–59, pour une autre approche de l'exceptionnalité par les identités polynomiales, et *A Taste of Jordan Algebras*
(2004) de McCrimmon pour l'histoire et la théorie générale.
