[English](../13_jordan_algebras.md) · [Español](../es/13_jordan_algebras.md) · **Français** · [Deutsch](../de/13_jordan_algebras.md) · [中文](../zh/13_jordan_algebras.md) · [日本語](../ja/13_jordan_algebras.md) · [Русский](../ru/13_jordan_algebras.md)

# Tutoriel 13. Algèbres de Jordan : pourquoi la mécanique quantique utilise des *-algèbres

Ce tutoriel accompagne le [chapitre 13](../../lecture_notes/fr/13_jordan_algebras.md) des notes de cours, qui étudie le produit de Jordan
$x\circ y=\frac12(xy+yx)$ et la représentation quadratique $P(y)x=yxy$, les facteurs de spin $V_m$ et le cône de Lorentz,
et montre que la structure de Jordan ne détermine le commutateur qu'au signe près, ce que les systèmes composés
détectent ([proposition 13.12](../../lecture_notes/fr/13_jordan_algebras.md#r-13-12)). Nous utilisons ici le module `torch_vn_algebra.jordan` et la classe `SpinFactor` pour vérifier
numériquement ces résultats.

<a id="s-13-1"></a>

## 13.1 Observables et produit de Jordan

*Notes de cours :* [section 13.1](../../lecture_notes/fr/13_jordan_algebras.md#s-13-1) : [lemme 13.1](../../lecture_notes/fr/13_jordan_algebras.md#r-13-1) (l'identité de Jordan) et [proposition 13.2](../../lecture_notes/fr/13_jordan_algebras.md#r-13-2) (la représentation quadratique).

Le module `torch_vn_algebra.jordan` implémente le produit de Jordan $x\circ y$, les puissances de Jordan et la représentation
quadratique $P(y)=2L(y)^2-L(y^2)$ sur les opérateurs d'une `TypeIAlgebra`, de façon paresseuse et par lots sur tous les secteurs.
Le code vérifie le [lemme 13.1](../../lecture_notes/fr/13_jordan_algebras.md#r-13-1) et la [proposition 13.2](../../lecture_notes/fr/13_jordan_algebras.md#r-13-2) sur des éléments autoadjoints aléatoires de $M_2\oplus M_3$ :

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

<a id="s-13-2"></a>

## 13.2 Facteurs de spin et cône de Lorentz

*Notes de cours :* [section 13.2](../../lecture_notes/fr/13_jordan_algebras.md#s-13-2) : [lemme 13.3](../../lecture_notes/fr/13_jordan_algebras.md#r-13-3) (symétries anticommutantes), [théorème 13.4](../../lecture_notes/fr/13_jordan_algebras.md#r-13-4) (théorie spectrale de $V_m$) et [proposition 13.5](../../lecture_notes/fr/13_jordan_algebras.md#r-13-5).

La classe `SpinFactor(m)` stocke les éléments $x=(a,u)$ de $V_m$ sous forme de tenseurs de forme `(batch, m+1)`, avec le produit,
la décomposition spectrale et le calcul fonctionnel du [théorème 13.4](../../lecture_notes/fr/13_jordan_algebras.md#r-13-4) ; `to_matrix` est l'application $\phi$ du
[lemme 13.3](../../lecture_notes/fr/13_jordan_algebras.md#r-13-3), construite à partir des $s_i$ de `jordan.anticommuting_symmetries(m)`. Le code vérifie l'identité de Jordan
et le théorème 13.4 dans $V_5$, la [proposition 13.5](../../lecture_notes/fr/13_jordan_algebras.md#r-13-5)(b) pour un $O$ orthogonal aléatoire et (c) pour $m\le6$, ainsi que les
opérateurs de Majorana de la section 9.4 pour $p=3$ ([notes de cours, section 13.2](../../lecture_notes/fr/13_jordan_algebras.md#s-13-2)) :

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

Pour $m=1$, $s_1=Z$ engendre les matrices diagonales, de dimension $2=m+1$, et $V_1\cong\mathbb R\oplus\mathbb R$
n'est pas simple. Les dimensions affichées sont $3,4$ pour $m=2$ ; $5,16$ pour $m=4$ ; $7,64$ pour $m=6$.

<a id="s-13-3"></a>

## 13.3 Sous-algèbres de Jordan qui ne sont pas des *-sous-algèbres

*Notes de cours :* [section 13.3](../../lecture_notes/fr/13_jordan_algebras.md#s-13-3) : [proposition 13.6](../../lecture_notes/fr/13_jordan_algebras.md#r-13-6).

`jordan.compare_with_star_algebra` calcule les dimensions de l'algèbre de Jordan engendrée par des matrices hermitiennes
données et $1$, et de la partie autoadjointe de la $*$-algèbre qu'elles engendrent. Le code vérifie la
[proposition 13.6](../../lecture_notes/fr/13_jordan_algebras.md#r-13-6)(a) pour un couple aléatoire de matrices symétriques réelles, $n=2,\dots,5$, et (b) : pour deux matrices
hermitiennes $2\times2$ aléatoires, les dimensions sont $3$ et $4$, et l'unitaire $U$ de la démonstration les rend toutes deux réelles.

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

Trois matrices hermitiennes $2\times2$ génériques engendrent $\mathrm{Herm}_2(\mathbb C)$ tout entier :
dimension $4$ dans les deux colonnes. Dans un bloc $M_2$ de $\mathcal M$, deux observables ne révèlent donc jamais la
structure complexe ; une troisième le peut.

<a id="s-13-5"></a>

## 13.5 Applications de Jordan et *-applications

*Notes de cours :* [section 13.5](../../lecture_notes/fr/13_jordan_algebras.md#s-13-5) : [proposition 13.9](../../lecture_notes/fr/13_jordan_algebras.md#r-13-9) (la transposée), [proposition 13.11](../../lecture_notes/fr/13_jordan_algebras.md#r-13-11), [proposition 13.12](../../lecture_notes/fr/13_jordan_algebras.md#r-13-12) (systèmes composés) et [proposition 13.13](../../lecture_notes/fr/13_jordan_algebras.md#r-13-13) (dynamique).

`jordan.transpose_map` est la transposée $T$, appliquée bloc par bloc. Le code vérifie la [proposition 13.9](../../lecture_notes/fr/13_jordan_algebras.md#r-13-9) (avec la
matrice de Choi de $T$ sur $M_3$), la matrice $O_T=\mathrm{diag}(1,-1,1)$ de la transposée sur les vecteurs de Bloch, tirée de la
démonstration de la [proposition 13.11](../../lecture_notes/fr/13_jordan_algebras.md#r-13-11), la formule du produit et $\circ_\theta$ de la [proposition 13.12](../../lecture_notes/fr/13_jordan_algebras.md#r-13-12) avec le
spectre de la transposée partielle de l'état de Bell, et la [proposition 13.13](../../lecture_notes/fr/13_jordan_algebras.md#r-13-13)(a), (b) :

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

<a id="s-13-6"></a>

## 13.6 Retour aux inégalités de traces

*Notes de cours :* [section 13.6](../../lecture_notes/fr/13_jordan_algebras.md#s-13-6) : [proposition 13.14](../../lecture_notes/fr/13_jordan_algebras.md#r-13-14).

D'après [notes de cours, proposition 13.14](../../lecture_notes/fr/13_jordan_algebras.md#r-13-14), $\mathrm{tr}|P(y)x|\le\mathrm{tr}\,P(y)|x|$ dans $V_m$, avec égalité si $x\in Q$ ou
$-x\in Q$. Via $\phi$, qui vérifie $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ et $\phi(|x|)=|\phi(x)|$, c'est aussi
l'inégalité $\mathrm{Tr}|YXY|\le\mathrm{Tr}(Y|X|Y)$ de la proposition 11.9 pour $Y=Y^*$ ([notes de cours, section 13.6](../../lecture_notes/fr/13_jordan_algebras.md#s-13-6)).
Le code vérifie les deux pour $m=2,4,7$ :

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

L'écart $\mathrm{tr}\,P(y)|x|-\mathrm{tr}|P(y)x|$ s'annule pour $x$ défini (proposition 13.14). Observation
numérique (pas un théorème) : pour ces échantillons gaussiens, les écarts maximaux affichés sont de l'ordre de
$10^1$–$10^2$, donc l'inégalité est loin d'être une égalité pour $x$ indéfini.

<a id="s-13-7"></a>

## 13.7 Points de comparaison

À notre connaissance, il n'existe pas de bibliothèque numérique courante pour les algèbres de Jordan ; les opérations
de ce chapitre sont d'ordinaire écrites à la main en termes de matrices. Le calcul formel couvre le côté algébrique :
`JordanAlgebra` de SageMath construit des algèbres de Jordan spéciales, des algèbres de Jordan de formes bilinéaires
symétriques (qui comprennent les facteurs de spin) et, à partir d'une algèbre d'octonions, l'algèbre de Jordan
exceptionnelle de dimension 27 ; Magma dispose lui aussi de constructeurs d'algèbres de Jordan, y compris
l'exceptionnelle. Ce sont des outils algébriques exacts plutôt que des outils numériques par lots. Les cônes de
Lorentz sont standard en optimisation : les solveurs SOCP (ECOS, SCS, Clarabel, MOSEK, accessibles par la contrainte
`SOC` de cvxpy) utilisent la décomposition spectrale de $V_m$ (pour projeter sur le cône) ou sa représentation
quadratique (pour la mise à l'échelle de Nesterov–Todd des méthodes de points intérieurs), et les solveurs
semi-définis font de même pour $\mathrm{Herm}_n$. La bibliothèque ajoute les produits de Jordan, la
représentation quadratique et la transposée comme opérateurs par lots et différentiables sur $\bigoplus_cM_{k_c}$ ;
les facteurs de spin avec leur calcul fonctionnel et leur représentation matricielle ; et un calcul randomisé de
l'algèbre de Jordan engendrée par des matrices données, comparée à la $*$-algèbre qu'elles engendrent.

<a id="x-13"></a>

## Exercices

[Exercices](../../lecture_notes/fr/13_jordan_algebras.md#x-13) 1–3 des notes de cours, numériquement :

- *Exercice 1* (homogénéité du cône de Lorentz) : $P(x^{1/2})1=x$, $P(y^{-1})P(y)=\mathrm{id}$ et $P(y)Q\subseteq Q$
  dans $V_4$.
- *Exercice 2* (états isotropes) : la plus petite valeur propre de $\theta(\rho_p)$ est $\frac{1-3p}4$.
- *Exercice 3* (produits de points du bord) : $(1,w)\circ(1,w')\notin Q$ pour des vecteurs unitaires aléatoires $w,w'$, et
  $(1,w)\circ(1,w)\in Q$.

<details>
<summary>Solutions</summary>

*Exercice 1.*

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

*Exercice 2.*

```python
for pp in (0.2, 1 / 3, 0.5):
    rho_p = pp * omega + (1 - pp) * torch.eye(4, dtype=dt) / 4
    lam_min = torch.linalg.eigvalsh(ptranspose(rho_p, 2, 2))[0].item()
    assert abs(lam_min - (1 - 3 * pp) / 4) < 1e-12
```

*Exercice 3.*

```python
V3c = SpinFactor(3)
wv = torch.nn.functional.normalize(torch.randn(2, 1000, 3, dtype=torch.float64), dim=-1)
prod = V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[1]))
assert not V3c.in_cone(prod).any()
assert V3c.in_cone(V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[0]))).all()
```

</details>

**Suite :** [Tutoriel, annexe A. Résultats d'analyse et de probabilités](A_analysis.md)
