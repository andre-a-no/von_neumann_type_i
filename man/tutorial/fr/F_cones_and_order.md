[English](../F_cones_and_order.md) · [Español](../es/F_cones_and_order.md) · **Français** · [Deutsch](../de/F_cones_and_order.md) · [中文](../zh/F_cones_and_order.md) · [日本語](../ja/F_cones_and_order.md) · [Русский](../ru/F_cones_and_order.md)

# Tutoriel, annexe F. Cônes symétriques, isomorphismes d'ordre et applications de Jordan

Ce tutoriel accompagne l'[annexe F](../../lecture_notes/fr/F_cones_and_order.md) des notes de cours, qui démontre que le cône des carrés d'une algèbre de
Jordan euclidienne est symétrique, le théorème de Koecher–Vinberg ([théorème F.12](../../lecture_notes/fr/F_cones_and_order.md#r-F-12)), ainsi que le théorème de Kadison et
la forme des $*$-isomorphismes de Jordan pour les algèbres $\bigoplus_cM_{k_c}(\mathbb C)$
([théorèmes F.15](../../lecture_notes/fr/F_cones_and_order.md#r-F-15) et [F.18](../../lecture_notes/fr/F_cones_and_order.md#r-F-18)). Nous vérifions numériquement ces résultats sur des matrices, sur le facteur
de spin $V_4$ de la bibliothèque et sur le cône des matrices définies positives.

<a id="s-F-2"></a>

## F.2 Le cône des carrés est symétrique

*Notes de cours :* [section F.2](../../lecture_notes/fr/F_cones_and_order.md#s-F-2) ; [théorème F.5](../../lecture_notes/fr/F_cones_and_order.md#r-F-5), [théorème F.6](../../lecture_notes/fr/F_cones_and_order.md#r-F-6), [corollaire F.7](../../lecture_notes/fr/F_cones_and_order.md#r-F-7).

Dans $V_4$, le code utilise le `SpinFactor` de la bibliothèque ([tutoriel, section 13.2](../../tutorial/fr/13_jordan_algebras.md#s-13-2)) : `random` (avec `cone=True` pour des éléments de $Q$),
`inner`, `spectral_decomposition`, `in_cone`, `product`, `quadratic_representation`, `exp`, `log` et `apply_function`.

Le code vérifie le théorème F.5 pour $\mathrm{Herm}_4(\mathbb C)$ et pour le cône de Lorentz $Q\subseteq V_4$ (le
produit scalaire négatif est le témoin $c_i$ de (1), $\langle y,c_i\rangle=\lambda_i\|c_i\|^2<0$), le théorème F.6 dans $V_4$, où, dans les
coordonnées $(a,u)$, les opérateurs $L(b)$ et $P(y)$ sont des matrices symétriques $5\times5$, et le corollaire F.7 dans
$\mathrm{Herm}_5(\mathbb R)$ :

```python
import math
import torch
from torch_vn_algebra import SpinFactor

torch.manual_seed(0)
dt = torch.complex128
f64 = torch.float64
close = lambda A, B, tol=1e-10: (A - B).abs().max().item() < tol

G = torch.randn(2, 500, 4, 4, dtype=dt)
X, Y = G[0] @ G[0].mH, G[1] @ G[1].mH                             # squares in Herm_4(C)
assert ((X @ Y).diagonal(dim1=-2, dim2=-1).sum(-1).real > -1e-10).all()   # <Q, Q> >= 0
Z = torch.randn(500, 4, 4, dtype=dt)
Z = Z + Z.mH
w, U = torch.linalg.eigh(Z)
c1 = U[:, :, :1] @ U[:, :, :1].mH                                  # eigenprojection of the smallest eigenvalue
pair = (Z @ c1).diagonal(dim1=-2, dim2=-1).sum(-1).real            # <Z, c1> = lambda_min ||c1||^2
assert close(pair, w[:, 0]) and (pair < 0).float().mean() > 0.9   # witnesses of Z not in Q

V = SpinFactor(4)                                                  # <x, y> = tr(x o y) = V.inner
x, y = V.random(2000, cone=True), V.random(2000, cone=True)
assert (V.inner(x, y) > -1e-12).all()
z = V.random(2000)
lam, idem = V.spectral_decomposition(z)                            # (lambda-, lambda+), (c-, c+); ||c-||^2 = 1
assert close(V.inner(z, idem[:, 0]), lam[:, 0])
assert ((lam[:, 0] < 0) == ~V.in_cone(z, tol=0.0)).all()

I5 = torch.eye(5, dtype=f64)
Lmat = lambda b: V.product(b.expand(5, 5), I5).T                   # columns b o e_k
Pmat = lambda b: V.quadratic_representation(b.expand(5, 5), I5).T
for _ in range(20):
    b, yi, yg = V.random(), V.random(cone=True), V.random()
    E = torch.linalg.matrix_exp(Lmat(b))
    assert close(E @ V.unit(), V.exp(b))                           # Theorem F.6(1)
    assert V.in_cone((E @ x.T).T, tol=1e-9).all()
    assert close(Pmat(yi), torch.linalg.matrix_exp(2 * Lmat(V.log(yi))))      # Theorem F.6(3)
    lam_g, idem_g = V.spectral_decomposition(yg)
    s = (torch.where(lam_g >= 0, 1.0, -1.0)[:, None] * idem_g).sum(0)
    assert close(Pmat(yg), Pmat(s) @ Pmat(V.apply_function(yg, torch.abs)))  # P(y) = P(s) P(|y|)
    assert V.in_cone((Pmat(yg) @ x.T).T, tol=1e-9).all()                     # Theorem F.6(4)

xs = torch.randn(300, 5, 5, dtype=f64)
xs, ys = xs + xs.mT, torch.randn(300, 5, 5, dtype=f64)
ys = ys + ys.mT                                                    # Corollary F.7 in Herm_5(R), tau = Tr
absm = lambda A: (lambda w, U: U @ torch.diag_embed(w.abs()) @ U.mT)(*torch.linalg.eigh(A))
lhs = absm(ys @ xs @ ys).diagonal(dim1=-2, dim2=-1).sum(-1)
rhs = (ys @ absm(xs) @ ys).diagonal(dim1=-2, dim2=-1).sum(-1)
assert (rhs - lhs > -1e-9 * rhs).all()
```

<a id="s-F-3"></a>

## F.3 Le théorème de Koecher–Vinberg

*Notes de cours :* [section F.3](../../lecture_notes/fr/F_cones_and_order.md#s-F-3) ; [proposition F.11](../../lecture_notes/fr/F_cones_and_order.md#r-F-11), [théorème F.12](../../lecture_notes/fr/F_cones_and_order.md#r-F-12), [exemple F.13](../../lecture_notes/fr/F_cones_and_order.md#r-F-13).

Le code calcule $D^3f$ en $e$ par différentiation automatique de $f(x)=-\frac{n+1}2\log\det x$ en coordonnées
orthonormées de $\mathrm{Herm}_3(\mathbb R)$ et vérifie les énoncés de la proposition F.11 et du théorème F.12 :

```python
from torch.autograd.functional import jacobian, hessian

n = 3
pairs = [(i, j) for i in range(n) for j in range(i, n)]
Nd = len(pairs)                                                    # N = 6
Bs = torch.zeros(Nd, n, n, dtype=f64)                              # orthonormal basis for Tr(xy)
for k, (i, j) in enumerate(pairs):
    Bs[k, i, j] = Bs[k, j, i] = 1.0 if i == j else 1 / math.sqrt(2)
mat = lambda v: torch.einsum('k,kij->ij', v, Bs)
vec = lambda A: torch.einsum('kij,ij->k', Bs, A)
f = lambda v: -(n + 1) / 2 * torch.logdet(mat(v))                  # log phi up to a constant (Example F.13)
kappa = math.sqrt((n + 1) / 2)
e = vec(kappa * torch.eye(n, dtype=f64))
assert close(-jacobian(f, e), e)                                   # e* = e
assert close(hessian(f, e), torch.eye(Nd, dtype=f64))             # H(e) = I (Proposition F.11(2))
T = jacobian(lambda v: hessian(f, v, create_graph=True), e)        # T[a, b, c] = D^3 f(e)[B_a, B_b, B_c]
prod = lambda u, v: torch.einsum('a,b,abc->c', u, v, -0.5 * T)     # Theorem F.12
u, v = torch.randn(2, Nd, dtype=f64)
assert close(prod(e, v), v)                                        # unit
assert close(prod(u, v), vec(mat(u) @ mat(v) + mat(v) @ mat(u)) / (2 * kappa))
u2 = prod(u, u)
assert close(prod(u, prod(u2, v)), prod(u2, prod(u, v)))           # Jordan identity
```

<a id="s-F-4"></a>

## F.4 Isomorphismes d'ordre : le théorème de Kadison en dimension finie

*Notes de cours :* [section F.4](../../lecture_notes/fr/F_cones_and_order.md#s-F-4) ; [théorème F.15](../../lecture_notes/fr/F_cones_and_order.md#r-F-15), [théorème F.18](../../lecture_notes/fr/F_cones_and_order.md#r-F-18).

D'après [notes de cours, théorème F.15](../../lecture_notes/fr/F_cones_and_order.md#r-F-15), une bijection linéaire unifère $\Phi$ telle que $\Phi$ et $\Phi^{-1}$ soient
positives est un $*$-isomorphisme de Jordan. Le code construit une telle application de $M_2\oplus M_2\oplus M_3$, de la forme du
[théorème F.18](../../lecture_notes/fr/F_cones_and_order.md#r-F-18) : elle échange les deux blocs $M_2$, transpose deux blocs et conjugue par des unitaires. Il
vérifie que $\Phi$ est unifère et positive dans les deux sens, qu'elle préserve le produit de Jordan et qu'elle n'est pas
multiplicative :

```python
def blocks_rand(ks, psd=False):
    out = []
    for k in ks:
        A = torch.randn(k, k, dtype=dt)
        out.append(A @ A.mH if psd else A + A.mH)
    return out

ks = (2, 2, 3)                                                     # M = M_2 (+) M_2 (+) M_3
Us = [torch.linalg.qr(torch.randn(k, k, dtype=dt))[0] for k in ks]
def Phi(x):                                                        # swap the M_2 blocks, transpose two blocks
    x1, x2, x3 = x
    return [Us[0] @ x2.mT @ Us[0].mH, Us[1] @ x1 @ Us[1].mH, Us[2] @ x3.mT @ Us[2].mH]
def Phi_inv(y):
    y1, y2, y3 = y
    return [Us[1].mH @ y2 @ Us[1], (Us[0].mH @ y1 @ Us[0]).mT, (Us[2].mH @ y3 @ Us[2]).mT]
is_psd = lambda x: all(torch.linalg.eigvalsh(b)[0] > -1e-10 for b in x)
one = [torch.eye(k, dtype=dt) for k in ks]
assert all(close(a, b) for a, b in zip(Phi(one), one))            # unital
for _ in range(50):
    p = blocks_rand(ks, psd=True)
    assert is_psd(Phi(p)) and is_psd(Phi_inv(p))                   # positive in both directions
x, y = blocks_rand(ks), blocks_rand(ks)
assert all(close(a, b) for a, b in zip(Phi_inv(Phi(x)), x))
jordan = lambda a, b: [0.5 * (s @ t + t @ s) for s, t in zip(a, b)]
assert all(close(a, b) for a, b in zip(Phi(jordan(x, y)), jordan(Phi(x), Phi(y))))   # Theorem F.15
assert not close(Phi([s @ t for s, t in zip(x, y)])[0], (Phi(x)[0] @ Phi(y)[0]))      # not multiplicative
```

<a id="s-F-5"></a>

## F.5 $*$-isomorphismes de Jordan de $\bigoplus_cM_{k_c}$

*Notes de cours :* [section F.5](../../lecture_notes/fr/F_cones_and_order.md#s-F-5) ; [théorème F.17](../../lecture_notes/fr/F_cones_and_order.md#r-F-17).

Le code construit l'automorphisme de Jordan $x\mapsto Ux^TU^*$ de $M_3$ et retrouve $U$ (à une phase près) en suivant les
étapes de la démonstration ([notes de cours, théorème F.17](../../lecture_notes/fr/F_cones_and_order.md#r-F-17)) :

```python
k = 3
E = lambda i, j: torch.zeros(k, k, dtype=dt).index_put_((torch.tensor(i), torch.tensor(j)), torch.tensor(1.0, dtype=dt))
U0 = torch.linalg.qr(torch.randn(k, k, dtype=dt))[0]
Psi = lambda A: U0 @ A.mT @ U0.mH                                  # the unknown Jordan *-automorphism

W = torch.stack([torch.linalg.eigh(Psi(E(i, i)))[1][:, -1] for i in range(k)], dim=1)   # Psi(E_ii) = xi_i xi_i^*
Psi1 = lambda A: W.mH @ Psi(A) @ W
assert all(close(Psi1(E(i, i)), E(i, i)) for i in range(k))
transposed = Psi1(E(0, 1))[1, 0].abs().item() > 0.5               # epsilon = - : Psi1(E_01) is a multiple of E_10
Psi2 = (lambda A: Psi1(A.mT)) if transposed else Psi1
alpha = torch.stack([Psi2(E(i, 0))[i, 0] for i in range(k)])       # u_i = alpha_{i1}
for i in range(k):
    for j in range(k):
        assert close(Psi2(E(i, j)), alpha[i] * alpha[j].conj() * E(i, j))
Ufound = W @ torch.diag(alpha)
X = torch.randn(k, k, dtype=dt)
assert transposed and close(Psi(X), Ufound @ X.mT @ Ufound.mH)     # Theorem F.17, second form
phase = (U0.mH @ Ufound)[0, 0]
assert close(U0.mH @ Ufound, phase * torch.eye(k, dtype=dt)) and abs(abs(phase) - 1) < 1e-12
print("U recovered up to the phase", complex(phase))
```

<a id="x-F"></a>

## Exercices

[Exercice 2](../../lecture_notes/fr/F_cones_and_order.md#x-F) des notes de cours : pour $0<t<1$, $\Psi_t(x)=(1-t)x+\frac tn\mathrm{Tr}(x)1$ sur $M_n$ est une
bijection positive unifère dont l'inverse n'est pas positive, et elle n'envoie pas les projecteurs sur des projecteurs. Le
vérifier numériquement pour $n=3$, $t=\frac14$ et un projecteur de rang un $p$.

<details>
<summary>Solutions</summary>

```python
nn, t = 3, 0.25
p = E(0, 0)
Psi_t = lambda A: (1 - t) * A + t / nn * torch.trace(A) * torch.eye(nn, dtype=dt)
Psi_t_inv = lambda A: (A - t / nn * torch.trace(A) * torch.eye(nn, dtype=dt)) / (1 - t)
assert close(Psi_t_inv(Psi_t(X)), X)
assert abs(torch.linalg.eigvalsh(Psi_t_inv(p))[0].item() + t / (nn * (1 - t))) < 1e-12
assert not close(Psi_t(p) @ Psi_t(p), Psi_t(p))                   # not a projection
```

</details>
