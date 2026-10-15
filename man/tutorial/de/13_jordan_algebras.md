[English](../13_jordan_algebras.md) · [Español](../es/13_jordan_algebras.md) · [Français](../fr/13_jordan_algebras.md) · **Deutsch** · [中文](../zh/13_jordan_algebras.md) · [日本語](../ja/13_jordan_algebras.md) · [Русский](../ru/13_jordan_algebras.md)

# Tutorial 13. Jordan-Algebren: warum die Quantenmechanik *-Algebren verwendet

Dieses Tutorial begleitet [Kapitel 13](../../lecture_notes/de/13_jordan_algebras.md) der Vorlesungsnotizen, das das Jordan-Produkt
$`x\circ y=\frac12(xy+yx)`$ und die quadratische Darstellung $`P(y)x=yxy`$, die Spinfaktoren $`V_m`$ und den Lorentz-Kegel untersucht
und zeigt, dass die Jordan-Struktur den Kommutator nur bis auf das Vorzeichen bestimmt, was zusammengesetzte Systeme
erkennen ([Proposition 13.12](../../lecture_notes/de/13_jordan_algebras.md#r-13-12)). Hier verwenden wir das Modul `torch_vn_algebra.jordan` und die Klasse `SpinFactor`, um
diese Ergebnisse numerisch zu überprüfen.

<a id="s-13-1"></a>

## 13.1 Observablen und das Jordan-Produkt

*Vorlesungsnotizen:* [Abschnitt 13.1](../../lecture_notes/de/13_jordan_algebras.md#s-13-1): [Lemma 13.1](../../lecture_notes/de/13_jordan_algebras.md#r-13-1) (die Jordan-Identität) und [Proposition 13.2](../../lecture_notes/de/13_jordan_algebras.md#r-13-2) (die quadratische Darstellung).

Das Modul `torch_vn_algebra.jordan` implementiert das Jordan-Produkt $`x\circ y`$, die Jordan-Potenzen und die quadratische
Darstellung $`P(y)=2L(y)^2-L(y^2)`$ auf Operatoren einer `TypeIAlgebra`, verzögert ausgewertet und gebündelt über alle Sektoren.
Der Code prüft [Lemma 13.1](../../lecture_notes/de/13_jordan_algebras.md#r-13-1) und [Proposition 13.2](../../lecture_notes/de/13_jordan_algebras.md#r-13-2) an zufälligen selbstadjungierten Elementen von $`M_2\oplus M_3`$:

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

## 13.2 Spinfaktoren und der Lorentz-Kegel

*Vorlesungsnotizen:* [Abschnitt 13.2](../../lecture_notes/de/13_jordan_algebras.md#s-13-2): [Lemma 13.3](../../lecture_notes/de/13_jordan_algebras.md#r-13-3) (antikommutierende Symmetrien), [Satz 13.4](../../lecture_notes/de/13_jordan_algebras.md#r-13-4) (Spektraltheorie von $`V_m`$) und [Proposition 13.5](../../lecture_notes/de/13_jordan_algebras.md#r-13-5).

Die Klasse `SpinFactor(m)` speichert Elemente $`x=(a,u)`$ von $`V_m`$ als Tensoren der Form `(batch, m+1)`, mit dem Produkt,
der Spektralzerlegung und dem Funktionalkalkül aus [Satz 13.4](../../lecture_notes/de/13_jordan_algebras.md#r-13-4); `to_matrix` ist die Abbildung $`\phi`$ aus
[Lemma 13.3](../../lecture_notes/de/13_jordan_algebras.md#r-13-3), aufgebaut aus den $`s_i`$ von `jordan.anticommuting_symmetries(m)`. Der Code prüft die Jordan-Identität
und Satz 13.4 in $`V_5`$, [Proposition 13.5](../../lecture_notes/de/13_jordan_algebras.md#r-13-5)(b) für ein zufälliges orthogonales $`O`$ und (c) für $`m\le6`$ sowie die
Majorana-Operatoren aus Abschnitt 9.4 für $`p=3`$ ([Vorlesungsnotizen, Abschnitt 13.2](../../lecture_notes/de/13_jordan_algebras.md#s-13-2)):

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

Für $`m=1`$ erzeugt $`s_1=Z`$ die Diagonalmatrizen, der Dimension $`2=m+1`$, und $`V_1\cong\mathbb R\oplus\mathbb R`$ ist
nicht einfach. Die ausgegebenen Dimensionen sind $`3,4`$ für $`m=2`$; $`5,16`$ für $`m=4`$; $`7,64`$ für $`m=6`$.

<a id="s-13-3"></a>

## 13.3 Jordan-Unteralgebren, die keine *-Unteralgebren sind

*Vorlesungsnotizen:* [Abschnitt 13.3](../../lecture_notes/de/13_jordan_algebras.md#s-13-3): [Proposition 13.6](../../lecture_notes/de/13_jordan_algebras.md#r-13-6).

`jordan.compare_with_star_algebra` berechnet die Dimensionen der von gegebenen hermiteschen Matrizen und $`1`$ erzeugten
Jordan-Algebra und des selbstadjungierten Teils der von ihnen erzeugten $`*`$-Algebra. Der Code prüft
[Proposition 13.6](../../lecture_notes/de/13_jordan_algebras.md#r-13-6)(a) für ein zufälliges Paar reeller symmetrischer Matrizen, $`n=2,\dots,5`$, und (b): Für zwei zufällige
hermitesche $`2\times2`$-Matrizen sind die Dimensionen $`3`$ und $`4`$, und die unitäre Matrix $`U`$ aus dem Beweis macht beide reell.

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

Drei generische hermitesche $`2\times2`$-Matrizen erzeugen ganz $`\mathrm{Herm}_2(\mathbb C)`$: Dimension $`4`$ in
beiden Spalten. In einem Block $`M_2`$ von $`\mathcal M`$ verraten zwei Observablen die komplexe Struktur also nie; eine
dritte kann es.

<a id="s-13-5"></a>

## 13.5 Jordan-Abbildungen und *-Abbildungen

*Vorlesungsnotizen:* [Abschnitt 13.5](../../lecture_notes/de/13_jordan_algebras.md#s-13-5): [Proposition 13.9](../../lecture_notes/de/13_jordan_algebras.md#r-13-9) (die Transposition), [Proposition 13.11](../../lecture_notes/de/13_jordan_algebras.md#r-13-11), [Proposition 13.12](../../lecture_notes/de/13_jordan_algebras.md#r-13-12) (zusammengesetzte Systeme) und [Proposition 13.13](../../lecture_notes/de/13_jordan_algebras.md#r-13-13) (Dynamik).

`jordan.transpose_map` ist die Transposition $`T`$, blockweise angewandt. Der Code prüft [Proposition 13.9](../../lecture_notes/de/13_jordan_algebras.md#r-13-9) (mit der
Choi-Matrix von $`T`$ auf $`M_3`$), die Matrix $`O_T=\mathrm{diag}(1,-1,1)`$ der Transposition auf Bloch-Vektoren aus dem
Beweis von [Proposition 13.11](../../lecture_notes/de/13_jordan_algebras.md#r-13-11), die Produktformel und $`\circ_\theta`$ aus [Proposition 13.12](../../lecture_notes/de/13_jordan_algebras.md#r-13-12) mit dem
Spektrum der partiellen Transposition des Bell-Zustands sowie [Proposition 13.13](../../lecture_notes/de/13_jordan_algebras.md#r-13-13)(a), (b):

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

## 13.6 Zurück zu den Spurungleichungen

*Vorlesungsnotizen:* [Abschnitt 13.6](../../lecture_notes/de/13_jordan_algebras.md#s-13-6): [Proposition 13.14](../../lecture_notes/de/13_jordan_algebras.md#r-13-14).

Nach [Vorlesungsnotizen, Proposition 13.14](../../lecture_notes/de/13_jordan_algebras.md#r-13-14) gilt $`\mathrm{tr}|P(y)x|\le\mathrm{tr}\,P(y)|x|`$ in $`V_m`$, mit Gleichheit, falls $`x\in Q`$ oder
$`-x\in Q`$. Über $`\phi`$, das $`\phi(P(y)x)=\phi(y)\phi(x)\phi(y)`$ und $`\phi(|x|)=|\phi(x)|`$ erfüllt, ist das auch die
Ungleichung $`\mathrm{Tr}|YXY|\le\mathrm{Tr}(Y|X|Y)`$ aus Proposition 11.9 für $`Y=Y^*`$ ([Vorlesungsnotizen, Abschnitt 13.6](../../lecture_notes/de/13_jordan_algebras.md#s-13-6)).
Der Code prüft beides für $`m=2,4,7`$:

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

Die Lücke $`\mathrm{tr}\,P(y)|x|-\mathrm{tr}|P(y)x|`$ verschwindet für definite $`x`$ (Proposition 13.14).
Numerische Beobachtung (kein Satz): Für diese Gaußschen Stichproben liegen die ausgegebenen maximalen Lücken in der
Größenordnung $`10^1`$–$`10^2`$; die Ungleichung ist für indefinite $`x`$ also weit von einer Gleichheit entfernt.

<a id="s-13-7"></a>

## 13.7 Womit man vergleichen sollte

Unseres Wissens gibt es keine verbreitete numerische Bibliothek für Jordan-Algebren; die Operationen dieses Kapitels
werden gewöhnlich von Hand mit Matrizen geschrieben. Computeralgebra deckt die algebraische Seite ab: SageMaths
`JordanAlgebra` konstruiert spezielle Jordan-Algebren, Jordan-Algebren symmetrischer Bilinearformen (zu denen die
Spinfaktoren gehören) und, ausgehend von einer Oktonionenalgebra, die 27-dimensionale exzeptionelle Jordan-Algebra;
auch Magma hat Konstruktoren für Jordan-Algebren, einschließlich der exzeptionellen. Das sind exakte algebraische
Werkzeuge, keine gebündelten numerischen. Die Lorentz-Kegel sind Standard in der Optimierung: SOCP-Löser (ECOS, SCS,
Clarabel, MOSEK, zugänglich über die `SOC`-Nebenbedingung von cvxpy) verwenden die Spektralzerlegung von $`V_m`$ (zur
Projektion auf den Kegel) oder seine quadratische Darstellung (für die Nesterov-Todd-Skalierung von
Innere-Punkte-Verfahren), und semidefinite Löser tun dasselbe für $`\mathrm{Herm}_n`$. Die Bibliothek fügt
Jordan-Produkte, die quadratische Darstellung und die Transposition als gebündelte, differenzierbare Operatoren auf
$`\bigoplus_cM_{k_c}`$ hinzu; Spinfaktoren mit ihrem Funktionalkalkül und ihrer Matrixdarstellung; und eine
randomisierte Berechnung der von gegebenen Matrizen erzeugten Jordan-Algebra, verglichen mit der von ihnen erzeugten
$`*`$-Algebra.

<a id="x-13"></a>

## Aufgaben

[Aufgaben](../../lecture_notes/de/13_jordan_algebras.md#x-13) 1–3 der Vorlesungsnotizen, numerisch:

- *Aufgabe 1* (Homogenität des Lorentz-Kegels): $`P(x^{1/2})1=x`$, $`P(y^{-1})P(y)=\mathrm{id}`$ und $`P(y)Q\subseteq Q`$
  in $`V_4`$.
- *Aufgabe 2* (isotrope Zustände): Der kleinste Eigenwert von $`\theta(\rho_p)`$ ist $`\frac{1-3p}4`$.
- *Aufgabe 3* (Produkte von Randpunkten): $`(1,w)\circ(1,w')\notin Q`$ für zufällige Einheitsvektoren $`w,w'`$, und
  $`(1,w)\circ(1,w)\in Q`$.

<details>
<summary>Lösungen</summary>

*Aufgabe 1.*

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

*Aufgabe 2.*

```python
for pp in (0.2, 1 / 3, 0.5):
    rho_p = pp * omega + (1 - pp) * torch.eye(4, dtype=dt) / 4
    lam_min = torch.linalg.eigvalsh(ptranspose(rho_p, 2, 2))[0].item()
    assert abs(lam_min - (1 - 3 * pp) / 4) < 1e-12
```

*Aufgabe 3.*

```python
V3c = SpinFactor(3)
wv = torch.nn.functional.normalize(torch.randn(2, 1000, 3, dtype=torch.float64), dim=-1)
prod = V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[1]))
assert not V3c.in_cone(prod).any()
assert V3c.in_cone(V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[0]))).all()
```

</details>

**Weiter:** [Tutorial, Anhang A. Tatsachen aus Analysis und Wahrscheinlichkeitstheorie](A_analysis.md)
