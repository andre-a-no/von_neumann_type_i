[English](../13_jordan_algebras.md) · **Español** · [Français](../fr/13_jordan_algebras.md) · [Deutsch](../de/13_jordan_algebras.md) · [中文](../zh/13_jordan_algebras.md) · [日本語](../ja/13_jordan_algebras.md) · [Русский](../ru/13_jordan_algebras.md)

# Tutorial 13. Álgebras de Jordan: por qué la mecánica cuántica usa *-álgebras

Este tutorial acompaña al [capítulo 13](../../lecture_notes/es/13_jordan_algebras.md) de los apuntes de clase, que estudia el producto de Jordan
$x\circ y=\frac12(xy+yx)$ y la representación cuadrática $P(y)x=yxy$, los factores de espín $V_m$ y el cono de Lorentz,
y muestra que la estructura de Jordan determina el conmutador solo salvo signo, algo que los sistemas compuestos
detectan ([proposición 13.12](../../lecture_notes/es/13_jordan_algebras.md#r-13-12)). Aquí usamos el módulo `torch_vn_algebra.jordan` y la clase `SpinFactor` para comprobar
numéricamente estos resultados.

<a id="s-13-1"></a>

## 13.1 Los observables y el producto de Jordan

*Apuntes:* [sección 13.1](../../lecture_notes/es/13_jordan_algebras.md#s-13-1): [lema 13.1](../../lecture_notes/es/13_jordan_algebras.md#r-13-1) (la identidad de Jordan) y [proposición 13.2](../../lecture_notes/es/13_jordan_algebras.md#r-13-2) (la representación cuadrática).

El módulo `torch_vn_algebra.jordan` implementa el producto de Jordan $x\circ y$, las potencias de Jordan y la representación
cuadrática $P(y)=2L(y)^2-L(y^2)$ sobre operadores de un `TypeIAlgebra`, de forma perezosa y en lotes sobre todos los sectores.
El código comprueba el [lema 13.1](../../lecture_notes/es/13_jordan_algebras.md#r-13-1) y la [proposición 13.2](../../lecture_notes/es/13_jordan_algebras.md#r-13-2) sobre elementos autoadjuntos aleatorios de $M_2\oplus M_3$:

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

## 13.2 Factores de espín y el cono de Lorentz

*Apuntes:* [sección 13.2](../../lecture_notes/es/13_jordan_algebras.md#s-13-2): [lema 13.3](../../lecture_notes/es/13_jordan_algebras.md#r-13-3) (simetrías que anticonmutan), [teorema 13.4](../../lecture_notes/es/13_jordan_algebras.md#r-13-4) (teoría espectral de $V_m$) y [proposición 13.5](../../lecture_notes/es/13_jordan_algebras.md#r-13-5).

La clase `SpinFactor(m)` almacena los elementos $x=(a,u)$ de $V_m$ como tensores de forma `(batch, m+1)`, con el producto,
la descomposición espectral y el cálculo funcional del [teorema 13.4](../../lecture_notes/es/13_jordan_algebras.md#r-13-4); `to_matrix` es la aplicación $\phi$ del
[lema 13.3](../../lecture_notes/es/13_jordan_algebras.md#r-13-3), construida a partir de las $s_i$ de `jordan.anticommuting_symmetries(m)`. El código comprueba la identidad de Jordan
y el teorema 13.4 en $V_5$, la [proposición 13.5](../../lecture_notes/es/13_jordan_algebras.md#r-13-5)(b) para una $O$ ortogonal aleatoria y (c) para $m\le6$, y los
operadores de Majorana de la sección 9.4 para $p=3$ ([apuntes de clase, sección 13.2](../../lecture_notes/es/13_jordan_algebras.md#s-13-2)):

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

Para $m=1$, $s_1=Z$ genera las matrices diagonales, de dimensión $2=m+1$, y $V_1\cong\mathbb R\oplus\mathbb R$ no es
simple. Las dimensiones impresas son $3,4$ para $m=2$; $5,16$ para $m=4$; $7,64$ para $m=6$.

<a id="s-13-3"></a>

## 13.3 Subálgebras de Jordan que no son *-subálgebras

*Apuntes:* [sección 13.3](../../lecture_notes/es/13_jordan_algebras.md#s-13-3): [proposición 13.6](../../lecture_notes/es/13_jordan_algebras.md#r-13-6).

`jordan.compare_with_star_algebra` calcula las dimensiones del álgebra de Jordan generada por unas matrices hermíticas
dadas y $1$, y de la parte autoadjunta de la $*$-álgebra que generan. El código comprueba la
[proposición 13.6](../../lecture_notes/es/13_jordan_algebras.md#r-13-6)(a) para un par aleatorio de matrices simétricas reales, $n=2,\dots,5$, y (b): para dos matrices
hermíticas $2\times2$ aleatorias, las dimensiones son $3$ y $4$, y la unitaria $U$ de la demostración hace reales a ambas.

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

Tres matrices hermíticas $2\times2$ genéricas generan todo $\mathrm{Herm}_2(\mathbb C)$: dimensión $4$ en ambas
columnas. Por tanto, en un bloque $M_2$ de $\mathcal M$, dos observables nunca revelan la estructura compleja; un
tercero puede hacerlo.

<a id="s-13-5"></a>

## 13.5 Aplicaciones de Jordan frente a *-aplicaciones

*Apuntes:* [sección 13.5](../../lecture_notes/es/13_jordan_algebras.md#s-13-5): [proposición 13.9](../../lecture_notes/es/13_jordan_algebras.md#r-13-9) (la traspuesta), [proposición 13.11](../../lecture_notes/es/13_jordan_algebras.md#r-13-11), [proposición 13.12](../../lecture_notes/es/13_jordan_algebras.md#r-13-12) (sistemas compuestos) y [proposición 13.13](../../lecture_notes/es/13_jordan_algebras.md#r-13-13) (dinámica).

`jordan.transpose_map` es la traspuesta $T$, aplicada bloque a bloque. El código comprueba la [proposición 13.9](../../lecture_notes/es/13_jordan_algebras.md#r-13-9) (con la
matriz de Choi de $T$ sobre $M_3$), la matriz $O_T=\mathrm{diag}(1,-1,1)$ de la traspuesta sobre los vectores de Bloch,
de la demostración de la [proposición 13.11](../../lecture_notes/es/13_jordan_algebras.md#r-13-11), la fórmula del producto y $\circ_\theta$ de la [proposición 13.12](../../lecture_notes/es/13_jordan_algebras.md#r-13-12) con el
espectro de la traspuesta parcial del estado de Bell, y la [proposición 13.13](../../lecture_notes/es/13_jordan_algebras.md#r-13-13)(a), (b):

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

## 13.6 De vuelta a las desigualdades de trazas

*Apuntes:* [sección 13.6](../../lecture_notes/es/13_jordan_algebras.md#s-13-6): [proposición 13.14](../../lecture_notes/es/13_jordan_algebras.md#r-13-14).

Según los [apuntes de clase, proposición 13.14](../../lecture_notes/es/13_jordan_algebras.md#r-13-14), $\mathrm{tr}|P(y)x|\le\mathrm{tr}\,P(y)|x|$ en $V_m$, con igualdad si $x\in Q$ o
$-x\in Q$. A través de $\phi$, que cumple $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ y $\phi(|x|)=|\phi(x)|$, es también la
desigualdad $\mathrm{Tr}|YXY|\le\mathrm{Tr}(Y|X|Y)$ de la proposición 11.9 para $Y=Y^*$ ([apuntes de clase, sección 13.6](../../lecture_notes/es/13_jordan_algebras.md#s-13-6)).
El código comprueba ambas para $m=2,4,7$:

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

La diferencia $\mathrm{tr}\,P(y)|x|-\mathrm{tr}|P(y)x|$ se anula para $x$ definido (proposición 13.14).
Observación numérica (no un teorema): para estas muestras gaussianas, las diferencias máximas impresas son del orden de
$10^1$–$10^2$, así que la desigualdad está lejos de ser una igualdad para $x$ indefinido.

<a id="s-13-7"></a>

## 13.7 Con qué comparar

Hasta donde sabemos, no existe ninguna biblioteca numérica de uso extendido para álgebras de Jordan; las operaciones
de este capítulo suelen escribirse a mano en términos de matrices. El álgebra computacional cubre el lado algebraico:
`JordanAlgebra` de SageMath construye álgebras de Jordan especiales, álgebras de Jordan de formas bilineales simétricas
(que incluyen los factores de espín) y, a partir de un álgebra de octoniones, el álgebra de Jordan excepcional de
dimensión 27; Magma también tiene constructores de álgebras de Jordan, incluida la excepcional. Son herramientas
algebraicas exactas, más que numéricas y en lotes. Los conos de Lorentz son estándar en optimización: los solvers de
SOCP (ECOS, SCS, Clarabel, MOSEK, accesibles mediante la restricción `SOC` de cvxpy) usan la descomposición espectral de
$V_m$ (para proyectar sobre el cono) o su representación cuadrática (para el escalado de Nesterov–Todd de los métodos de
punto interior), y los solvers semidefinidos hacen lo mismo para $\mathrm{Herm}_n$. La biblioteca añade productos
de Jordan, la representación cuadrática y la traspuesta como operadores en lotes y diferenciables sobre
$\bigoplus_cM_{k_c}$; factores de espín con su cálculo funcional y su representación matricial; y un cálculo
aleatorizado del álgebra de Jordan generada por unas matrices dadas, comparada con la $*$-álgebra que generan.

<a id="x-13"></a>

## Ejercicios

[Ejercicios](../../lecture_notes/es/13_jordan_algebras.md#x-13) 1–3 de los apuntes de clase, numéricamente:

- *Ejercicio 1* (homogeneidad del cono de Lorentz): $P(x^{1/2})1=x$, $P(y^{-1})P(y)=\mathrm{id}$ y $P(y)Q\subseteq Q$
  en $V_4$.
- *Ejercicio 2* (estados isótropos): el menor valor propio de $\theta(\rho_p)$ es $\frac{1-3p}4$.
- *Ejercicio 3* (productos de puntos de la frontera): $(1,w)\circ(1,w')\notin Q$ para vectores unitarios aleatorios
  $w,w'$, y $(1,w)\circ(1,w)\in Q$.

<details>
<summary>Soluciones</summary>

*Ejercicio 1.*

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

*Ejercicio 2.*

```python
for pp in (0.2, 1 / 3, 0.5):
    rho_p = pp * omega + (1 - pp) * torch.eye(4, dtype=dt) / 4
    lam_min = torch.linalg.eigvalsh(ptranspose(rho_p, 2, 2))[0].item()
    assert abs(lam_min - (1 - 3 * pp) / 4) < 1e-12
```

*Ejercicio 3.*

```python
V3c = SpinFactor(3)
wv = torch.nn.functional.normalize(torch.randn(2, 1000, 3, dtype=torch.float64), dim=-1)
prod = V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[1]))
assert not V3c.in_cone(prod).any()
assert V3c.in_cone(V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[0]))).all()
```

</details>

**Siguiente:** [Tutorial, apéndice A. Hechos de análisis y de probabilidad](A_analysis.md)
