[English](../05_random_matrices.md) · **Español** · [Français](../fr/05_random_matrices.md) · [Deutsch](../de/05_random_matrices.md) · [中文](../zh/05_random_matrices.md) · [日本語](../ja/05_random_matrices.md) · [Русский](../ru/05_random_matrices.md)

# Tutorial 5. Matrices aleatorias en un álgebra

Este tutorial acompaña al [capítulo 5](../../lecture_notes/es/05_random_matrices.md) de los apuntes de clase, que construye la medida de Haar
sobre el grupo unitario de $`\mathcal M=\bigoplus_cM_{k_c}`$ ([corolario 5.4](../../lecture_notes/es/05_random_matrices.md#r-5-4)), calcula sus momentos de primer y
segundo orden ([proposición 5.7](../../lecture_notes/es/05_random_matrices.md#r-5-7), [teorema 5.10](../../lecture_notes/es/05_random_matrices.md#r-5-10)) y muestra por qué la estadística de niveles debe
calcularse sector por sector ([teorema 5.15](../../lecture_notes/es/05_random_matrices.md#r-5-15)). Generamos unitarias de Haar y conjuntos gaussianos con
la biblioteca y comprobamos estos resultados mediante Monte Carlo: una media debe coincidir con el valor exacto
dentro de $`5`$ errores estándar. Los resultados que no se apoyan en un teorema se señalan como observaciones numéricas.

En todo el capítulo, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ con $`D=\sum_ck_c`$, $`P_c`$ son las
proyecciones centrales minimales y $`E`$ es la esperanza condicional sobre el centro de la
proposición 2.8. $`U(k)`$ es el grupo de las matrices unitarias $`k\times k`$ y
$`U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}`$ es el grupo unitario de $`\mathcal M`$. Las variables
aleatorias viven en un espacio de probabilidad con esperanza $`\mathbb E`$; la esperanza de una matriz
aleatoria se toma entrada a entrada. Una variable *gaussiana compleja estándar* tiene partes real e
imaginaria independientes $`\sim N(0,\tfrac12)`$, de modo que $`\mathbb E|z|^2=1`$ (esto es lo que
genera `torch.randn` para un dtype complejo, según su documentación; el primer bloque de código
comprueba los segundos momentos); una *matriz de Ginibre* tiene entradas independientes gaussianas
complejas estándar.

<a id="s-5-1"></a>

## 5.1 Medida de Haar

*Apuntes:* [sección 5.1](../../lecture_notes/es/05_random_matrices.md#s-5-1); [teorema 5.3](../../lecture_notes/es/05_random_matrices.md#r-5-3) (el factor unitario de la descomposición QR de una matriz de Ginibre, con las fases fijadas por $`\Lambda=\mathrm{diag}(R'_{ii}/|R'_{ii}|)`$, está distribuido según la medida de Haar), [corolario 5.4](../../lecture_notes/es/05_random_matrices.md#r-5-4) (los bloques de una unitaria de Haar de $`U(\mathcal M)`$ son independientes y están distribuidos según la medida de Haar) y [corolario 5.5](../../lecture_notes/es/05_random_matrices.md#r-5-5) ($`\mathbb EU_{ij}=0`$, $`\mathbb E|U_{ij}|^2=1/k`$). La medida de Haar de $`SU(k)`$ es la ley de $`e^{-i\theta/k}U`$, $`\theta=\arg\det U`$.

Sin la corrección $`\Lambda`$ del [teorema 5.3](../../lecture_notes/es/05_random_matrices.md#r-5-3), la ley de $`Q'`$ depende del convenio de signos de la
rutina QR. Para la rutina a la que llama `torch.linalg.qr` (LAPACK), el código siguiente encuentra
$`\mathbb E\mathrm{Re}\,\mathrm{Tr}\,Q'`$ a muchos errores estándar de $`0`$, mientras que una
unitaria de Haar tiene $`\mathbb EU=0`$ (corolario 5.5): así que este $`Q'`$ no está distribuido según
la medida de Haar. Esta es una afirmación sobre esa rutina, comprobada por el código, no un teorema.

La biblioteca genera unitarias de Haar de $`U(\mathcal M)`$ con `random_unitary_operator`, mediante
una única descomposición QR por lotes por cada tamaño de bloque, para todos los sectores y todas
las muestras a la vez (en una GPU, mediante Gram–Schmidt, que da la misma $`Q`$ por la unicidad del
lema 5.2). La medida `'haar'` es el teorema 5.3 en cada bloque; para un álgebra compleja,
`'haar_su'` multiplica por $`e^{-i\theta/k}`$ como antes; `'coe'` es el conjunto $`S=W^TW`$ del
ejercicio 3 (las demás medidas de la biblioteca no se usan en este curso). Comprobamos los
corolarios 5.4 y 5.5, y que el muestreador de la biblioteca es exactamente la QR corregida del
teorema 5.3:

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

z = torch.randn(B, dtype=dt)                                         # torch.randn draws standard complex Gaussians:
mc(z.real ** 2, 0.5, "E (Re z)^2")                                   # Re z, Im z ~ N(0, 1/2), uncorrelated
mc(z.imag ** 2, 0.5, "E (Im z)^2")
mc(z.real * z.imag, 0.0, "E Re z Im z")

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

*Observación numérica.* Sin la corrección, el código encuentra
$`\mathbb E\mathrm{Re}\,\mathrm{Tr}\,Q'\approx-0.84`$ para $`k=2`$, mientras que una unitaria de
Haar tiene $`\mathbb EU=0`$.

<a id="s-5-2"></a>

## 5.2 Momentos y el cálculo de Weingarten

*Apuntes:* [sección 5.2](../../lecture_notes/es/05_random_matrices.md#s-5-2); [proposición 5.7](../../lecture_notes/es/05_random_matrices.md#r-5-7) ($`\mathbb E[UXU^*]=E(X)`$), [lema 5.8](../../lecture_notes/es/05_random_matrices.md#r-5-8) (el conmutante de todos los $`V\otimes V`$ es $`\mathrm{span}\{1,F\}`$, con $`F`$ el intercambio), [teorema 5.10](../../lecture_notes/es/05_random_matrices.md#r-5-10) (la fórmula de Weingarten de orden dos) y [corolario 5.11](../../lecture_notes/es/05_random_matrices.md#r-5-11).

Comprobamos la proposición 5.7 frente a `channels.center_expectation`, el lema 5.8 mediante álgebra
lineal (el conmutante de tres $`V\otimes V`$ muestreados, calculado como un núcleo, contiene a $`F`$ y
tiene dimensión $`2`$, así que, por el lema 5.8, ya es $`\mathrm{span}\{1,F\}`$), y el
teorema 5.10 y el corolario 5.11 mediante Monte Carlo con el bloque $`U(3)`$ de las muestras
anteriores:

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

<a id="s-5-3"></a>

## 5.3 Estadística de niveles

*Apuntes:* [sección 5.3](../../lecture_notes/es/05_random_matrices.md#s-5-3): el GOE $`H=(G+G^T)/2`$, el GUE $`H=(G+G^*)/\sqrt2`$ y el cociente de espaciamientos $`r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})`$ de los espaciamientos $`s_i=\lambda_{i+1}-\lambda_i`$; [proposición 5.13](../../lecture_notes/es/05_random_matrices.md#r-5-13) (la conjetura de Wigner para matrices $`2\times2`$) y [proposición 5.14](../../lecture_notes/es/05_random_matrices.md#r-5-14) (el cociente medio es $`2\ln2-1`$ para niveles independientes, y $`4-2\sqrt3`$ y $`\frac{2\sqrt3}\pi-\frac12`$ para el GOE y el GUE $`3\times3`$, donde $`\mathbb Er=I_1/I_0`$ con $`I_j=\int_0^1x^jf_\beta(x)\,dx`$).

La función `SpinChain.level_spacing_ratio` (en `chains`) devuelve el cociente medio de cada fila de
niveles. Comprobamos las proposiciones 5.13 y 5.14 y registramos el valor del GOE para $`n=64`$:

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
mc(ratio(levels), 2 * math.log(2) - 1, "<r>, independent uniform levels")      # Proposition 5.14(a)

n3 = 100_000                                                         # Proposition 5.14(b): 3 x 3
G = torch.randn(n3, 1, 3, 3, dtype=torch.float64)
Gc = torch.randn(n3, 1, 3, 3, dtype=dt)
goe3 = TypeIAlgebra([3], [3], complex_valued=False, precision='double').operator((G + G.mT) / 2, is_self_adjoint=True)
gue3 = TypeIAlgebra([3], [3], precision='double').operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
h = 1e-4
xs = (torch.arange(int(1 / h), dtype=torch.float64) + 0.5) * h      # midpoint rule on [0, 1]
for beta, H, exact in ((1, goe3, 4 - 2 * math.sqrt(3)), (2, gue3, 2 * math.sqrt(3) / math.pi - 0.5)):
    f = (xs + xs ** 2) ** beta / (1 + xs + xs ** 2) ** (1 + 1.5 * beta)
    assert abs((xs * f).sum() / f.sum() - exact) < 1e-7              # E r = I_1 / I_0
    r3 = ratio(H.eigenvalues()[0])                                   # one ratio per matrix
    mc(r3, exact, f"beta={beta}: <r>, 3 x 3")
    mc((r3 <= 0.5).double(), (f * (xs <= 0.5)).sum().item() / f.sum().item(), f"beta={beta}: P(r <= 1/2), 3 x 3")

n = 64                                                               # numerical observation, no theorem
alg64 = TypeIAlgebra([n], [n], complex_valued=False, precision='double')
G = torch.randn(800, 1, n, n, dtype=torch.float64)
lam = alg64.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()[0]
r_goe = ratio(lam[:, n // 4: 3 * n // 4])                            # 800 samples of the mean ratio
print(f"<r>, GOE(64), middle half of the spectrum: {r_goe.mean():.4f} +- {r_goe.std() / math.sqrt(800):.4f}")

def same_mean(a, b, name):
    """Two independent Monte Carlo samples of the same law: means agree within 5 combined standard errors."""
    ea, eb = a.std().item() / math.sqrt(a.numel()), b.std().item() / math.sqrt(b.numel())
    print(f"{name}: {a.mean():.4f} vs {b.mean():.4f} +- {math.hypot(ea, eb):.4f}")
    assert abs(a.mean().item() - b.mean().item()) <= 5 * math.hypot(ea, eb)
```

*Observación numérica.* El código encuentra $`\langle r\rangle\approx0.53`$ para la mitad central del
espectro del GOE con $`n=64`$, por debajo del valor exacto $`3\times3`$, $`4-2\sqrt3\approx0.536`$, y muy
por encima del valor $`2\ln2-1\approx0.386`$ de los niveles independientes. Este número es un
resultado del cálculo, no un teorema; más abajo sirve de referencia para los bloques GOE del mismo
tamaño $`64`$, cuyos cocientes de espaciamientos tienen exactamente la misma ley.

<a id="s-5-4"></a>

## 5.4 Por qué no deben mezclarse los sectores

*Apuntes:* [sección 5.4](../../lecture_notes/es/05_random_matrices.md#s-5-4); por el [teorema 5.15](../../lecture_notes/es/05_random_matrices.md#r-5-15), los niveles de bloques independientes distintos no se repelen: para dos bloques GOE, el menor espaciamiento $`\delta`$ del espectro fusionado tiene $`P(\delta<\varepsilon)`$ al menos del orden de $`\varepsilon`$, mientras que dentro de un bloque GOE es a lo sumo del orden de $`\varepsilon^2`$.

Con más sectores el efecto es más fuerte. No demostramos ningún teorema límite para muchos
sectores; el experimento siguiente es una observación numérica. En la biblioteca, los bloques están
separados desde el principio: `eigenvalues()` devuelve un espectro por sector, y mezclarlos requiere
un `torch.cat` explícito. Ponemos $`C`$ bloques GOE independientes de tamaño $`64`$ en una misma álgebra
y comparamos el cociente medio por sector con el del espectro fusionado. Cada sector es una matriz
GOE de tamaño $`64`$, así que sus cocientes tienen exactamente la ley de `r_goe` de más arriba, y el
código lo comprueba; el comportamiento del espectro fusionado solo se observa:

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
    same_mean(per_sector.flatten(), r_goe, f"C = {C:2d}: sectors vs GOE(64)")   # every sector is GOE(64)
assert merged_r[1] > merged_r[2] > merged_r[4] > merged_r[16]       # numerical observation
assert abs(merged_r[16] - (2 * math.log(2) - 1)) < 0.01              # numerical observation
```

*Observación numérica.* Cada sector tiene $`\langle r\rangle\approx0.53`$, el valor GOE(64) de la
sección 5.3 (esta coincidencia es la comprobación de un hecho: los sectores son matrices GOE(64)).
El espectro fusionado baja a $`0.42`$ con dos sectores y a aproximadamente $`0.39`$, a menos de $`0.01`$
del valor $`2\ln2-1`$ de los niveles independientes (proposición 5.14(a)), con $`8`$ y $`16`$ sectores.
Que el espectro fusionado se acerque a niveles independientes cuando $`C`$ crece es lo que observa el
código; no se demuestra en este curso.

Lo mismo ocurre en un modelo físico. La cadena de Heisenberg con campos aleatorios (capítulo 9)
conserva la magnetización, así que su hamiltoniano está en $`\bigoplus_{N=0}^LM_{\binom LN}`$
(proposición 9.1). Sus bloques *no* son independientes —todos dependen de los mismos campos
aleatorios—, así que el teorema 5.15 no se aplica tal como está enunciado; el experimento siguiente
muestra el mismo efecto, como observación numérica.

```python
chain = SpinChain(10, complex_valued=False)                         # sectors N = 0, ..., 10
lam = chain.random_field_heisenberg(W=0.5, batch_size=20).eigenvalues()
k5 = chain.sector_dim(5)                                            # 252
r_sector = ratio(lam[5][:, k5 // 4: 3 * k5 // 4]).mean().item()
spectrum = torch.cat(lam, dim=-1).sort(dim=-1).values               # all 1024 levels mixed
r_all = ratio(spectrum[:, 256:768]).mean().item()
print(f"Heisenberg chain, L = 10: <r> in sector N = 5: {r_sector:.4f}, all sectors mixed: {r_all:.4f}")
assert r_sector > 0.49 and r_all < 0.42                             # numerical observation
```

*Observación numérica.* Dentro del sector $`N=5`$ el código encuentra $`\langle r\rangle\approx0.53`$,
cerca del valor GOE(64) de más arriba, y el espectro con todos los sectores mezclados da $`0.39`$,
cerca de $`2\ln2-1`$.

<a id="s-5-5"></a>

## 5.5 Con qué comparar

Unitarias de Haar: `scipy.stats.unitary_group` y `ortho_group`, las unitarias aleatorias de Haar
de QuTiP, y cualquier QR de `numpy`/`torch` seguida de la corrección de fases. Cálculo de
Weingarten: los paquetes simbólicos IntU (Mathematica; Z. Puchała, J. A. Miszczak), RTNI
(Mathematica y Python; M. Fukuda, R. König, I. Nechita, J. Phys. A 52 (2019), 425303) y el paquete
de Python Haarpy calculan momentos de cualquier orden de forma simbólica; el teorema 5.10 es su
caso de orden dos, y el teorema 5.12 (demostrado en el apéndice B), la fórmula general que
implementan. La estadística de niveles suele calcularse con unas pocas líneas de `numpy`, y los
códigos de diagonalización exacta como QuSpin proporcionan los espectros sector a sector. Lo que
añade la biblioteca es el grupo unitario de un álgebra como un único objeto: bloques de Haar
independientes de todos los tamaños generados en llamadas por lotes, el twirl identificado con la
esperanza condicional `channels.center_expectation`, y espectros devueltos por sector, de modo que
mezclar sectores nunca es la opción por defecto.

<a id="x-5"></a>

## Ejercicios

Las partes numéricas de los [ejercicios](../../lecture_notes/es/05_random_matrices.md#x-5) de los apuntes de clase:
1. Para $`U`$ de Haar en $`U(\mathcal M)`$, $`\mathcal M=M_1\oplus M_2\oplus M_3`$, compruebe que $`\mathbb E|\mathrm{Tr}\,U|^2=C=3`$ y
   $`\mathbb E|\mathrm{Tr}\,U|^4=17`$.
2. Calcule numéricamente la matriz de Weingarten de orden $`3`$ para $`k=3`$ ([teorema 5.12](../../lecture_notes/es/05_random_matrices.md#r-5-12)) y compruebe que
   $`\mathbb E|U_{11}|^6=\frac1{10}`$.
3. El **conjunto ortogonal circular** (`measure='coe'`) es $`S=W^TW`$ con $`W`$ de Haar en $`U(k)`$. Compruebe que el
   COE de la biblioteca es $`W^TW`$, que $`S`$ es unitaria y simétrica, que $`\mathbb E\mathrm{Tr}\,S^m=0`$ para $`m=1,2,3`$
   y que la densidad media de las fases propias es uniforme. Calcule $`\langle r\rangle`$ para las fases propias del
   COE con $`k=32`$ y compárelo con el valor del GOE de la sección 5.3 (una observación numérica).

<details>
<summary>Soluciones</summary>

*Ejercicio 1.*

```python
tr = U.trace                                                        # Tr U = sum_c Tr U_c, shape (B,)
mc(tr.abs() ** 2, 3.0, "E|Tr U|^2")
mc(tr.abs() ** 4, 17.0, "E|Tr U|^4")
```

*Ejercicio 2.*

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

*Ejercicio 3.*

   *Observación numérica.* El código encuentra $`\langle r\rangle\approx0.53`$ para el COE con
   $`k=32`$, a menos de $`0.015`$ del valor GOE(64) de la sección 5.3. Esta coincidencia es un resultado
   del cálculo; ningún teorema de este curso relaciona la estadística de espaciamientos del COE con
   la del GOE (sobre los conjuntos circulares, véase F. J. Dyson, J. Math. Phys. 3 (1962), 140–156).

```python
alg_s = TypeIAlgebra([32], [32], precision='double')
torch.manual_seed(2)
S = alg_s.random_unitary_operator(400, measure='coe').matrix[:, 0]
torch.manual_seed(2)
W = alg_s.random_unitary(32, batch_size=400)                         # same random numbers
assert close(S, W.mT @ W)                                            # the library's COE is W^T W
assert close(S, S.mT) and close(S @ S.mH, torch.eye(32, dtype=dt))
for m in (1, 2, 3):                                                  # E Tr S^m = 0
    trm = torch.linalg.matrix_power(S, m).diagonal(dim1=-2, dim2=-1).sum(-1)
    mc(trm.real, 0.0, f"E Re Tr S^{m}")
    mc(trm.imag, 0.0, f"E Im Tr S^{m}")
phases = torch.linalg.eigvals(S).angle().sort(dim=-1).values         # in (-pi, pi]
mc(((phases > 0) & (phases <= 1)).sum(-1).double(), 32 / (2 * math.pi), "E #eigenphases in (0, 1]")
s = torch.diff(torch.cat([phases, phases[:, :1] + 2 * math.pi], dim=-1), dim=-1)   # around the circle
r = torch.minimum(s, s.roll(-1, -1)) / torch.maximum(s, s.roll(-1, -1))
r_coe = r.mean(-1)
print(f"<r>, COE(32) eigenphases: {r_coe.mean():.4f}; GOE(64), Section 5.3: {r_goe.mean():.4f}")
assert abs(r_coe.mean() - r_goe.mean()) < 0.015                      # numerical observation
```

</details>

**Siguiente:** [Tutorial 6](06_quantum_channels.md): los canales cuánticos en la biblioteca.
