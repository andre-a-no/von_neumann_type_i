[English](../01_matrix_algebras.md) · **Español** · [Français](../fr/01_matrix_algebras.md) · [Deutsch](../de/01_matrix_algebras.md) · [中文](../zh/01_matrix_algebras.md) · [日本語](../ja/01_matrix_algebras.md) · [Русский](../ru/01_matrix_algebras.md)

# Tutorial 1. Álgebras de matrices y cálculo funcional

Este tutorial acompaña al [capítulo 1](../../lecture_notes/es/01_matrix_algebras.md) de los apuntes de clase, que demuestra el teorema
espectral, el cálculo funcional, la descomposición polar y la dualidad de la norma traza en $`M_n(\mathbb C)`$.
Aquí se comprueba cada uno de estos resultados con la biblioteca.

**Aprenderá** a construir operadores con un espectro prescrito, a aplicar una función a un
operador, cómo se almacenan los lotes de operadores y cómo contrastar un cálculo con un
resultado exacto.

<a id="s-1-5"></a>

## 1.5 Operadores en la biblioteca

*Apuntes:* [sección 1.2](../../lecture_notes/es/01_matrix_algebras.md#s-1-2) y [sección 1.3](../../lecture_notes/es/01_matrix_algebras.md#s-1-3): [teorema 1.3](../../lecture_notes/es/01_matrix_algebras.md#r-1-3), [corolario 1.4](../../lecture_notes/es/01_matrix_algebras.md#r-1-4), [teorema 1.5](../../lecture_notes/es/01_matrix_algebras.md#r-1-5), [teorema 1.6](../../lecture_notes/es/01_matrix_algebras.md#r-1-6), [teorema 1.11](../../lecture_notes/es/01_matrix_algebras.md#r-1-11), [proposición 1.12](../../lecture_notes/es/01_matrix_algebras.md#r-1-12), [proposición 1.13](../../lecture_notes/es/01_matrix_algebras.md#r-1-13).

La biblioteca trabaja con álgebras $`\bigoplus_cM_{k_c}(\mathbb C)`$ (capítulo 2). Un único
bloque corresponde al caso $`C=1`$: `TypeIAlgebra([n], [n])` es $`M_n(\mathbb C)`$. (Las dos
listas son los tamaños de bloque nominales y los reales; por ahora coinciden.)

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
M3 = TypeIAlgebra([3], [3], precision='double')       # the algebra M_3(C), complex128

# a positive operator with eigenvalues 1, 2, 4 and Haar-random eigenvectors
A = M3.operator_from_eigenvalues(lambda d: torch.tensor([[1.0, 2.0, 4.0]]),
                                 force_positive=True)
print(A.matrix.shape)          # (batch, C, k_max, k_max) = (1, 1, 3, 3)
print(A.eigenvalues()[0])      # one tensor per block, eigenvalues in ascending order
```

Cada operador se almacena como un único tensor de forma `(batch, C, k_max, k_max)`: un lote
de operadores, cada uno de ellos una lista de `C` bloques. Aquí `batch = C = 1`. La
construcción es exactamente el teorema 1.3 leído al revés: $`A=U\Lambda U^*`$ con $`U`$
distribuida según la medida de Haar (la biblioteca genera $`U`$ mediante la construcción QR del
teorema 5.3). El indicador `force_positive=True` se verifica sobre los valores propios prescritos
(la biblioteca lanza un error si alguno de ellos es menor que $`-10^{-12}`$); por el corolario 1.4(3),
el resultado es entonces positivo (de forma exacta cuando los valores propios prescritos son $`\ge0`$,
como aquí; en general, salvo esta tolerancia), y la biblioteca lo registra en lugar de comprobarlo
de nuevo. Por la misma razón, `A.eigenvalues()` devuelve más arriba el espectro prescrito,
ordenado, sin calcular nada; no es una comprobación. La comprobación de
$`\mathrm{spec}\,f(A)=f(\mathrm{spec}\,A)`$ de más abajo sí calcula: `apply_function`
diagonaliza la matriz de $`A`$, y los valores propios de la matriz resultante $`e^A`$ se calculan
desde cero.

El cálculo funcional para operadores autoadjuntos es `apply_function` (diagonaliza cada
bloque con `eigh`). La función $`f`$ debe actuar elemento a elemento: recibe los valores propios de
todos los sectores de un mismo tamaño $`k`$ como un único tensor de forma `(batch, n·k)`, donde $`n`$
es el número de tales sectores, y devuelve un tensor de la misma forma. Las funciones habituales
tienen sus propios métodos. Comprobamos los teoremas 1.5 y 1.6:

```python
def close(X, Y, tol=1e-10):
    """Maximal entry of |X - Y| for two operators (or tensors) is below tol."""
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

R = A.sqrt()
assert close(R @ R, A)                                  # Theorem 1.6
assert close(A.log().expm(), A)                         # exp(log A) = A, Theorem 1.5(e)
assert close(A.apply_function(lambda w: w ** 3), A @ A @ A)   # Theorem 1.5(b)
assert close(A.power(-1), A.inverse())                  # f(t) = 1/t
spec = A.apply_function(torch.exp).eigenvalues()[0]     # spec f(A) = f(spec A)
assert torch.allclose(spec, torch.exp(torch.tensor([[1.0, 2.0, 4.0]], dtype=spec.dtype)))
```

`R @ R` todavía no calcula nada: los productos, las sumas y las funciones construyen una
receta, que se evalúa cuando se solicita `.matrix` o un número, como una traza.

**Derivadas.** Por el [teorema 1.11](../../lecture_notes/es/01_matrix_algebras.md#r-1-11) (Daleckii–Krein), la derivada de
$`f(A)`$ en una dirección autoadjunta $`E`$ es $`U(L\circ U^*EU)U^*`$, con la matriz de Loewner
$`L_{ij}=f^{[1]}(\lambda_i,\lambda_j)`$ de diferencias divididas ($`f'(\lambda_i)`$ para valores propios iguales).
`functional.divided_differences` calcula $`L`$ a partir de los valores propios y de los valores de $`f`$ y $`f'`$, y
autograd deriva `apply_function` con ella. Comprobamos la fórmula en un espectro degenerado frente a un cociente de
diferencias centrado de la exponencial densa, cuyo error es $`O(s^2)`$ más un redondeo del orden de $`\varepsilon/s`$:

```python
from torch_vn_algebra import functional

lam = torch.tensor([1.0, 1.0, 2.0], dtype=torch.float64)          # a degenerate spectrum
Q = torch.linalg.qr(torch.randn(3, 3, dtype=torch.complex128))[0]  # a unitary
A0 = Q @ torch.diag(lam).to(torch.complex128) @ Q.mH
E = torch.randn(3, 3, dtype=torch.complex128)
E = (E + E.mH) / 2                                                 # a self-adjoint direction
L = functional.divided_differences(lam, torch.exp(lam), torch.exp(lam))
D = Q @ (L * (Q.mH @ E @ Q)) @ Q.mH                                # Theorem 1.11 for f = exp
s = 1e-6
fd = (torch.linalg.matrix_exp(A0 + s * E) - torch.linalg.matrix_exp(A0 - s * E)) / (2 * s)
print(f"{(D - fd).abs().max().item():.1e}")
assert (D - fd).abs().max() < 1e-8
```

Para valores propios próximos $`a,b`$ a distancia $`h`$, el cociente $`\hat q`$ pierde precisión por cancelación, y el
valor del punto medio $`\hat m=(f'(a)+f'(b))/2`$ yerra a lo sumo en $`T=h^2\max|f'''|/12`$ más el redondeo
([proposición 1.12](../../lecture_notes/es/01_matrix_algebras.md#r-1-12)). `divided_differences` toma $`\hat m`$ si y solo si
$`|\hat q-\hat m|\le2R`$, con la cota de redondeo $`R`$ de la [proposición 1.13](../../lecture_notes/es/01_matrix_algebras.md#r-1-13)
y el error de evaluación $`\varepsilon(|\hat f_a|+|a\hat f'_a|)`$; su error es entonces a lo sumo $`\min(T,3R)`$. Lo
comprobamos para $`f(x)=e^{x-c}`$ cerca de $`c=10^4`$, donde $`f^{[1]}(a,b)=e^{a-c}(e^h-1)/h`$ se calcula con precisión con
`expm1`. La última distancia, $`h=0.06`$, está por debajo del umbral $`\varepsilon^{1/3}\max(|a|,|b|)`$ de versiones
anteriores, que allí tomaban el valor del punto medio:

```python
eps = torch.finfo(torch.float64).eps
c = 1e4
f = lambda x: torch.exp(x - c)
for h in (1e-12, 1e-9, 1e-6, 1e-4, 1e-2, 0.06):
    w = torch.tensor([c, c + h], dtype=torch.float64)
    h = (w[1] - w[0]).item()                                       # the distance actually stored
    fw = f(w)                                                      # f' = f
    exact = math.expm1(h) / h                                      # f[a, b], as e^(a - c) = 1
    q = ((fw[1] - fw[0]) / h).item()                               # quotient
    m = (fw.sum() / 2).item()                                      # midpoint value
    y = functional.divided_differences(w, fw, fw)[0, 1].item()     # the library's choice
    R = eps * (fw.abs() + (w * fw).abs()).sum().item() / h + 4 * eps * abs(q)
    T = h ** 2 * math.exp(h) / 12 + 2 * eps                        # Proposition 1.12, plus rounding
    print(f"h = {h:.0e}: quotient {abs(q - exact):.0e}, midpoint {abs(m - exact):.0e}, library {abs(y - exact):.0e}")
    assert abs(y - exact) <= min(T, 3 * R) + 2 * eps              # Proposition 1.13(b)
```

*Observación numérica:* el cociente yerra hasta en $`5\cdot10^{-10}`$ para $`h\le10^{-6}`$ y en $`10^{-15}`$ en
$`h=0.06`$; el valor del punto medio es exacto salvo redondeo para $`h\le10^{-9}`$ y yerra en $`3\cdot10^{-4}`$ en $`h=0.06`$. La
biblioteca tiene el menor de los dos errores en todas las filas salvo $`h=10^{-4}`$, donde toma el valor del punto
medio ($`8\cdot10^{-10}`$) aunque el cociente es mejor ($`10^{-12}`$): $`R`$ admite un argumento perturbado en
$`\varepsilon|a|\approx2\cdot10^{-12}`$, lo cual es pesimista para `torch.exp` pero necesario para funciones como
$`e^{-i\tau x}`$, que redondean $`\tau x`$. La garantía $`\min(T,3R)`$ se cumple en todas las filas.

<a id="s-1-6"></a>

## 1.6 Normas y la dualidad de la norma traza

*Apuntes:* [sección 1.4](../../lecture_notes/es/01_matrix_algebras.md#s-1-4): [teorema 1.8](../../lecture_notes/es/01_matrix_algebras.md#r-1-8), las tres normas, [teorema 1.9](../../lecture_notes/es/01_matrix_algebras.md#r-1-9).

Una matriz general (no normal) se introduce mediante `alg.operator(tensor)`. Comparamos las
tres normas con los valores singulares calculados directamente por PyTorch, y comprobamos la
descomposición polar (teorema 1.8) y la dualidad (teorema 1.9):

```python
G = torch.randn(1, 1, 3, 3, dtype=torch.complex128)    # Ginibre matrix: no symmetry at all
X = M3.operator(G)
s = torch.linalg.svdvals(G[0, 0])                       # singular values, descending

assert abs(X.operator_norm().item() - s[0].item()) < 1e-12
assert abs(X.frobenius_norm().item() - s.norm().item()) < 1e-12
assert abs(X.trace_norm().item() - s.sum().item()) < 1e-12

absX = X.abs()
W = X @ absX.inverse()                                  # |X| is invertible almost surely, see below
assert close(W.adjoint() @ W, M3.identity())            # W is unitary
assert abs((X @ W.adjoint()).trace.item() - s.sum().item()) < 1e-10    # Tr(X W*) = ||X||_1
assert abs(X.trace.item()) <= X.trace_norm().item()     # |Tr X| <= ||X||_1
```

`torch.randn` con un dtype complejo genera entradas independientes con partes real e imaginaria
independientes $`\sim N(0,\tfrac12)`$, de modo que `G` es una matriz de Ginibre, que tiene densidad
(proposición A.16(b)). $`|X|`$ es invertible casi seguramente:
$`(\det|X|)^2=\det(|X|^2)=\det(X^*X)=|\det X|^2`$, y $`|\det X|^2`$ es un polinomio en las partes real e
imaginaria de las entradas que vale $`1`$ en $`X=1`$, así que se anula con probabilidad cero
(corolario A.7(b); este es también el primer paso del teorema 5.3).

<a id="s-1-7"></a>

## 1.7 Lotes: muchos operadores a la vez

*Apuntes:* [sección 1.7](../../lecture_notes/es/01_matrix_algebras.md#s-1-7), [proposición 1.10](../../lecture_notes/es/01_matrix_algebras.md#r-1-10).

El eje de lote es la razón de ser de la biblioteca: la teoría de matrices aleatorias, los
métodos de Monte Carlo y la optimización con múltiples arranques necesitan miles de
operadores, y una única llamada por lotes sustituye a un bucle de Python sobre ellos. Como
ejemplo, comprobamos un momento exacto del **conjunto gaussiano unitario** (GUE).

Según los [apuntes de clase, proposición 1.10](../../lecture_notes/es/01_matrix_algebras.md#r-1-10), si $`G`$ tiene entradas independientes con partes
real e imaginaria independientes $`\sim N(0,\tfrac12)`$ y $`X=(G+G^*)/\sqrt{2n}`$, entonces $`\mathbb E\,\tfrac1n\mathrm{Tr}\,X^4=2+\tfrac1{n^2}`$.

```python
n, B = 4, 200_000
G = torch.complex(torch.randn(B, 1, n, n), torch.randn(B, 1, n, n)).to(torch.complex128) / math.sqrt(2)
M4 = TypeIAlgebra([n], [n], precision='double')
X = M4.operator((G + G.conj().transpose(-2, -1)) / math.sqrt(2 * n), is_self_adjoint=True)

X2 = X @ X
m4 = (X2 @ X2).trace.real / n                           # shape (B,): one value per sample
mean, err = m4.mean().item(), m4.std().item() / math.sqrt(B)
print(f"E tr X^4 / n = {mean:.4f} +- {err:.4f}   (exact {2 + 1 / n**2:.4f})")
assert abs(mean - (2 + 1 / n**2)) < 5 * err
```

Una comprobación de Monte Carlo es un contraste estadístico: la estimación debe coincidir
con el valor exacto salvo unos pocos errores típicos. La variable $`\tfrac1n\mathrm{Tr}\,X^4`$ es un
polinomio en variables gaussianas, así que tiene varianza finita $`\sigma^2`$ (lema A.23: todos los
momentos gaussianos son finitos; los productos se tratan con $`|y_1\cdots y_r|\le\sum_i|y_i|^r`$), y por
la desigualdad de Chebyshev la media de $`B`$ muestras independientes se aparta del valor exacto en
más de $`5\sigma/\sqrt B`$ con probabilidad a lo sumo $`1/25`$. El código usa la desviación típica
muestral en lugar de $`\sigma`$, y la semilla fija hace que la comprobación sea reproducible. Este es
el esquema de toda comprobación frente a resultados demostrados en este tutorial: un enunciado
exacto, una estimación numérica y una tolerancia explícita.

<a id="s-1-8"></a>

## 1.8 Con qué comparar

Para matrices individuales, `numpy.linalg` y `scipy.linalg` (`sqrtm`, `expm`, `logm`,
`polar`) son las implementaciones de referencia, y las comprobaciones anteriores pueden
repetirse con ellas. La biblioteca añade tres cosas que importarán más adelante: el
procesamiento por lotes (sección 1.7), operadores dotados de una estructura de bloques
(capítulo 2) y la diferenciación automática (el autograd de PyTorch) a través de sus
operaciones.

<a id="x-1"></a>

## Ejercicios

[Ejercicios 1 y 3](../../lecture_notes/es/01_matrix_algebras.md#x-1) de los apuntes de clase, numéricamente:
- Compruebe que $`(fg)(A)=f(A)g(A)`$ para $`f=\sin`$, $`g=\exp`$ y una $`A`$ autoadjunta aleatoria (ejercicio 1).
- La **media geométrica** de $`A,B>0`$ es $`A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}`$, la única solución
  positiva $`X`$ de $`XA^{-1}X=B`$ (ejercicio 3). Compruebe numéricamente esta ecuación y que $`A\#B=B\#A`$.

<details>
<summary>Soluciones</summary>

```python
Hs = M3.operator_from_eigenvalues(lambda d: torch.randn(1, d), force_self_adjoint=True)
assert close(Hs.apply_function(lambda w: torch.sin(w) * torch.exp(w)),
             Hs.apply_function(torch.sin) @ Hs.apply_function(torch.exp))
```

```python
def geometric_mean(A, B):
    a, ai = A.sqrt(), A.power(-0.5)
    return a @ (ai @ B @ ai).sqrt() @ a

P = M3.operator_from_eigenvalues(lambda d: 0.5 + torch.rand(1, d), force_positive=True)
Q = M3.operator_from_eigenvalues(lambda d: 0.5 + torch.rand(1, d), force_positive=True)
GM = geometric_mean(P, Q)
assert close(GM, geometric_mean(Q, P), 1e-9)
assert close(GM @ P.inverse() @ GM, Q, 1e-9)
```

</details>

**Siguiente:** [Tutorial 2](02_von_neumann_algebras.md): conmutantes, el centro y el álgebra $`\bigoplus_cM_{k_c}(\mathbb C)`$ en la biblioteca.
