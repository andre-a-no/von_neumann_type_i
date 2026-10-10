[English](../01_matrix_algebras.md) · **Español** · [Français](../fr/01_matrix_algebras.md) · [Deutsch](../de/01_matrix_algebras.md) · [中文](../zh/01_matrix_algebras.md) · [日本語](../ja/01_matrix_algebras.md) · [Русский](../ru/01_matrix_algebras.md)

# 1. Álgebras de matrices y cálculo funcional

Este capítulo fija la notación y demuestra los hechos sobre $M_n(\mathbb C)$ que se usan en
todos los capítulos posteriores: el teorema espectral, el cálculo funcional, la
descomposición polar y las normas definidas por los valores singulares. Cada resultado va
seguido del código correspondiente de la biblioteca y de una comprobación numérica.

**Aprenderá** a construir operadores con un espectro prescrito, a aplicar una función a un
operador, cómo se almacenan los lotes de operadores y cómo contrastar un cálculo con un
resultado exacto.

**Convenios.** Los productos escalares $\langle x,y\rangle=\sum_i\bar x_iy_i$ son lineales en
el segundo argumento. $A^*$ es la traspuesta conjugada, de modo que $\langle x,Ay\rangle=\langle A^*x,y\rangle$.
$e_1,\dots,e_n$ es la base canónica, $1$ la matriz identidad y $\operatorname{Tr}$ la traza
usual.

## 1.1 El álgebra $M_n(\mathbb C)$

$M_n(\mathbb C)$ es el álgebra de las matrices complejas $n\times n$ con la involución
$A\mapsto A^*$ y la norma de operador $\|A\|=\sup_{\|x\|=1}\|Ax\|$. Un elemento $A$ es

| nombre | definición |
|---|---|
| normal | $A^*A=AA^*$ |
| autoadjunto | $A=A^*$ |
| positivo, $A\ge0$ | $\langle x,Ax\rangle\ge0$ para todo $x$ |
| unitario | $A^*A=AA^*=1$ |
| proyección | $A=A^*=A^2$ |

Escribimos $A\le B$ si $B-A\ge0$, y $A>0$ si $A\ge0$ y $A$ es invertible.

**Lema 1.1 (identidad C\*).** $\|A^*\|=\|A\|$ y $\|A^*A\|=\|A\|^2$.

*Demostración.* $\|Ax\|=\sup_{\|y\|=1}|\langle y,Ax\rangle|=\sup_{\|y\|=1}|\langle A^*y,x\rangle|$,
luego $\|A\|=\sup_{\|x\|=\|y\|=1}|\langle A^*y,x\rangle|=\|A^*\|$. Además,
$\|Ax\|^2=\langle x,A^*Ax\rangle\le\|A^*A\|$ para $\|x\|=1$, de donde $\|A\|^2\le\|A^*A\|\le\|A^*\|\,\|A\|=\|A\|^2$. $\square$

> **Significado físico.** En mecánica cuántica, un sistema con $n$ configuraciones
> distinguibles se describe mediante $\mathbb C^n$. Las matrices autoadjuntas son los
> *observables* (energía, espín, posición en una red); sus valores propios son los
> resultados posibles de una medición. Las proyecciones son preguntas de respuesta sí/no.
> Las matrices unitarias son las simetrías y las evoluciones temporales, $e^{-itH}$ para un
> hamiltoniano $H$. Las matrices positivas de traza uno son los *estados* (matrices
> densidad, capítulo 3).

## 1.2 El teorema espectral

**Lema 1.2 (Schur).** Para toda $A\in M_n(\mathbb C)$ existe una unitaria $U$ tal que
$U^*AU$ es triangular superior.

*Demostración.* Inducción sobre $n$; el caso $n=1$ es trivial. El polinomio característico
de $A$ tiene una raíz $\lambda$ (teorema fundamental del álgebra), luego $Av=\lambda v$ para
algún vector unitario $v$. Completamos $v$ hasta una base ortonormal $v,w_2,\dots,w_n$ y
sea $V$ la unitaria cuyas columnas son estos vectores. Entonces
$V^*AV=\begin{pmatrix}\lambda&*\\0&A'\end{pmatrix}$ con $A'\in M_{n-1}(\mathbb C)$.
Por inducción, $U'^*A'U'$ es triangular superior para cierta unitaria $U'$, y
$U=V\,(1\oplus U')$ cumple lo pedido. $\square$

**Teorema 1.3 (teorema espectral).** $A$ es normal si y solo si
$A=U\,\mathrm{diag}(\lambda_1,\dots,\lambda_n)\,U^*$ para alguna unitaria $U$. Los
$\lambda_i$ son los valores propios de $A$ contados con multiplicidad.

*Demostración.* Si $A=U\Lambda U^*$ con $\Lambda$ diagonal, entonces
$A^*A=U\bar\Lambda\Lambda U^*=U\Lambda\bar\Lambda U^*=AA^*$. Recíprocamente, sea $A$
normal y $T=U^*AU$ triangular superior (lema 1.2); $T$ es normal. Comparando las
entradas $(1,1)$: $(T^*T)_{11}=|t_{11}|^2$ y $(TT^*)_{11}=\sum_j|t_{1j}|^2$, luego
$t_{1j}=0$ para $j>1$. Así $T=t_{11}\oplus T'$ con $T'$ normal y triangular superior, y
por inducción $T$ es diagonal. La diagonal de $U^*AU$ está formada por las raíces de
$\det(A-\lambda)=\det(T-\lambda)$. $\square$

**Corolario 1.4.** Sea $A$ normal con valores propios $\lambda_i$. Entonces

1. $\|A\|=\max_i|\lambda_i|$;
2. $A$ es autoadjunta si y solo si todos los $\lambda_i\in\mathbb R$; unitaria si y solo si
   todos los $|\lambda_i|=1$; una proyección si y solo si todos los $\lambda_i\in\{0,1\}$;
3. $A\ge0$ si y solo si $A$ es autoadjunta y todos los $\lambda_i\ge0$.

*Demostración.* Con $y=U^*x$, se tiene $\|Ax\|^2=\sum_i|\lambda_i|^2|y_i|^2$ y $\|y\|=\|x\|$,
lo que da (1), con igualdad para $x=Ue_i$ con $|\lambda_i|$ máximo. En (2), $A^*=U\bar\Lambda U^*$,
$A^*A=U|\Lambda|^2U^*$ y $A^2=U\Lambda^2U^*$, y las condiciones se cumplen si y solo si se
cumplen entrada a entrada para $\Lambda$. Para (3), en primer lugar: toda $A$ positiva es
autoadjunta, ya que
$\langle x,(A-A^*)x\rangle=\langle x,Ax\rangle-\overline{\langle x,Ax\rangle}=0$ para todo $x$,
y una matriz $D$ con $\langle x,Dx\rangle=0$ para todo $x$ es nula (polarización:
$4\langle x,Dy\rangle=\sum_{k=0}^3 i^{-k}\langle x+i^ky,D(x+i^ky)\rangle$). Entonces
$\langle x,Ax\rangle=\sum_i\lambda_i|y_i|^2$, que es $\ge0$ para todo $x$ si y solo si todos los $\lambda_i\ge0$. $\square$

## 1.3 El cálculo funcional

Sea $A=U\Lambda U^*$ normal y $f$ una función compleja definida en $\operatorname{spec}A$. Se define
$$
f(A)=U\,\mathrm{diag}\big(f(\lambda_1),\dots,f(\lambda_n)\big)\,U^*.
$$

**Teorema 1.5.** (a) $f(A)$ no depende de la elección de $U$. (b) Para un polinomio $p$,
$p(A)$ es el polinomio usual en $A$. (c) $f\mapsto f(A)$ es lineal y multiplicativa,
$\bar f(A)=f(A)^*$, y $\operatorname{spec}f(A)=f(\operatorname{spec}A)$.
(d) Si $B$ conmuta con $A$, entonces $B$ conmuta con $f(A)$. (e) Si $g$ está definida en
$f(\operatorname{spec}A)$, entonces $g(f(A))=(g\circ f)(A)$.

*Demostración.* (b) $A^k=U\Lambda^kU^*$, luego $p(A)=Up(\Lambda)U^*$ por linealidad. (a) El
espectro es finito, así que la interpolación de Lagrange proporciona un polinomio $p$ con
$p=f$ en $\operatorname{spec}A$. Entonces $U\,f(\Lambda)\,U^*=U\,p(\Lambda)\,U^*=p(A)$ para
*toda* unitaria $U$ que diagonalice $A$, y $p(A)$ está definido sin recurrir a $U$. (c) se
cumple para las matrices diagonales $f(\Lambda)$ y se conserva al conjugar por $U$; el
espectro de $Uf(\Lambda)U^*$ es la diagonal de $f(\Lambda)$. (d) $B$ conmuta con todo
polinomio en $A$, en particular con $p(A)=f(A)$. (e) $f(A)=Uf(\Lambda)U^*$ es normal y la
diagonaliza la misma $U$, luego por (a) $g(f(A))=U\,g(f(\Lambda))\,U^*=(g\circ f)(A)$. $\square$

**Teorema 1.6 (raíz cuadrada).** Para $A\ge0$ existe exactamente una $B\ge0$ con
$B^2=A$, a saber, $B=\sqrt A$ definida mediante el cálculo funcional.

*Demostración.* $\sqrt A\ge0$ y $(\sqrt A)^2=A$ por el teorema 1.5(c) y el corolario 1.4(3).
Sea $B\ge0$ con $B^2=A$, $B=V\,\mathrm{diag}(\mu_j)\,V^*$, $\mu_j\ge0$. Entonces
$A=V\,\mathrm{diag}(\mu_j^2)\,V^*$, luego $\operatorname{spec}A=\{\mu_j^2\}$. Sea $p$ un
polinomio con $p(t)=\sqrt t$ en $\operatorname{spec}A$. Entonces
$p(A)=\sqrt A$ y $p(A)=p(B^2)=V\,\mathrm{diag}(p(\mu_j^2))\,V^*=V\,\mathrm{diag}(\mu_j)\,V^*=B$. $\square$

**Corolario 1.7.** $A\ge0$ si y solo si $A=C^*C$ para alguna $C$. (Si $A\ge0$, tómese
$C=\sqrt A$; recíprocamente, $\langle x,C^*Cx\rangle=\|Cx\|^2\ge0$.)

## 1.4 Descomposición polar y valores singulares

Para cualquier $A$, la matriz $A^*A$ es positiva (corolario 1.7); su raíz cuadrada
$|A|=(A^*A)^{1/2}$ es el **valor absoluto** de $A$.

**Teorema 1.8 (descomposición polar).** Toda $A\in M_n(\mathbb C)$ puede escribirse como
$A=W|A|$ con $W$ unitaria.

*Demostración.* Para todo $x$, $\||A|x\|^2=\langle x,|A|^2x\rangle=\langle x,A^*Ax\rangle=\|Ax\|^2$.
Por tanto $W_0:|A|x\mapsto Ax$ es una isometría lineal bien definida de la imagen
$R_1$ de $|A|$ sobre la imagen $R_2$ de $A$: si $|A|x=|A|x'$, entonces
$\|A(x-x')\|=\||A|(x-x')\|=0$. En particular $\dim R_1=\dim R_2$, de modo que los
complementos ortogonales también tienen la misma dimensión; sea $W_1:R_1^\perp\to R_2^\perp$
una isometría cualquiera entre ellos (que envíe una base ortonormal en una base
ortonormal). Entonces $W=W_0\oplus W_1$ es unitaria y $W|A|x=W_0|A|x=Ax$. $\square$

Los valores propios $s_1\ge\dots\ge s_n\ge0$ de $|A|$ son los **valores singulares** de
$A$. A lo largo del curso se usan tres normas:
$$
\|A\|=s_1,\qquad \|A\|_2=\Big(\sum_is_i^2\Big)^{1/2}=(\operatorname{Tr}A^*A)^{1/2},
\qquad\|A\|_1=\sum_is_i=\operatorname{Tr}|A|.
$$
La primera identidad se cumple porque $\|Ax\|=\||A|x\|$ para todo $x$, luego
$\|A\|=\||A|\|=s_1$ por el corolario 1.4(1); la segunda, porque
$\operatorname{Tr}A^*A=\operatorname{Tr}|A|^2=\sum_is_i^2$.

**Teorema 1.9 (dualidad de la norma traza).**
$\|A\|_1=\max_{U\text{ unitary}}|\operatorname{Tr}(AU)|$, y el máximo se alcanza en
$U=W^*$, donde $A=W|A|$. En particular $|\operatorname{Tr}A|\le\|A\|_1$.

*Demostración.* $\operatorname{Tr}(AW^*)=\operatorname{Tr}(W|A|W^*)=\operatorname{Tr}|A|$. Para cualquier unitaria $U$,
por la ciclicidad de la traza y la desigualdad de Cauchy–Schwarz para $\langle X,Y\rangle_{HS}=\operatorname{Tr}(X^*Y)$,
$$
|\operatorname{Tr}(AU)|=|\operatorname{Tr}(|A|^{1/2}\,UW\,|A|^{1/2})|
=|\langle|A|^{1/2},\,UW|A|^{1/2}\rangle_{HS}|
\le\||A|^{1/2}\|_2\,\|UW|A|^{1/2}\|_2=\operatorname{Tr}|A|,
$$
ya que $\|UW X\|_2=\|X\|_2$ por ser $UW$ unitaria, y $\||A|^{1/2}\|_2^2=\operatorname{Tr}|A|$. $\square$

En particular, las tres cantidades anteriores son normas. $\|\cdot\|$ es una norma de
operador; $\|\cdot\|_2$ es la norma asociada al producto escalar $\langle X,Y\rangle_{HS}$; y
para $\|\cdot\|_1$, el teorema 1.9 da
$\|A+B\|_1=\max_U|\operatorname{Tr}(AU)+\operatorname{Tr}(BU)|\le\|A\|_1+\|B\|_1$, mientras que
$\|A\|_1=0$ obliga a que $|A|=0$ y $A=W|A|=0$. Para $P\ge0$, $|P|=(P^2)^{1/2}=P$ por el
teorema 1.6, luego $\|P\|_1=\operatorname{Tr}P$.

## 1.5 Operadores en la biblioteca

La biblioteca trabaja con álgebras $\bigoplus_cM_{k_c}(\mathbb C)$ (capítulo 2). Un único
bloque corresponde al caso $C=1$: `TypeIAlgebra([n], [n])` es $M_n(\mathbb C)$. (Las dos
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
construcción es exactamente el teorema 1.3 leído al revés: $A=U\Lambda U^*$ con $U$
distribuida según la medida de Haar. El indicador `force_positive=True` se verifica sobre
los valores propios (corolario 1.4(3)), de modo que la biblioteca sabe que `A` es positivo
sin tener que comprobarlo de nuevo.

El cálculo funcional para operadores autoadjuntos es `apply_function` (diagonaliza cada
bloque con `eigh`); las funciones habituales tienen sus propios métodos. Comprobamos los
teoremas 1.5 y 1.6:

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

## 1.6 Normas y la dualidad de la norma traza

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
W = X @ absX.inverse()                                  # |X| is invertible (almost surely, for Ginibre)
assert close(W.adjoint() @ W, M3.identity())            # W is unitary
assert abs((X @ W.adjoint()).trace.item() - s.sum().item()) < 1e-10    # Tr(X W*) = ||X||_1
assert abs(X.trace.item()) <= X.trace_norm().item()     # |Tr X| <= ||X||_1
```

## 1.7 Lotes: muchos operadores a la vez

El eje de lote es la razón de ser de la biblioteca: la teoría de matrices aleatorias, los
métodos de Monte Carlo y la optimización con múltiples arranques necesitan miles de
operadores, y una única llamada por lotes es mucho más rápida que un bucle de Python. Como
ejemplo, comprobamos un momento exacto del **conjunto gaussiano unitario** (GUE).

Sea $G$ una matriz de entradas independientes, con partes real e imaginaria independientes
$\sim N(0,\tfrac12)$, de modo que $\mathbb E|G_{ij}|^2=1$, y sea $H=(G+G^*)/\sqrt2$.

**Proposición 1.10.** $\mathbb E\,\operatorname{Tr}H^4=2n^3+n$. Por tanto, para
$X=H/\sqrt n$, $\ \mathbb E\,\tfrac1n\operatorname{Tr}X^4=2+\tfrac1{n^2}$.

*Demostración.* *Fórmula de Isserlis (de Wick).* Sean $Y_1,\dots,Y_4$ variables aleatorias
reales centradas conjuntamente gaussianas (se admiten repeticiones) con
$c_{ij}=\mathbb E[Y_iY_j]$. Su función generatriz de momentos conjunta es
$\mathbb E\,e^{\sum_it_iY_i}=e^{q(t)/2}$ con $q(t)=\sum_{ij}c_{ij}t_it_j$ (la función
generatriz de momentos gaussiana, un hecho estándar de probabilidad). Apliquemos
$\partial^4/\partial t_1\cdots\partial t_4$ en $t=0$: a la izquierda se obtiene
$\mathbb E[Y_1Y_2Y_3Y_4]$; a la derecha, solo el término $q^2/8$ de la serie exponencial
contiene $t_1t_2t_3t_4$, con coeficiente $c_{12}c_{34}+c_{13}c_{24}+c_{14}c_{23}$.
Por tanto $\mathbb E[Y_1Y_2Y_3Y_4]$ es la suma, sobre los tres emparejamientos, de los
productos de las esperanzas por pares. Por multilinealidad, lo mismo vale para variables
complejas que son combinaciones lineales complejas de variables reales conjuntamente
gaussianas, como las entradas de $H$.

*El momento.* Las esperanzas por pares son
$\mathbb E[H_{ij}H_{kl}]=\tfrac12\mathbb E[(G_{ij}+\bar G_{ji})(G_{kl}+\bar G_{lk})]=\delta_{il}\delta_{jk}$,
porque $\mathbb E[G_{ij}G_{kl}]=0$ y $\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}$.
En $\operatorname{Tr}H^4=\sum_{ijkl}H_{ij}H_{jk}H_{kl}H_{li}$ los emparejamientos dan:
$\{12\}\{34\}$: $\delta_{ik}\delta_{ki}$, cuya suma es $n^3$ (índices libres $i=k$, $j$, $l$);
$\{14\}\{23\}$: $\delta_{jl}\delta_{jl}$, cuya suma es $n^3$;
$\{13\}\{24\}$: $\delta_{il}\delta_{jk}\delta_{ji}\delta_{kl}$, que obliga a $i=j=k=l$ y cuya suma es $n$.
Dividiendo entre $n\cdot n^2$ se obtiene la segunda afirmación. $\square$

Los dos emparejamientos *planos* dan el número de Catalan $2$ de la ley del semicírculo; el
emparejamiento cruzado da la corrección $1/n^2$.

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
con el valor exacto salvo unos pocos errores típicos. Este es el esquema de toda
comprobación frente a resultados conocidos en este tutorial: un enunciado exacto, una
estimación numérica y una tolerancia explícita.

## 1.8 Con qué comparar

Para matrices individuales, `numpy.linalg` y `scipy.linalg` (`sqrtm`, `expm`, `logm`,
`polar`) son las implementaciones de referencia, y las comprobaciones anteriores pueden
repetirse con ellas. La biblioteca añade tres cosas que importarán más adelante: el
procesamiento por lotes (sección 1.7), operadores dotados de una estructura de bloques
(capítulo 2) y la diferenciación automática a través de todas las operaciones.

## Ejercicios

1. Compruebe que $(fg)(A)=f(A)g(A)$ para $f=\sin$, $g=\exp$ y una $A$ autoadjunta aleatoria.
2. Demuestre que $\|A\|_1=\|A^*\|_1$ y que $\|AB\|_1\le\|A\|\,\|B\|_1$, usando el teorema 1.9.
3. La **media geométrica** de $A,B>0$ es $A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$.
   Demuestre que $A\#B$ es la única solución positiva $X$ de $XA^{-1}X=B$, deduzca que
   $A\#B=B\#A$ y compruebe ambas cosas numéricamente.

<details>
<summary>Soluciones</summary>

1. Por el teorema 1.5(c).

```python
Hs = M3.operator_from_eigenvalues(lambda d: torch.randn(1, d), force_self_adjoint=True)
assert close(Hs.apply_function(lambda w: torch.sin(w) * torch.exp(w)),
             Hs.apply_function(torch.sin) @ Hs.apply_function(torch.exp))
```

2. $|\operatorname{Tr}(A^*U)|=|\overline{\operatorname{Tr}(U^*A)}|=|\operatorname{Tr}(AU^*)|$, y
   $U\mapsto U^*$ es una biyección del grupo unitario, luego los máximos del teorema 1.9
   coinciden. Para la segunda afirmación, sea $B=V|B|$. Para $U$ unitaria,
   $|\operatorname{Tr}(ABU)|=|\operatorname{Tr}(|B|^{1/2}UAV|B|^{1/2})|
   \le\||B|^{1/2}\|_2\,\|UAV|B|^{1/2}\|_2\le\|A\|\operatorname{Tr}|B|$, usando
   $\|YX\|_2\le\|Y\|\,\|X\|_2$ (aplíquese $\|Yx\|\le\|Y\|\,\|x\|$ a cada columna de $X$) y
   $\|UAV\|=\|A\|$ para unitarias $U,V$ (conservan la longitud de los vectores).

3. Escribamos $T=A^{-1/2}XA^{-1/2}$. Entonces $XA^{-1}X=B$ si y solo si $T^2=A^{-1/2}BA^{-1/2}$, y
   $X\ge0$ si y solo si $T\ge0$, ya que $\langle x,Xx\rangle=\langle A^{1/2}x,TA^{1/2}x\rangle$ y
   $A^{1/2}$ es invertible. Por el teorema 1.6, $T=(A^{-1/2}BA^{-1/2})^{1/2}$ es la única
   solución positiva, es decir, $X=A\#B$. Toda solución $X$ de $XA^{-1}X=B$ es invertible,
   porque $\det(X)^2\det(A)^{-1}=\det B\neq0$. Así, $XA^{-1}X=B$ equivale a
   $XB^{-1}X=A$ (de $A^{-1}=X^{-1}BX^{-1}$, invirtiendo ambos miembros), y la unicidad da
   $A\#B=B\#A$.

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

**Siguiente:** el [capítulo 2](02_von_neumann_algebras.md) sustituye $M_n$ por una suma
directa de álgebras de matrices y demuestra que toda álgebra de von Neumann de dimensión
finita es de esta forma.
