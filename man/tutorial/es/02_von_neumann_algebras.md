[English](../02_von_neumann_algebras.md) · **Español** · [Français](../fr/02_von_neumann_algebras.md) · [Deutsch](../de/02_von_neumann_algebras.md) · [中文](../zh/02_von_neumann_algebras.md) · [日本語](../ja/02_von_neumann_algebras.md) · [Русский](../ru/02_von_neumann_algebras.md)

# Tutorial 2. Álgebras de von Neumann de dimensión finita

Este tutorial acompaña al [capítulo 2](../../lecture_notes/es/02_von_neumann_algebras.md) de los apuntes de clase, que demuestra
que toda álgebra de von Neumann sobre $\mathbb C^n$ es, salvo un cambio de base unitario,
$\bigoplus_cM_{k_c}(\mathbb C)\otimes1_{m_c}$ ([teorema 2.6](../../lecture_notes/es/02_von_neumann_algebras.md#r-2-6)).
Calculamos conmutantes numéricamente en el ejemplo de un grupo de simetría y vemos cómo la biblioteca
almacena un álgebra de este tipo, su centro y la esperanza condicional sobre el centro.

<a id="s-2-4"></a>

## 2.4 Ejemplo: un grupo de simetría y su conmutante

*Apuntes:* [sección 2.4](../../lecture_notes/es/02_von_neumann_algebras.md#s-2-4); [teorema 2.2](../../lecture_notes/es/02_von_neumann_algebras.md#r-2-2), [teorema 2.6](../../lecture_notes/es/02_von_neumann_algebras.md#r-2-6), [proposición 2.7](../../lecture_notes/es/02_von_neumann_algebras.md#r-2-7).

La [sección 2.4](../../lecture_notes/es/02_von_neumann_algebras.md#s-2-4) de los apuntes de clase descompone la representación $\pi$ de
$S_3$ por permutación de los tres factores de $(\mathbb C^2)^{\otimes3}$: $\pi(S_3)'\cong M_4\oplus M_2$ tiene dimensión
$20$ y $\pi(S_3)''\cong M_1\oplus M_2$ tiene dimensión $5$. El código de abajo confirma numéricamente las dos dimensiones,
y que el conmutante de $\pi(S_3)'$ es de nuevo $\pi(S_3)''$
([teorema 2.2](../../lecture_notes/es/02_von_neumann_algebras.md#r-2-2)): el conmutante de un conjunto $S$ es el
núcleo de la aplicación lineal $X\mapsto(AX-XA)_{A\in S}$, cuya dimensión se calcula aquí en coma
flotante, contando los valores singulares por encima de una tolerancia.

```python
import itertools
import torch

torch.manual_seed(0)
dt = torch.complex128

def permutation_operator(perm, n_qubits=3):
    """Unitary on (C^2)^{(x)3} that permutes the tensor factors."""
    dim = 2 ** n_qubits
    P = torch.zeros(dim, dim, dtype=dt)
    for bits in itertools.product([0, 1], repeat=n_qubits):
        src = int(''.join(map(str, bits)), 2)
        dst = int(''.join(str(bits[perm[i]]) for i in range(n_qubits)), 2)
        P[dst, src] = 1
    return P

group = [permutation_operator(p) for p in itertools.permutations(range(3))]

def commutant_basis(ops, tol=1e-10):
    """Orthonormal basis (Hilbert-Schmidt) of {X : X A = A X for all A in ops}."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    # vec(AX - XA) = (A (x) 1 - 1 (x) A^T) vec(X) for row-major vec
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    _, s, Vh = torch.linalg.svd(L)
    null = Vh[(s > tol).sum():].conj().resolve_conj()   # rows span the null space
    return [v.reshape(d, d) for v in null]

def span_dim(ops, tol=1e-10):
    return int((torch.linalg.svdvals(torch.stack([A.reshape(-1) for A in ops])) > tol).sum())

comm = commutant_basis(group)
print(len(comm), span_dim(group), len(commutant_basis(comm)))
assert len(comm) == 20                                   # Step 1: pi(S3)'  = M_4 (+) M_2
assert span_dim(group) == 5                              # Step 2: pi(S3)'' = M_1 (+) M_2
assert len(commutant_basis(comm)) == 5                   # double commutant theorem
```

Las multiplicidades se observan en el espectro ([apuntes de clase, sección 2.4](../../lecture_notes/es/02_von_neumann_algebras.md#s-2-4), a partir del
[teorema 2.6](../../lecture_notes/es/02_von_neumann_algebras.md#r-2-6)):
todo $K\in\pi(S_3)'$ autoadjunto tiene valores propios de multiplicidad par sobre el complemento del subespacio
simétrico $F$, y para $K=E(Z+Z^*)$, con $Z$ una matriz de Ginibre y $E$ el promedio sobre el grupo, con probabilidad uno
$K$ tiene $4$ valores propios simples sobre $F$ y $2$ doblemente degenerados sobre su complemento. El código comprueba los
pares degenerados (que se dan para todo $K$) y, para el $K$ muestreado, que los dos pares y los cuatro
valores propios sobre $F$ son distintos, usando la proyección ortogonal $P_F=\frac1{|G|}\sum_g\pi(g)$ sobre $F$.

```python
def group_average(X):
    return sum(g @ X @ g.conj().T for g in group) / len(group)

H0 = torch.randn(8, 8, dtype=dt)                        # Ginibre: Re, Im independent N(0, 1/2)
Hsym = group_average(H0 + H0.conj().T)                  # K = E(Z + Z^*), a random symmetric Hamiltonian
assert all(torch.allclose(g @ Hsym, Hsym @ g) for g in group)

P_sym = sum(group) / len(group)                         # P_F, the projection onto the symmetric subspace
w, Q = torch.linalg.eigh(P_sym)
Q_rest = Q[:, w.real < 0.5]                             # orthonormal basis of the complement (dim 4)
e_rest = torch.linalg.eigvalsh(Q_rest.conj().T @ Hsym @ Q_rest)
print(e_rest)
assert Q_rest.shape[1] == 4
assert torch.allclose(e_rest[0], e_rest[1]) and torch.allclose(e_rest[2], e_rest[3])   # even multiplicities
e_sym = torch.linalg.eigvalsh(Q[:, w.real > 0.5].conj().T @ Hsym @ Q[:, w.real > 0.5])   # eigenvalues on F
assert e_sym.shape[0] == 4 and e_rest[2] - e_rest[1] > 1e-6                           # this sample: generic
assert (e_sym[1:] - e_sym[:-1]).min() > 1e-6                                           # 4 simple eigenvalues on F
assert (e_sym[:, None] - e_rest[None, :]).abs().min() > 1e-6                           # all six distinct
```

<a id="s-2-5"></a>

## 2.5 El álgebra en la biblioteca

*Apuntes:* [sección 2.3](../../lecture_notes/es/02_von_neumann_algebras.md#s-2-3), [teorema 2.6](../../lecture_notes/es/02_von_neumann_algebras.md#r-2-6); el centro: [sección 2.2](../../lecture_notes/es/02_von_neumann_algebras.md#s-2-2), [lema 2.3](../../lecture_notes/es/02_von_neumann_algebras.md#r-2-3).

La biblioteca almacena el álgebra abstracta $\bigoplus_cM_{k_c}(\mathbb C)$ del teorema 2.6
sobre $\bigoplus_c\mathbb C^{k_c}$, sin las multiplicidades: un operador es la lista de sus
bloques $A_c$, rellenados hasta un tamaño común. El conmutante del ejemplo es
`TypeIAlgebra([4, 2], [4, 2])`.

```python
from torch_vn_algebra import TypeIAlgebra, channels

alg = TypeIAlgebra([4, 2], [4, 2], precision='double')    # M_4 (+) M_2
print(alg.k_factors, alg.C)

A = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
print(A.matrix.shape)              # (1, 2, 4, 4): block 1 occupies the top-left 2 x 2 corner
assert A.matrix[0, 1, 2:, :].abs().max() == 0 and A.matrix[0, 1, :, 2:].abs().max() == 0
```

El centro está formado por los elementos $\sum_cz_cP_c$; `alg.central` los construye, y
conmutan con todo:

```python
def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

P0, P1 = alg.central([1.0, 0.0]), alg.central([0.0, 1.0])   # minimal central projections
Z = alg.central([2.0, -1.0])
assert P0.is_projection and P1.is_projection
assert close(P0 + P1, alg.identity())
assert close(Z @ A, A @ Z)
assert close(A @ P0 + A @ P1, A)                         # A = sum_c A P_c
```

<a id="s-2-6"></a>

## 2.6 Esperanzas condicionales

*Apuntes:* [sección 2.6](../../lecture_notes/es/02_von_neumann_algebras.md#s-2-6).

Según [apuntes de clase, proposición 2.8](../../lecture_notes/es/02_von_neumann_algebras.md#r-2-8), $E(A)=\sum_c\frac{\mathrm{Tr}\,A_c}{k_c}\,P_c$
es una esperanza condicional de $\bigoplus_cM_{k_c}$ sobre su centro: es positiva, es la identidad sobre el centro,
es un morfismo de bimódulos sobre el centro y conserva toda traza $\sum_cw_c\mathrm{Tr}\,A_c$.

En la biblioteca, esta aplicación es el canal `channels.center_expectation(alg)`. Al
aplicarla a operadores se obtiene su parte clásica (promediada en cada sector):

```python
E = channels.center_expectation(alg)
EA = E(A)
assert close(E(EA), EA)                                  # idempotent
assert close(E(Z @ A @ P1), Z @ EA @ P1)                 # bimodule property
A2 = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
assert close(EA @ A2, A2 @ EA)                           # E(A) commutes with a fresh element
for c, k in enumerate(alg.k_factors):                    # E(A)_c = (Tr A_c / k_c) 1, so Tr E(A)_c = Tr A_c
    trA = A.matrix[0, c].diagonal().sum()
    assert close(EA.matrix[0, c, :k, :k], trA / k * torch.eye(k, dtype=dt))
    assert abs(EA.matrix[0, c].diagonal().sum() - trA) < 1e-10
```

<a id="s-2-7"></a>

## 2.7 Con qué comparar

Los ejemplos de la sección 2.4 pueden comprobarse con cualquier paquete de álgebra lineal;
el software de teoría de representaciones (GAP, SageMath) da la descomposición
$\pi\cong\bigoplus\rho_\lambda\otimes1_{m_\lambda}$ de forma simbólica. Los códigos de
diagonalización exacta en física (QuSpin, los módulos de simetría de QuTiP) usan la misma
estructura de bloques para diagonalizar hamiltonianos sector a sector. La biblioteca toma
la estructura como punto de partida: una vez conocidos los bloques, todos los cálculos
posteriores —funciones, trazas, estados, canales, dinámica— permanecen dentro del álgebra.

<a id="x-2"></a>

## Ejercicios

[Ejercicio 2](../../lecture_notes/es/02_von_neumann_algebras.md#x-2) de los apuntes de clase: para el desplazamiento cíclico $S$ en $\mathbb C^3$,
$Se_i=e_{i+1\bmod3}$, tanto $\{S\}'$ como $\{S\}''$ son iguales a $\mathrm{span}\{1,S,S^2\}$, las matrices
circulantes. Compruebe numéricamente las dimensiones.

<details>
<summary>Soluciones</summary>

```python
S = torch.roll(torch.eye(3, dtype=dt), 1, dims=0)
assert len(commutant_basis([S])) == 3
assert len(commutant_basis(commutant_basis([S]))) == 3
assert span_dim([torch.eye(3, dtype=dt), S, S @ S]) == 3
```

</details>

**Siguiente:** [Tutorial 3. Trazas, estados y entropía](03_traces_and_states.md)
