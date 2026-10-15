[English](../04_symmetries.md) · **Español** · [Français](../fr/04_symmetries.md) · [Deutsch](../de/04_symmetries.md) · [中文](../zh/04_symmetries.md) · [日本語](../ja/04_symmetries.md) · [Русский](../ru/04_symmetries.md)

# Tutorial 4. Simetrías y representaciones

Este tutorial acompaña al [capítulo 4](../../lecture_notes/es/04_symmetries.md) de los apuntes de clase, que desarrolla la teoría de
representaciones que hay detrás de las álgebras de simetrías, la convierte en un algoritmo aleatorizado que descompone el
álgebra de von Neumann generada por cualquier lista de matrices ([proposición 4.9](../../lecture_notes/es/04_symmetries.md#r-4-9)), la aplica a $`U(1)`$ y
$`SU(2)`$, clasifica los sectores de un álgebra generada por matrices reales en tipo real, complejo y cuaterniónico
([teorema 4.20](../../lecture_notes/es/04_symmetries.md#r-4-20)) y analiza la descomposición de generadores con ruido
([teorema 4.25](../../lecture_notes/es/04_symmetries.md#r-4-25)). Calculamos caracteres, ejecutamos la descomposición sobre grupos y sobre generadores,
construimos cadenas de espines sector por sector con `SpinChain` y `su2.SU2Chain`, calculamos la estructura real de
una descomposición y certificamos una descomposición de generadores con ruido.

<a id="s-4-2"></a>

## 4.2 La descomposición isotípica y el conmutante

*Apuntes:* [sección 4.2](../../lecture_notes/es/04_symmetries.md#s-4-2); [teorema 4.3](../../lecture_notes/es/04_symmetries.md#r-4-3), [teorema 4.4](../../lecture_notes/es/04_symmetries.md#r-4-4), [corolario 4.5](../../lecture_notes/es/04_symmetries.md#r-4-5).

La [sección 4.2](../../lecture_notes/es/04_symmetries.md#s-4-2) de los apuntes de clase calcula a mano los caracteres de tres representaciones
irreducibles de $`S_3`$ —la trivial, el signo y la representación estándar sobre $`(1,1,1)^\perp\subseteq\mathbb C^3`$— y de la
representación $`\pi`$ que permuta los tres factores de $`(\mathbb C^2)^{\otimes3}`$, $`\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}`$
(tabla de abajo), y obtiene del [corolario 4.5](../../lecture_notes/es/04_symmetries.md#r-4-5) las multiplicidades $`m_{\rm triv}=4`$, $`m_{\rm sgn}=0`$,
$`m_{\rm std}=2`$ y $`\dim\pi(S_3)'=20`$. El código recalcula la tabla: obtiene $`P_\sigma`$ como la restricción de $`\pi(\sigma)`$
al subespacio invariante generado por $`|100\rangle,|010\rangle,|001\rangle`$ (sobre el que $`\pi(\sigma)`$ permuta las
posiciones del único $`1`$, lema B.1(a)), y $`\#\text{cycles}(\sigma)=\dim\ker(P_\sigma-1)`$ (los vectores constantes en los ciclos).

| | identidad ($`1`$) | transposiciones ($`3`$) | $`3`$-ciclos ($`2`$) |
|---|---|---|---|
| $`\chi_{\rm triv}`$ | $`1`$ | $`1`$ | $`1`$ |
| $`\mathrm{sgn}`$ | $`1`$ | $`-1`$ | $`1`$ |
| $`\chi_{\rm std}`$ | $`2`$ | $`0`$ | $`-1`$ |
| $`\chi_\pi`$ | $`8`$ | $`4`$ | $`2`$ |

```python
import math
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    return (torch.as_tensor(X) - torch.as_tensor(Y)).abs().max().item() < tol

G = dc.permutation_representation(3, 2)                  # (6, 8, 8): pi(sigma) permutes three qubits
chi_pi = G.diagonal(dim1=-2, dim2=-1).sum(-1).real       # characters Tr pi(sigma)
P3 = G[:, [4, 2, 1]][:, :, [4, 2, 1]].real              # on |100>, |010>, |001>: S_3 permuting C^3
cycles = 3 - torch.linalg.matrix_rank(P3 - torch.eye(3, dtype=torch.float64)).double()   # dim ker(P - 1)
assert close(chi_pi, 2.0 ** cycles)
irreps = {'trivial': torch.ones(6, dtype=torch.float64), 'sign': (-1.0) ** (3 - cycles),
          'standard': P3.diagonal(dim1=-2, dim2=-1).sum(-1) - 1}                 # fixed points - 1
chis = torch.stack(list(irreps.values()))
assert close(chis @ chis.T / 6, torch.eye(3, dtype=torch.float64))    # Corollary 4.5(a): all irreducible
mult = chis @ chi_pi / 6                                               # Corollary 4.5(b)
print(dict(zip(irreps, mult.tolist())), (chi_pi ** 2).mean().item())
assert close(mult, torch.tensor([4.0, 0.0, 2.0], dtype=torch.float64))
assert dc.commutant_dimension(G) == 20                                # Corollary 4.5(c): M_4 (+) M_2
```

Es la descomposición de la sección 2.4, obtenida ahora a partir de una tabla de seis números por representación.

<a id="s-4-3"></a>

## 4.3 De los generadores al álgebra: un algoritmo aleatorizado

*Apuntes:* [sección 4.3](../../lecture_notes/es/04_symmetries.md#s-4-3); [lema 4.6](../../lecture_notes/es/04_symmetries.md#r-4-6), [proposición 4.7](../../lecture_notes/es/04_symmetries.md#r-4-7), [lema 4.8](../../lecture_notes/es/04_symmetries.md#r-4-8), [proposición 4.9](../../lecture_notes/es/04_symmetries.md#r-4-9).

El algoritmo aleatorizado de los apuntes de clase es `decompose.wedderburn`, con las interfaces `from_group` (todos los
elementos de un grupo finito; el conmutante se muestrea con el twirl `decompose.twirl`) y `from_generators` (cualquier
lista de matrices; el conmutante es un núcleo, `decompose.commutant_basis`). Divide $`\mathbb C^n`$ en los subespacios
propios de un $`K\in\mathcal A'_h`$ aleatorio, los agrupa y los alinea con los bloques de un $`Y\in\mathcal A'`$ aleatorio, y
devuelve una `Decomposition` con la unitaria `V`, los $`k_c`$ (`k_factors`) y los $`m_c`$ (`multiplicities`). En coma
flotante trata como iguales los valores propios más cercanos que una tolerancia y divide los subespacios resultantes
comprimiendo en ellos nuevos elementos aleatorios del conmutante ([lema 4.8(a)](../../lecture_notes/es/04_symmetries.md#r-4-8)).

`Decomposition.verify()` informa de las distancias relativas de los generadores (y de nuevos elementos del
conmutante) a $`\mathcal B`$ (resp. $`\mathcal B'`$), y compara $`\sum_cm_c^2`$ con $`\dim\mathcal A'`$ (a partir de los
caracteres para un grupo, corolario 4.5(c); a partir del núcleo para generadores). Lo ejecutamos sobre $`S_3`$ y tres
cúbits, y luego miramos dentro: los subespacios propios de un $`K`$ aleatorio y los operadores de entrelazamiento de
un $`Y`$ aleatorio.

```python
D = dc.from_group(G)                                     # the twirl samples pi(S_3)'
print(D)
assert D.k_factors == (1, 2) and D.multiplicities == (4, 2)            # Section 2.4: k = (1, 2), m = (4, 2)
assert D.verify()['ok']
rho = D.to_algebra(G).matrix                             # (6, 2, 2, 2): the blocks rho_c(sigma), padded
assert close(rho[:, 0, 0, 0], torch.ones(6))                           # sector 1: trivial representation
assert close(rho[:, 1].diagonal(dim1=-2, dim2=-1).sum(-1), chis[2])    # sector 2: the standard one
schur = torch.einsum('gij,gkl->ijkl', rho[:, 1], rho[:, 1].conj()) / 6
assert close(schur.reshape(4, 4), torch.eye(4) / 2)                    # Theorem 4.4 with k = 2
assert close(D.from_algebra(D.to_algebra(G)), G)

Z = torch.randn(8, 8, dtype=dt)
K = dc.twirl(G, Z + Z.mH)                                # a random element of pi(S_3)'_h
w, Q = torch.linalg.eigh(K)
cuts = [0] + [i + 1 for i in range(7) if w[i + 1] - w[i] > 1e-8] + [8]
spaces = [Q[:, a:b] for a, b in zip(cuts[:-1], cuts[1:])]
assert sorted(W.shape[1] for W in spaces) == [1, 1, 1, 1, 2, 2]      # Lemma 4.6: m_c pieces of size k_c
assert all(close(G @ W, W @ (W.mH @ G @ W)) for W in spaces)           # invariant subspaces
Y = dc.twirl(G, torch.randn(8, 8, dtype=dt))
Wa, Wb = [W for W in spaces if W.shape[1] == 2]
s = torch.linalg.svdvals(Wa.mH @ Y @ Wb)                               # Lemma 4.8(b): a multiple of a unitary
assert s[1] > 1e-3 and abs(s[0] - s[1]) < 1e-10
assert all(close(W.mH @ Y @ Wa, 0) for W in spaces if W.shape[1] == 1)   # different summands: zero
```

El conmutante domina el coste: una muestra del twirl cuesta $`|G|`$ productos de matrices $`n\times n`$. Para generadores,
la biblioteca diagonaliza primero un elemento hermítico $`R=U\mathrm{diag}(\lambda)U^*`$ de $`\mathcal A`$ (una
combinación real aleatoria de las partes hermíticas de los generadores más un múltiplo de un producto simetrizado de
dos de tales combinaciones). Todo $`X\in\mathcal A'`$ conmuta con $`R`$, y por tanto con sus proyecciones espectrales
(teorema 1.5(d)), así que $`U^*XU`$ es diagonal por bloques respecto de los subespacios propios de $`R`$, de dimensiones
$`d_\lambda`$; esto se impone de antemano, y el problema del núcleo tiene $`N_R=\sum_\lambda d_\lambda^2`$ incógnitas en
lugar de $`n^2`$. Agrupar varios valores propios solo agranda los bloques, así que el problema reducido tiene el mismo
conjunto de soluciones $`\mathcal A'`$. En la descomposición, $`V^*RV=\bigoplus_cr_c\otimes1_{m_c}`$ con $`r_c\in M_{k_c}`$
hermíticas; si $`\lambda`$ tiene multiplicidad $`\mu_c(\lambda)`$ en $`r_c`$, entonces $`d_\lambda=\sum_c\mu_c(\lambda)m_c`$ y,
como $`\mu^2\ge\mu`$ para enteros $`\mu\ge0`$,
```math
N_R=\sum_\lambda\Big(\sum_c\mu_c(\lambda)m_c\Big)^2\ge\sum_{\lambda,c}\mu_c(\lambda)^2m_c^2\ge\sum_{\lambda,c}\mu_c(\lambda)m_c^2=\sum_ck_cm_c^2,
```
con igualdad si y solo si los valores propios de todas las $`r_c`$ son simples y distintos dos a dos. Para el espín
total de $`L=3,4,5`$ espines, el código de abajo encuentra la igualdad (observación numérica). Incluso el mínimo es
grande: para ocho espines $`\frac12`$ ($`n=256`$), $`\sum_ck_cm_c^2=4900`$, y el sistema reducido, una matriz
$`(n^2+N_R)\times N_R`$, tiene $`3.5\cdot10^8`$ entradas complejas. Por eso, para cadenas de espines, la sección 4.6 usa
el módulo `su2`, que conoce la respuesta de antemano.

```python
for n in (3, 4, 5):
    Dn = dc.from_generators(dc.total_spin(n))
    _, info = dc.commutant_basis(dc.total_spin(n), return_info=True)
    N_min = sum(k * m * m for k, m in zip(Dn.k_factors, Dn.multiplicities))
    assert info['unknowns'] >= N_min                                  # N_R >= sum_c k_c m_c^2
    print(n, info['unknowns'], N_min)                                # observed: equality
```

<a id="s-4-4"></a>

## 4.4 $`U(1)`$ y sectores de carga

*Apuntes:* [sección 4.4](../../lecture_notes/es/04_symmetries.md#s-4-4); [proposición 4.10](../../lecture_notes/es/04_symmetries.md#r-4-10).

Para $`n`$ cúbits, la carga es el operador número $`\hat N=\sum_jn_j`$, con $`n_j=|1\rangle\langle1|`$ sobre el cúbit $`j`$, y los
operadores que la conservan forman $`U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C)`$ ([apuntes de clase, sección 4.4](../../lecture_notes/es/04_symmetries.md#s-4-4)),
que es `chains.SpinChain(n).algebra` (bit $`1`$ = espín hacia arriba, $`\hat N=S^z+\frac n2`$; en `decompose.total_spin`, en
cambio, el espín hacia arriba es el *primer* vector de la base de $`\mathbb C^2`$, lo que reetiqueta la base pero no los
sectores). Todo hamiltoniano XXZ está en $`U(1)'`$. El código encuentra los sectores a partir únicamente de $`\hat N`$ y
diagonaliza sector por sector un anillo XXZ de cuatro espines.

```python
import functools
from torch_vn_algebra import SpinChain

L = 4
Sx, Sy, Sz = dc.total_spin(L)                              # (16, 16) each, site 0 most significant
N_op = Sz + L / 2 * torch.eye(2 ** L, dtype=dt)            # number of up spins
DU = dc.from_generators(N_op)                              # the algebra generated by U(1)
print(DU)
chain = SpinChain(L, boundary='periodic')                  # U(1)' = (+)_w M_binom(L, w)
assert set(DU.k_factors) == {1} and sorted(DU.multiplicities) == sorted(chain.algebra.k_factors)
assert sorted(chain.algebra.k_factors) == sorted(math.comb(L, w) for w in range(L + 1))

sx, sy, sz = dc.spin_matrices(0.5)
def two_site(a, b, i, j, n=L):                             # a on site i, b on site j of n sites
    ops = [a if k == i else b if k == j else torch.eye(2, dtype=dt) for k in range(n)]
    return functools.reduce(torch.kron, ops)

H = sum(two_site(sx, sx, i, j) + two_site(sy, sy, i, j) + 0.5 * two_site(sz, sz, i, j)
        for i, j in chain.bonds)                           # XXZ ring, Delta = 0.5
assert close(H @ N_op, N_op @ H)                           # [H, N] = 0 ...
assert close(DU.from_commutant(DU.to_commutant(H)), H)     # ... so H is block diagonal (Proposition 4.10)
Hc = chain.xxz(J=1.0, Delta=0.5)                           # the same H, built sector by sector
E_sec = torch.cat([torch.linalg.eigvalsh(Hc.matrix[0, N, :k, :k]) for N, k in enumerate(chain.algebra.k_factors)])
assert close(torch.sort(E_sec).values, torch.linalg.eigvalsh(H))
```

<a id="s-4-5"></a>

## 4.5 $`SU(2)`$: el espín

*Apuntes:* [sección 4.5](../../lecture_notes/es/04_symmetries.md#s-4-5); [teorema 4.11](../../lecture_notes/es/04_symmetries.md#r-4-11), [proposición 4.13](../../lecture_notes/es/04_symmetries.md#r-4-13), [teorema 4.14](../../lecture_notes/es/04_symmetries.md#r-4-14). El espín $`s`$ es `decompose.spin_matrices(s)`, y el espín total de $`L`$ espines $`\frac12`$ es `decompose.total_spin(L)`.

Para $`L=4`$, según los [apuntes de clase, proposición 4.13](../../lecture_notes/es/04_symmetries.md#r-4-13), $`\{S_x,S_y,S_z\}''\cong M_5\oplus M_3\oplus M_1`$
con multiplicidades $`1,3,2`$, de dimensión $`35`$, y el conmutante es $`\{S_x,S_y,S_z\}'\cong M_1\oplus M_3\oplus M_2`$, de
dimensión $`14`$. La biblioteca encuentra lo mismo a partir únicamente de los tres generadores (sectores ordenados por
$`(k_c,m_c)`$), y $`\mathbf S^2=S(S+1)`$ sobre el espín $`S`$ ([teorema 4.11](../../lecture_notes/es/04_symmetries.md#r-4-11)):

```python
S4 = dc.total_spin(4)
assert close(S4[0] @ S4[1] - S4[1] @ S4[0], 1j * S4[2])                # [S_x, S_y] = i S_z
D4 = dc.from_generators(S4)
print(D4)
assert D4.k_factors == (1, 3, 5) and D4.multiplicities == (2, 3, 1)     # M_1, M_3, M_5 with m = 2, 3, 1
assert D4.verify()['ok'] and D4.dims == (35, 14)
m_S = {S: math.comb(4, 2 - S) - (math.comb(4, 1 - S) if S < 2 else 0) for S in (0, 1, 2)}
assert [m_S[(k - 1) // 2] for k in D4.k_factors] == list(D4.multiplicities)    # Proposition 4.13
casimir = D4.to_algebra((S4 @ S4).sum(0))                                  # S^2 = S(S+1) on spin S
for c, k in enumerate(D4.k_factors):
    assert close(casimir.matrix[0, c, :k, :k], (k - 1) / 2 * (k + 1) / 2 * torch.eye(k))
```

Según los [apuntes de clase, teorema 4.14](../../lecture_notes/es/04_symmetries.md#r-4-14) (dualidad de Schur–Weyl), las permutaciones de los
sitios generan el conmutante del espín total: $`\mathrm{span}\,\pi(S_L)=\{S_x,S_y,S_z\}'`$. Para $`L\le5`$ el código compara
numéricamente los dos espacios: por la parte (a) del teorema, $`\mathrm{span}\,\pi(S_L)\subseteq\{S_x,S_y,S_z\}'`$, basta con
que tengan la misma dimensión.

```python
for n in range(2, 6):
    comm = dc.commutant_basis(dc.total_spin(n)).reshape(-1, 4 ** n)          # a basis of {S_a}'
    perm = dc.permutation_representation(n, 2).reshape(-1, 4 ** n)          # all n! permutations
    dims = [torch.linalg.matrix_rank(A).item() for A in (comm, perm, torch.cat([comm, perm]))]
    print(n, dims)
    assert dims[0] == dims[1] == dims[2]                                     # span pi(S_n) = {S_a}'
Dp = dc.from_group(dc.permutation_representation(4, 2))                   # pi(S_4)'' = {S_a}'
print(Dp)
assert sorted(zip(Dp.multiplicities, Dp.k_factors)) == sorted(zip(D4.k_factors, D4.multiplicities))
```

Las dimensiones impresas son $`2,5,14,42`$, los números de Catalan $`\binom{2L}L/(L+1)`$; el corolario B.10 demuestra
$`\dim\{S_a\}'=\sum_Sm_S^2=\binom{2L}L/(L+1)`$ para todo $`L`$. La última línea es la
dualidad en acción: $`\pi(S_4)''`$ tiene bloques $`k=(1,2,3)`$ con multiplicidades $`(5,1,3)`$, los pares $`(k_c,m_c)`$ del
espín total con los papeles intercambiados. En general (corolario B.9), el espacio de multiplicidad del espín $`S`$
porta una representación irreducible de $`S_L`$, de dimensión $`m_S`$, y estas son no equivalentes dos a dos; aquí las
dimensiones son $`1,3,2`$ para $`S=2,1,0`$.

<a id="s-4-6"></a>

## 4.6 Sistemas de espines

*Apuntes:* [sección 4.6](../../lecture_notes/es/04_symmetries.md#s-4-6); [proposición 4.15](../../lecture_notes/es/04_symmetries.md#r-4-15), [teorema 4.16](../../lecture_notes/es/04_symmetries.md#r-4-16), [corolario 4.17](../../lecture_notes/es/04_symmetries.md#r-4-17).

Según los [apuntes de clase, corolario 4.17](../../lecture_notes/es/04_symmetries.md#r-4-17), las bases ortonormales de los vectores de peso máximo de
$`L`$ espines $`\frac12`$ pueden construirse sitio a sitio con los coeficientes de Clebsch–Gordan del [teorema 4.16](../../lecture_notes/es/04_symmetries.md#r-4-16).
Etiquetar cada vector de la base por su sucesión de espines intermedios $`S_1=\frac12,S_2,\dots,S_L=S`$ da la base de
Yamanouchi que construye `su2.SU2Chain` (su código fuente usa exactamente estos coeficientes; el código de abajo solo
comprueba las conclusiones del corolario 4.17), dentro del sector $`U(1)`$ $`\hat N=\frac L2+S`$ de `SpinChain` (cuya
ordenación de los sitios es un reetiquetado de los factores tensoriales, que conmuta con el espín total por el lema
B.1(b)). El código comprueba el teorema 4.16 frente a matrices densas y frente a `su2.clebsch_gordan` de la biblioteca, y
después las conclusiones del corolario 4.17 para las bases de ocho espines. Comprobamos también la
[proposición 4.15(a)](../../lecture_notes/es/04_symmetries.md#r-4-15) sobre un anillo de cuatro espines con `D4`, y comparamos los tamaños de los bloques
para ocho espines.

```python
from torch_vn_algebra import su2

def ket(j, m):                                            # |j, m> of Theorem 4.11 (zero if |m| > j)
    v = torch.zeros(int(2 * j + 1), dtype=dt)
    if abs(m) <= j:
        v[int(j - m)] = 1
    return v

up, dn = ket(0.5, 0.5), ket(0.5, -0.5)
for j in (0.5, 1.0, 1.5, 2.0):                            # Theorem 4.16
    Sj = dc.spin_matrices(j)
    sp_j = [math.sqrt(j * (j + 1) - (j - k) * (j - k + 1)) for k in range(1, int(2 * j) + 1)]
    assert close(Sj[0] + 1j * Sj[1], torch.diag(torch.tensor(sp_j, dtype=dt), 1))     # S_+ of Theorem 4.11
    assert close(Sj[2], torch.diag(torch.tensor([j - k for k in range(int(2 * j) + 1)], dtype=dt)))
    Sc = torch.stack([torch.kron(a, torch.eye(2, dtype=dt)) + torch.kron(torch.eye(len(a), dtype=dt), b)
                      for a, b in zip(Sj, dc.spin_matrices(0.5))])
    Sp, Sm = Sc[0] + 1j * Sc[1], Sc[0] - 1j * Sc[1]
    cols = []
    for J in (j + 0.5, j - 0.5):
        for M in [J - k for k in range(int(2 * J + 1))]:
            al, be = math.sqrt((j + M + 0.5) / (2 * j + 1)), math.sqrt((j - M + 0.5) / (2 * j + 1))
            p, q = (al, be) if J > j else (-be, al)
            x = p * torch.kron(ket(j, M - 0.5), up) + q * torch.kron(ket(j, M + 0.5), dn)
            assert close(Sc[2] @ x, M * x)
            if M == J:
                assert close(Sp @ x, 0)                               # highest weight
            else:
                assert close(Sm @ cols[-1], math.sqrt((J + M + 1) * (J - M)) * x)   # lowering, (c)
            for s2, e in ((0.5, up), (-0.5, dn)):                     # the library's Racah formula agrees
                cg = su2.clebsch_gordan(j, M - s2, 0.5, s2, J, M) if abs(M - s2) <= j else 0.0
                assert abs(cg - torch.vdot(torch.kron(ket(j, M - s2), e), x).real.item()) < 1e-12
            cols.append(x)
    U = torch.stack(cols, 1)
    assert U.shape == (2 * len(Sj[0]), 2 * len(Sj[0])) and close(U.mH @ U, torch.eye(U.shape[1]))   # (a)

def swap(i, j, n):
    p = list(range(n))
    p[i], p[j] = j, i
    return dc.permutation_representation(n, 2, [p])[0]

ring = [(i, (i + 1) % 4) for i in range(4)]
H4 = sum(swap(i, j, 4) / 2 - torch.eye(16, dtype=dt) / 4 for i, j in ring)        # Dirac's identity
assert close(H4, sum(two_site(a, a, i, j) for a in (sx, sy, sz) for i, j in ring))
h = D4.to_commutant(H4)                                   # blocks H_S of sizes m_S = 2, 3, 1
assert close(D4.from_commutant(h), H4)                    # H4 lies in {S_a}'
levels = torch.cat([torch.linalg.eigvalsh(h.matrix[0, c, :m, :m]).repeat_interleave(k)
                    for c, (k, m) in enumerate(zip(D4.k_factors, D4.multiplicities))])
assert close(torch.sort(levels).values, torch.linalg.eigvalsh(H4))   # level of H_S repeated 2S+1 times

u1, sym = SpinChain(8, boundary='periodic'), su2.SU2Chain(8, boundary='periodic')
for S, m in zip(sym.spins, sym.multiplicities):           # Corollary 4.17: orthonormal highest-weight bases
    B, N = sym.hw_basis(S), int(4 + S)                     # columns in the U(1) sector N = L/2 + S
    assert B.shape == (math.comb(8, N), m) and close(B.T @ B, torch.eye(m, dtype=B.dtype))
    assert N == 8 or close(sym.raise_(B, N), 0)            # S_+ B = 0 (S_z = S holds in sector N)
H_u1, H_su2 = u1.xxz(J=1.0, Delta=1.0), sym.heisenberg(J=1.0)
E_u1 = torch.sort(torch.cat([torch.linalg.eigvalsh(H_u1.matrix[0, N, :k, :k])
                             for N, k in enumerate(u1.algebra.k_factors)])).values
E_su2, S_lev, _ = sym.full_spectrum(H_su2, return_quantum_numbers=True)   # level of block S, 2S+1 times
E_su2, S_lev = E_su2[0], S_lev[0]
assert close(E_u1, E_su2)
assert S_lev[0] == 0 and (S_lev[1:4] == 1).all() and E_su2[4] - E_su2[3] > 1e-3     # singlet, then triplet
assert E_su2[1] - E_su2[0] > 1e-3 and close(E_su2[1:4], E_su2[1].expand(3))       # simple singlet (Theorem 9.9)
assert close(sym.restrict(H_u1).matrix, H_su2.matrix)      # Proposition 4.15(b)
flops = lambda ks: sum(k ** 3 for k in ks)
print(u1.algebra.k_factors, sym.multiplicities, flops(u1.algebra.k_factors), flops(sym.multiplicities))
print(f"ground state energy {E_su2[0].item():.6f}, first excited {E_su2[1].item():.6f} (x3)")
```

Para el anillo de ocho espines, los bloques de $`U(1)`$ tienen tamaños $`\binom8w=1,8,28,56,70,56,28,8,1`$ (sección 4.4), y
los bloques de $`SU(2)`$ tamaños $`m_S=14,28,20,7,1`$ para $`S=0,\dots,4`$ (proposición 4.13). El coste de diagonalizar,
$`\sum k^3`$, baja de $`739162`$ a $`33040`$, en un factor $`22`$; la dimensión del álgebra, $`\sum k^2`$, de
$`\sum_w\binom8w\binom8{8-w}=\binom{16}8=12870`$ (elegir $`8`$ de $`16`$ elementos, $`w`$ de ellos entre los $`8`$ primeros) a
$`\sum_Sm_S^2=1430=\binom{16}8/9`$ (corolario B.10). Por el teorema 9.9 (el anillo tiene longitud par y $`J>0`$), el estado
fundamental es un único singlete y $`E_0(0)<E_0(1)`$. *Observación numérica* (calculada por el código de
arriba): la energía fundamental es $`-3.651093`$, y el primer nivel excitado es un triplete en $`-3.128419`$.

<a id="s-4-7"></a>

## 4.7 Álgebras reales: tipo real, complejo y cuaterniónico

*Apuntes:* [sección 4.7](../../lecture_notes/es/04_symmetries.md#s-4-7). Para una $`*`$-álgebra $`\mathcal A`$ cerrada para la conjugación compleja,
con una descomposición $`V`$ cuyos bloques de columnas $`V_c`$ pertenecen a los sectores, el [teorema 4.20](../../lecture_notes/es/04_symmetries.md#r-4-20) da la
permutación $`\sigma`$ de los sectores con $`\bar P_c=P_{\sigma(c)}`$, unitarias $`J_c`$, $`M_c`$ con $`\overline{V_c}=V_{\sigma(c)}(J_c\otimes M_c)`$
y, para $`\sigma(c)=c`$, el signo $`\varepsilon_c`$ con $`J_c\bar J_c=\varepsilon_c1`$; los tipos y el álgebra real $`\mathcal A_{\mathbb R}`$
están en el [teorema 4.21](../../lecture_notes/es/04_symmetries.md#r-4-21), el lema de Schur real es el [teorema 4.22](../../lecture_notes/es/04_symmetries.md#r-4-22), y los
elementos reales del conmutante los describe la [proposición 4.23](../../lecture_notes/es/04_symmetries.md#r-4-23).

`Decomposition.real_structure()` calcula $`\sigma`$, $`J_c`$, $`M_c`$ y los tipos para una descomposición devuelta por
`wedderburn` (`from_generators`, `from_group`). Comprueba que los conjugados de los generadores almacenados están en
$`V(\bigoplus_cM_{k_c}\otimes1)V^*`$ (residuo relativo `'conjugation'`) y forma la unitaria $`P=V^*\bar V`$, cuyo bloque
de columnas $`c`$ es $`V^*\overline{V_c}`$: por el teorema 4.20(a), está soportado en las filas del sector $`\sigma(c)`$,
donde es igual a $`J_c\otimes M_c`$. La biblioteca toma como $`\sigma(c)`$ el sector que lleva la mayor parte del bloque,
exige $`(k_{\sigma(c)},m_{\sigma(c)})=(k_c,m_c)`$ y acota el resto (`'leak'`); reordena el bloque en la matriz
$`k_c^2\times m_c^2`$ de entradas $`T_{(a,b),(j,l)}=T_{(a,j),(b,l)}`$, que para $`T=J\otimes M`$ es la matriz de rango uno
$`\mathrm{vec}\,J\,(\mathrm{vec}\,M)^T`$, lee $`J_c`$ y $`M_c`$ en sus vectores singulares principales (con
$`\|J_c\|_F^2=k_c`$), y registra la distancia del bloque a $`J_c\otimes M_c`$ (`'kron'`). Si uno de estos tres residuos
supera la tolerancia, los tamaños no coinciden o los compañeros no forman una involución, lanza `ValueError` (el
álgebra no es cerrada para la conjugación). Para $`\sigma(c)=c`$ toma como $`\varepsilon_c`$ el signo de
$`\mathrm{Re}\,\mathrm{Tr}(J_c\bar J_c)/k_c`$ —exactamente $`\pm1`$ en aritmética exacta, de modo que los dos
candidatos distan $`2`$— y registra como `'sign'` el máximo sobre estos sectores de
$`\|J_c\bar J_c-\varepsilon_c\|_F/\sqrt{k_c}`$ y $`\|M_c\bar M_c-\varepsilon_c\|_F/\sqrt{m_c}`$, elevado al menos a $`1`$ si
algún $`\varepsilon_c=-1`$ con $`k_c`$ o $`m_c`$ impar (imposible por el teorema 4.20(b)). El resultado tiene `types`,
`partner` ($`\sigma`$), `J`, `M`, `residuals` y `ok`, que indica que todos los residuos son a lo sumo la tolerancia (por
defecto, el umbral de aceptación de la descomposición); entonces
$`|\mathrm{Tr}(J_c\bar J_c)/k_c-\varepsilon_c|\le\|J_c\bar J_c-\varepsilon_c\|_F/\sqrt{k_c}`$ (Cauchy–Schwarz)
también es a lo sumo la tolerancia.

`Decomposition.real_basis()` devuelve una matriz ortogonal real $`O`$ de tamaño $`n\times n`$ y una lista que describe sus
bloques de columnas (tipo, sectores, rango de columnas): para un sector de tipo real, las columnas
$`O_c=V_c(W_J\otimes W_M)`$ del teorema 4.21(a), de modo que $`O_c^TaO_c=\alpha_c(a)\otimes1_{m_c}`$ con $`\alpha_c(a)`$
real; para un sector de tipo cuaterniónico, las columnas $`V_cW`$ con $`J_c\otimes M_c=WW^T`$ y $`(J_c\otimes M_c)\bar W=W`$
(lema 4.19(c) con $`\varepsilon=1`$, pues $`(J_c\otimes M_c)\overline{(J_c\otimes M_c)}=\varepsilon_c^2=1`$), que son
reales por el cálculo del teorema 4.21(a); y para un par $`c,\sigma(c)`$ (con $`c<\sigma(c)`$), las columnas
$`\sqrt2\mathrm{Re}\,V_c`$ seguidas de $`\sqrt2\mathrm{Im}\,V_c`$. Son ortonormales:
$`V_c^*\overline{V_c}=V_c^*V_{\sigma(c)}(J_c\otimes M_c)=0`$ y $`V_c^T\overline{V_c}=\overline{V_c^*V_c}=1`$, junto con
sus conjugadas, dan
$`(\mathrm{Re}\,V_c)^T\mathrm{Re}\,V_c=(\mathrm{Im}\,V_c)^T\mathrm{Im}\,V_c=\frac12`$ y
$`(\mathrm{Re}\,V_c)^T\mathrm{Im}\,V_c=0`$; las columnas de bloques distintos están en las imágenes de las
proyecciones ortogonales $`P_{\mathcal O}=\sum_{e\in\mathcal O}P_e`$ de órbitas distintas $`\mathcal O`$ de $`\sigma`$ (pues
$`\overline{V_c}=V_{\sigma(c)}(J_c\otimes M_c)`$), así que $`O`$ es ortogonal. Tomando partes real e imaginaria de
$`aV_c=V_cx`$, $`x=\pi_c(a)\otimes1_{m_c}`$, para $`a`$ real,
```math
O^TaO=\begin{pmatrix}\mathrm{Re}\,x&\mathrm{Im}\,x\\-\mathrm{Im}\,x&\mathrm{Re}\,x\end{pmatrix}\quad\text{on these columns.}
```
Los factores de Takagi $`W`$ con $`WW^T=J`$ de una unitaria simétrica $`J=X+iY`$ se calculan a partir de una diagonalización
ortogonal de $`X+tY`$ para un $`t`$ aleatorio ($`X`$ e $`Y`$ son matrices reales simétricas que conmutan, pues la parte
imaginaria de $`J\bar J=JJ^*=1`$ es $`YX-XY`$); la biblioteca acepta el resultado según sus residuos
$`\|WW^T-J\|_F/\|J\|_F`$ y $`\|\mathrm{Im}(V_cW)\|_F/\sqrt{k_cm_c}`$ y emite un aviso (`DecompositionWarning`) si
uno de ellos supera $`\max(\text{tolerancia},10^3\,\epsilon_{\rm mach}\,n)`$.

La [sección 4.7](../../lecture_notes/es/04_symmetries.md#s-4-7) de los apuntes de clase estudia tres álgebras generadas por matrices reales: (i)
$`L_a\otimes1_2`$ sobre $`\mathbb R^3\otimes\mathbb R^2`$, con los generadores $`(L_a)_{bc}=-\epsilon_{abc}`$ de las rotaciones de
$`\mathbb R^3`$: un sector de tipo real, $`k=3`$, $`m=2`$; (ii) $`\Omega_2\otimes1_2`$ sobre $`\mathbb R^2\otimes\mathbb R^2`$: un par de
tipo complejo, $`k=1`$, $`m=2`$; (iii) las multiplicaciones por la izquierda $`1_2\otimes L_q`$ por los cuaterniones sobre
$`\mathbb H^2=\mathbb R^8`$: un sector de tipo cuaterniónico, $`k=2`$, $`m=4`$. El código calcula la descomposición y la
estructura real de cada una, y comprueba el teorema 4.20(a),(b) y la forma por bloques de $`O^TaO`$ para el tipo real
(teorema 4.21(a)) y para el par.

```python
def qmul(p, q):                                           # H = CD(C), Section E.6: (a, b)(c, d)
    (a1, b1), (a2, b2) = p, q
    return (a1 * a2 - b2.conjugate() * b1, b2 * a1 + b1 * a2.conjugate())

def q_left(q):                                            # left multiplication by q on H = R^4
    cols = []                                             # coordinates (Re a, Im a, Re b, Im b)
    for e in torch.eye(4, dtype=torch.float64).tolist():
        u, v = qmul(q, (complex(e[0], e[1]), complex(e[2], e[3])))
        cols.append([u.real, u.imag, v.real, v.imag])
    return torch.tensor(cols, dtype=torch.float64).T

L_so3 = torch.zeros(3, 3, 3, dtype=torch.float64)         # (L_a)_bc = -epsilon_abc
for i, (b1, b2) in enumerate([(1, 2), (2, 0), (0, 1)]):
    L_so3[i, b1, b2], L_so3[i, b2, b1] = -1.0, 1.0
R_rot = torch.tensor([[[0.0, -1.0], [1.0, 0.0]]])         # Omega_2
L_H = torch.stack([q_left(u) for u in ((1j, 0), (0, 1), (0, 1j))])            # L_I, L_J, L_K
assert close(L_H[0] @ L_H[1], L_H[2]) and close(q_left((2 - 1j, 3j)).T, q_left((2 + 1j, -3j)))   # IJ = K, L_qbar = L_q^T
eye2 = torch.eye(2, dtype=torch.float64)
real_examples = {'R': (torch.kron(L_so3, eye2[None]), (3,), (2,)),            # L_a (x) 1_2 on R^3 (x) R^2
                 'C': (torch.kron(R_rot, eye2[None]), (1, 1), (2, 2)),        # Omega_2 (x) 1_2
                 'H': (torch.kron(eye2[None], L_H), (2,), (4,))}              # L_q on H^2 = R^8
dec_real = {}
for t, (gens, ks, ms) in real_examples.items():
    D_t = dc.from_generators(gens)
    rs_t, (O_t, blocks_t) = D_t.real_structure(), D_t.real_basis()
    print(t, D_t, rs_t)
    assert D_t.verify()['ok'] and rs_t.ok and (D_t.k_factors, D_t.multiplicities) == (ks, ms)
    assert set(rs_t.types) == {t} and close(O_t.T @ O_t, torch.eye(D_t.n))     # real orthogonal basis
    for c, d in enumerate(rs_t.partner):                                      # Theorem 4.20
        o_c, o_d, w = D_t.offsets[c], D_t.offsets[d], ks[c] * ms[c]
        assert close(D_t.V[:, o_c:o_c + w].conj(), D_t.V[:, o_d:o_d + w] @ torch.kron(rs_t.J[c], rs_t.M[c]))
        if d == c:                                                            # J conj(J) = M conj(M) = +-1
            eps = 1.0 if t == 'R' else -1.0
            assert close(rs_t.J[c] @ rs_t.J[c].conj(), eps * torch.eye(ks[c]))
            assert close(rs_t.M[c] @ rs_t.M[c].conj(), eps * torch.eye(ms[c]))
    dec_real[t] = D_t, rs_t, O_t

D_t, rs_t, O_t = dec_real['R']
for g in real_examples['R'][0]:                           # Theorem 4.21(a): O^T a O = alpha (x) 1_2, alpha real
    G_t = O_t.T @ g @ O_t
    assert close(G_t, torch.kron(G_t[0::2, 0::2], eye2))
D_t, rs_t, O_t = dec_real['C']                            # the pair: [[Re x, Im x], [-Im x, Re x]]
x_t = torch.kron(D_t.to_algebra(real_examples['C'][0]).matrix[0, 0, :1, :1], torch.eye(2, dtype=dt))
assert close(O_t.T @ real_examples['C'][0][0] @ O_t,
             torch.cat([torch.cat([x_t.real, x_t.imag], 1), torch.cat([-x_t.imag, x_t.real], 1)]))
```

El código comprueba el teorema 4.22 sobre los módulos reales irreducibles $`\mathbb R^3`$, $`\mathbb R^2`$, $`\mathbb R^4`$
de los tres ejemplos (los tipos y multiplicidades que encuentra la biblioteca —real con $`m=1`$, un par complejo con
$`m=1`$, cuaterniónico con $`m=2`$— los hacen irreducibles por el teorema 4.22(c)), calculando el conmutante real como el
espacio nulo de las ecuaciones lineales reales $`Xs=sX`$, $`Xs^T=s^TX`$, y después la proposición 4.23 para un elemento
real simétrico aleatorio del conmutante en cada ejemplo.

```python
def real_commutant_dim(gens):                             # dim_R {X real : X s = s X, X s^T = s^T X}
    gens = torch.cat([gens, gens.mT])
    eye_n = torch.eye(gens.shape[-1], dtype=torch.float64)
    eqs = torch.cat([torch.kron(eye_n, s) - torch.kron(s.T.contiguous(), eye_n) for s in gens])
    return gens.shape[-1] ** 2 - torch.linalg.matrix_rank(eqs).item()

for gens, types, ms, dim_D in ((L_so3, ('R',), (1,), 1), (R_rot, ('C', 'C'), (1, 1), 2), (L_H, ('H',), (2,), 4)):
    D_irr = dc.from_generators(gens)                      # Theorem 4.22(c): R^3, R^2, R^4 are irreducible
    assert D_irr.real_structure().types == types and D_irr.multiplicities == ms
    assert real_commutant_dim(gens) == dim_D              # Theorem 4.22(b): D_U = R, C, H

gen_X = torch.Generator().manual_seed(1)
for t, (gens, ks, ms) in real_examples.items():           # Proposition 4.23
    D_t, rs_t, O_t = dec_real[t]
    B_t = dc.commutant_basis(gens)
    X_t = torch.einsum('r,rij->ij', torch.randn(len(B_t), dtype=dt, generator=gen_X), B_t).real
    X_t = X_t + X_t.T                                     # real symmetric, in A' (Lemma 4.18(b))
    assert close(X_t @ gens, gens @ X_t)
    x_t = D_t.to_commutant(X_t).matrix[0]
    ev = [torch.linalg.eigvalsh(x_t[c, :m, :m]) for c, m in enumerate(ms)]
    if t == 'R':                                          # (a): O^T X O = 1_3 (x) y, y real symmetric
        Y_t = O_t.T @ X_t @ O_t
        assert close(Y_t, torch.kron(torch.eye(3, dtype=torch.float64), Y_t[:2, :2]))
    if t == 'C':                                          # (c): the two sectors have the same spectrum
        assert close(ev[0], ev[1])
    if t == 'H':                                          # (b): Kramers pairs, each level 2 k_c = 4 times in X
        assert close(ev[0][0::2], ev[0][1::2])
        assert close(torch.linalg.eigvalsh(X_t), ev[0].repeat_interleave(2))
        print('quaternionic block x_c:', [round(v, 4) for v in ev[0].tolist()])
```

<a id="s-4-8"></a>

## 4.8 Con qué comparar

Para grupos finitos, GAP y SageMath calculan tablas de caracteres de forma exacta, y el corolario 4.5 da entonces
las multiplicidades; RepLAB (MATLAB; D. Rosset, F. Montealegre-Mora, J.-D. Bancal) descompone numéricamente
representaciones de grupos finitos y compactos, también con elementos aleatorios del conmutante, y
SymbolicWedderburn.jl (Julia) calcula descomposiciones de Wedderburn para la reducción por simetría de programas
semidefinidos. En física, QuSpin y otros códigos de diagonalización exacta diagonalizan por bloques las simetrías
$`U(1)`$ y de red, ITensor y TeNPy manejan cargas abelianas en redes tensoriales, y TensorKit.jl trata $`SU(2)`$ con
árboles de fusión. El módulo `decompose` de la biblioteca es un método denso para matrices pequeñas ($`n`$ de hasta
unos pocos cientos); acepta cualquier lista de matrices (elementos de un grupo o generadores), certifica su
resultado (proposición 4.9) y devuelve el álgebra y su conmutante como objetos `TypeIAlgebra` con `to_algebra`,
`from_algebra`, `to_commutant`, `from_commutant`, de modo que los capítulos 3 y 5–8 se aplican de inmediato. Para
cadenas de espín $`\frac12`$, `su2.SU2Chain` da directamente los bloques de $`SU(2)`$ (corolario 4.17); acepta hasta
$`L=18`$ espines, un límite fijado por la memoria de las bases.

<a id="s-4-9"></a>

## 4.9 Estabilidad de la descomposición numérica

*Apuntes:* [sección 4.9](../../lecture_notes/es/04_symmetries.md#s-4-9); [lema 4.24](../../lecture_notes/es/04_symmetries.md#r-4-24), [teorema 4.25](../../lecture_notes/es/04_symmetries.md#r-4-25), [proposición 4.26](../../lecture_notes/es/04_symmetries.md#r-4-26), [proposición 4.27](../../lecture_notes/es/04_symmetries.md#r-4-27), [teorema 4.28](../../lecture_notes/es/04_symmetries.md#r-4-28).

El conmutante de generadores con ruido se decide mediante una tolerancia: el [teorema 4.25](../../lecture_notes/es/04_symmetries.md#r-4-25) dice que un umbral
$`\theta`$ sobre los valores singulares de la aplicación conmutador $`X\mapsto([H_j,X])_j`$ da la dimensión exacta $`d`$ del conmutante si
$`2\eta<\theta<\gamma-2\eta`$, con $`\eta`$ el ruido y $`\gamma`$ el hueco de la aplicación exacta, y que la cota $`2\eta`$ sobre los valores
singulares nulos perturbados se alcanza. El primer bloque calcula los valores singulares del ejemplo del teorema 4.25(c),
$`\tilde H_1=\mathrm{diag}(1,1,0)+\eta(E_{12}+E_{21})`$, y decide el conmutante con `commutant_basis` por debajo, dentro y por encima de la
ventana (la biblioteca normaliza los generadores a norma de Frobenius unidad; de ahí el `tol` reescalado).

```python
import warnings

def commutator_svals(Hs):
    """Singular values, in increasing order, of X -> ([H_j, X])_j with Hilbert-Schmidt norms (row-major vec)."""
    n = Hs.shape[-1]
    eye = torch.eye(n, dtype=dt)
    L = torch.cat([torch.kron(H, eye) - torch.kron(eye, H.T.contiguous()) for H in Hs])
    return torch.linalg.svdvals(L).flip(0)

eta = 1e-3
H1 = torch.diag(torch.tensor([1.0, 1.0, 0.0], dtype=dt))
H1t = H1.clone()
H1t[0, 1] = H1t[1, 0] = eta                               # a symmetry broken at the level eta
s = commutator_svals(H1t[None])
expected = torch.tensor([0, 0, 0, 2 * eta, 2 * eta, 1 - eta, 1 - eta, 1 + eta, 1 + eta], dtype=torch.float64)
assert close(s, expected, 1e-12)                          # Theorem 4.25(a) and (c): the bound 2 eta is attained
assert close(commutator_svals(H1[None])[4:6], [0, 1])      # exact generator: d = 5, gap 1
scale = torch.linalg.matrix_norm(H1t).item()              # theta = tol for unit-norm generators
for theta, dim in ((0.5 * eta, 3), (20 * eta, 5), (2.0, 9)):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', dc.DecompositionWarning)
        assert dc.commutant_basis(H1t, tol=theta / scale).shape[0] == dim
print('theta < 2 eta: dim 3 (symmetry lost); in the window: 5; above the gap: 9')
```

A continuación, los tres operadores de espín de tres espines $`\frac12`$ ($`k=(2,4)`$, $`m=(2,1)`$, $`d=5`$, hueco $`\gamma=1/\sqrt3`$ tras la
normalización), perturbados por ruido hermítico de tamaño creciente $`\eta`$. Para cada nivel el código comprueba el teorema 4.25(a) en los
valores singulares $`\sigma_d`$ y $`\sigma_{d+1}`$, coloca `tol` en el centro de la ventana observada (en escala logarítmica), comprueba la
dimensión y la cota $`\|\tilde\Pi-\Pi\|\le4\eta/\gamma`$ del teorema 4.25(b), y descompone. Con el `tol` por defecto
$`=\epsilon_{\rm mach}^{1/2}`$ el ruido se toma en serio en cuanto supera la tolerancia: los generadores con ruido generan $`M_8`$, que es la
respuesta correcta para ellos, pero no la simetría que se busca. (El `DecompositionWarning` filtrado aquí indica que un valor singular está a
menos de un factor $`\epsilon_{\rm mach}^{-1/8}`$ del corte.)

```python
S = dc.total_spin(3)
unit = lambda S: S / torch.linalg.matrix_norm(S)[:, None, None]     # the Hermitian generators the library uses
sv = commutator_svals(unit(S))
d, gamma = 5, sv[5].item()
assert sv[4] < 1e-12 and abs(gamma - 1 / math.sqrt(3)) < 1e-12
B0 = dc.commutant_basis(S).reshape(d, -1)
Pi0 = B0.T @ B0.conj()                                    # orthogonal projection onto vec(A')
gen = torch.Generator().manual_seed(1)
print(' eta       sigma_d   sigma_d+1  tol       (k, m)            default tol')
for eta in (1e-8, 1e-6, 1e-4, 1e-3):
    E = torch.randn(S.shape, dtype=dt, generator=gen)
    E = E + E.mH
    St = S + eta * E / torch.linalg.matrix_norm(E, ord=2)[:, None, None]
    eta_n = torch.linalg.matrix_norm(unit(St) - unit(S), ord=2).norm().item()   # eta of Theorem 4.25
    st = commutator_svals(unit(St))
    assert st[d - 1] <= 2 * eta_n and st[d] >= gamma - 2 * eta_n                # Theorem 4.25(a)
    tol = math.sqrt(st[d - 1] * st[d]) / max(1.0, st[-1].item())
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', dc.DecompositionWarning)
        Bt = dc.commutant_basis(St, tol=tol).reshape(d, -1)                     # Theorem 4.25(b): dimension d
        Dt = dc.from_generators(St, tol=tol)
        Dd = dc.from_generators(St)
    assert torch.linalg.matrix_norm(Bt.T @ Bt.conj() - Pi0, ord=2) <= 4 * eta_n / gamma
    assert (Dt.k_factors, Dt.multiplicities) == ((2, 4), (2, 1))
    print(f'{eta_n:.1e}  {st[d - 1]:.1e}   {st[d]:.3f}      {tol:.1e}   {Dt.k_factors}, {Dt.multiplicities}'
          f'   {Dd.k_factors}, {Dd.multiplicities}')
```

`Decomposition.certificate()` evalúa la [proposición 4.26](../../lecture_notes/es/04_symmetries.md#r-4-26) sobre la salida: el error hacia atrás $`r`$ (la distancia de los
generadores normalizados al álgebra calculada $`\tilde{\mathcal B}`$), el hueco $`\gamma_{\mathcal B}`$ de la aplicación conmutador de sus
proyecciones $`H_j=E(\tilde H_j)`$ —calculado a partir de los bloques y comprobado aquí frente a los valores singulares densos— y la ventana
$`(2r,\gamma_{\mathcal B}-2r)`$ de umbrales que certifican $`\dim\mathcal A'=\sum_cm_c^2`$. El resultado es `trusted` si el umbral efectivamente usado
está en esa ventana; si no, el certificado avisa. Una tolerancia por debajo del ruido da $`M_8`$, certificada para estos generadores, pero con
un hueco del orden del ruido: un modelo con más simetría puede estar a distancia $`\gamma_{\mathcal B}/2`$. Una tolerancia por encima del ruido
pero por debajo de $`2r`$ da la estructura correcta, que el certificado no puede garantizar.

```python
c = dc.from_generators(S).certificate()
assert c.certified and c.trusted and c.residuals['backward'] < 1e-12 and abs(c.gap - gamma) < 1e-10
E = torch.randn(S.shape, dtype=dt, generator=gen)
St = S + 1e-4 * (E + E.mH)
with warnings.catch_warnings():
    warnings.simplefilter('ignore', dc.DecompositionWarning)
    Dt = dc.from_generators(St, tol=1e-2)
    D_low = dc.from_generators(St, tol=1e-12)                          # below the noise
    D_mid = dc.from_generators(St, tol=1e-3)                           # above the noise, below 2 r
ct = Dt.certificate()
print(ct, 'window', [f'{w:.1e}' for w in ct.window], 'threshold', f"{ct.rank['threshold']:.1e}")
assert ct.trusted and 1e-5 < ct.residuals['backward'] < 1e-2 and ct.window[0] < ct.rank['threshold'] < ct.window[1]
sb = commutator_svals(Dt.conditional_expectation(unit(St)))           # the projected generators H_j
assert sb[d - 1] < 1e-12 and abs(sb[d] - ct.gap) < 1e-10              # Proposition 4.26(a), (b)
c_low = D_low.certificate()
assert D_low.k_factors == (8,) and c_low.certified and c_low.gap < 1e-3
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    c_mid = D_mid.certificate()
assert D_mid.k_factors == (2, 4) and not c_mid.trusted and caught
print('tol 1e-12:', D_low.k_factors, f'gap {c_low.gap:.1e};  tol 1e-3:', D_mid.k_factors, 'trusted:', c_mid.trusted)
```

Por último, la [proposición 4.27](../../lecture_notes/es/04_symmetries.md#r-4-27): para un vector gaussiano estándar $`G`$ de $`\mathcal A'_h`$ los bloques son
$`\kappa_c=k_c^{-1/2}H_c`$ con matrices GUE independientes, y la probabilidad de que autovalores de sectores distintos se acerquen a menos de
$`x`$ es a lo sumo de orden $`x`$ (aquí la cota es $`\sqrt{2/\pi}\,x\,m_1m_2(\frac1{k_1m_1}+\frac1{k_2m_2})^{-1/2}`$), mientras que dentro del bloque
$`2\times2`$ el espaciamiento es de orden $`x^3`$. Una comprobación de Monte Carlo con $`2\cdot10^4`$ muestras, tomadas en la base ortonormal
`D.basis('commutant')` de $`\mathcal A'`$:

```python
D = dc.from_generators(S)
coef = torch.randn(20000, d, dtype=dt, generator=gen)
Z = torch.einsum('br,rij->bij', coef, D.basis('commutant'))         # standard complex Gaussian vectors of A'
x_c = D.to_commutant((Z + Z.mH) / math.sqrt(2)).matrix              # blocks of G: kappa_1 (2 x 2), kappa_2 (1 x 1)
assert abs((x_c[:, 0, :2, :2].abs() ** 2).sum((-2, -1)).mean().item() - 4 / 2) < 0.05   # E Tr kappa_1^2 = m_1^2 / k_1
ev = torch.linalg.eigvalsh(x_c[:, 0, :2, :2])
cross = (ev - x_c[:, 1, :1, 0].real).abs().min(1).values           # delta_x
(k1, k2), (m1, m2) = D.k_factors, D.multiplicities
for x in (0.01, 0.03, 0.1):
    bound = math.sqrt(2 / math.pi) * x * m1 * m2 / math.sqrt(1 / (k1 * m1) + 1 / (k2 * m2))
    p = (cross < x).double().mean().item()
    print(f'x = {x}:  P(delta_x < x) = {p:.4f} <= {bound:.4f}')
    assert p <= bound
p_in = (ev[:, 1] - ev[:, 0] < 0.2).double().mean().item()
assert p_in <= k1 ** 1.5 * 0.2 ** 3 / (6 * math.sqrt(math.pi))      # Proposition 4.27(c), m = 2
```

<a id="x-4"></a>

## Ejercicios

Los [ejercicios](../../lecture_notes/es/04_symmetries.md#x-4) de los apuntes de clase, comprobados numéricamente (las soluciones de abajo siguen este orden):
1. $`\dim\pi(S_L)'=\binom{d^2+L-1}L`$ para $`d=2`$ y $`L=3,4`$ (los valores $`20`$ y $`35`$);
2. $`1\otimes\frac12\cong\frac12\oplus\frac32`$ y $`1\otimes1\cong0\oplus1\oplus2`$, con `from_generators`;
3. el conmutante de $`e^{i\theta\hat N}`$ sobre cuatro cúbits para $`\theta=1`$ (las cargas $`0,\dots,4`$ dan $`e^{iq}`$ distintos, pues
   $`0<|q-q'|\le4<2\pi`$) y para $`\theta=\pi`$;
4. el desdoblamiento Zeeman: el espectro de $`H+hS_z`$ para el anillo de Heisenberg de ocho espines, con `su2.SU2Chain` y en los
   sectores $`U(1)`$;
5. los sectores de la representación regular de $`Q_8`$, con `real_structure`.

<details>
<summary>Soluciones</summary>

```python
for n, value in ((3, 20), (4, 35)):
    assert dc.commutant_dimension(dc.permutation_representation(n, 2)) == math.comb(n + 3, n) == value
```

```python
def coupled(a, b):                                        # S^(a) (x) 1 + 1 (x) S^(b)
    Sa, Sb = dc.spin_matrices(a), dc.spin_matrices(b)
    ea, eb = torch.eye(Sa.shape[-1], dtype=dt), torch.eye(Sb.shape[-1], dtype=dt)
    return torch.stack([torch.kron(x, eb) + torch.kron(ea, y) for x, y in zip(Sa, Sb)])
for a, b, ks in ((1.0, 0.5, (2, 4)), (1.0, 1.0, (1, 3, 5))):  # 1 (x) 1/2 = 1/2 (+) 3/2, 1 (x) 1 = 0 (+) 1 (+) 2
    D_cg = dc.from_generators(coupled(a, b))
    assert D_cg.k_factors == ks and set(D_cg.multiplicities) == {1}
```

```python
for theta, m in ((1.0, (1, 1, 4, 4, 6)), (math.pi, (8, 8))):   # theta = 1: e^{iq}, q = 0..4, distinct as 4 < 2 pi
    D_theta = dc.from_generators(torch.linalg.matrix_exp(1j * theta * N_op))
    assert D_theta.multiplicities == m and set(D_theta.k_factors) == {1}
```

```python
assert not close(H4 @ Sx - Sx @ H4 + 0.3 * (Sz @ Sx - Sx @ Sz), 0)          # [H4 + 0.3 S_z, S_x] != 0
E_field = sym.full_spectrum(sym.heisenberg(J=1.0, h=0.3))[0]                # E + h M on every multiplet
H_field = u1.xxz(J=1.0, Delta=1.0, h=0.3)
E_u1_field = torch.cat([torch.linalg.eigvalsh(H_field.matrix[0, N, :k, :k])
                        for N, k in enumerate(u1.algebra.k_factors)])
assert close(torch.sort(E_u1_field).values, E_field)
```

```python
Q8 = [(s * a, s * b) for a, b in ((1, 0), (1j, 0), (0, 1), (0, 1j)) for s in (1, -1)]
key = lambda q: tuple(round(v) for z in q for v in (complex(z).real, complex(z).imag))
pos = {key(q): i for i, q in enumerate(Q8)}
reg = torch.zeros(8, 8, 8, dtype=torch.float64)           # pi(g) e_h = e_{gh}
for i, g in enumerate(Q8):
    for j, h in enumerate(Q8):
        reg[i, pos[key(qmul(g, h))], j] = 1.0
D_Q8 = dc.from_group(reg)
rs_Q8 = D_Q8.real_structure()
assert rs_Q8.ok and sorted(zip(rs_Q8.types, D_Q8.k_factors, D_Q8.multiplicities)) == [('H', 2, 2)] + [('R', 1, 1)] * 4
assert D_Q8.dims[0] == 8                                  # dim_R A_R = dim_C A = |Q_8| (Theorem 4.21)
```

</details>

**Siguiente:** [Tutorial 5](05_random_matrices.md): matrices aleatorias en un álgebra (unitarias de Haar, cálculo de Weingarten y estadística de niveles) en la biblioteca.
