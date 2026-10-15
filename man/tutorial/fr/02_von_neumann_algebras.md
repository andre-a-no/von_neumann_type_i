[English](../02_von_neumann_algebras.md) · [Español](../es/02_von_neumann_algebras.md) · **Français** · [Deutsch](../de/02_von_neumann_algebras.md) · [中文](../zh/02_von_neumann_algebras.md) · [日本語](../ja/02_von_neumann_algebras.md) · [Русский](../ru/02_von_neumann_algebras.md)

# Tutoriel 2. Algèbres de von Neumann de dimension finie

Ce tutoriel accompagne le [chapitre 2](../../lecture_notes/fr/02_von_neumann_algebras.md) des notes de cours, qui démontre
que toute algèbre de von Neumann sur $`\mathbb C^n`$ est, à un changement de base unitaire près,
$`\bigoplus_cM_{k_c}(\mathbb C)\otimes1_{m_c}`$ ([théorème 2.6](../../lecture_notes/fr/02_von_neumann_algebras.md#r-2-6)).
Nous calculons numériquement des commutants sur l'exemple d'un groupe de symétrie, et voyons comment la bibliothèque
stocke une telle algèbre, son centre et l'espérance conditionnelle sur le centre.

<a id="s-2-4"></a>

## 2.4 Exemple : un groupe de symétrie et son commutant

*Notes de cours :* [section 2.4](../../lecture_notes/fr/02_von_neumann_algebras.md#s-2-4) ; [théorème 2.2](../../lecture_notes/fr/02_von_neumann_algebras.md#r-2-2), [théorème 2.6](../../lecture_notes/fr/02_von_neumann_algebras.md#r-2-6), [proposition 2.7](../../lecture_notes/fr/02_von_neumann_algebras.md#r-2-7).

La [section 2.4](../../lecture_notes/fr/02_von_neumann_algebras.md#s-2-4) des notes de cours décompose la représentation $`\pi`$ de
$`S_3`$ par permutation des trois facteurs de $`(\mathbb C^2)^{\otimes3}`$ : $`\pi(S_3)'\cong M_4\oplus M_2`$ est de dimension
$`20`$ et $`\pi(S_3)''\cong M_1\oplus M_2`$ est de dimension $`5`$. Le code ci-dessous confirme numériquement les deux dimensions,
ainsi que le fait que le commutant de $`\pi(S_3)'`$ est à nouveau $`\pi(S_3)''`$
([théorème 2.2](../../lecture_notes/fr/02_von_neumann_algebras.md#r-2-2)) : le commutant d'une partie $`S`$ est le
noyau de l'application linéaire $`X\mapsto(AX-XA)_{A\in S}`$, dont la dimension est calculée ici en virgule
flottante, en comptant les valeurs singulières supérieures à une tolérance. La tolérance doit séparer les valeurs singulières qui
s'annulent en arithmétique exacte des autres ; pour des matrices bruitées, le [théorème 4.25](../../lecture_notes/fr/04_symmetries.md#r-4-25)
dit où elle peut se trouver ([tutoriel 4, section 4.9](04_symmetries.md#s-4-9)).

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

Les multiplicités se lisent sur le spectre ([notes de cours, section 2.4](../../lecture_notes/fr/02_von_neumann_algebras.md#s-2-4), d'après le
[théorème 2.6](../../lecture_notes/fr/02_von_neumann_algebras.md#r-2-6)) :
tout $`K\in\pi(S_3)'`$ autoadjoint a des valeurs propres de multiplicité paire sur le supplémentaire du sous-espace
symétrique $`F`$, et pour $`K=E(Z+Z^*)`$, où $`Z`$ est une matrice de Ginibre et $`E`$ la moyenne sur le groupe, presque sûrement
$`K`$ a $`4`$ valeurs propres simples sur $`F`$ et $`2`$ valeurs propres doubles sur son supplémentaire. Le code vérifie les
paires dégénérées (valables pour tout $`K`$) et, pour le $`K`$ tiré, que les deux paires et les quatre
valeurs propres sur $`F`$ sont distinctes, en utilisant le projecteur orthogonal $`P_F=\frac1{|G|}\sum_g\pi(g)`$ sur $`F`$.

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

## 2.5 L'algèbre dans la bibliothèque

*Notes de cours :* [section 2.3](../../lecture_notes/fr/02_von_neumann_algebras.md#s-2-3), [théorème 2.6](../../lecture_notes/fr/02_von_neumann_algebras.md#r-2-6) ; le centre : [section 2.2](../../lecture_notes/fr/02_von_neumann_algebras.md#s-2-2), [lemme 2.3](../../lecture_notes/fr/02_von_neumann_algebras.md#r-2-3).

La bibliothèque stocke l'algèbre abstraite $`\bigoplus_cM_{k_c}(\mathbb C)`$ du théorème 2.6
sur $`\bigoplus_c\mathbb C^{k_c}`$, sans les multiplicités : un opérateur est la liste de ses
blocs $`A_c`$, complétés par des zéros jusqu'à une taille commune. Le commutant de l'exemple
est `TypeIAlgebra([4, 2], [4, 2])`.

```python
from torch_vn_algebra import TypeIAlgebra, channels

alg = TypeIAlgebra([4, 2], [4, 2], precision='double')    # M_4 (+) M_2
print(alg.k_factors, alg.C)

A = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
print(A.matrix.shape)              # (1, 2, 4, 4): block 1 occupies the top-left 2 x 2 corner
assert A.matrix[0, 1, 2:, :].abs().max() == 0 and A.matrix[0, 1, :, 2:].abs().max() == 0
```

Le centre est formé des éléments $`\sum_cz_cP_c`$ ; `alg.central` les construit, et ils
commutent avec tout :

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

## 2.6 Espérances conditionnelles

*Notes de cours :* [section 2.6](../../lecture_notes/fr/02_von_neumann_algebras.md#s-2-6).

D'après [notes de cours, proposition 2.8](../../lecture_notes/fr/02_von_neumann_algebras.md#r-2-8), $`E(A)=\sum_c\frac{\mathrm{Tr}\,A_c}{k_c}\,P_c`$
est une espérance conditionnelle de $`\bigoplus_cM_{k_c}`$ sur son centre : elle est positive, elle est l'identité sur le centre,
c'est un morphisme de bimodules sur le centre, et elle préserve toute trace $`\sum_cw_c\mathrm{Tr}\,A_c`$.

Dans la bibliothèque, cette application est le canal `channels.center_expectation(alg)`.
Appliquée aux opérateurs, elle en donne la partie classique (moyennée sur chaque secteur) :

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

## 2.7 Points de comparaison

Les exemples de la section 2.4 peuvent être vérifiés avec n'importe quelle bibliothèque
d'algèbre linéaire ; les logiciels de théorie des représentations (GAP, SageMath)
fournissent la décomposition $`\pi\cong\bigoplus\rho_\lambda\otimes1_{m_\lambda}`$ de façon
symbolique. Les codes de diagonalisation exacte utilisés en physique (QuSpin, les modules
de symétrie de QuTiP) exploitent la même structure par blocs pour diagonaliser les
hamiltoniens secteur par secteur. La bibliothèque prend cette structure comme point de
départ : une fois les blocs connus, tous les calculs ultérieurs — fonctions, traces,
états, canaux, dynamique — restent à l'intérieur de l'algèbre.

<a id="x-2"></a>

## Exercices

[Exercice 2](../../lecture_notes/fr/02_von_neumann_algebras.md#x-2) des notes de cours : pour la permutation circulaire $`S`$ de $`\mathbb C^3`$,
$`Se_i=e_{i+1\bmod3}`$, $`\{S\}'`$ et $`\{S\}''`$ sont tous deux égaux à $`\mathrm{span}\{1,S,S^2\}`$, les matrices
circulantes. Vérifier numériquement les dimensions.

<details>
<summary>Solutions</summary>

```python
S = torch.roll(torch.eye(3, dtype=dt), 1, dims=0)
assert len(commutant_basis([S])) == 3
assert len(commutant_basis(commutant_basis([S]))) == 3
assert span_dim([torch.eye(3, dtype=dt), S, S @ S]) == 3
```

</details>

**Suite :** [Tutoriel 3. Traces, états et entropie](03_traces_and_states.md)
