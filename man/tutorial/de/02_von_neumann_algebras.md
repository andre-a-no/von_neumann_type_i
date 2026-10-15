[English](../02_von_neumann_algebras.md) · [Español](../es/02_von_neumann_algebras.md) · [Français](../fr/02_von_neumann_algebras.md) · **Deutsch** · [中文](../zh/02_von_neumann_algebras.md) · [日本語](../ja/02_von_neumann_algebras.md) · [Русский](../ru/02_von_neumann_algebras.md)

# Tutorial 2. Endlichdimensionale von-Neumann-Algebren

Dieses Tutorial begleitet [Kapitel 2](../../lecture_notes/de/02_von_neumann_algebras.md) der Vorlesungsnotizen, das beweist,
dass jede von-Neumann-Algebra auf $`\mathbb C^n`$ bis auf einen unitären Basiswechsel gleich
$`\bigoplus_cM_{k_c}(\mathbb C)\otimes1_{m_c}`$ ist ([Satz 2.6](../../lecture_notes/de/02_von_neumann_algebras.md#r-2-6)).
Wir berechnen Kommutanten numerisch am Beispiel einer Symmetriegruppe und sehen, wie die Bibliothek
eine solche Algebra, ihr Zentrum und die bedingte Erwartung auf das Zentrum speichert.

<a id="s-2-4"></a>

## 2.4 Beispiel: eine Symmetriegruppe und ihre Kommutante

*Vorlesungsnotizen:* [Abschnitt 2.4](../../lecture_notes/de/02_von_neumann_algebras.md#s-2-4); [Satz 2.2](../../lecture_notes/de/02_von_neumann_algebras.md#r-2-2), [Satz 2.6](../../lecture_notes/de/02_von_neumann_algebras.md#r-2-6), [Proposition 2.7](../../lecture_notes/de/02_von_neumann_algebras.md#r-2-7).

[Abschnitt 2.4](../../lecture_notes/de/02_von_neumann_algebras.md#s-2-4) der Vorlesungsnotizen zerlegt die Darstellung $`\pi`$ von
$`S_3`$, die die drei Faktoren von $`(\mathbb C^2)^{\otimes3}`$ permutiert: $`\pi(S_3)'\cong M_4\oplus M_2`$ hat die Dimension
$`20`$ und $`\pi(S_3)''\cong M_1\oplus M_2`$ die Dimension $`5`$. Der Code unten bestätigt die beiden Dimensionen numerisch
sowie, dass die Kommutante von $`\pi(S_3)'`$ wieder $`\pi(S_3)''`$ ist
([Satz 2.2](../../lecture_notes/de/02_von_neumann_algebras.md#r-2-2)): Die Kommutante einer Menge $`S`$ ist der
Kern der linearen Abbildung $`X\mapsto(AX-XA)_{A\in S}`$, dessen Dimension hier in Gleitkommaarithmetik
berechnet wird, indem die Singulärwerte oberhalb einer Toleranz gezählt werden. Die Toleranz muss die Singulärwerte, die in
exakter Arithmetik verschwinden, von den übrigen trennen; für verrauschte Matrizen sagt [Satz 4.25](../../lecture_notes/de/04_symmetries.md#r-4-25),
wo sie liegen darf ([Tutorial 4, Abschnitt 4.9](04_symmetries.md#s-4-9)).

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

Die Vielfachheiten sind im Spektrum sichtbar ([Vorlesungsnotizen, Abschnitt 2.4](../../lecture_notes/de/02_von_neumann_algebras.md#s-2-4), nach
[Satz 2.6](../../lecture_notes/de/02_von_neumann_algebras.md#r-2-6)):
Jedes selbstadjungierte $`K\in\pi(S_3)'`$ hat auf dem Komplement des symmetrischen Unterraums $`F`$ Eigenwerte
gerader Vielfachheit, und für $`K=E(Z+Z^*)`$, mit einer Ginibre-Matrix $`Z`$ und der Gruppenmittelung $`E`$, hat $`K`$ mit
Wahrscheinlichkeit eins $`4`$ einfache Eigenwerte auf $`F`$ und $`2`$ zweifach entartete auf dessen Komplement. Der Code überprüft die
entarteten Paare (die für jedes $`K`$ gelten) und, für das gezogene $`K`$, dass die beiden Paare und die vier
Eigenwerte auf $`F`$ verschieden sind; dabei verwendet er die Orthogonalprojektion $`P_F=\frac1{|G|}\sum_g\pi(g)`$ auf $`F`$.

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

## 2.5 Die Algebra in der Bibliothek

*Vorlesungsnotizen:* [Abschnitt 2.3](../../lecture_notes/de/02_von_neumann_algebras.md#s-2-3), [Satz 2.6](../../lecture_notes/de/02_von_neumann_algebras.md#r-2-6); das Zentrum: [Abschnitt 2.2](../../lecture_notes/de/02_von_neumann_algebras.md#s-2-2), [Lemma 2.3](../../lecture_notes/de/02_von_neumann_algebras.md#r-2-3).

Die Bibliothek speichert die abstrakte Algebra $`\bigoplus_cM_{k_c}(\mathbb C)`$ aus Satz 2.6
auf $`\bigoplus_c\mathbb C^{k_c}`$, ohne die Vielfachheiten: Ein Operator ist die Liste seiner
Blöcke $`A_c`$, aufgefüllt auf eine gemeinsame Größe. Die Kommutante aus dem Beispiel ist
`TypeIAlgebra([4, 2], [4, 2])`.

```python
from torch_vn_algebra import TypeIAlgebra, channels

alg = TypeIAlgebra([4, 2], [4, 2], precision='double')    # M_4 (+) M_2
print(alg.k_factors, alg.C)

A = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
print(A.matrix.shape)              # (1, 2, 4, 4): block 1 occupies the top-left 2 x 2 corner
assert A.matrix[0, 1, 2:, :].abs().max() == 0 and A.matrix[0, 1, :, 2:].abs().max() == 0
```

Das Zentrum besteht aus den Elementen $`\sum_cz_cP_c`$; `alg.central` konstruiert sie, und
sie kommutieren mit allem:

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

## 2.6 Bedingte Erwartungen

*Vorlesungsnotizen:* [Abschnitt 2.6](../../lecture_notes/de/02_von_neumann_algebras.md#s-2-6).

Nach [Vorlesungsnotizen, Proposition 2.8](../../lecture_notes/de/02_von_neumann_algebras.md#r-2-8) ist $`E(A)=\sum_c\frac{\mathrm{Tr}\,A_c}{k_c}\,P_c`$
eine bedingte Erwartung von $`\bigoplus_cM_{k_c}`$ auf sein Zentrum: Sie ist positiv, auf dem Zentrum die Identität,
ein Bimodulhomomorphismus über dem Zentrum, und sie erhält jede Spur $`\sum_cw_c\mathrm{Tr}\,A_c`$.

In der Bibliothek ist diese Abbildung der Kanal `channels.center_expectation(alg)`. Wendet
man ihn auf Operatoren an, erhält man den klassischen (über die Sektoren gemittelten) Teil:

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

## 2.7 Womit man vergleichen sollte

Die Beispiele aus Abschnitt 2.4 lassen sich mit jedem Paket für lineare Algebra
überprüfen; Software für Darstellungstheorie (GAP, SageMath) liefert die Zerlegung
$`\pi\cong\bigoplus\rho_\lambda\otimes1_{m_\lambda}`$ symbolisch. Programme zur exakten
Diagonalisierung in der Physik (QuSpin, die Symmetriemodule von QuTiP) nutzen dieselbe
Blockstruktur, um Hamilton-Operatoren Sektor für Sektor zu diagonalisieren. Die Bibliothek
nimmt die Struktur als Ausgangspunkt: Sind die Blöcke einmal bekannt, so bleiben alle
weiteren Rechnungen – Funktionen, Spuren, Zustände, Kanäle, Dynamik – innerhalb der Algebra.

<a id="x-2"></a>

## Aufgaben

[Aufgabe 2](../../lecture_notes/de/02_von_neumann_algebras.md#x-2) der Vorlesungsnotizen: Für die zyklische Verschiebung $`S`$ auf $`\mathbb C^3`$,
$`Se_i=e_{i+1\bmod3}`$, sind $`\{S\}'`$ und $`\{S\}''`$ beide gleich $`\mathrm{span}\{1,S,S^2\}`$, den zirkulanten
Matrizen. Man überprüfe die Dimensionen numerisch.

<details>
<summary>Lösungen</summary>

```python
S = torch.roll(torch.eye(3, dtype=dt), 1, dims=0)
assert len(commutant_basis([S])) == 3
assert len(commutant_basis(commutant_basis([S]))) == 3
assert span_dim([torch.eye(3, dtype=dt), S, S @ S]) == 3
```

</details>

**Weiter:** [Tutorial 3. Spuren, Zustände und Entropie](03_traces_and_states.md)
