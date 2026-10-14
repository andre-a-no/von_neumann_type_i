**English** · [Español](es/02_von_neumann_algebras.md) · [Français](fr/02_von_neumann_algebras.md) · [Deutsch](de/02_von_neumann_algebras.md) · [中文](zh/02_von_neumann_algebras.md) · [日本語](ja/02_von_neumann_algebras.md) · [Русский](ru/02_von_neumann_algebras.md)

# Tutorial 2. Finite-dimensional von Neumann algebras

This tutorial accompanies [Chapter 2](../lecture_notes/02_von_neumann_algebras.md) of the lecture notes, which proves
that every von Neumann algebra on $\mathbb C^n$ is, up to a unitary change of basis,
$\bigoplus_cM_{k_c}(\mathbb C)\otimes1_{m_c}$ ([Theorem 2.6](../lecture_notes/02_von_neumann_algebras.md#r-2-6)).
We compute commutants numerically on the example of a symmetry group, and see how the library
stores such an algebra, its centre and the conditional expectation onto the centre.

<a id="s-2-4"></a>

## 2.4 Example: a symmetry group and its commutant

*Lecture notes:* [Section 2.4](../lecture_notes/02_von_neumann_algebras.md#s-2-4); [Theorem 2.2](../lecture_notes/02_von_neumann_algebras.md#r-2-2), [Theorem 2.6](../lecture_notes/02_von_neumann_algebras.md#r-2-6), [Proposition 2.7](../lecture_notes/02_von_neumann_algebras.md#r-2-7).

[Section 2.4](../lecture_notes/02_von_neumann_algebras.md#s-2-4) of the lecture notes decomposes the representation $\pi$ of
$S_3$ permuting the three factors of $(\mathbb C^2)^{\otimes3}$: $\pi(S_3)'\cong M_4\oplus M_2$ has dimension
$20$ and $\pi(S_3)''\cong M_1\oplus M_2$ has dimension $5$. The code below confirms the two dimensions
numerically, and that the commutant of $\pi(S_3)'$ is again $\pi(S_3)''$
([Theorem 2.2](../lecture_notes/02_von_neumann_algebras.md#r-2-2)): the commutant of a set $S$ is the null
space of the linear map $X\mapsto(AX-XA)_{A\in S}$, whose dimension is computed here in floating
point, by counting the singular values above a tolerance.

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

The multiplicities are visible in the spectrum ([lecture notes, Section 2.4](../lecture_notes/02_von_neumann_algebras.md#s-2-4), from
[Theorem 2.6](../lecture_notes/02_von_neumann_algebras.md#r-2-6)):
every self-adjoint $K\in\pi(S_3)'$ has eigenvalues of even multiplicity on the complement of the symmetric
subspace $F$, and for $K=E(Z+Z^*)$, with $Z$ a Ginibre matrix and $E$ the group average, with probability one
$K$ has $4$ simple eigenvalues on $F$ and $2$ doubly degenerate ones on its complement. The code checks the
degenerate pairs (which hold for every $K$) and, for the sampled $K$, that the two pairs and the four
eigenvalues on $F$ are distinct, using the orthogonal projection $P_F=\frac1{|G|}\sum_g\pi(g)$ onto $F$.

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

## 2.5 The algebra in the library

*Lecture notes:* [Section 2.3](../lecture_notes/02_von_neumann_algebras.md#s-2-3), [Theorem 2.6](../lecture_notes/02_von_neumann_algebras.md#r-2-6); the centre: [Section 2.2](../lecture_notes/02_von_neumann_algebras.md#s-2-2), [Lemma 2.3](../lecture_notes/02_von_neumann_algebras.md#r-2-3).

The library stores the abstract algebra $\bigoplus_cM_{k_c}(\mathbb C)$ of Theorem 2.6 on
$\bigoplus_c\mathbb C^{k_c}$, without the multiplicities: an operator is the list of its
blocks $A_c$, padded to a common size. The commutant of the example is
`TypeIAlgebra([4, 2], [4, 2])`.

```python
from torch_vn_algebra import TypeIAlgebra, channels

alg = TypeIAlgebra([4, 2], [4, 2], precision='double')    # M_4 (+) M_2
print(alg.k_factors, alg.C)

A = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
print(A.matrix.shape)              # (1, 2, 4, 4): block 1 occupies the top-left 2 x 2 corner
assert A.matrix[0, 1, 2:, :].abs().max() == 0 and A.matrix[0, 1, :, 2:].abs().max() == 0
```

The centre consists of the elements $\sum_cz_cP_c$; `alg.central` builds them, and
they commute with everything:

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

## 2.6 Conditional expectations

*Lecture notes:* [Section 2.6](../lecture_notes/02_von_neumann_algebras.md#s-2-6).

By [lecture notes, Proposition 2.8](../lecture_notes/02_von_neumann_algebras.md#r-2-8), $E(A)=\sum_c\frac{\mathrm{Tr}\,A_c}{k_c}\,P_c$
is a conditional expectation of $\bigoplus_cM_{k_c}$ onto its centre: it is positive, it is the identity on the centre,
it is a bimodule map over the centre, and it preserves every trace $\sum_cw_c\mathrm{Tr}\,A_c$.

In the library this map is the channel `channels.center_expectation(alg)`. Applying it
to operators gives the classical (sector-averaged) part:

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

## 2.7 What to compare with

The examples of Section 2.4 can be checked with any linear-algebra package;
representation-theory software (GAP, SageMath) gives the decomposition
$\pi\cong\bigoplus\rho_\lambda\otimes1_{m_\lambda}$ symbolically. Exact-diagonalisation
codes in physics (QuSpin, the symmetry modules of QuTiP) use the same block structure to
diagonalise Hamiltonians sector by sector. The library takes the structure as the
starting point: once the blocks are known, all further computations — functions,
traces, states, channels, dynamics — stay inside the algebra.

<a id="x-2"></a>

## Exercises

[Exercise 2](../lecture_notes/02_von_neumann_algebras.md#x-2) of the lecture notes: for the cyclic shift $S$ on $\mathbb C^3$,
$Se_i=e_{i+1\bmod3}$, both $\{S\}'$ and $\{S\}''$ equal $\mathrm{span}\{1,S,S^2\}$, the circulant
matrices. Check the dimensions numerically.

<details>
<summary>Solutions</summary>

```python
S = torch.roll(torch.eye(3, dtype=dt), 1, dims=0)
assert len(commutant_basis([S])) == 3
assert len(commutant_basis(commutant_basis([S]))) == 3
assert span_dim([torch.eye(3, dtype=dt), S, S @ S]) == 3
```

</details>

**Next:** [Tutorial 3. Traces, states and entropy](03_traces_and_states.md)
