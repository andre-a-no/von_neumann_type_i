[English](../07_tensor_products.md) · [Español](../es/07_tensor_products.md) · **Français** · [Deutsch](../de/07_tensor_products.md) · [中文](../zh/07_tensor_products.md) · [日本語](../ja/07_tensor_products.md) · [Русский](../ru/07_tensor_products.md)

# Tutoriel 7. Produits tensoriels et intrication

Ce tutoriel accompagne le [chapitre 7](../../lecture_notes/fr/07_tensor_products.md) des notes de cours, qui établit la structure du produit
tensoriel simple et du produit tensoriel fusionné de deux algèbres à secteurs, définit la trace partielle
comme la duale d'un plongement et étudie l'intrication et les états aléatoires. Nous construisons ici ces
produits avec la bibliothèque, calculons des traces partielles et des entropies d'intrication, et
vérifions numériquement les théorèmes.

Dans tout ce chapitre, $`\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)`$ agit sur
$`H_1=\bigoplus_c\mathbb C^{k_c}`$, avec les projecteurs centraux minimaux $`P_c`$, et
$`\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)`$ agit sur $`H_2=\bigoplus_d\mathbb C^{l_d}`$, avec
les projecteurs centraux minimaux $`Q_d`$. $`\mathrm{Tr}`$ désigne la trace de l'espace sur lequel
agit un opérateur ; sur $`\mathcal M_i`$, c'est la trace brute du chapitre 3. Pour $`X`$ sur
$`\mathbb C^k\otimes\mathbb C^l`$, on note $`X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle`$,
et l'on identifie $`\mathbb C^k\otimes\mathbb C^l`$ à $`\mathbb C^{kl}`$ par $`e_a\otimes e_b\mapsto e_{(a-1)l+b}`$
(la convention de Kronecker), de sorte que $`(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}`$.

<a id="s-7-1"></a>

## 7.1 Le produit tensoriel de deux algèbres

*Notes de cours :* [section 7.1](../../lecture_notes/fr/07_tensor_products.md#s-7-1) : [théorème 7.1](../../lecture_notes/fr/07_tensor_products.md#r-7-1) (la structure de $`\mathcal M_1\otimes\mathcal M_2`$).

`tensor_product(alg1, alg2)` construit cette algèbre. Le secteur $`(c,d)`$ a pour indice $`s=cC_2+d`$
(indices à partir de $`0`$ dans le code) et pour charge $`(q_c,r_d)`$, et, à l'intérieur de ce
secteur, $`e_a\otimes e_b`$ a pour indice $`a\,l_d+b`$. `kron(A, B, alg12)` calcule $`A\otimes B`$. L'ordre
des secteurs diffère de l'ordre de Kronecker de $`H_1\otimes H_2`$ par une permutation de la base, que
`kron_rows` construit explicitement ci-dessous ; elle permet de comparer avec des produits de
Kronecker de matrices pleines.

```python
import math
import numpy as np
import torch
from torch_vn_algebra import (TypeIAlgebra, DensityMatrix, states, composite,
                              tensor_product, fused_tensor_product, kron, partial_trace)

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def full(X):
    """An operator of (+)_c M_{k_c} as a block-diagonal matrix on (+)_c C^{k_c}, shape (batch, D, D)."""
    return torch.stack([torch.block_diag(*[x[c, :k, :k] for c, k in enumerate(X.algebra.k_factors)])
                        for x in X.matrix])

def kron_rows(alg12):
    """For every sector of a product algebra, the positions of its basis vectors in H1 (x) H2 (Kronecker order)."""
    k1, k2 = alg12.tensor_factors[0].k_factors, alg12.tensor_factors[1].k_factors
    return [torch.tensor([(sum(k1[:c]) + a) * sum(k2) + sum(k2[:d]) + b
                          for c, d, off in pairs for a in range(k1[c]) for b in range(k2[d])])
            for pairs in alg12.product_blocks]              # pair blocks (c, d, offset), by offset

def full12(X):
    """An operator of a (plain or fused) product algebra as a matrix on H1 (x) H2 in Kronecker order."""
    f1, f2 = X.algebra.tensor_factors
    N = sum(f1.k_factors) * sum(f2.k_factors)
    out = torch.zeros(X.matrix.shape[0], N, N, dtype=X.matrix.dtype)
    for s, r in enumerate(kron_rows(X.algebra)):
        out[:, r[:, None], r[None, :]] = X.matrix[:, s, :len(r), :len(r)]
    return out

def random_operator(alg):
    return alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])

alg1 = TypeIAlgebra([1, 2], [1, 2], precision='double')            # C (+) M_2
alg2 = TypeIAlgebra([2, 3], [2, 3], precision='double')            # M_2 (+) M_3
alg12 = tensor_product(alg1, alg2)
print(alg12.k_factors, alg12.charges, sum(k * k for k in alg12.k_factors))     # 4 sectors, dimension 65
assert sum(k * k for k in alg12.k_factors) == 5 * 13                           # Theorem 7.1(b)

A, A2, B, B2 = random_operator(alg1), random_operator(alg1), random_operator(alg2), random_operator(alg2)
AB = kron(A, B, alg12)
assert close(full12(AB)[0], torch.kron(full(A)[0], full(B)[0]))               # Theorem 7.1(a)
assert close(AB @ kron(A2, B2, alg12), kron(A @ A2, B @ B2, alg12))
for tr in ('Tr_blunt', 'Tr_norm', 'tau_vN'):                                   # Theorem 7.1(c)
    assert (getattr(AB, tr)() - getattr(A, tr)() * getattr(B, tr)()).abs().max() < 1e-10
```

<a id="s-7-2"></a>

## 7.2 La trace partielle

*Notes de cours :* [section 7.2](../../lecture_notes/fr/07_tensor_products.md#s-7-2) : la définition de $`\mathrm{Tr}_2`$, [proposition 7.2](../../lecture_notes/fr/07_tensor_products.md#r-7-2), [théorème 7.3](../../lecture_notes/fr/07_tensor_products.md#r-7-3) (la trace partielle est un canal) et [proposition 7.4](../../lecture_notes/fr/07_tensor_products.md#r-7-4) (une espérance conditionnelle d'indice $`m^2`$).

`partial_trace(X, keep=1)` remet chaque bloc $`\rho_{cd}`$ sous la forme $`(k_c,l_d,k_c,l_d)`$ et prend
la trace ; `composite.partial_trace_channel(alg12, keep=1)` est le canal du théorème 7.3, avec les
opérateurs de Kraus $`1_{k_c}\otimes\langle j|`$ du secteur $`(c,d)`$ vers le secteur $`c`$. On vérifie les
deux par rapport à la définition, à un `numpy.einsum` indépendant et à la dualité :

```python
rho = states.random_density_matrix(alg12, batch_size=3)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
assert isinstance(rho1, DensityMatrix)                                      # states give states
R = full12(rho).numpy().reshape(3, 3, 5, 3, 5)                              # indices (batch, i, j, i', j')
assert np.abs(full(rho1).numpy() - np.einsum('xijkj->xik', R)).max() < 1e-12
assert np.abs(full(rho2).numpy() - np.einsum('xijil->xjl', R)).max() < 1e-12
one2 = alg2.identity()
assert (states.expectation(rho, kron(A, one2, alg12)) - states.expectation(rho1, A)).abs().max() < 1e-12
assert close(partial_trace(AB, keep=1), A * B.trace.item())                 # Proposition 7.2(b)

Ptr = composite.partial_trace_channel(alg12, keep=1)
print(Ptr)
assert Ptr.is_trace_preserving()
assert close(Ptr(rho), rho1)
assert close(Ptr.adjoint('blunt')(A), kron(A, one2, alg12))                 # Theorem 7.3: the dual
```

Le module `torch_vn_algebra.inclusions` implémente les inclusions, leurs espérances conditionnelles et leurs
indices ; `from_tensor_product(alg1, m)` construit l'inclusion $`\mathcal M_1\otimes1\subset\mathcal M_1\otimes M_m`$ de la
[proposition 7.4](../../lecture_notes/fr/07_tensor_products.md#r-7-4), de matrice d'inclusion $`\Lambda=m1_{C_1}`$ et d'indice $`m^2`$
([notes de cours, section 7.2](../../lecture_notes/fr/07_tensor_products.md#s-7-2)) :

```python
from torch_vn_algebra import inclusions
m = 3
inc = inclusions.from_tensor_product(alg1, m)                # alg1 (x) 1 c alg1 (x) M_3
X = random_operator(inc.M)
for trace in ('Tr_blunt', 'Tr_norm', 'tau_vN'):
    E = inc.conditional_expectation(trace.replace('Tr_', ''))
    assert close(E(X), partial_trace(X, keep=1) * (1 / m))   # Proposition 7.4(a): the same E for every trace
    assert abs(getattr(inc.embed(E(X)), trace)() - getattr(X, trace)()).max() < 1e-10
    assert close(E.index(), inc.M.identity() * m**2)         # Proposition 7.4(c)

def unit(a, b):
    e = torch.zeros(m, m, dtype=dt)
    e[a, b] = 1
    return math.sqrt(m) * kron(alg1.identity(), inc.M.tensor_factors[1].from_blocks([e]), inc.M)

u = [unit(a, b) for a in range(m) for b in range(m)]         # the quasi-basis of Proposition 7.4(b)
assert close(sum((v @ inc.embed(E(v.adjoint() @ X))).matrix for v in u), X)
assert close(sum((v @ v.adjoint()).matrix for v in u), inc.M.identity() * m**2)
print(inc.matrix.tolist(), inc.component_indices())           # Lambda = 3 * 1_2, index 9 per component
```

<a id="s-7-3"></a>

## 7.3 Fusion des charges

*Notes de cours :* [section 7.3](../../lecture_notes/fr/07_tensor_products.md#s-7-3) : [proposition 7.5](../../lecture_notes/fr/07_tensor_products.md#r-7-5) (le produit fusionné) et [théorème 7.6](../../lecture_notes/fr/07_tensor_products.md#r-7-6) (l'état réduit d'un produit fusionné).

`fused_tensor_product(alg1, alg2, fuse)` (par défaut, `fuse` est l'addition) ordonne ses secteurs
selon la charge fusionnée, et les couples $`(c,d)`$ à l'intérieur d'un secteur dans l'ordre
lexicographique, chacun avec la base $`e_a\otimes e_b`$. Exemple : les bosons de cœur dur
([notes de cours, section 7.3](../../lecture_notes/fr/07_tensor_products.md#s-7-3)). Les opérateurs qui conservent le nombre forment $`M_1\oplus M_2\oplus M_1`$
(charges $`0,1,2`$) sur deux sites et $`M_1\oplus M_3\oplus M_3\oplus M_1`$ sur trois sites, et, d'après la
[proposition 7.5(c)](../../lecture_notes/fr/07_tensor_products.md#r-7-5), leur produit fusionné est $`\{N\}'`$ pour le nombre d'occupation total de cinq sites,
avec $`n_q=\binom5q`$.

```python
bos1 = TypeIAlgebra([1, 2, 1], [1, 2, 1], precision='double', charges=[0, 1, 2])
bos2 = TypeIAlgebra([1, 3, 3, 1], [1, 3, 3, 1], precision='double', charges=[0, 1, 2, 3])
F = fused_tensor_product(bos1, bos2)                         # U(1): total particle number
P = tensor_product(bos1, bos2)                               # particle numbers conserved separately
print(F.charges, F.k_factors, sum(k * k for k in F.k_factors), sum(k * k for k in P.k_factors))
assert F.k_factors == [math.comb(5, q) for q in range(6)]
assert sum(k * k for k in F.k_factors) == 252 and sum(k * k for k in P.k_factors) == 120

N1, N2 = full(bos1.central([0.0, 1.0, 2.0]))[0], full(bos2.central([0.0, 1.0, 2.0, 3.0]))[0]
N = torch.kron(N1, torch.eye(8, dtype=dt)) + torch.kron(torch.eye(4, dtype=dt), N2)
X = full12(F.operator(torch.randn(1, F.C, F.k_max, F.k_max, dtype=dt)))[0]
assert close(X @ N, N @ X)                                   # Proposition 7.5(c)
A, B = random_operator(bos1), random_operator(bos2)
assert close(full12(kron(A, B, F)), full12(kron(A, B, P)))   # Proposition 7.5(b)
```

L'algèbre fusionnée est de dimension $`\sum_q\binom5q^2=\binom{10}5=252`$ (encore Vandermonde), le
produit simple de dimension $`(1+4+1)(1+9+9+1)=120`$ (théorème 7.1(b)). D'après le théorème 7.1(a) et
la proposition 7.5(b), le produit simple est formé des éléments du produit fusionné qui sont
diagonaux par blocs pour les couples $`(c,d)`$ ; les $`132`$ dimensions restantes sont engendrées par
les opérateurs qui envoient un bloc de couple $`(c,d)`$ dans un autre bloc de couple $`(c',d')`$ du même
$`H_q`$, $`c'\ne c`$ (donc $`d'\ne d`$) : ils changent la charge de chaque partie en conservant la charge
totale, c'est-à-dire qu'ici ils déplacent des particules d'une partie à l'autre.

D'après la formule par blocs du théorème 7.6, $`(\rho_1)_c=\sum_d\mathrm{tr}_{l_d}\,\rho_{(c,d),(c,d)}`$. Les
deux implémentations de la bibliothèque sont conçues pour calculer cette application pour toute
règle de fusion : `partial_trace` ne prend la trace que des blocs diagonaux de couples
$`\rho_{(c,d),(c,d)}`$, et `partial_trace_channel` utilise des indices de Kraus distincts pour les
couples d'un même secteur ayant le même $`c`$, de sorte que sa duale est de nouveau
$`A\mapsto A\otimes1`$. Le code vérifie les deux, pour les deux règles de fusion ci-dessous et les
deux choix du facteur conservé, contre $`E_1(\mathrm{Tr}_{H_2}\,\rho)`$ calculé indépendamment avec
`einsum`, et la duale contre $`A\otimes1`$. Avec
l'addition, $`E_1`$ n'a aucun effet dans notre exemple ; avec la règle de parité
$`f(q,r)=q+r\bmod2`$, les charges $`0`$ et $`2`$ de `bos1` fusionnent de la même manière avec tout $`r`$, et
$`E_1`$ supprime de véritables cohérences :

```python
G = fused_tensor_product(bos1, bos2, fuse=lambda q, r: (q + r) % 2)       # only the parity is conserved
print(G.k_factors)

def block_mask(alg):
    """Ones in the diagonal blocks of (+)_c C^{k_c}: X * block_mask(alg) = sum_c P_c X P_c."""
    return full(alg.operator(torch.ones(1, alg.C, alg.k_max, alg.k_max, dtype=dt)))

for alg, coherent in ((F, False), (G, True)):
    rho = states.random_density_matrix(alg, batch_size=2)
    R = full12(rho).reshape(2, 4, 8, 4, 8)
    for keep, X_H, alg_k in ((1, torch.einsum('xijkj->xik', R), bos1), (2, torch.einsum('xijil->xjl', R), bos2)):
        E_X = X_H * block_mask(alg_k)                                     # E(Tr_H rho), Theorem 7.6
        assert close(full(partial_trace(rho, keep=keep)), E_X)
        assert close(full(composite.partial_trace_channel(alg, keep=keep)(rho)), E_X)
        lost = (X_H - E_X).abs().max().item()                             # coherences removed by E
        print(f"keep={keep}: largest coherence removed: {lost:.3e}")
        assert (lost > 1e-3) if coherent else (lost == 0)
    assert close(composite.partial_trace_channel(alg, keep=1).adjoint('blunt')(A), kron(A, bos2.identity(), alg))
```

<a id="s-7-4"></a>

## 7.4 Intrication des états purs

*Notes de cours :* [section 7.4](../../lecture_notes/fr/07_tensor_products.md#s-7-4) : [lemme 7.7](../../lecture_notes/fr/07_tensor_products.md#r-7-7) (décomposition de Schmidt) et [proposition 7.8](../../lecture_notes/fr/07_tensor_products.md#r-7-8) (intrication résolue en symétrie) : $`E(\psi)=H(p)+\sum_cp_cE(\hat\psi_c)`$, où $`p_c`$ est la probabilité du secteur $`c`$ de la première partie.

$`H(p)`$ est l'**intrication de nombre** et $`\sum_cp_cE(\hat\psi_c)`$ l'**intrication
configurationnelle**. Pour un état aléatoire de $`2`$ particules sur les cinq sites, on calcule
$`E(\psi)`$ à partir des traces partielles des deux côtés et des valeurs singulières de la matrice
$`\Psi`$ de taille $`4\times8`$ :

```python
s = F.charges.index(2)                                      # the sector of two particles, n_q = 10
n = F.k_factors[s]
psi = torch.randn(n, dtype=dt)
psi = psi / psi.norm()
mat = torch.zeros(1, F.C, F.k_max, F.k_max, dtype=dt)
mat[0, s, :n, :n] = torch.outer(psi, psi.conj())
rho = DensityMatrix(F, matrix=mat)
r1, r2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)

psi_H = torch.zeros(4 * 8, dtype=dt)
psi_H[kron_rows(F)[s]] = psi                                # psi as a vector of H1 (x) H2
lam = torch.linalg.svdvals(psi_H.reshape(4, 8)) ** 2        # squared Schmidt coefficients
lam = lam[lam > 1e-14]
E_psi = -(lam * lam.log()).sum().item()
assert abs(r1.von_neumann_entropy().item() - E_psi) < 1e-10
assert abs(r2.von_neumann_entropy().item() - E_psi) < 1e-10  # Lemma 7.7

def entropy(sigma):
    w = torch.linalg.eigvalsh(sigma).clamp(min=1e-300)
    return -(w * w.log()).sum().item()

p = r1.sector_probabilities()[0]                            # c = 0, 1, 2 particles in the first part
assert torch.allclose(r2.sector_probabilities()[0, [2, 1, 0]], p)   # d(c) = 2 - c
H_num = -(p * p.log()).sum().item()
S_conf = sum(p[c].item() * entropy(r1.matrix[0, c, :k, :k] / p[c]) for c, k in enumerate(bos1.k_factors))
print(f"E = {E_psi:.4f} = H(p) {H_num:.4f} + configurational {S_conf:.4f}")
assert abs(E_psi - H_num - S_conf) < 1e-10                  # Proposition 7.8
```

<a id="s-7-5"></a>

## 7.5 États aléatoires et formule de Page

*Notes de cours :* [section 7.5](../../lecture_notes/fr/07_tensor_products.md#s-7-5) : [théorème 7.11](../../lecture_notes/fr/07_tensor_products.md#r-7-11) (Lubkin), [corollaire 7.12](../../lecture_notes/fr/07_tensor_products.md#r-7-12), [théorème 7.13](../../lecture_notes/fr/07_tensor_products.md#r-7-13) (formule de Page) et [proposition 7.14](../../lecture_notes/fr/07_tensor_products.md#r-7-14) (la mesure induite).

On vérifie les théorèmes 7.11 et 7.13 sur $`10^4`$ états purs de $`\mathbb C^3\otimes\mathbb C^6`$ réduits
avec `partial_trace` et sur la mesure induite, ainsi que la proposition 7.14 pour
$`M_1\oplus M_2\oplus M_3`$, $`r=2`$ :

```python
a, b, n_samples = 3, 6, 10_000
g = torch.randn(n_samples, a * b, dtype=dt)                 # complex Gaussian: E|g_i|^2 = 1
psi = g / g.norm(dim=-1, keepdim=True)
AB_alg = tensor_product(TypeIAlgebra([a], [a], precision='double'), TypeIAlgebra([b], [b], precision='double'))
pure = DensityMatrix(AB_alg, matrix=(psi[:, :, None] * psi.conj()[:, None, :])[:, None],
                     validate=False)                        # pure by construction
rho_A = partial_trace(pure, keep=1)
purity_exact = (a + b) / (a * b + 1)
page = sum(1 / k for k in range(b + 1, a * b + 1)) - (a - 1) / (2 * b)

def mean_se(x):
    return x.mean().item(), x.std().item() / math.sqrt(x.numel())

for name, rho_s in (('reduced pure states', rho_A),
                    ('induced measure', states.random_density_matrix(rho_A.algebra, 40_000, rank=b))):
    (m_p, e_p), (m_S, e_S) = mean_se(rho_s.purity()), mean_se(rho_s.von_neumann_entropy())
    print(f"{name}: E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {purity_exact:.4f}), "
          f"E S = {m_S:.4f} +- {e_S:.4f} (Page {page:.4f})")
    assert abs(m_p - purity_exact) < 5 * e_p and abs(m_S - page) < 5 * e_S     # Theorems 7.11, 7.13
    assert math.log((a * b + 1) / (a + b)) < m_S < math.log(a)                # Corollary 7.12
assert math.log(a) - (a - 1) / (a * b + 1) - (a - 1) / (2 * b) <= page <= math.log(a) - (a - 1) / (2 * b)
assert abs(sum(1 / k for k in (3, 4)) - 1 / 4 - 1 / 3) < 1e-15                # two qubits: E S = 1/3

alg3 = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
m_p, e_p = mean_se(states.random_density_matrix(alg3, 40_000, rank=2).purity())
print(f"E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {26 / 78:.4f})")
assert abs(m_p - 26 / 78) < 5 * e_p                                          # Proposition 7.14
```

<a id="s-7-6"></a>

## 7.6 Points de comparaison

Sur un seul espace de Hilbert, QuTiP (`Qobj.ptrace`, `entropy_vn`, `rand_ket`, `rand_dm` avec la
distribution de Ginibre et un rang) et `numpy` (un `reshape` et un `einsum`, comme ci-dessus)
calculent traces partielles, entropies et états induits pour les produits simples (les
vérifications par `einsum` ci-dessus comparent un tel calcul avec la bibliothèque). Les bibliothèques de réseaux de tenseurs à charges conservées — ITensor
(tenseurs creux par blocs indexés par des nombres quantiques), TeNPy (pattes chargées, spectres
d'intrication par charge) et, pour les symétries non abéliennes, TensorKit.jl (arbres de fusion) —
fusionnent les secteurs de la même manière pour de grands systèmes et donnent l'intrication
résolue en symétrie. La bibliothèque ajoute la couche algébrique : les produits sont des algèbres
de von Neumann avec secteurs et traces, la règle de fusion est une application quelconque des
étiquettes, la trace partielle est un canal avec sa duale pour chaque trace et la restriction
correcte pour les produits fusionnés (théorèmes 7.3 et 7.6), et tout cela fonctionne par lots. Elle
n'a pas de réseaux de tenseurs, et `fused_tensor_product` prend une règle avec un seul résultat par
couple de charges, comme pour les symétries abéliennes ; une règle non abélienne, comme pour
$`SU(2)`$, où deux spins $`j_1,j_2`$ se combinent en les différents spins totaux
$`|j_1-j_2|,\dots,j_1+j_2`$ (chapitre 4, exercice 2), n'est pas ici un produit fusionné (la bibliothèque traite $`SU(2)`$ séparément pour les chaînes de spins $`\frac12`$,
chapitre 4).

<a id="x-7"></a>

## Exercices

[Exercices](../../lecture_notes/fr/07_tensor_products.md#x-7) des notes de cours, numériquement :

1. Vérifier que la duale de $`\mathrm{Tr}_2`$ relativement à $`\tau_{\mathrm{vN}}`$,
   `partial_trace_channel(alg12).adjoint('tau_vN')`, est $`A\mapsto A\otimes Z`$ avec $`Z=\bigoplus_dC_2l_d1_{l_d}`$.
2. (Sous-additivité.) Vérifier sur des états aléatoires de `alg12` que l'information mutuelle
   $`S(\rho_1)+S(\rho_2)-S(\rho)`$ est positive et égale à $`S(\rho\|\rho_1\otimes\rho_2)`$.
3. Pour une particule partagée entre deux sites, $`\psi=(|01\rangle+|10\rangle)/\sqrt2`$, vérifier que
   $`E(\psi)=H(p)=\log2`$.

<details>
<summary>Solutions</summary>

*Exercice 1.*

```python
dual_tau = composite.partial_trace_channel(alg12, keep=1).adjoint('tau_vN')
A = random_operator(alg1)
Z = alg2.central([alg2.C * l for l in alg2.k_factors])
assert close(dual_tau(A), kron(A, Z, alg12))
```

*Exercice 2.*

```python
rho = states.random_density_matrix(alg12, batch_size=4)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
mutual = rho1.von_neumann_entropy() + rho2.von_neumann_entropy() - rho.von_neumann_entropy()
assert (mutual > 0).all()
assert torch.allclose(mutual, states.relative_entropy(rho, kron(rho1, rho2, alg12)))
```

*Exercice 3.*

```python
site = TypeIAlgebra([1, 1], [1, 1], precision='double', charges=[0, 1])
S2 = fused_tensor_product(site, site)
assert S2.k_factors == [1, 2, 1]
mat = torch.zeros(1, 3, 2, 2, dtype=dt)
mat[0, 1] = 0.5                                       # |psi><psi| with psi = (|01> + |10>)/sqrt 2
r1 = partial_trace(DensityMatrix(S2, matrix=mat), keep=1)
assert torch.allclose(r1.sector_probabilities()[0], torch.tensor([0.5, 0.5], dtype=torch.float64))
assert abs(r1.von_neumann_entropy().item() - math.log(2)) < 1e-12
```

</details>

**Suite :** [Tutoriel 8](08_dynamics.md) : les équations de Schrödinger et de Lindblad dans la bibliothèque.
