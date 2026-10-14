[English](../04_symmetries.md) · [Español](../es/04_symmetries.md) · **Français** · [Deutsch](../de/04_symmetries.md) · [中文](../zh/04_symmetries.md) · [日本語](../ja/04_symmetries.md) · [Русский](../ru/04_symmetries.md)

# Tutoriel 4. Symétries et représentations

Ce tutoriel accompagne le [chapitre 4](../../lecture_notes/fr/04_symmetries.md) des notes de cours, qui développe la théorie des représentations
sous-jacente aux algèbres de symétries, en tire un algorithme randomisé qui décompose l'algèbre de von Neumann engendrée par
une liste quelconque de matrices ([proposition 4.9](../../lecture_notes/fr/04_symmetries.md#r-4-9)), l'applique à $`U(1)`$ et $`SU(2)`$, et répartit les secteurs
d'une algèbre engendrée par des matrices réelles en type réel, complexe et quaternionique
([théorème 4.20](../../lecture_notes/fr/04_symmetries.md#r-4-20)). Nous calculons des caractères, appliquons la décomposition à des groupes et à des
générateurs, construisons des chaînes de spins secteur par secteur avec `SpinChain` et `su2.SU2Chain`, et calculons la
structure réelle d'une décomposition.

<a id="s-4-2"></a>

## 4.2 La décomposition isotypique et le commutant

*Notes de cours :* [section 4.2](../../lecture_notes/fr/04_symmetries.md#s-4-2) ; [théorème 4.3](../../lecture_notes/fr/04_symmetries.md#r-4-3), [théorème 4.4](../../lecture_notes/fr/04_symmetries.md#r-4-4), [corollaire 4.5](../../lecture_notes/fr/04_symmetries.md#r-4-5).

La [section 4.2](../../lecture_notes/fr/04_symmetries.md#s-4-2) des notes de cours calcule à la main les caractères de trois représentations
irréductibles de $`S_3`$ — la triviale, la signature et la représentation standard sur $`(1,1,1)^\perp\subseteq\mathbb C^3`$ — et de
la représentation $`\pi`$ qui permute les trois facteurs de $`(\mathbb C^2)^{\otimes3}`$, $`\chi_\pi(\sigma)=2^{\#\text{cycles}(\sigma)}`$
(table ci-dessous), et tire du [corollaire 4.5](../../lecture_notes/fr/04_symmetries.md#r-4-5) les multiplicités $`m_{\rm triv}=4`$, $`m_{\rm sgn}=0`$,
$`m_{\rm std}=2`$ et $`\dim\pi(S_3)'=20`$. Le code recalcule la table : il obtient $`P_\sigma`$ comme restriction de $`\pi(\sigma)`$
au sous-espace invariant engendré par $`|100\rangle,|010\rangle,|001\rangle`$ (sur lequel $`\pi(\sigma)`$ permute les positions
de l'unique $`1`$, lemme B.1(a)), et $`\#\text{cycles}(\sigma)=\dim\ker(P_\sigma-1)`$ (les vecteurs constants sur les cycles).

| | identité ($`1`$) | transpositions ($`3`$) | $`3`$-cycles ($`2`$) |
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

C'est la décomposition de la section 2.4, obtenue cette fois à partir d'une table de six nombres par représentation.

<a id="s-4-3"></a>

## 4.3 Des générateurs à l'algèbre : un algorithme randomisé

*Notes de cours :* [section 4.3](../../lecture_notes/fr/04_symmetries.md#s-4-3) ; [lemme 4.6](../../lecture_notes/fr/04_symmetries.md#r-4-6), [proposition 4.7](../../lecture_notes/fr/04_symmetries.md#r-4-7), [lemme 4.8](../../lecture_notes/fr/04_symmetries.md#r-4-8), [proposition 4.9](../../lecture_notes/fr/04_symmetries.md#r-4-9).

L'algorithme randomisé des notes de cours est `decompose.wedderburn`, avec les interfaces `from_group` (tous les éléments
d'un groupe fini ; le commutant est échantillonné avec le twirl `decompose.twirl`) et `from_generators` (une liste
quelconque de matrices ; le commutant est un noyau, `decompose.commutant_basis`). Il scinde $`\mathbb C^n`$ en les
sous-espaces propres d'un $`K\in\mathcal A'_h`$ aléatoire, les regroupe et les aligne à l'aide des blocs d'un
$`Y\in\mathcal A'`$ aléatoire, et renvoie une `Decomposition` avec l'unitaire `V`, les $`k_c`$ (`k_factors`) et les $`m_c`$
(`multiplicities`). En virgule flottante, il traite comme égales des valeurs propres plus proches qu'une tolérance et
scinde davantage les sous-espaces obtenus en y comprimant de nouveaux éléments aléatoires du commutant
([lemme 4.8(a)](../../lecture_notes/fr/04_symmetries.md#r-4-8)).

`Decomposition.verify()` renvoie les distances relatives des générateurs (et de nouveaux éléments du commutant) à
$`\mathcal B`$ (resp. $`\mathcal B'`$), et compare $`\sum_cm_c^2`$ à $`\dim\mathcal A'`$ (obtenue à partir des caractères pour
un groupe, corollaire 4.5(c) ; à partir du noyau pour des générateurs). Nous l'exécutons sur $`S_3`$ et trois qubits,
puis nous regardons à l'intérieur : les sous-espaces propres d'un $`K`$ aléatoire et les opérateurs d'entrelacement
d'un $`Y`$ aléatoire.

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

Le commutant domine le coût : un échantillon du twirl coûte $`|G|`$ produits de matrices $`n\times n`$. Pour des
générateurs, la bibliothèque diagonalise d'abord un élément hermitien $`R=U\mathrm{diag}(\lambda)U^*`$ de $`\mathcal A`$
(une combinaison réelle aléatoire des parties hermitiennes des générateurs, plus un multiple d'un produit symétrisé de
deux telles combinaisons). Tout $`X\in\mathcal A'`$ commute avec $`R`$, donc avec ses projecteurs spectraux (théorème
1.5(d)), si bien que $`U^*XU`$ est diagonale par blocs relativement aux sous-espaces propres de $`R`$, de dimensions
$`d_\lambda`$ ; on l'impose d'avance, et le problème de noyau a $`N_R=\sum_\lambda d_\lambda^2`$ inconnues au lieu de $`n^2`$.
Regrouper plusieurs valeurs propres ne fait qu'agrandir les blocs, donc le problème réduit a le même ensemble de
solutions $`\mathcal A'`$. Dans la décomposition, $`V^*RV=\bigoplus_cr_c\otimes1_{m_c}`$ avec des $`r_c\in M_{k_c}`$
hermitiennes ; si $`\lambda`$ a la multiplicité $`\mu_c(\lambda)`$ dans $`r_c`$, alors $`d_\lambda=\sum_c\mu_c(\lambda)m_c`$, et comme
$`\mu^2\ge\mu`$ pour les entiers $`\mu\ge0`$,
```math
N_R=\sum_\lambda\Big(\sum_c\mu_c(\lambda)m_c\Big)^2\ge\sum_{\lambda,c}\mu_c(\lambda)^2m_c^2\ge\sum_{\lambda,c}\mu_c(\lambda)m_c^2=\sum_ck_cm_c^2,
```
avec égalité si et seulement si les valeurs propres de toutes les $`r_c`$ sont simples et deux à deux distinctes. Pour le
spin total de $`L=3,4,5`$ spins, le code ci-dessous trouve l'égalité (observation numérique). Même le minimum est
grand : pour huit spins $`\frac12`$ ($`n=256`$), $`\sum_ck_cm_c^2=4900`$, et le système réduit, une matrice
$`(n^2+N_R)\times N_R`$, a $`3.5\cdot10^8`$ coefficients complexes. C'est pourquoi, pour les chaînes de spins, la
section 4.6 utilise le module `su2`, qui connaît la réponse à l'avance.

```python
for n in (3, 4, 5):
    Dn = dc.from_generators(dc.total_spin(n))
    _, info = dc.commutant_basis(dc.total_spin(n), return_info=True)
    N_min = sum(k * m * m for k, m in zip(Dn.k_factors, Dn.multiplicities))
    assert info['unknowns'] >= N_min                                  # N_R >= sum_c k_c m_c^2
    print(n, info['unknowns'], N_min)                                # observed: equality
```

<a id="s-4-4"></a>

## 4.4 $`U(1)`$ et secteurs de charge

*Notes de cours :* [section 4.4](../../lecture_notes/fr/04_symmetries.md#s-4-4) ; [proposition 4.10](../../lecture_notes/fr/04_symmetries.md#r-4-10).

Pour $`n`$ qubits, la charge est l'opérateur nombre $`\hat N=\sum_jn_j`$, où $`n_j=|1\rangle\langle1|`$ agit sur le qubit $`j`$, et
les opérateurs qui la conservent forment $`U(1)'\cong\bigoplus_{w=0}^nM_{\binom nw}(\mathbb C)`$
([notes de cours, section 4.4](../../lecture_notes/fr/04_symmetries.md#s-4-4)), qui est `chains.SpinChain(n).algebra` (bit $`1`$ = spin haut,
$`\hat N=S^z+\frac n2`$ ; dans `decompose.total_spin`, le spin haut est au contraire le *premier* vecteur de base de
$`\mathbb C^2`$, ce qui renomme la base mais pas les secteurs). Tout hamiltonien XXZ appartient à $`U(1)'`$. Le code trouve
les secteurs à partir de $`\hat N`$ seul et diagonalise un anneau XXZ de quatre spins secteur par secteur.

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

## 4.5 $`SU(2)`$ : le spin

*Notes de cours :* [section 4.5](../../lecture_notes/fr/04_symmetries.md#s-4-5) ; [théorème 4.11](../../lecture_notes/fr/04_symmetries.md#r-4-11), [proposition 4.13](../../lecture_notes/fr/04_symmetries.md#r-4-13), [théorème 4.14](../../lecture_notes/fr/04_symmetries.md#r-4-14). Le spin $`s`$ est `decompose.spin_matrices(s)`, et le spin total de $`L`$ spins $`\frac12`$ est `decompose.total_spin(L)`.

Pour $`L=4`$, [notes de cours, proposition 4.13](../../lecture_notes/fr/04_symmetries.md#r-4-13) donne $`\{S_x,S_y,S_z\}''\cong M_5\oplus M_3\oplus M_1`$
avec les multiplicités $`1,3,2`$, de dimension $`35`$, et le commutant $`\{S_x,S_y,S_z\}'\cong M_1\oplus M_3\oplus M_2`$, de
dimension $`14`$. La bibliothèque retrouve la même chose à partir des trois seuls générateurs (secteurs triés selon
$`(k_c,m_c)`$), ainsi que $`\mathbf S^2=S(S+1)`$ sur le spin $`S`$ ([théorème 4.11](../../lecture_notes/fr/04_symmetries.md#r-4-11)) :

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

D'après [notes de cours, théorème 4.14](../../lecture_notes/fr/04_symmetries.md#r-4-14) (dualité de Schur–Weyl), les permutations des sites
engendrent le commutant du spin total : $`\mathrm{span}\,\pi(S_L)=\{S_x,S_y,S_z\}'`$. Pour $`L\le5`$, le code compare
numériquement les deux espaces : d'après la partie (a) du théorème, $`\mathrm{span}\,\pi(S_L)\subseteq\{S_x,S_y,S_z\}'`$, il
suffit qu'ils aient la même dimension.

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

Les dimensions affichées sont $`2,5,14,42`$, les nombres de Catalan $`\binom{2L}L/(L+1)`$ ; le corollaire B.10 démontre
$`\dim\{S_a\}'=\sum_Sm_S^2=\binom{2L}L/(L+1)`$ pour tout $`L`$. La dernière ligne montre la
dualité à l'œuvre : $`\pi(S_4)''`$ a des blocs $`k=(1,2,3)`$ de multiplicités $`(5,1,3)`$ — les couples $`(k_c,m_c)`$ du
spin total, rôles échangés. En général (corollaire B.9), l'espace de multiplicité du spin $`S`$ porte une
représentation irréductible de $`S_L`$, de dimension $`m_S`$, et ces représentations sont deux à deux non équivalentes ;
ici les dimensions sont $`1,3,2`$ pour $`S=2,1,0`$.

<a id="s-4-6"></a>

## 4.6 Systèmes de spins

*Notes de cours :* [section 4.6](../../lecture_notes/fr/04_symmetries.md#s-4-6) ; [proposition 4.15](../../lecture_notes/fr/04_symmetries.md#r-4-15), [théorème 4.16](../../lecture_notes/fr/04_symmetries.md#r-4-16), [corollaire 4.17](../../lecture_notes/fr/04_symmetries.md#r-4-17).

D'après [notes de cours, corollaire 4.17](../../lecture_notes/fr/04_symmetries.md#r-4-17), des bases orthonormées des vecteurs de plus haut poids de
$`L`$ spins $`\frac12`$ se construisent site par site avec les coefficients de Clebsch–Gordan du [théorème 4.16](../../lecture_notes/fr/04_symmetries.md#r-4-16).
Étiqueter chaque vecteur de base par sa suite de spins intermédiaires $`S_1=\frac12,S_2,\dots,S_L=S`$ donne la base de
Yamanouchi que construit `su2.SU2Chain` (son code source utilise exactement ces coefficients ; le code ci-dessous ne
vérifie que les conclusions du corollaire 4.17), à l'intérieur du secteur $`U(1)`$ $`\hat N=\frac L2+S`$ de `SpinChain`
(dont l'ordre des sites est un réétiquetage des facteurs tensoriels, qui commute avec le spin total d'après le
lemme B.1(b)). Le code vérifie le théorème 4.16 sur des matrices denses et par comparaison avec `su2.clebsch_gordan` de
la bibliothèque, puis les conclusions du corollaire 4.17 pour les bases de huit spins. Nous vérifions aussi la
[proposition 4.15(a)](../../lecture_notes/fr/04_symmetries.md#r-4-15) sur un anneau de quatre spins avec `D4`, et comparons les tailles des blocs
pour huit spins.

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

Pour l'anneau de huit spins, les blocs de $`U(1)`$ ont pour tailles $`\binom8w=1,8,28,56,70,56,28,8,1`$ (section 4.4), les
blocs de $`SU(2)`$ les tailles $`m_S=14,28,20,7,1`$ pour $`S=0,\dots,4`$ (proposition 4.13). Le coût de la diagonalisation,
$`\sum k^3`$, passe de $`739162`$ à $`33040`$, soit un facteur $`22`$ ; la dimension de l'algèbre, $`\sum k^2`$, passe de
$`\sum_w\binom8w\binom8{8-w}=\binom{16}8=12870`$ (choisir $`8`$ éléments parmi $`16`$, dont $`w`$ parmi les $`8`$ premiers) à
$`\sum_Sm_S^2=1430=\binom{16}8/9`$ (corollaire B.10). D'après le théorème 9.9 (l'anneau est de longueur paire et $`J>0`$),
l'état fondamental est un singulet unique et $`E_0(0)<E_0(1)`$. *Observation numérique* (calculée par le code
ci-dessus) : l'énergie fondamentale vaut $`-3.651093`$, et le premier niveau excité est un triplet à $`-3.128419`$.

<a id="s-4-7"></a>

## 4.7 Algèbres réelles : type réel, complexe et quaternionique

*Notes de cours :* [section 4.7](../../lecture_notes/fr/04_symmetries.md#s-4-7). Pour une $`*`$-algèbre $`\mathcal A`$ stable par conjugaison complexe,
munie d'une décomposition $`V`$ dont les blocs de colonnes $`V_c`$ appartiennent aux secteurs, le
[théorème 4.20](../../lecture_notes/fr/04_symmetries.md#r-4-20) donne la permutation $`\sigma`$ des secteurs telle que $`\bar P_c=P_{\sigma(c)}`$, des unitaires
$`J_c`$, $`M_c`$ tels que $`\overline{V_c}=V_{\sigma(c)}(J_c\otimes M_c)`$ et, pour $`\sigma(c)=c`$, le signe $`\varepsilon_c`$ tel que
$`J_c\bar J_c=\varepsilon_c1`$ ; les types et l'algèbre réelle $`\mathcal A_{\mathbb R}`$ sont donnés par le
[théorème 4.21](../../lecture_notes/fr/04_symmetries.md#r-4-21), le lemme de Schur réel est le [théorème 4.22](../../lecture_notes/fr/04_symmetries.md#r-4-22), et les éléments
réels du commutant sont décrits par la [proposition 4.23](../../lecture_notes/fr/04_symmetries.md#r-4-23).

`Decomposition.real_structure()` calcule $`\sigma`$, $`J_c`$, $`M_c`$ et les types pour une décomposition renvoyée par
`wedderburn` (`from_generators`, `from_group`). Elle vérifie que les conjugués des générateurs stockés sont dans
$`V(\bigoplus_cM_{k_c}\otimes1)V^*`$ (résidu relatif `'conjugation'`) et forme la matrice unitaire $`P=V^*\bar V`$, dont
le bloc de colonnes $`c`$ est $`V^*\overline{V_c}`$ : d'après le théorème 4.20(a), il est porté par les lignes du secteur
$`\sigma(c)`$, où il vaut $`J_c\otimes M_c`$. La bibliothèque prend pour $`\sigma(c)`$ le secteur qui porte la plus grande
partie du bloc, exige $`(k_{\sigma(c)},m_{\sigma(c)})=(k_c,m_c)`$ et majore le reste (`'leak'`) ; elle réarrange le bloc
en la matrice $`k_c^2\times m_c^2`$ de coefficients $`T_{(a,b),(j,l)}=T_{(a,j),(b,l)}`$, qui pour $`T=J\otimes M`$ est la
matrice de rang un $`\mathrm{vec}\,J\,(\mathrm{vec}\,M)^T`$, lit $`J_c`$ et $`M_c`$ sur ses vecteurs singuliers
dominants (avec $`\|J_c\|_F^2=k_c`$), et enregistre la distance du bloc à $`J_c\otimes M_c`$ (`'kron'`). Si l'un de ces
trois résidus dépasse la tolérance, si les tailles diffèrent ou si les partenaires ne forment pas une involution, elle
lève `ValueError` (l'algèbre n'est pas stable par conjugaison). Pour $`\sigma(c)=c`$, elle prend pour $`\varepsilon_c`$ le
signe de $`\mathrm{Re}\,\mathrm{Tr}(J_c\bar J_c)/k_c`$ — exactement $`\pm1`$ en arithmétique exacte, de sorte
que les deux candidats sont distants de $`2`$ — et enregistre comme `'sign'` le maximum sur ces secteurs de
$`\|J_c\bar J_c-\varepsilon_c\|_F/\sqrt{k_c}`$ et de $`\|M_c\bar M_c-\varepsilon_c\|_F/\sqrt{m_c}`$, porté à au moins $`1`$
si un $`\varepsilon_c=-1`$ va avec $`k_c`$ ou $`m_c`$ impair (impossible d'après le théorème 4.20(b)). Le résultat a les
attributs `types`, `partner` ($`\sigma`$), `J`, `M`, `residuals` et `ok`, qui indique que tous les résidus sont au plus
égaux à la tolérance (par défaut le seuil d'acceptation de la décomposition) ; alors
$`|\mathrm{Tr}(J_c\bar J_c)/k_c-\varepsilon_c|\le\|J_c\bar J_c-\varepsilon_c\|_F/\sqrt{k_c}`$ (Cauchy–Schwarz) est
lui aussi au plus égal à la tolérance.

`Decomposition.real_basis()` renvoie une matrice orthogonale réelle $`O`$ de taille $`n\times n`$ et une liste qui décrit
ses blocs de colonnes (type, secteurs, plage de colonnes) : pour un secteur de type réel, les colonnes
$`O_c=V_c(W_J\otimes W_M)`$ du théorème 4.21(a), de sorte que $`O_c^TaO_c=\alpha_c(a)\otimes1_{m_c}`$ avec $`\alpha_c(a)`$
réelle ; pour un secteur de type quaternionique, les colonnes $`V_cW`$ avec $`J_c\otimes M_c=WW^T`$ et
$`(J_c\otimes M_c)\bar W=W`$ (lemme 4.19(c) avec $`\varepsilon=1`$, car
$`(J_c\otimes M_c)\overline{(J_c\otimes M_c)}=\varepsilon_c^2=1`$), qui sont réelles par le calcul du théorème 4.21(a) ;
et pour une paire $`c,\sigma(c)`$ (avec $`c<\sigma(c)`$), les colonnes $`\sqrt2\mathrm{Re}\,V_c`$ suivies de
$`\sqrt2\mathrm{Im}\,V_c`$. Elles sont orthonormées : $`V_c^*\overline{V_c}=V_c^*V_{\sigma(c)}(J_c\otimes M_c)=0`$ et
$`V_c^T\overline{V_c}=\overline{V_c^*V_c}=1`$, avec leurs conjuguées, donnent
$`(\mathrm{Re}\,V_c)^T\mathrm{Re}\,V_c=(\mathrm{Im}\,V_c)^T\mathrm{Im}\,V_c=\frac12`$ et
$`(\mathrm{Re}\,V_c)^T\mathrm{Im}\,V_c=0`$ ; les colonnes de blocs différents sont dans les images des
projecteurs orthogonaux $`P_{\mathcal O}=\sum_{e\in\mathcal O}P_e`$ d'orbites $`\mathcal O`$ différentes de $`\sigma`$ (car
$`\overline{V_c}=V_{\sigma(c)}(J_c\otimes M_c)`$), donc $`O`$ est orthogonale. En prenant les parties réelle et imaginaire
de $`aV_c=V_cx`$, $`x=\pi_c(a)\otimes1_{m_c}`$, pour $`a`$ réelle,
```math
O^TaO=\begin{pmatrix}\mathrm{Re}\,x&\mathrm{Im}\,x\\-\mathrm{Im}\,x&\mathrm{Re}\,x\end{pmatrix}\quad\text{on these columns.}
```
Les facteurs de Takagi $`W`$ avec $`WW^T=J`$ d'une matrice unitaire symétrique $`J=X+iY`$ sont calculés à partir d'une
diagonalisation orthogonale de $`X+tY`$ pour un $`t`$ aléatoire ($`X`$ et $`Y`$ sont des matrices réelles symétriques qui
commutent, car la partie imaginaire de $`J\bar J=JJ^*=1`$ est $`YX-XY`$) ; la bibliothèque accepte le résultat au vu de
ses résidus $`\|WW^T-J\|_F/\|J\|_F`$ et $`\|\mathrm{Im}(V_cW)\|_F/\sqrt{k_cm_c}`$ et émet un avertissement
(`DecompositionWarning`) si l'un d'eux dépasse $`\max(\text{tolérance},10^3\,\epsilon_{\rm mach}\,n)`$.

La [section 4.7](../../lecture_notes/fr/04_symmetries.md#s-4-7) des notes de cours traite trois algèbres engendrées par des matrices réelles :
(i) $`L_a\otimes1_2`$ sur $`\mathbb R^3\otimes\mathbb R^2`$, avec les générateurs $`(L_a)_{bc}=-\epsilon_{abc}`$ des rotations de
$`\mathbb R^3`$ : un secteur de type réel, $`k=3`$, $`m=2`$ ; (ii) $`\Omega_2\otimes1_2`$ sur $`\mathbb R^2\otimes\mathbb R^2`$ : une
paire de type complexe, $`k=1`$, $`m=2`$ ; (iii) les multiplications à gauche $`1_2\otimes L_q`$ par les quaternions sur
$`\mathbb H^2=\mathbb R^8`$ : un secteur de type quaternionique, $`k=2`$, $`m=4`$. Le code calcule la décomposition et la
structure réelle de chacune, et vérifie le théorème 4.20(a),(b) et la forme par blocs de $`O^TaO`$ pour le type réel
(théorème 4.21(a)) et pour la paire.

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

Le code vérifie le théorème 4.22 sur les modules réels irréductibles $`\mathbb R^3`$, $`\mathbb R^2`$, $`\mathbb R^4`$ des
trois exemples (les types et les multiplicités trouvés par la bibliothèque — réel avec $`m=1`$, une paire complexe avec
$`m=1`$, quaternionique avec $`m=2`$ — les rendent irréductibles d'après le théorème 4.22(c)), en calculant le commutant
réel comme le noyau des équations linéaires réelles $`Xs=sX`$, $`Xs^T=s^TX`$, puis la proposition 4.23 pour un élément
réel symétrique aléatoire du commutant dans chaque exemple.

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

## 4.8 Points de comparaison

Pour les groupes finis, GAP et SageMath calculent exactement les tables de caractères, et le corollaire 4.5 donne
alors les multiplicités ; RepLAB (MATLAB ; D. Rosset, F. Montealegre-Mora, J.-D. Bancal) décompose numériquement les
représentations des groupes finis et compacts, lui aussi à l'aide d'éléments aléatoires du commutant, et
SymbolicWedderburn.jl (Julia) calcule des décompositions de Wedderburn pour la réduction par symétrie de programmes
semi-définis. En physique, QuSpin et d'autres codes de diagonalisation exacte diagonalisent par blocs les symétries
$`U(1)`$ et de réseau, ITensor et TeNPy gèrent des charges abéliennes dans les réseaux de tenseurs, et TensorKit.jl
traite $`SU(2)`$ avec des arbres de fusion. Le module `decompose` de la bibliothèque est une méthode dense pour de
petites matrices ($`n`$ jusqu'à quelques centaines) ; il accepte n'importe quelle liste de matrices (éléments de groupe
ou générateurs), certifie son résultat (proposition 4.9) et renvoie l'algèbre et son commutant sous forme d'objets
`TypeIAlgebra` munis de `to_algebra`, `from_algebra`, `to_commutant`, `from_commutant`, de sorte que les chapitres 3
et 5–8 s'appliquent immédiatement. Pour les chaînes de spins $`\frac12`$, `su2.SU2Chain` fournit directement les blocs
de $`SU(2)`$ (corollaire 4.17) ; il accepte jusqu'à $`L=18`$ spins, une limite fixée par la mémoire des bases.

<a id="x-4"></a>

## Exercices

Les [exercices](../../lecture_notes/fr/04_symmetries.md#x-4) des notes de cours, vérifiés numériquement (les solutions ci-dessous suivent cet ordre) :
1. $`\dim\pi(S_L)'=\binom{d^2+L-1}L`$ pour $`d=2`$ et $`L=3,4`$ (les valeurs $`20`$ et $`35`$) ;
2. $`1\otimes\frac12\cong\frac12\oplus\frac32`$ et $`1\otimes1\cong0\oplus1\oplus2`$, avec `from_generators` ;
3. le commutant de $`e^{i\theta\hat N}`$ sur quatre qubits pour $`\theta=1`$ (les charges $`0,\dots,4`$ donnent des $`e^{iq}`$
   distincts, car $`0<|q-q'|\le4<2\pi`$) et pour $`\theta=\pi`$ ;
4. l'effet Zeeman : le spectre de $`H+hS_z`$ pour l'anneau de Heisenberg de huit spins, avec `su2.SU2Chain` et dans les
   secteurs $`U(1)`$ ;
5. les secteurs de la représentation régulière de $`Q_8`$, avec `real_structure`.

<details>
<summary>Solutions</summary>

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

**Suite :** [Tutoriel 5](05_random_matrices.md) : les matrices aléatoires dans une algèbre — unitaires de Haar, calcul de Weingarten et statistiques de niveaux dans la bibliothèque.
