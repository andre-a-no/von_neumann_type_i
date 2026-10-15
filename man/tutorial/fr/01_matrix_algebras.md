[English](../01_matrix_algebras.md) · [Español](../es/01_matrix_algebras.md) · **Français** · [Deutsch](../de/01_matrix_algebras.md) · [中文](../zh/01_matrix_algebras.md) · [日本語](../ja/01_matrix_algebras.md) · [Русский](../ru/01_matrix_algebras.md)

# Tutoriel 1. Algèbres de matrices et calcul fonctionnel

Ce tutoriel accompagne le [chapitre 1](../../lecture_notes/fr/01_matrix_algebras.md) des notes de cours, qui démontre le théorème
spectral, le calcul fonctionnel, la décomposition polaire et la dualité de la norme trace dans $`M_n(\mathbb C)`$.
Chacun de ces résultats est ici vérifié avec la bibliothèque.

**Vous apprendrez** à construire des opérateurs de spectre prescrit, à appliquer une
fonction à un opérateur, à manipuler des lots d'opérateurs, et à confronter un calcul à un
résultat exact.

<a id="s-1-5"></a>

## 1.5 Les opérateurs dans la bibliothèque

*Notes de cours :* [section 1.2](../../lecture_notes/fr/01_matrix_algebras.md#s-1-2) et [section 1.3](../../lecture_notes/fr/01_matrix_algebras.md#s-1-3) : [théorème 1.3](../../lecture_notes/fr/01_matrix_algebras.md#r-1-3), [corollaire 1.4](../../lecture_notes/fr/01_matrix_algebras.md#r-1-4), [théorème 1.5](../../lecture_notes/fr/01_matrix_algebras.md#r-1-5), [théorème 1.6](../../lecture_notes/fr/01_matrix_algebras.md#r-1-6), [théorème 1.11](../../lecture_notes/fr/01_matrix_algebras.md#r-1-11), [proposition 1.12](../../lecture_notes/fr/01_matrix_algebras.md#r-1-12), [proposition 1.13](../../lecture_notes/fr/01_matrix_algebras.md#r-1-13).

La bibliothèque travaille avec des algèbres $`\bigoplus_cM_{k_c}(\mathbb C)`$ (chapitre 2). Un
bloc unique correspond au cas $`C=1`$ : `TypeIAlgebra([n], [n])` est $`M_n(\mathbb C)`$. (Les
deux listes sont les tailles nominales et les tailles effectives des blocs ; pour l'instant
elles coïncident.)

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

Chaque opérateur est stocké sous la forme d'un unique tenseur de forme
`(batch, C, k_max, k_max)` : un lot d'opérateurs, chacun étant une liste de `C` blocs. Ici
`batch = C = 1`. La construction n'est autre que le théorème 1.3 lu à l'envers :
$`A=U\Lambda U^*`$ avec $`U`$ distribuée selon la mesure de Haar (la bibliothèque tire $`U`$ par la
construction QR du théorème 5.3). L'option `force_positive=True` est contrôlée sur les valeurs
propres prescrites (la bibliothèque lève une erreur si l'une d'elles est inférieure à
$`-10^{-12}`$) ; d'après le corollaire 1.4(3), le résultat est alors positif (exactement lorsque les
valeurs propres prescrites sont $`\ge0`$, comme ici ; à cette tolérance près en général), et la
bibliothèque l'enregistre au lieu de le tester de nouveau. Pour la même raison,
`A.eigenvalues()` ci-dessus renvoie le spectre prescrit, trié, sans rien calculer ; ce n'est pas
une vérification. La vérification de $`\mathrm{spec}\,f(A)=f(\mathrm{spec}\,A)`$ ci-dessous,
elle, calcule : `apply_function` diagonalise la matrice de $`A`$, et les valeurs propres de la
matrice $`e^A`$ obtenue sont calculées à partir de zéro.

Le calcul fonctionnel pour les opérateurs autoadjoints est assuré par `apply_function` (qui
diagonalise chaque bloc au moyen de `eigh`). La fonction $`f`$ doit agir coefficient par
coefficient : elle reçoit les valeurs propres de tous les secteurs d'une même taille $`k`$ sous la
forme d'un seul tenseur de forme `(batch, n·k)`, où $`n`$ est le nombre de ces secteurs, et renvoie
un tenseur de même forme. Les fonctions usuelles disposent de leurs propres méthodes. Vérifions
les théorèmes 1.5 et 1.6 :

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

`R @ R` ne calcule encore rien : produits, sommes et fonctions construisent une recette,
qui n'est évaluée que lorsqu'on demande `.matrix` ou un nombre tel qu'une trace.

**Dérivées.** D'après le [théorème 1.11](../../lecture_notes/fr/01_matrix_algebras.md#r-1-11) (Daleckii–Krein), la dérivée
de $`f(A)`$ dans une direction autoadjointe $`E`$ est $`U(L\circ U^*EU)U^*`$, avec la matrice de Loewner
$`L_{ij}=f^{[1]}(\lambda_i,\lambda_j)`$ des différences divisées ($`f'(\lambda_i)`$ pour des valeurs propres égales).
`functional.divided_differences` calcule $`L`$ à partir des valeurs propres et des valeurs de $`f`$ et $`f'`$, et
autograd dérive `apply_function` avec elle. Nous vérifions la formule en un spectre dégénéré, contre un quotient
différentiel centré de l'exponentielle dense, dont l'erreur est $`O(s^2)`$ plus un arrondi de l'ordre de $`\varepsilon/s`$ :

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

Pour des valeurs propres proches $`a,b`$ à distance $`h`$, le quotient $`\hat q`$ perd en précision par compensation, et
la valeur au point milieu $`\hat m=(f'(a)+f'(b))/2`$ se trompe d'au plus $`T=h^2\max|f'''|/12`$ plus l'arrondi
([proposition 1.12](../../lecture_notes/fr/01_matrix_algebras.md#r-1-12)). `divided_differences` prend $`\hat m`$ si et
seulement si $`|\hat q-\hat m|\le2R`$, avec la borne d'arrondi $`R`$ de la
[proposition 1.13](../../lecture_notes/fr/01_matrix_algebras.md#r-1-13) et l'erreur d'évaluation
$`\varepsilon(|\hat f_a|+|a\hat f'_a|)`$ ; son erreur est alors au plus $`\min(T,3R)`$. Nous le vérifions pour
$`f(x)=e^{x-c}`$ près de $`c=10^4`$, où $`f^{[1]}(a,b)=e^{a-c}(e^h-1)/h`$ se calcule précisément avec `expm1`. La dernière
distance, $`h=0.06`$, est sous le seuil $`\varepsilon^{1/3}\max(|a|,|b|)`$ des versions antérieures, qui y prenaient la
valeur au point milieu :

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

*Observation numérique :* le quotient se trompe jusqu'à $`5\cdot10^{-10}`$ pour $`h\le10^{-6}`$ et de $`10^{-15}`$ en
$`h=0.06`$ ; la valeur au point milieu est exacte à l'arrondi près pour $`h\le10^{-9}`$ et se trompe de $`3\cdot10^{-4}`$ en
$`h=0.06`$. La bibliothèque a la plus petite des deux erreurs sur chaque ligne sauf $`h=10^{-4}`$, où elle prend la
valeur au point milieu ($`8\cdot10^{-10}`$) bien que le quotient soit meilleur ($`10^{-12}`$) : $`R`$ admet un argument
perturbé de $`\varepsilon|a|\approx2\cdot10^{-12}`$, ce qui est pessimiste pour `torch.exp` mais nécessaire pour des
fonctions comme $`e^{-i\tau x}`$, qui arrondissent $`\tau x`$. La garantie $`\min(T,3R)`$ est respectée sur chaque ligne.

<a id="s-1-6"></a>

## 1.6 Normes et dualité de la norme trace

*Notes de cours :* [section 1.4](../../lecture_notes/fr/01_matrix_algebras.md#s-1-4) : [théorème 1.8](../../lecture_notes/fr/01_matrix_algebras.md#r-1-8), les trois normes, [théorème 1.9](../../lecture_notes/fr/01_matrix_algebras.md#r-1-9).

Une matrice quelconque (non normale) s'introduit par `alg.operator(tensor)`. Nous
comparons les trois normes aux valeurs singulières calculées directement par PyTorch, puis
nous vérifions la décomposition polaire (théorème 1.8) et la dualité (théorème 1.9) :

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

Avec un dtype complexe, `torch.randn` tire des coefficients indépendants dont les parties réelles
et imaginaires sont indépendantes et de loi $`N(0,\tfrac12)`$ ; `G` est donc une matrice de Ginibre,
qui admet une densité (proposition A.16(b)). $`|X|`$ est presque sûrement inversible :
$`(\det|X|)^2=\det(|X|^2)=\det(X^*X)=|\det X|^2`$, et $`|\det X|^2`$ est un polynôme en les parties
réelles et imaginaires des coefficients qui vaut $`1`$ en $`X=1`$, donc il s'annule avec probabilité
nulle (corollaire A.7(b) ; c'est aussi la première étape du théorème 5.3).

<a id="s-1-7"></a>

## 1.7 Lots : de nombreux opérateurs à la fois

*Notes de cours :* [section 1.7](../../lecture_notes/fr/01_matrix_algebras.md#s-1-7), [proposition 1.10](../../lecture_notes/fr/01_matrix_algebras.md#r-1-10).

L'axe des lots est la raison d'être de la bibliothèque : la théorie des matrices
aléatoires, les méthodes de Monte-Carlo et l'optimisation à départs multiples requièrent
des milliers d'opérateurs, et un seul appel vectorisé remplace une boucle Python sur
ceux-ci. À titre d'exemple, vérifions un moment exact de l'**ensemble gaussien unitaire**
(GUE).

D'après [notes de cours, proposition 1.10](../../lecture_notes/fr/01_matrix_algebras.md#r-1-10), si $`G`$ a des coefficients indépendants dont les
parties réelles et imaginaires sont indépendantes et $`\sim N(0,\tfrac12)`$, et si $`X=(G+G^*)/\sqrt{2n}`$, alors
$`\mathbb E\,\tfrac1n\mathrm{Tr}\,X^4=2+\tfrac1{n^2}`$.

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

Une vérification par Monte-Carlo est un test statistique : l'estimation doit coïncider avec
la valeur exacte à quelques écarts-types près. La variable $`\tfrac1n\mathrm{Tr}\,X^4`$ est un
polynôme en des variables gaussiennes ; elle a donc une variance finie $`\sigma^2`$ (lemme A.23 : tous
les moments gaussiens sont finis ; les produits se traitent par $`|y_1\cdots y_r|\le\sum_i|y_i|^r`$), et
d'après l'inégalité de Bienaymé–Tchebychev la moyenne de $`B`$ échantillons indépendants s'écarte de
la valeur exacte de plus de $`5\sigma/\sqrt B`$ avec une probabilité d'au plus $`1/25`$. Le code utilise
l'écart-type empirique à la place de $`\sigma`$, et la graine fixée rend la vérification
reproductible. C'est le schéma de toutes les vérifications par rapport à des résultats démontrés
dans ce cours : un énoncé exact, une estimation numérique et une tolérance explicite.

<a id="s-1-8"></a>

## 1.8 Points de comparaison

Pour des matrices isolées, `numpy.linalg` et `scipy.linalg` (`sqrtm`, `expm`, `logm`,
`polar`) sont les implémentations de référence, et les vérifications ci-dessus peuvent être
refaites avec elles. La bibliothèque apporte trois éléments qui compteront par la suite :
le traitement par lots (section 1.7), des opérateurs porteurs d'une structure par blocs
(chapitre 2), et la différentiation automatique (l'autograd de PyTorch) à travers ses
opérations.

<a id="x-1"></a>

## Exercices

[Exercices 1 et 3](../../lecture_notes/fr/01_matrix_algebras.md#x-1) des notes de cours, numériquement :
- Vérifier que $`(fg)(A)=f(A)g(A)`$ pour $`f=\sin`$, $`g=\exp`$ et $`A`$ autoadjointe aléatoire (exercice 1).
- La **moyenne géométrique** de $`A,B>0`$ est $`A\#B=A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}`$, l'unique solution
  positive $`X`$ de $`XA^{-1}X=B`$ (exercice 3). Vérifier numériquement cette équation et $`A\#B=B\#A`$.

<details>
<summary>Solutions</summary>

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

**Suite :** [Tutoriel 2](02_von_neumann_algebras.md) : commutants, centre et algèbre $`\bigoplus_cM_{k_c}(\mathbb C)`$ dans la bibliothèque.
