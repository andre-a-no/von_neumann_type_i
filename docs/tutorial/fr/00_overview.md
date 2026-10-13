[English](../00_overview.md) · [Español](../es/00_overview.md) · **Français** · [Deutsch](../de/00_overview.md) · [中文](../zh/00_overview.md) · [日本語](../ja/00_overview.md) · [Русский](../ru/00_overview.md)

# 0. L'objet de ce cours

Les algèbres de von Neumann s'étudient presque exclusivement avec un crayon et du papier. En
dimension finie, chacune d'elles est une somme directe d'algèbres de matrices pleines
(chapitre 2) ; tout est donc calculable en principe, mais en pratique on aplatit le plus souvent
l'algèbre en une seule grande matrice, et sa structure (les secteurs, le centre, les nombreuses
traces) disparaît du calcul. Ce cours et la bibliothèque `torch_vn_algebra` adoptent la démarche
inverse : **l'algèbre elle-même est le type de données**, et les calculs y sont menés
efficacement.

## 0.1 L'observation

Une algèbre de von Neumann de dimension finie est, à un changement de base unitaire près,
$$
\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c}
$$
(théorème 2.6). En tant qu'algèbre abstraite, c'est $\bigoplus_cM_{k_c}(\mathbb C)$ : un élément est
une liste de $C$ blocs carrés. C'est exactement la forme d'un tenseur muni d'un *axe des
secteurs* : un élément de $\mathcal M$ est un tableau de forme $(C,k,k)$ (les blocs étant
complétés jusqu'à $k=\max_ck_c$), et $B$ éléments (échantillons de Monte-Carlo, points de départ
d'une optimisation, pas de temps, réalisations du désordre) forment un tableau de forme
$(B,C,k,k)$. Les bibliothèques modernes de tableaux (PyTorch sur CPU et GPU) sont conçues
précisément pour de tels « lots de petites matrices », dont les réseaux de neurones ont besoin.

Grâce à cette identification, les opérations de l'algèbre deviennent des opérations par lots sur
des tableaux, exactement et sans aucune approximation :

| dans l'algèbre $\mathcal M$ | sur le tableau de forme $(B,C,k,k)$ |
|---|---|
| produit $AB$, somme, adjoint | produit matriciel par lots, somme, transposée conjuguée sur $(B,C)$ |
| calcul fonctionnel $f(A)$ (chapitre 1) | diagonalisation par lots des blocs |
| une trace $\sum_cw_c\operatorname{Tr}A_c$ (chapitre 3) | somme diagonale de chaque bloc, somme pondérée le long de l'axe des secteurs |
| le centre, les projecteurs centraux (chapitre 2) | un nombre par secteur |
| espérance conditionnelle sur le centre | trace de chaque bloc divisée par $k_c$, rediffusée sur le bloc |
| mesure de Haar sur le groupe unitaire $\prod_cU(k_c)$ | décomposition QR par lots de blocs gaussiens (chapitre 5) |
| états, canaux, dynamique (chapitres 3, 6, 8) | les mêmes opérations, composées |

Deux conséquences en font davantage qu'une simple affaire de comptabilité.

- **Coût.** Une décomposition spectrale dans $\mathcal M$ coûte $\sum_ck_c^3$ opérations (les blocs
  sont regroupés par taille), contre $D^3$ pour le plongement dans $M_D(\mathbb C)$, $D=\sum_ck_c$ ;
  pour $C$ blocs égaux, cela représente un facteur $C^2$. Les produits et les sommes s'effectuent
  sur le tenseur complété par des zéros et coûtent $C\,k_{\max}^3$, ce qui vaut $\sum_ck_c^3$
  lorsque tous les blocs ont la même taille. Les produits, sommes et traces traitent les
  $B\cdot C$ blocs en un seul appel ; les décompositions (valeurs propres, valeurs singulières)
  s'exécutent en un appel par lots par taille de bloc distincte — un seul appel lorsque tous les
  secteurs ont la même taille — et, sur un CPU multicœur, la bibliothèque découpe chacun de ces
  appels en un morceau par thread. Cela occupe pleinement un GPU, ou tous les cœurs d'un CPU.
- **Correction par construction.** Tout ce qui est calculé reste dans $\mathcal M$ : les unitaires
  aléatoires suivent la loi de Haar sur $\prod_cU(k_c)$ (un ensemble négligeable dans $U(D)$, que
  l'on n'atteint jamais en échantillonnant $U(D)$), les pas de gradient ne sortent pas de
  l'algèbre, et la trace, les densités et les espérances conditionnelles se rapportent à l'algèbre
  et non à un plongement fortuit.

## 0.2 Un premier calcul

Nous comparons le même calcul effectué dans $\mathcal M=\bigoplus_{c=1}^{64}M_{16}(\mathbb C)$ et dans
son plongement dans $M_{1024}(\mathbb C)$ : la racine carrée d'un élément positif aléatoire, et sa
trace. `TypeIAlgebra(n_factors, k_factors)` prend deux listes de tailles de blocs : les tailles
nominales $n_c$ et les tailles $k_c\le n_c$ des blocs effectivement utilisés ; dans ce cours, les
deux listes sont toujours égales (section 1.5).

```python
import time
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
C, k = 64, 16
alg = TypeIAlgebra([k] * C, [k] * C, precision='double')       # M = (+)_{c=1}^{64} M_16(C)

A = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(1, d), force_positive=True)

t0 = time.perf_counter()
root_blocks = A.sqrt().matrix                                    # (1, 64, 16, 16): 64 small eigenproblems
t_blocks = time.perf_counter() - t0

big = torch.block_diag(*A.matrix[0])                             # the same element in M_1024(C)
t0 = time.perf_counter()
w, V = torch.linalg.eigh(big)                                    # one large eigenproblem
root_big = (V * w.sqrt()) @ V.conj().T
t_big = time.perf_counter() - t0

assert torch.allclose(torch.block_diag(*root_blocks[0]), root_big, atol=1e-10)   # same result
assert abs(A.sqrt().Tr_blunt().item() - root_big.trace().item()) < 1e-9
print(f"in the algebra: {t_blocks * 1e3:.1f} ms, as one 1024 x 1024 matrix: {t_big * 1e3:.1f} ms")
```

Les résultats coïncident aux erreurs d'arrondi près. La version dans l'algèbre est des dizaines de
fois plus rapide sur un CPU (environ 10 ms contre 400 ms sur une machine à 4 cœurs ; les chiffres
exacts dépendent du matériel — le code les affiche). Les 64 blocs passent au solveur de valeurs
propres en un seul appel par lots, que la bibliothèque découpe, sur un CPU multicœur, en un
morceau par thread. Le gain reste bien inférieur à ce que laisse penser le nombre d'opérations
$C^2=4096$, car un seul appel à un grand solveur de valeurs propres exploite mieux le matériel que
de nombreux petits. L'écart croît avec le nombre de secteurs, et le gain de mémoire est d'un
facteur $C$. Avec un lot de $B$ éléments, la version dans l'algèbre reste un seul appel par lots ;
la version en grande matrice devient $B$ grands problèmes aux valeurs propres.

## 0.3 Ce qui existe, et ce qui est nouveau

Nous avons cherché des logiciels qui calculent numériquement *à l'intérieur* des algèbres de von
Neumann de dimension finie. Ce que nous avons trouvé se répartit en quatre groupes, tous utiles
et tous conçus pour d'autres fins :

- **Les logiciels matriciels généraux** (NumPy/SciPy, QuTiP et d'autres boîtes à outils
  quantiques) travaillent sur une matrice à la fois. La structure par blocs, s'il y en a une,
  doit être gérée à la main.
- **Les bibliothèques de réseaux de tenseurs symétriques** stockent les tenseurs sous forme de
  tableaux creux par blocs indexés par les secteurs de charge d'une symétrie abélienne, certaines
  avec un backend PyTorch pour les GPU et la différentiation automatique (par exemple YASTN [1] ;
  aussi TeNPy, ITensor, TensorKit). Techniquement, ce sont les plus proches parents. Leur objet est
  un réseau de tenseurs à contracter, non une algèbre d'opérateurs : elles ne fournissent ni la
  famille des traces, états et densités, ni le centre et les espérances conditionnelles, ni les
  mesures de Haar sur $\prod_cU(k_c)$, ni les canaux entre secteurs.
- **La diagonalisation numérique par blocs** des $*$-algèbres de matrices [2, 3, 4] détermine la
  décomposition $\mathcal M\cong\bigoplus_cM_{k_c}\otimes1_{m_c}$ d'une algèbre donnée, typiquement
  pour réduire la taille de programmes semi-définis présentant une symétrie. C'est
  complémentaire : ces algorithmes trouvent les blocs ; cette bibliothèque calcule dedans. Le
  calcul formel dispose de décompositions symboliques exactes des algèbres de groupe (par exemple
  le paquet Wedderga de GAP).
- **L'optimisation polynomiale non commutative** (NCSOStools, ncpol2sdpa) majore les valeurs
  propres et les traces de polynômes en variables non commutatives, pour toutes les tailles de
  matrices à la fois, au moyen de relaxations semi-définies [5]. Elle répond à une autre question
  (une borne valable pour toutes les algèbres) par d'autres moyens (relaxation convexe).

À notre connaissance, il n'existe pas de bibliothèque générale dans laquelle une algèbre de von
Neumann de dimension finie soit le type de données de base pour un calcul numérique par lots et
différentiable. La nouveauté n'est donc pas mathématique mais réside dans un **calcul
efficace** : l'identification du §0.1, implémentée de sorte que la structure de l'algèbre soit
préservée et que chaque opération s'exécute en noyaux par lots (un par taille de bloc distincte pour les décompositions, découpé en un morceau par
thread sur un CPU multicœur). Si vous connaissez des
logiciels apparentés, nous serons heureux d'en entendre parler.

## 0.4 Comment lire ce cours

| chapitres | contenu |
|---|---|
| 1–3 | les mathématiques : algèbres de matrices, théorème de structure, traces, états et entropie, avec des démonstrations complètes |
| 4–8 | d'où viennent ces algèbres et ce que l'on y calcule : symétries, matrices aléatoires, canaux, produits tensoriels, dynamique |
| 9–12 | applications : chaînes de spins, calcul quantique, optimisation, apprentissage automatique équivariant |

Chaque chapitre énonce des résultats avec leurs démonstrations, présente le code de la
bibliothèque, le confronte à un résultat exact et indique les outils standard auxquels le
comparer.

## Références

1. M. M. Rams, G. Wójtowicz, A. Sinha, J. Hasik, *YASTN: Yet another symmetric tensor networks; A
   Python library for Abelian symmetric tensor network calculations*, SciPost Phys. Codebases 52
   (2025).
2. K. Murota, Y. Kanno, M. Kojima, S. Kojima, *A numerical algorithm for block-diagonal
   decomposition of matrix $*$-algebras*, Technical Report METR 2007-52, University of Tokyo (2007).
3. T. Maehara, K. Murota, *A numerical algorithm for block-diagonal decomposition of matrix
   $*$-algebras, Part II: general algorithm*, Optimization Online (2009).
4. E. de Klerk, C. Dobre, D. V. Pasechnik, *Numerical block diagonalization of matrix
   $*$-algebras with application to semidefinite programming*, Math. Program. 129 (2011) 91–111.
5. I. Klep, J. Povh, *Semidefinite programming and sums of hermitian squares of noncommutative
   polynomials*, J. Pure Appl. Algebra 214 (2010) 740–749.

**Suite :** [chapitre 1](01_matrix_algebras.md).
