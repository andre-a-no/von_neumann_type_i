[English](../00_overview.md) · [Español](../es/00_overview.md) · **Français** · [Deutsch](../de/00_overview.md) · [中文](../zh/00_overview.md) · [日本語](../ja/00_overview.md) · [Русский](../ru/00_overview.md)

# Tutoriel 0. L'objet de ce cours

Ce tutoriel accompagne le [chapitre 0](../../lecture_notes/fr/00_overview.md) des notes de cours, qui identifie un élément d'une
algèbre de von Neumann $\bigoplus_cM_{k_c}(\mathbb C)$ à un tableau de forme $(C,k,k)$ et démontre que ses
opérations sont des opérations par lots sur des tableaux. Nous voyons ici ce que cela donne dans la bibliothèque
`torch_vn_algebra` : ce que cela coûte, un premier calcul, et la place de la bibliothèque parmi les logiciels existants.

Les algèbres de von Neumann s'étudient surtout avec un crayon et du papier. En dimension finie,
chacune d'elles est $*$-isomorphe à une somme directe d'algèbres de matrices pleines
(théorème 2.6) ; tout est donc calculable en principe, mais lorsqu'on traite numériquement une
telle algèbre, on l'aplatit souvent en une seule grande matrice, et sa structure (les secteurs,
le centre, les nombreuses traces) disparaît du calcul. Ce cours et la bibliothèque
`torch_vn_algebra` adoptent la démarche inverse : **l'algèbre elle-même est le type de
données**, et les calculs y sont menés efficacement.

<a id="s-0-1"></a>

## 0.1 L'observation

*Notes de cours :* [section 0.1](../../lecture_notes/fr/00_overview.md#s-0-1), qui démontre, à partir du [théorème 2.6](../../lecture_notes/fr/02_von_neumann_algebras.md#r-2-6),
du [théorème 3.1](../../lecture_notes/fr/03_traces_and_states.md#r-3-1), du [théorème 5.3](../../lecture_notes/fr/05_random_matrices.md#r-5-3) et du
[corollaire 5.4](../../lecture_notes/fr/05_random_matrices.md#r-5-4), que les six premières lignes du tableau suivant sont des identités :

| dans l'algèbre $\mathcal M$ | sur le tableau de forme $(B,C,k,k)$ |
|---|---|
| produit $AB$, somme, adjoint | produit matriciel par lots, somme, transposée conjuguée sur $(B,C)$ |
| calcul fonctionnel $f(A)$ (chapitre 1) | diagonalisation par lots des blocs |
| une trace $\sum_cw_c\mathrm{Tr}\,A_c$ (chapitre 3) | somme diagonale de chaque bloc, somme pondérée le long de l'axe des secteurs |
| le centre, les projecteurs centraux (chapitre 2) | un nombre par secteur |
| espérance conditionnelle sur le centre | trace de chaque bloc divisée par $k_c$, rediffusée sur le bloc |
| mesure de Haar sur le groupe unitaire $\prod_cU(k_c)$ | décomposition QR par lots de blocs gaussiens, avec la correction de phase du théorème 5.3 |
| états, canaux, dynamique (chapitres 3, 6, 8) | construits à partir des opérations ci-dessus |

Deux conséquences en font davantage qu'une simple affaire de comptabilité.

- **Coût.** Une décomposition spectrale dans $\mathcal M$ se calcule bloc par bloc (ligne 2 du
  tableau). Les solveurs denses standard de valeurs propres (LAPACK, que PyTorch appelle sur un
  CPU) effectuent un nombre d'opérations en virgule flottante qui croît comme le cube de la
  taille de la matrice ; c'est une description de ces algorithmes, non un résultat de ce cours.
  Dans le modèle où un problème aux valeurs propres $k\times k$ coûte une constante fois $k^3$, une
  décomposition spectrale dans $\mathcal M$ coûte $\sum_ck_c^3$, contre $D^3$ pour le plongement
  dans $M_D(\mathbb C)$, $D=\sum_ck_c$ ; pour $C$ blocs égaux de taille $k$, le rapport vaut
  $(Ck)^3/(Ck^3)=C^2$. Les produits s'effectuent sur le tenseur complété par des zéros ; le produit
  naïf de deux matrices $k\times k$ utilise $k^3$ multiplications ($k^2$ coefficients, chacun somme
  de $k$ produits), donc un produit dans $\mathcal M$ en utilise $C\,k_{\max}^3$, ce qui vaut
  $\sum_ck_c^3$ lorsque tous les blocs ont la même taille. Un élément est stocké sous forme de
  $C\,k_{\max}^2$ nombres, contre $D^2$ ; pour des blocs égaux, le rapport vaut $(Ck)^2/(Ck^2)=C$.
  Ce que fait le code : produits, sommes et traces sont des opérations tensorielles uniques sur
  l'ensemble du tableau $(B,C,k_{\max},k_{\max})$ ; les diagonalisations hermitiennes, les valeurs
  singulières et les décompositions QR s'appliquent à la pile des $B\cdot n_k$ blocs de chaque
  taille distincte $k$ ($n_k$ étant le nombre de secteurs de cette taille) en un seul appel par
  lots — un seul appel au total lorsque tous les secteurs ont la même taille. Sur un CPU avec
  $t\ge2$ threads PyTorch, la bibliothèque découpe un tel appel en $t$ morceaux exécutés dans des
  threads parallèles, pourvu que la pile compte au moins $2t$ blocs et que (nombre de blocs)$\,\cdot k^3\ge2^{18}$
  (module `cost`, fonction `_run`) ; les résultats sont les mêmes appels LAPACK sur les mêmes
  matrices. Le but est d'occuper pleinement un GPU, ou tous les cœurs d'un CPU.
- **Correction par construction.** Tout tableau de forme $(C,k_{\max},k_{\max})$ complété par des
  zéros est un élément de $\mathcal M$, de sorte que tout ce qui est calculé à partir des blocs
  reste dans $\mathcal M$. Les unitaires aléatoires de la bibliothèque sont tirés par la
  construction du théorème 5.3 indépendamment dans chaque bloc (`random_unitary_operator` appelle
  l'échantillonneur par bloc `random_unitary` une fois par taille de bloc, et la section 5.1
  vérifie `random_unitary` par rapport à cette construction), ce qui donne la mesure de Haar de
  $\prod_cU(k_c)$ d'après le corollaire 5.4 (pour des nombres aléatoires gaussiens idéaux) ; un
  unitaire de Haar de $U(D)$ appartient à ce sous-groupe avec probabilité nulle lorsque $C\ge2$
  (section 5.1, après le corollaire 5.4). Un pas de gradient modifie les blocs et ne peut donc pas
  sortir de l'algèbre. Enfin, les traces sont les fonctionnelles $\sum_cw_c\mathrm{Tr}\,A_c$ du
  théorème 3.1 avec des poids choisis explicitement (et les densités et espérances
  conditionnelles sont prises par rapport à elles), et non la trace d'un $\mathbb C^n$ sur lequel
  $\mathcal M$ se trouve agir, qui imposerait $w_c=m_c$ (section 2.3).

<a id="s-0-2"></a>

## 0.2 Un premier calcul

*Notes de cours :* [section 0.1](../../lecture_notes/fr/00_overview.md#s-0-1), ligne 2 du tableau : $f(A)=\bigoplus_cf(A_c)$, une diagonalisation par bloc.

Nous comparons le même calcul effectué dans $\mathcal M=\bigoplus_{c=1}^{64}M_{16}(\mathbb C)$ et dans
son plongement dans $M_{1024}(\mathbb C)$ : la racine carrée d'un élément positif aléatoire, et sa
trace. `TypeIAlgebra(n_factors, k_factors)` prend deux listes de tailles de blocs : les tailles
nominales $n_c$ et les tailles $k_c\le n_c$ des blocs effectivement utilisés (section 1.5). Toute
algèbre que le cours construit avec `TypeIAlgebra` a deux listes égales.

```python
import time
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)
C, k = 64, 16
alg = TypeIAlgebra([k] * C, [k] * C, precision='double')       # M = (+)_{c=1}^{64} M_16(C)

A = alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(1, d), force_positive=True)
blocks = A.matrix                                                # build A before timing

t0 = time.perf_counter()
root_blocks = A.sqrt().matrix                                    # (1, 64, 16, 16): 64 small eigenproblems
t_blocks = time.perf_counter() - t0

big = torch.block_diag(*blocks[0])                               # the same element in M_1024(C)
t0 = time.perf_counter()
w, V = torch.linalg.eigh(big)                                    # one large eigenproblem
root_big = (V * w.sqrt()) @ V.conj().T
t_big = time.perf_counter() - t0

assert torch.allclose(torch.block_diag(*root_blocks[0]), root_big, atol=1e-10)   # row 2: f(A) = (+)_c f(A_c)
assert abs(A.sqrt().Tr_blunt().item() - root_big.trace().item()) < 1e-9
assert big.numel() == C * blocks.numel()                          # storage: a factor C = 64
print(f"in the algebra: {t_blocks * 1e3:.1f} ms, as one 1024 x 1024 matrix: {t_big * 1e3:.1f} ms")
```

Les résultats coïncident aux erreurs d'arrondi près (ligne 2 du tableau), et le rapport de
stockage vaut $C$, comme calculé plus haut. Les 64 blocs passent au solveur de valeurs propres en
un seul appel par lots sur une pile de forme $(64,16,16)$ ; avec $2\le t\le32$ threads PyTorch, la
bibliothèque le découpe en $t$ morceaux ($64\cdot16^3=2^{18}$ atteint exactement le seuil du module
`cost`), avec un seul thread elle ne le découpe pas. *Observation numérique* (les durées
dépendent du matériel ; le code les affiche) : sur deux CPU à 4 cœurs, la version dans l'algèbre
a pris 3–10 ms et la version en grande matrice 400–760 ms, soit un facteur compris entre
plusieurs dizaines et quelques centaines (sur la seconde machine, avec 4 threads comme avec un
seul). Le facteur mesuré est bien inférieur au rapport $C^2=4096$ du modèle de coût cubique : un
décompte d'opérations est un modèle du temps d'exécution, non le temps d'exécution lui-même, et
sur ces machines les petits problèmes aux valeurs propres se sont exécutés à un débit
d'opérations par seconde bien plus faible que le grand. Dans le modèle de coût, le rapport $C^2$
croît avec le nombre de secteurs. Avec un lot de $B$ éléments, la version dans l'algèbre reste un
seul appel par lots ; la version en grande matrice doit résoudre $B$ problèmes aux valeurs propres
de taille $1024$.

<a id="s-0-3"></a>

## 0.3 Ce qui existe, et ce qui est nouveau

Nous avons cherché des logiciels qui calculent numériquement *à l'intérieur* des algèbres de von
Neumann de dimension finie. Ce que nous avons trouvé se répartit en quatre groupes, tous utiles
et tous conçus pour d'autres fins :

- **Les logiciels matriciels généraux** (NumPy/SciPy, QuTiP et d'autres boîtes à outils
  quantiques) travaillent sur une matrice à la fois. La structure par blocs, s'il y en a une,
  doit être gérée à la main.
- **Les bibliothèques de réseaux de tenseurs symétriques** stockent les tenseurs sous forme de
  tableaux creux par blocs indexés par les secteurs de charge d'une symétrie (abélienne dans
  YASTN [1], TeNPy et ITensor ; TensorKit traite aussi les symétries non abéliennes) ; YASTN
  dispose d'un backend PyTorch pour les GPU et la différentiation automatique. Techniquement, ce
  sont les plus proches parents. Leur objet est un réseau de tenseurs à contracter, non une
  algèbre d'opérateurs : pour autant que nous ayons pu le constater dans leur documentation,
  elles ne fournissent ni la famille des traces, états et densités, ni le centre et les
  espérances conditionnelles, ni les mesures de Haar sur $\prod_cU(k_c)$, ni les canaux entre
  secteurs.
- **La diagonalisation numérique par blocs** des $*$-algèbres de matrices [2, 3, 4] détermine la
  décomposition $\mathcal M\cong\bigoplus_cM_{k_c}\otimes1_{m_c}$ d'une algèbre donnée, typiquement
  pour réduire la taille de programmes semi-définis présentant une symétrie. C'est
  complémentaire : ces algorithmes trouvent les blocs ; cette bibliothèque calcule dedans. Son
  étape d'entrée (module `decompose`, section 4.3) reprend l'idée de [2, 3] selon laquelle les
  sous-espaces propres d'un élément aléatoire du commutant font apparaître les blocs ; la
  section 4.3 démontre qu'en arithmétique exacte cela réussit avec probabilité un (lemme 4.6,
  proposition 4.7, lemme 4.8), et la proposition 4.9 est le critère que `Decomposition.verify()`
  vérifie numériquement. Le calcul formel dispose de décompositions symboliques exactes des
  algèbres de groupe (par exemple le paquet Wedderga de GAP).
- **L'optimisation polynomiale non commutative** (NCSOStools, ncpol2sdpa) majore les valeurs
  propres et les traces de polynômes en variables non commutatives, pour toutes les tailles de
  matrices à la fois, au moyen de relaxations semi-définies [5]. Elle répond à une autre question
  (une borne valable pour toutes les tailles de matrices) par d'autres moyens (relaxation
  convexe).

*Notes de cours :* l'algorithme de décomposition randomisé et son critère — [section 4.3](../../lecture_notes/fr/04_symmetries.md#s-4-3), [proposition 4.7](../../lecture_notes/fr/04_symmetries.md#r-4-7), [proposition 4.9](../../lecture_notes/fr/04_symmetries.md#r-4-9).

Ce panorama décrit les logiciels tels que nous les avons trouvés ; ce n'est pas un théorème. À
notre connaissance, il n'existe pas de bibliothèque générale dans laquelle une algèbre de von
Neumann de dimension finie soit le type de données de base pour un calcul numérique par lots et
différentiable. La nouveauté n'est donc pas mathématique mais réside dans un **calcul
efficace** : l'identification du §0.1, implémentée de sorte que la structure de l'algèbre soit
préservée et que les opérations s'exécutent en noyaux par lots (un appel par taille de bloc
distincte pour les décompositions, découpé en un morceau par thread sur un CPU multicœur au-delà
du seuil de taille du §0.1). Si vous connaissez des logiciels apparentés, nous serons heureux
d'en entendre parler.

<a id="s-0-4"></a>

## 0.4 Comment lire ce cours

| chapitres | contenu |
|---|---|
| 1–3 | les mathématiques : algèbres de matrices, théorème de structure, traces, états et entropie, avec des démonstrations complètes |
| 4–8 | d'où viennent ces algèbres et ce que l'on y calcule : symétries, matrices aléatoires, canaux, produits tensoriels, dynamique |
| 9–12 | applications : chaînes de spins, calcul quantique, optimisation, apprentissage automatique équivariant |
| 13 | fondements : algèbres de Jordan, et pourquoi la mécanique quantique utilise des *-algèbres |
| annexes | démonstrations complètes des résultats auxiliaires que les chapitres citent par leur numéro (par exemple le théorème A.24) |

Chaque tutoriel accompagne le chapitre des [notes de cours](../../lecture_notes/fr/00_overview.md) de même numéro : il présente
le code de la bibliothèque, le confronte numériquement à des résultats qui y sont démontrés et indique les outils
standard auxquels le comparer. Chaque section commence par une ligne *Notes de cours :* qui renvoie à la section
et aux résultats qu'elle utilise. Les blocs de code d'un tutoriel sont faits pour être exécutés dans l'ordre, dans
une même session Python. Les nombres que les notes de cours ne démontrent pas (les durées d'exécution, par
exemple) n'apparaissent que comme *observations numériques* du code du tutoriel lui-même.

## Références

1. M. M. Rams, G. Wójtowicz, A. Sinha, J. Hasik, *YASTN: Yet another symmetric tensor networks; A
   Python library for Abelian symmetric tensor network calculations*, SciPost Phys. Codebases 52
   (2025).
2. K. Murota, Y. Kanno, M. Kojima, S. Kojima, *A numerical algorithm for block-diagonal
   decomposition of matrix $*$-algebras with application to semidefinite programming*, Japan J.
   Indust. Appl. Math. 27 (2010) 125–160 (prépublication : Technical Report METR 2007-52,
   University of Tokyo, 2007).
3. T. Maehara, K. Murota, *A numerical algorithm for block-diagonal decomposition of matrix
   $*$-algebras with general irreducible components*, Japan J. Indust. Appl. Math. 27 (2010)
   263–293.
4. E. de Klerk, C. Dobre, D. V. Pasechnik, *Numerical block diagonalization of matrix
   $*$-algebras with application to semidefinite programming*, Math. Program. 129 (2011) 91–111.
5. I. Klep, J. Povh, *Semidefinite programming and sums of hermitian squares of noncommutative
   polynomials*, J. Pure Appl. Algebra 214 (2010) 740–749.

**Suite :** [Tutoriel 1](01_matrix_algebras.md) : opérateurs, fonctions d'opérateurs et lots dans la bibliothèque.
