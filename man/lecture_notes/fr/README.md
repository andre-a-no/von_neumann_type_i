[English](../README.md) · [Español](../es/README.md) · **Français** · [Deutsch](../de/README.md) · [中文](../zh/README.md) · [日本語](../ja/README.md) · [Русский](../ru/README.md)

# Algèbres de von Neumann de dimension finie : notes de cours

Ce cours présente la théorie des algèbres de von Neumann de dimension finie — les algèbres
$`\bigoplus_cM_{k_c}(\mathbb C)`$ de matrices diagonales par blocs — ainsi que leurs usages en physique quantique, en
information quantique, en théorie des matrices aléatoires et en apprentissage automatique : secteurs de symétrie et de
supersélection, théorème du bicommutant et théorème de structure de Wedderburn, traces, états et entropie, mesure de
Haar et calcul de Weingarten, canaux quantiques (Kraus, Stinespring, Choi), équation de Lindblad (GKSL), chaînes de
spins, correction d'erreurs quantiques, réseaux de neurones équivariants et algèbres de Jordan.

Il s'adresse simultanément à deux publics :

- **Les mathématiciens** y trouveront un texte autonome : chaque théorème, formule, constante
  et inégalité énoncés sont démontrés dans le cours — dans le chapitre lui-même ou dans l'une
  des annexes A–F, citée par son numéro — à partir des prérequis énumérés ci-dessous. Les
  nombres qui sont seulement calculés (par exemple la statistique des niveaux pour $`n`$ grand)
  sont signalés comme *observation numérique* et ne servent jamais de prémisses ni de cibles
  pour les tests. Les références relèvent de l'attribution historique ou de la lecture
  complémentaire, et ne remplacent jamais une démonstration.
- **Les physiciens et les ingénieurs** y trouveront la signification physique de chaque
  notion, dans des encadrés intitulés *Interprétation physique*. Ces encadrés relèvent de l'interprétation physique et ne sont
  pas des théorèmes du cours : rien de ce qui est démontré ailleurs n'en dépend. Du code
  fonctionnel qu'ils pourront adapter se trouve dans le [tutoriel](../../tutorial/fr/README.md).

Les chapitres 1 à 13 suivent le même plan : définitions et résultats avec démonstrations complètes,
exemples traités et exercices corrigés. Après chaque résultat que la bibliothèque `torch_vn_algebra` calcule ou
vérifie, une ligne *Calcul :* ou *Dans la bibliothèque :* renvoie à l'endroit du [tutoriel](../../tutorial/fr/README.md)
où cela est fait.

## Table des matières

| | chapitre | résultats principaux |
|---|---|---|
| 0 | [L'objet de ce cours](00_overview.md) | l'algèbre comme type de données : structure par blocs = tenseurs par lots ; coût ; ce qui existe et ce qui est nouveau |
| 1 | [Algèbres de matrices et calcul fonctionnel](01_matrix_algebras.md) | théorème spectral, calcul fonctionnel, décomposition polaire, dualité de la norme trace, un moment du GUE |
| 2 | [Algèbres de von Neumann de dimension finie](02_von_neumann_algebras.md) | théorème du bicommutant, théorème de structure, centre, groupes de symétrie et commutants, espérances conditionnelles |
| 3 | [Traces, états et entropie](03_traces_and_states.md) | classification des traces, densités, entropie des secteurs, inégalité de Klein, monotonie sous l'espérance conditionnelle sur le centre, contraste de Michelson |
| 4 | [Symétries et représentations](04_symmetries.md) | lemme de Schur, secteurs de charge, $`U(1)`$ et $`SU(2)`$, systèmes de spins, type réel, complexe et quaternionique |
| 5 | [Matrices aléatoires dans une algèbre](05_random_matrices.md) | mesure de Haar sur $`\prod_cU(k_c)`$, calcul de Weingarten, statistique des niveaux et pourquoi il ne faut pas mélanger les secteurs |
| 6 | [Canaux quantiques](06_quantum_channels.md) | Kraus et Stinespring, duaux relativement à une trace, applications entre secteurs, décohérence |
| 7 | [Produits tensoriels et intrication](07_tensor_products.md) | trace partielle, formule de Page, fusion des charges |
| 8 | [Dynamique](08_dynamics.md) | équations de Schrödinger et de Lindblad, sauts entre secteurs, post-sélection |
| 9 | [Chaînes de spins](09_spin_chains.md) | secteurs d'aimantation et d'impulsion, Lanczos, vérifications par l'ansatz de Bethe, désordre |
| 10 | [Calcul quantique](10_quantum_computing.md) | sous-espaces sans décohérence, correction d'erreurs par algèbres d'opérateurs, mélangeurs sous contraintes |
| 11 | [Optimisation et inégalités de traces](11_optimisation.md) | paramétrisations sous contraintes, échantillonnage ou optimisation |
| 12 | [Apprentissage automatique équivariant](12_equivariant_learning.md) | les couches équivariantes comme commutants |
| 13 | [Algèbres de Jordan : pourquoi la mécanique quantique utilise des *-algèbres](13_jordan_algebras.md) | produit de Jordan et représentation quadratique, facteurs de spin et cône de Lorentz, classification de Jordan–von Neumann–Wigner, transposée et transposée partielle, une inégalité de traces sous forme jordanienne |
| A | [Résultats d'analyse et de probabilités](A_analysis.md) | séries entières et ensembles négligeables, vecteurs gaussiens, mesures sphériques, inégalités sur les valeurs propres, Isserlis–Wick |
| B | [Dualité de Schur–Weyl et calcul de Weingarten](B_schur_weyl.md) | dualité de Schur–Weyl, fonctions de Weingarten, théorème de Wigner–Eckart, mesure de Haar sur $`O(k)`$ |
| C | [Densités de valeurs propres et formule de Page](C_eigenvalue_densities.md) | formule d'intégration de Weyl, GOE/GUE, densités de Wishart et densités induites, formule de Page |
| D | [Générateurs des semi-groupes dynamiques quantiques](D_gksl.md) | le théorème GKSL, générateurs sur les sommes directes, théorème de Frigerio, Bény–Kempf–Kribs |
| E | [Classification des algèbres de Jordan euclidiennes](E_jordan_classification.md) | décomposition de Peirce, théorèmes de Hurwitz et de Frobenius, coordinatisation, l'algèbre d'Albert, le théorème de Jordan–von Neumann–Wigner |
| F | [Cônes symétriques, isomorphismes d'ordre et applications de Jordan](F_cones_and_order.md) | cônes des carrés, le théorème de Koecher–Vinberg, le théorème de Kadison, automorphismes de Jordan de $`\bigoplus_cM_{k_c}`$ |

<a id="prerequisites"></a>

## Prérequis

Le cours suppose seulement ce qui suit, et démontre tout le reste.

- L'algèbre linéaire sur $`\mathbb R`$ et $`\mathbb C`$ (en dimension finie), et le théorème de
  d'Alembert-Gauss.
- L'analyse réelle et le calcul différentiel à plusieurs variables classiques : limites,
  continuité, compacité dans $`\mathbb R^n`$, différentiation, le théorème des accroissements
  finis, la formule de Taylor, les théorèmes d'inversion locale et des fonctions implicites.
- La théorie classique de la mesure et de l'intégration : mesure de Lebesgue, mesures
  produits et Fubini–Tonelli, convergence dominée et convergence monotone, changement de
  variables dans $`\mathbb R^n`$, mesures images, densités, indépendance.
- Les probabilités élémentaires : espérance, variance, indépendance, la loi gaussienne et sa
  fonction génératrice des moments, les inégalités de Markov et de Bienaymé-Tchebychev, la
  loi des grands nombres.

Aucune connaissance en physique ni en analyse fonctionnelle n'est supposée.

## Pour aller plus loin

M. Takesaki, *Theory of Operator Algebras I* (Springer, 2002) ; V. Jones,
*von Neumann algebras*, notes de cours ; D. Petz, *Quantum Information Theory and Quantum
Statistics* (Springer, 2008) ; J. Watrous, *The Theory of Quantum Information* (Cambridge,
2018) ; M. L. Mehta, *Random Matrices* (Elsevier, 2004).

La bibliothèque `torch_vn_algebra` et son usage dans le contexte de ces notes font l'objet du
[tutoriel](../../tutorial/fr/README.md). L'article qui accompagne la bibliothèque est
[`man/article/main.pdf`](https://github.com/andre-a-no/von_neumann_type_i/blob/main/man/article/main.pdf).
