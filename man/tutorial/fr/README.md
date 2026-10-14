[English](../README.md) · [Español](../es/README.md) · **Français** · [Deutsch](../de/README.md) · [中文](../zh/README.md) · [日本語](../ja/README.md) · [Русский](../ru/README.md)

# Tutoriel : torch_vn_algebra dans le contexte des notes de cours

Ce tutoriel montre comment calculer avec les algèbres de von Neumann de dimension finie $\bigoplus_cM_{k_c}(\mathbb C)$
dans la bibliothèque `torch_vn_algebra`, dans le contexte des [notes de cours](../../lecture_notes/fr/README.md). Chaque tutoriel
accompagne le chapitre des notes de cours de même numéro : il présente le code de la bibliothèque pour les résultats
qui y sont démontrés, vérifie ce code par rapport à un résultat exact (un théorème, une formule close ou un calcul
indépendant) et cite, dans une section *Points de comparaison*, les outils de référence pour la même tâche. Chaque
section commence par une ligne *Notes de cours :* qui renvoie à la section et aux résultats qu'elle utilise. Les
nombres qui sont seulement calculés (par exemple la statistique des niveaux pour $n$ grand) sont signalés comme
*observation numérique* et ne servent jamais de prémisses ni de cibles pour les tests. Outre les
[prérequis](../../lecture_notes/fr/README.md#prerequisites) des notes de cours, le tutoriel suppose Python.

Tout le code du tutoriel est exécuté par la suite de tests (`tests/test_man.py`), de
sorte que chaque exemple fonctionne avec la version actuelle de la bibliothèque.

## Exécuter le code

```bash
pip install -e .                # from the root of the repository
python -m pytest tests/test_man.py -q
```

(La première commande s'exécute depuis la racine du dépôt.) Chaque bloc de code d'un
chapitre peut être collé, dans l'ordre, dans une même session Python ou un même notebook
Jupyter. Les blocs d'un même chapitre partagent leurs variables.

## Table des matières

| | chapitre | résultats principaux |
|---|---|---|
| 0 | [L'objet de ce cours](00_overview.md) | l'algèbre comme type de données : structure par blocs = tenseurs par lots ; coût ; ce qui existe et ce qui est nouveau |
| 1 | [Algèbres de matrices et calcul fonctionnel](01_matrix_algebras.md) | théorème spectral, calcul fonctionnel, décomposition polaire, dualité de la norme trace, un moment du GUE |
| 2 | [Algèbres de von Neumann de dimension finie](02_von_neumann_algebras.md) | théorème du bicommutant, théorème de structure, centre, groupes de symétrie et commutants, espérances conditionnelles |
| 3 | [Traces, états et entropie](03_traces_and_states.md) | classification des traces, densités, entropie des secteurs, inégalité de Klein, monotonie sous l'espérance conditionnelle sur le centre, contraste de Michelson |
| 4 | [Symétries et représentations](04_symmetries.md) | lemme de Schur, secteurs de charge, $U(1)$ et $SU(2)$, systèmes de spins, type réel, complexe et quaternionique |
| 5 | [Matrices aléatoires dans une algèbre](05_random_matrices.md) | mesure de Haar sur $\prod_cU(k_c)$, calcul de Weingarten, statistique des niveaux et pourquoi il ne faut pas mélanger les secteurs |
| 6 | [Canaux quantiques](06_quantum_channels.md) | Kraus et Stinespring, duaux relativement à une trace, applications entre secteurs, décohérence |
| 7 | [Produits tensoriels et intrication](07_tensor_products.md) | trace partielle, formule de Page, fusion des charges |
| 8 | [Dynamique](08_dynamics.md) | équations de Schrödinger et de Lindblad, sauts entre secteurs, post-sélection |
| 9 | [Chaînes de spins](09_spin_chains.md) | secteurs d'aimantation et d'impulsion, Lanczos, vérifications par l'ansatz de Bethe, désordre |
| 10 | [Calcul quantique](10_quantum_computing.md) | sous-espaces sans décohérence, correction d'erreurs par algèbres d'opérateurs, mélangeurs sous contraintes |
| 11 | [Optimisation et inégalités de traces](11_optimisation.md) | paramétrisations sous contraintes, échantillonnage ou optimisation |
| 12 | [Apprentissage automatique équivariant](12_equivariant_learning.md) | les couches équivariantes comme commutants |
| 13 | [Algèbres de Jordan : pourquoi la mécanique quantique utilise des *-algèbres](13_jordan_algebras.md) | produit de Jordan et représentation quadratique, facteurs de spin et cône de Lorentz, classification de Jordan–von Neumann–Wigner, transposée et transposée partielle, une inégalité de traces sous forme jordanienne |
| A | [Résultats d'analyse et de probabilités](A_analysis.md) | séries entières et ensembles négligeables, vecteurs gaussiens, mesures sphériques, inégalités sur les valeurs propres, Isserlis–Wick |
| B | [Dualité de Schur–Weyl et calcul de Weingarten](B_schur_weyl.md) | dualité de Schur–Weyl, fonctions de Weingarten, théorème de Wigner–Eckart, mesure de Haar sur $O(k)$ |
| C | [Densités de valeurs propres et formule de Page](C_eigenvalue_densities.md) | formule d'intégration de Weyl, GOE/GUE, densités de Wishart et densités induites, formule de Page |
| D | [Générateurs des semi-groupes dynamiques quantiques](D_gksl.md) | le théorème GKSL, générateurs sur les sommes directes, théorème de Frigerio, Bény–Kempf–Kribs |
| E | [Classification des algèbres de Jordan euclidiennes](E_jordan_classification.md) | décomposition de Peirce, théorèmes de Hurwitz et de Frobenius, coordinatisation, l'algèbre d'Albert, le théorème de Jordan–von Neumann–Wigner |
| F | [Cônes symétriques, isomorphismes d'ordre et applications de Jordan](F_cones_and_order.md) | cônes des carrés, le théorème de Koecher–Vinberg, le théorème de Kadison, automorphismes de Jordan de $\bigoplus_cM_{k_c}$ |

Les mathématiques — définitions, théorèmes et démonstrations complètes — se trouvent dans les
[notes de cours](../../lecture_notes/fr/README.md), avec les prérequis, les lectures complémentaires et l'article qui accompagne
la bibliothèque.
