[English](../README.md) · [Español](../es/README.md) · **Français** · [Deutsch](../de/README.md) · [中文](../zh/README.md) · [日本語](../ja/README.md) · [Русский](../ru/README.md)

# Algèbres de von Neumann de dimension finie : un cours avec `torch_vn_algebra`

Ce cours présente la théorie des algèbres de von Neumann de dimension finie — les algèbres
$\bigoplus_cM_{k_c}(\mathbb C)$ — ainsi que leurs usages en physique quantique, en
information quantique, en théorie des matrices aléatoires et en apprentissage automatique,
et montre à chaque étape comment calculer avec elles dans `torch_vn_algebra`.

Il s'adresse simultanément à deux publics :

- **Les mathématiciens** y trouveront une démonstration complète de chaque résultat utilisé.
  Rien n'est admis sans preuve, hormis l'algèbre linéaire classique, le théorème de
  d'Alembert-Gauss et des faits élémentaires d'analyse et de probabilités (limites,
  convexité, fonction génératrice des moments d'une gaussienne) ; lorsqu'un résultat est
  énoncé avant d'être démontré, ou seulement appliqué, cela est dit explicitement.
- **Les physiciens et les ingénieurs** y trouveront la signification physique de chaque
  notion, dans des encadrés intitulés *Interprétation physique*, ainsi que du code
  fonctionnel qu'ils pourront adapter.

Les chapitres 1 à 3 suivent le même plan : la théorie avec ses démonstrations ; le code de la
bibliothèque correspondant à chaque résultat ; une vérification de ce code par rapport à un
résultat exact (un théorème, une formule close ou un calcul indépendant) ; une section
*Points de comparaison* qui cite les outils de référence pour la même tâche ; enfin des
exercices corrigés.

Tout le code du cours est exécuté par la suite de tests (`tests/test_tutorial.py`), de
sorte que chaque exemple fonctionne avec la version actuelle de la bibliothèque.

## Table des matières

| | chapitre | résultats principaux |
|---|---|---|
| 0 | [L'objet de ce cours](00_overview.md) | l'algèbre comme type de données : structure par blocs = tenseurs par lots ; coût ; ce qui existe et ce qui est nouveau |
| 1 | [Algèbres de matrices et calcul fonctionnel](01_matrix_algebras.md) | théorème spectral, calcul fonctionnel, décomposition polaire, dualité de la norme trace, un moment du GUE |
| 2 | [Algèbres de von Neumann de dimension finie](02_von_neumann_algebras.md) | théorème du bicommutant, théorème de structure, centre, groupes de symétrie et commutants, espérances conditionnelles |
| 3 | [Traces, états et entropie](03_traces_and_states.md) | classification des traces, densités, entropie des secteurs, inégalité de Klein, monotonie sous l'espérance conditionnelle sur le centre, contraste de Michelson |
| 4 | [Symétries et représentations](04_symmetries.md) | lemme de Schur, secteurs de charge, $U(1)$ et $SU(2)$, systèmes de spins |
| 5 | [Matrices aléatoires dans une algèbre](05_random_matrices.md) | mesure de Haar sur $\prod_cU(k_c)$, calcul de Weingarten, statistique des niveaux et pourquoi il ne faut pas mélanger les secteurs |
| 6 | [Canaux quantiques](06_quantum_channels.md) | Kraus et Stinespring, duaux relativement à une trace, applications entre secteurs, décohérence |
| 7 | [Produits tensoriels et intrication](07_tensor_products.md) | trace partielle, formule de Page, fusion des charges |
| 8 | [Dynamique](08_dynamics.md) | équations de Schrödinger et de Lindblad, sauts entre secteurs, post-sélection |
| 9 | Chaînes de spins *(en préparation)* | secteurs d'aimantation et d'impulsion, Lanczos, vérifications par l'ansatz de Bethe, désordre |
| 10 | [Calcul quantique](10_quantum_computing.md) | sous-espaces sans décohérence, correction d'erreurs par algèbres d'opérateurs, mélangeurs sous contraintes |
| 11 | [Optimisation et inégalités de traces](11_optimisation.md) | paramétrisations sous contraintes, échantillonnage ou optimisation |
| 12 | [Apprentissage automatique équivariant](12_equivariant_learning.md) | les couches équivariantes comme commutants |

## Prérequis

L'algèbre linéaire sur $\mathbb C$ (valeurs propres, produits scalaires, matrices
unitaires), des notions de base en probabilités, et Python. Aucune connaissance en physique
ni en analyse fonctionnelle n'est supposée.

## Exécuter le code

```bash
pip install -e .                # from the root of the repository
python -m pytest tests/test_tutorial.py -q
```

(La première commande s'exécute depuis la racine du dépôt.) Chaque bloc de code d'un
chapitre peut être collé, dans l'ordre, dans une même session Python ou un même notebook
Jupyter. Les blocs d'un même chapitre partagent leurs variables.

## Pour aller plus loin

M. Takesaki, *Theory of Operator Algebras I* (Springer, 2002) ; V. Jones,
*von Neumann algebras*, notes de cours ; D. Petz, *Quantum Information Theory and Quantum
Statistics* (Springer, 2008) ; J. Watrous, *The Theory of Quantum Information* (Cambridge,
2018) ; M. L. Mehta, *Random Matrices* (Elsevier, 2004). Pour le contexte physique dans un
langage mathématique, voir [La physique pour les mathématiciens](../../physics_for_mathematicians.fr.md).
L'article qui accompagne la bibliothèque est [`paper/main.pdf`](https://github.com/andre-a-no/von_neumann_type_i/blob/main/paper/main.pdf).
