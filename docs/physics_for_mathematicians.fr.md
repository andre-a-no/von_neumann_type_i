[English](physics_for_mathematicians.md) · [Español](physics_for_mathematicians.es.md) · **Français** · [Deutsch](physics_for_mathematicians.de.md) · [中文](physics_for_mathematicians.zh.md) · [日本語](physics_for_mathematicians.ja.md) · [Русский](physics_for_mathematicians.ru.md)

# La physique pour les mathématiciens : un dictionnaire pour `torch_vn_algebra`

*Traduction de la version anglaise ; en cas de divergence, la version anglaise fait foi.*

Ce guide s'adresse aux lecteurs venant des algèbres d'opérateurs qui souhaitent comprendre ce que
signifient les modèles physiques de la bibliothèque et de l'article. Tout est formulé en termes de
l'algèbre M = ⊕_c M_{k_c}(C) agissant sur H = ⊕_c C^{k_c}. Pour chaque notion physique, nous indiquons
ce qu'elle est dans M et où elle se trouve dans la bibliothèque.

---

## 0. L'idée principale

Une symétrie physique est une représentation unitaire π d'un groupe G sur H. Tout ce qui est compatible
avec la symétrie (le hamiltonien, l'évolution temporelle, les observables soumises à une règle de
supersélection) appartient au **commutant** π(G)'. D'après le théorème du bicommutant et le lemme de Schur,

    H ≅ ⊕_λ V_λ ⊗ C^{m_λ},     π(G)' ≅ ⊕_λ M_{m_λ}(C),

donc π(G)' est une algèbre de von Neumann de type I de dimension finie. Ses secteurs sont les
représentations irréductibles λ, les tailles des blocs sont les multiplicités m_λ, et le centre est formé
des fonctions des « charges ». Un physicien qui modélise un système doté d'une symétrie calcule donc dans
une algèbre de la forme (1), en général sans la désigner ainsi. Pour une seule quantité conservée Q (un
opérateur autoadjoint), les secteurs sont les valeurs propres de Q et les blocs sont ses sous-espaces propres.

Cela répond aussi à la question « pourquoi ne pas utiliser une seule grande matrice ? ». Le plongement
M ⊂ M_N (N = Σk_c) est exact, mais on y perd
- **le coût** : Σk_c³ au lieu de N³ ; pour C blocs égaux, c'est un facteur C² en temps ;
- **la mesure** : le groupe unitaire de M est Π U(k_c), un ensemble de mesure nulle dans U(N). Une matrice
  unitaire tirée selon la mesure de Haar sur U(N) n'appartient pas à M, et un pas de gradient dans M_N sort de M ;
- **les traces** : sur M_N, la trace est unique à un facteur près ; sur M, les traces fidèles forment la
  famille à C paramètres Σ w_c Tr_c. La « préservation de la trace » et les applications duales dépendent des poids ;
- **le centre** : les projecteurs de secteur constituent la partie classique du système (section 2).
  Les probabilités conditionnelles et la post-sélection sont des opérations sur le centre ; avec une seule
  matrice, il faut les reconstruire à partir des indices de blocs.

Le cas C = 1 est également pris en charge : il s'agit alors d'algèbre linéaire par lots dans M_n.

---

## 1. États, observables, mesures

| physique | mathématiques | dans la bibliothèque |
|---|---|---|
| observable (énergie, spin, nombre de particules) | A ∈ M autoadjoint | `Operator` |
| résultats de mesure possibles | spectre de A | `eigh` |
| état (mixte) | fonctionnelle positive normale ω, ω(1) = 1 ; densité ρ ≥ 0, Tr ρ = 1 | `DensityMatrix` |
| état pur, « fonction d'onde » | vecteur ψ ∈ H, ρ = \|ψ⟩⟨ψ\| (un projecteur minimal) | `vector_in_sector`, `basis_state` |
| valeur moyenne | ω(A) = Tr(ρA) | `expectation` |
| probabilité du résultat λ | ω(P_λ), où P_λ est le projecteur spectral | `apply_function` |
| état après une mesure (règle de Lüders) | ρ ↦ PρP / Tr(Pρ) | `lueders_update`, `condition_on` |
| entropie de von Neumann | S(ρ) = −Tr ρ log ρ | `entropy` |
| température, état thermique | ρ = e^{−βH}/Tr e^{−βH} (état de Gibbs = état KMS) | `gibbs_state` |
| température infinie | état tracial τ / τ(1) | `tracial_state` |

**Signification physique des trois traces.**
- Tr_blunt (poids 1) est la trace physique sur H. Normalisée, elle donne l'état à température infinie :
  les N états de base sont équiprobables, et le secteur c a la probabilité k_c/N.
- Tr_norm (poids 1/k_c) : l'état uniforme (microcanonique) à l'intérieur de chaque secteur pris séparément.
- τ_vN (poids 1/(C k_c)) : le mélange uniforme des états microcanoniques des secteurs. Chaque
  *secteur* (et non chaque état de base) a la probabilité 1/C.

La restriction d'un état quelconque au centre est une loi de probabilité classique sur les secteurs
(`sector_probabilities`). C'est ce que les physiciens appellent « la probabilité de la charge q ».

---

## 2. Supersélection, décohérence, classique et quantique

- **Règle de supersélection** (Wick–Wightman–Wigner) : les superpositions d'états de charges différentes
  (charge électrique, parité fermionique) ne sont pas observables. Mathématiquement, les observables
  admissibles forment le commutant de la charge, c'est-à-dire l'algèbre M et non B(H) tout entier.
- **Le centre de M** est formé des grandeurs que l'on peut mesurer sans perturber le système et qui
  commutent avec tout : ce sont les variables « classiques ». Si tous les k_c = 1, le système est purement
  classique (une algèbre commutative, c'est-à-dire un espace probabilisé à C points). Si C = 1, il est
  purement quantique. Le cas général est hybride.
- **Décohérence / einsélection** (Zurek) : l'interaction avec un environnement supprime les blocs
  hors diagonale ρ_{cd}, c ≠ d. C'est l'espérance conditionnelle E: B(H) → M, ρ ↦ Σ P_c ρ P_c, et
  l'algèbre effective des observables se réduit à M. Dans la bibliothèque, les opérateurs sont déjà dans M :
  le résultat de cette espérance conditionnelle est le point de départ des calculs.

---

## 3. Dynamique des systèmes fermés

| physique | mathématiques | dans la bibliothèque |
|---|---|---|
| hamiltonien | H ∈ M autoadjoint | — |
| équation de Schrödinger i dψ/dt = Hψ | ψ(t) = e^{−itH}ψ(0) | `dynamics.schrodinger`, `propagator` |
| équation de von Neumann dρ/dt = −i[H, ρ] | ρ(t) = U ρ U*, U = e^{−itH} | `von_neumann` |
| représentation de Heisenberg | A(t) = U* A U (un automorphisme de M) | `solve_operator_ode` |
| quantité conservée | Q tel que [H, Q] = 0 | secteurs |

Si H ∈ M, alors e^{−itH} ∈ M, et la dynamique ne quitte jamais les secteurs. Chaque secteur évolue
indépendamment, et l'exponentielle se calcule bloc par bloc.

---

## 4. Systèmes ouverts et canaux

- Un **canal quantique** est une application complètement positive qui préserve la trace (représentation
  de Schrödinger). Son adjoint (représentation de Heisenberg) est complètement positif et préserve l’unité. Forme de
  Kraus : Φ(ρ) = Σ_i K_i ρ K_i*. Dans la bibliothèque : `Channel`.
- Un **canal entre secteurs** (`InterSectorChannel`) décrit un processus qui modifie la charge :
  perte de particules, désintégration, sortie d'un électron par effet tunnel. Ses blocs de Kraus sont
  K^{dc}: C^{k_c} → C^{m_d}. La préservation de la trace et l'application duale dépendent des poids de
  trace en entrée et en sortie, d'où les facteurs w_out/w_in.
- L'**équation de Lindblad** (GKSL) engendre un semi-groupe de canaux e^{tL} :
  dρ/dt = −i[H, ρ] + Σ_j γ_j (L_j ρ L_j* − ½{L_j* L_j, ρ}).
  Les L_j sont les « opérateurs de saut » ; chacun décrit une manière dont l'environnement agit sur le système.
  - **T1** (relaxation d'énergie, amortissement d'amplitude) : L = σ⁻ = |0⟩⟨1|, un qubit retombe de l'état
    excité vers l'état fondamental. Le nombre d'excitations diminue d'une unité : c'est donc un saut **entre**
    secteurs. Exemple : `site_amplitude_damping`, `lowering`.
  - **T2** (déphasage) : L = σ^z. Il est diagonal et ne change pas le secteur ; il ne fait que détruire
    les phases. C'est un saut **à l'intérieur** d'un secteur.
  - T1 et T2 sont les temps de cohérence de référence indiqués pour tout ordinateur quantique.

---

## 5. Chaînes de spins

- **Spin ½ sur un site** : H_site = C², base |↑⟩, |↓⟩ (ou |1⟩, |0⟩). Matrices de Pauli σ^x, σ^y, σ^z ;
  S^a = σ^a/2 ; S^± = S^x ± iS^y sont les opérateurs de montée et de descente.
- **Une chaîne de L spins** : H = (C²)^{⊗L}, dim = 2^L. C'est cette croissance exponentielle qui oblige
  les physiciens à exploiter les symétries.
- **Hamiltonien XXZ** : H = Σ_⟨ij⟩ [J/2 (S_i⁺S_j⁻ + S_i⁻S_j⁺) + JΔ S_i^z S_j^z] + Σ_i h_i S_i^z.
  Le premier terme échange des ↑↓ voisins (un « saut de particule »), le deuxième décrit leur interaction,
  le troisième un champ extérieur. Δ = 1 correspond au modèle de Heisenberg (un aimant) ; Δ = 0 équivaut à
  des fermions libres (par la transformation de Jordan–Wigner), ce qui fournit des vérifications exactes.
- **Charge conservée** : S^z_tot, ou N = nombre de ↑ (le poids de Hamming d'une chaîne de bits). Chaque
  terme de H la conserve, donc H ∈ ⊕_{N=0}^{L} M_{C(L,N)} : les blocs ont des tailles binomiales.
- **Fusion des secteurs** (`fused_tensor_product`) : par produit tensoriel, les charges s'additionnent.
  C'est la règle de fusion de U(1) : (⊕_q M_{n_q}) ⊗ (⊕_r M_{m_r}) ⊂ ⊕_s M_{Σ_{q+r=s} n_q m_r}.
  `SpinChain` construit sa base de cette manière, site par site.
- **Impulsion** (`momentum_algebra`) : un anneau possède en outre la symétrie de translation T (le groupe
  Z_L). Les sous-espaces propres communs (N, k), où e^{ik} est une valeur propre de T, donnent une
  décomposition plus fine, avec des blocs environ L fois plus petits.
- L'**état fondamental** est le vecteur propre de H associé à la plus petite valeur propre ; il décrit le
  comportement à température nulle. Pour les grands blocs, on le calcule par la méthode de Lanczos (`krylov`).
- Ce que l'article vérifie, et pourquoi ce sont des tests classiques :
  - les énergies de l'état fondamental de l'anneau de Heisenberg (valeurs exactes connues, ansatz de Bethe) ;
  - la **règle des signes de Marshall** : pour le modèle de Heisenberg sur un réseau biparti, l'état
    fondamental a les signes (−1)^{nombre de ↑ sur un sous-réseau}. C'est un théorème, ce qui en fait un test sans ambiguïté ;
  - **emboîtement SU(2)** : le modèle de Heisenberg possède la symétrie SU(2) complète, donc pour N < L/2 le
    spectre du secteur N est contenu dans celui du secteur N+1 (multiplets) ;
  - **fonte d'une paroi de domaine** : |↑…↑↓…↓⟩ à Δ = 0 s'étale selon un profil connu ;
  - **formule de Peschel** : l'entropie d'intrication de fermions libres obtenue à partir de la matrice de corrélation.
- **Désordre et chaos** (`random_field_heisenberg`, `level_spacing_ratio`) : champs aléatoires
  h_i ∈ [−W, W]. Pour W petit, le système est « chaotique » et les niveaux se repoussent comme dans le GOE,
  ⟨r⟩ ≈ 0.53. Pour W grand, il se localise (localisation à N corps, MBL) et les niveaux sont indépendants
  comme pour un processus de Poisson, ⟨r⟩ ≈ 0.386. Les statistiques de niveaux n'ont de sens qu'**à
  l'intérieur d'un seul secteur** : mélanger les secteurs produit artificiellement une statistique de
  Poisson. C'est un bon argument physique pour travailler dans l'algèbre avec ses secteurs.

---

## 6. Intrication et sous-systèmes

- Un système composé A+B : H_A ⊗ H_B, algèbre M_A ⊗ M_B.
- L'**état réduit** (matrice densité réduite) ρ_A = Tr_B ρ est la restriction de la fonctionnelle ω à la
  sous-algèbre M_A ⊗ 1. En termes de densités, c'est une espérance conditionnelle (la trace partielle).
  Dans la bibliothèque : `partial_trace`, `reduced_state`.
- L'**entropie d'intrication** S(ρ_A) d'un état pur ρ mesure la corrélation quantique entre A et B.
- **Formule de Page** : l'entropie moyenne exacte d'un sous-système d'un état pur aléatoire. Elle permet de
  vérifier que l'échantillonneur suit bien la mesure de Haar.
- **Intrication résolue en symétrie** : si la charge est conservée, ρ_A est diagonale par blocs et son
  entropie se décompose en une partie « classique » (l'entropie de la distribution sur le centre) et une
  partie « quantique » (l'entropie moyenne des blocs). C'est une utilisation directe du centre.

---

## 7. Calcul quantique

- **Qubit** = C², n qubits = (C²)^{⊗n} ; une **porte** est un unitaire ; un **circuit** est un produit de
  portes ; une **mesure** dans la base de calcul utilise les projecteurs spectraux de σ^z.
- **QAOA** (*quantum approximate / alternating operator ansatz*) : on alterne e^{−iγ H_C} (H_C est
  diagonal et code la fonction objectif) et e^{−iβ H_M} (le « mélangeur », *mixer*). Pour les problèmes
  soumis à la contrainte « exactement k uns » (choisir k éléments parmi n), on utilise le **mélangeur XY**
  Σ(X_iX_j + Y_iY_j) = 2Σ(σ⁺σ⁻ + h.c.). Il conserve le poids de Hamming, donc tout le circuit reste dans
  ⊕_N M_{C(n,N)} et les solutions admissibles forment exactement le secteur N = k.
- **Bruit** : T1 fait sortir du secteur admissible (un 1 est perdu), T2 non.
- **Post-sélection** : on mesure le poids de Hamming et on rejette les exécutions erronées. C'est la règle
  de Lüders avec un projecteur **central** P_k : ρ ↦ P_kρP_k / Tr(P_kρ), avec la probabilité de succès
  ω(P_k). C'est ce que fait `ex9_qaoa_xy_mixer.py`. Pour un bruit purement T1, la post-sélection restitue
  exactement l'état idéal (F = 1), car chaque saut T1 fait sortir du secteur.
- Sous-espaces sans décohérence et correction d'erreurs par algèbres d'opérateurs : l'information protégée
  est codée dans un bloc (ou un facteur tensoriel d'un bloc) d'une algèbre que le bruit n'affecte pas. Là
  encore, c'est le langage des algèbres de type I (Knill–Laflamme–Viola, Bény–Kempf–Kribs).

---

## 8. Matrices aléatoires, réelles et complexes

Cette section est directement liée au paramètre `complex_valued`.
- Le renversement du temps T est une **symétrie antiunitaire**. Si T² = 1 (particules sans spin, pas de
  champ magnétique), H est **réel** dans une base adaptée : GOE, ensembles orthogonaux, COE. L'arithmétique
  réelle suffit ; elle divise la mémoire par deux et est plus rapide (voir le tableau des temps de calcul
  dans l'article).
- En l'absence d'une telle symétrie (champ magnétique, phases complexes, impulsion k ≠ 0, π), H est complexe :
  GUE, CUE. Les nombres complexes sont alors indispensables.
- Si T² = −1 (spin demi-entier avec couplage spin–orbite), on est dans le cas quaternionique : GSE, CSE.
- Le modèle XXZ avec des champs réels est réel. Les secteurs d'impulsion avec k ∉ {0, π} sont complexes. La
  dynamique e^{−itH} est toujours complexe, même pour H réel. C'est le tableau « ce qui suffit, et quand »
  de l'article.
- La double précision est nécessaire lorsque de petites différences comptent : spectres dégénérés,
  identités vérifiées à 1e−12 près, évolutions sur des temps longs, Lanczos avec de nombreuses itérations.
  La simple précision (et le TF32 sur les cœurs tensoriels) suffit pour des statistiques sur de nombreux
  échantillons.

---

## 9. Apprentissage automatique

Une couche linéaire équivariante est une application W: V → V qui commute avec une représentation d'un
groupe (les translations pour les convolutions, les permutations pour les réseaux de graphes, les rotations
pour les molécules). D'après le lemme de Schur, l'espace de ces W est ⊕_λ M_{m_λ}, avec les multiplicités
m_λ ; une couche équivariante est donc paramétrée par un élément d'une algèbre de la forme (1). La
bibliothèque fournit, pour de telles paramétrisations, une initialisation aléatoire selon la bonne mesure,
des contraintes spectrales et des gradients.

---

## 10. Périmètre : ce qu'est la bibliothèque, et ce qu'elle n'est pas

- La bibliothèque ne produit **pas** de physique nouvelle. C'est un outil validé pour le calcul dans des
  algèbres à secteurs. Les modèles physiques de l'article sont des systèmes connus, aux résultats connus,
  qui servent à vérifier l'outil.
- Pour un seul grand système (un unique état fondamental à L = 30, réseaux de tenseurs, simulation de
  grands circuits), les logiciels spécialisés (QuSpin, ITensor/TeNPy, Qiskit) sont plus appropriés.
- Son créneau, ce sont de nombreux opérateurs de taille modérée à structure de secteurs explicite :
  ensembles désordonnés, statistiques sur des états et des canaux aléatoires, opérateurs extrémaux sous
  contraintes spectrales, canaux entre secteurs avec des traces cohérentes, le tout avec différentiation
  automatique et calcul sur GPU.
- La question mathématique d'origine (inégalités de trace, contraste de Michelson) est un exemple de
  problème extrémal où l'échantillonnage Monte-Carlo échoue et où l'optimisation sous contraintes réussit.

---

## Pour aller plus loin

- M. A. Nielsen, I. L. Chuang, *Quantum Computation and Quantum Information* — chap. 2 (postulats),
  chap. 8 (canaux, T1/T2).
- M. M. Wilde, *Quantum Information Theory* — canaux et états, dans un esprit proche des algèbres d'opérateurs.
- A. W. Sandvik, "Computational studies of quantum spin systems", AIP Conf. Proc. 1297 (2010) —
  secteurs et diagonalisation exacte en pratique.
- H.-P. Breuer, F. Petruccione, *The Theory of Open Quantum Systems* — l'équation de Lindblad.
- F. Haake, *Quantum Signatures of Chaos* — matrices aléatoires, renversement du temps, GOE/GUE/GSE.
- S. D. Bartlett, T. Rudolph, R. W. Spekkens, Rev. Mod. Phys. 79, 555 (2007) — supersélection et
  référentiels, accessible aux mathématiciens.
