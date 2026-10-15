[English](../00_overview.md) · [Español](../es/00_overview.md) · **Français** · [Deutsch](../de/00_overview.md) · [中文](../zh/00_overview.md) · [日本語](../ja/00_overview.md) · [Русский](../ru/00_overview.md)

# 0. L'objet de ce cours

Les algèbres de von Neumann s'étudient surtout avec un crayon et du papier. En dimension finie,
chacune d'elles est $`*`$-isomorphe à une somme directe d'algèbres de matrices pleines
(théorème 2.6) ; tout est donc calculable en principe, mais lorsqu'on traite numériquement une
telle algèbre, on l'aplatit souvent en une seule grande matrice, et sa structure (les secteurs,
le centre, les nombreuses traces) disparaît du calcul. Ce cours et la bibliothèque
`torch_vn_algebra` adoptent la démarche inverse : **l'algèbre elle-même est le type de
données**, et les calculs y sont menés efficacement.

<a id="s-0-1"></a>

## 0.1 L'observation

Une algèbre de von Neumann $`\mathcal M`$ sur $`\mathbb C^n`$ est, à un changement de base unitaire près,
```math
\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c}
```
(théorème 2.6). En tant qu'algèbre abstraite, c'est $`\bigoplus_cM_{k_c}(\mathbb C)`$ : un élément est
une liste de $`C`$ blocs carrés. C'est exactement la forme d'un tenseur muni d'un *axe des
secteurs* : un élément de $`\mathcal M`$ est un tableau de forme $`(C,k,k)`$ (les blocs étant
complétés par des zéros jusqu'à $`k=\max_ck_c`$ ; c'est ainsi que la bibliothèque le stocke,
section 2.5), et $`B`$ éléments (échantillons de Monte-Carlo, points de départ d'une
optimisation, pas de temps, réalisations du désordre) forment un tableau de forme $`(B,C,k,k)`$.
Les bibliothèques modernes de tableaux (PyTorch sur CPU et GPU) fournissent des routines par
lots — produits matriciels, diagonalisations, décompositions QR — précisément pour de tels
« lots de petites matrices », dont se servent aussi les réseaux de neurones.

*Dans la bibliothèque :* un élément de $`\mathcal M`$ est un opérateur d'une `TypeIAlgebra`, stocké sous la forme d'un unique tenseur de forme $`(B,C,k_{\max},k_{\max})`$ — [tutoriel, section 1.5](../../tutorial/fr/01_matrix_algebras.md#s-1-5) et [tutoriel, section 2.5](../../tutorial/fr/02_von_neumann_algebras.md#s-2-5).

Grâce à cette identification, les opérations de l'algèbre deviennent des opérations par lots sur
des tableaux. La correspondance est exacte (chacune des six premières lignes du tableau est une
identité, démontrée ci-dessous à partir de résultats du cours) ; seules les routines numériques
qui exécutent les opérations travaillent en virgule flottante :

| dans l'algèbre $`\mathcal M`$ | sur le tableau de forme $`(B,C,k,k)`$ |
|---|---|
| produit $`AB`$, somme, adjoint | produit matriciel par lots, somme, transposée conjuguée sur $`(B,C)`$ |
| calcul fonctionnel $`f(A)`$ (chapitre 1) | diagonalisation par lots des blocs |
| une trace $`\sum_cw_c\mathrm{Tr}\,A_c`$ (chapitre 3) | somme diagonale de chaque bloc, somme pondérée le long de l'axe des secteurs |
| le centre, les projecteurs centraux (chapitre 2) | un nombre par secteur |
| espérance conditionnelle sur le centre | trace de chaque bloc divisée par $`k_c`$, rediffusée sur le bloc |
| mesure de Haar sur le groupe unitaire $`\prod_cU(k_c)`$ | décomposition QR par lots de blocs gaussiens, avec la correction de phase du théorème 5.3 |
| états, canaux, dynamique (chapitres 3, 6, 8) | construits à partir des opérations ci-dessus |

*Démonstrations des six premières lignes.* D'après le théorème 2.6, il existe un unitaire $`V`$ tel que
$`\mathcal M=V\big(\bigoplus_cM_{k_c}\otimes1_{m_c}\big)V^*`$, et
$`\varphi\big(V(\bigoplus_cA_c\otimes1_{m_c})V^*\big)=\bigoplus_cA_c`$ est un $`*`$-isomorphisme de $`\mathcal M`$ sur
$`\bigoplus_cM_{k_c}(\mathbb C)`$ (bien défini, car $`A_c\otimes1_{m_c}`$ détermine $`A_c`$ ; multiplicatif, car
$`(A_c\otimes1)(B_c\otimes1)=A_cB_c\otimes1`$). La bibliothèque stocke $`\varphi(A)`$ (section 2.5).
(1) Dans $`\bigoplus_cM_{k_c}`$, les sommes, produits et adjoints se calculent bloc par bloc.
(2) Soient $`A\in\mathcal M`$ normal, $`\varphi(A)=\bigoplus_cA_c`$ et $`f`$ une fonction sur $`\mathrm{spec}\,A`$. Chaque $`A_c`$ est
normal, puisque $`\varphi(A)^*\varphi(A)=\varphi(A^*A)=\varphi(AA^*)=\varphi(A)\varphi(A)^*`$ a lieu bloc par bloc. Si
$`A_cx=\lambda x`$ avec $`x\neq0`$, alors $`V(0\oplus\dots\oplus x\otimes e_1\oplus\dots\oplus0)`$ ($`x\otimes e_1`$ dans le
$`c`$-ième terme) est un vecteur propre de $`A`$ pour $`\lambda`$, donc $`\mathrm{spec}\,A_c\subseteq\mathrm{spec}\,A`$. Soit $`p`$ un
polynôme tel que $`p=f`$ sur $`\mathrm{spec}\,A`$ (interpolation de Lagrange, comme dans la démonstration du théorème 1.5). Pour tout
$`X`$ normal avec $`\mathrm{spec}\,X\subseteq\mathrm{spec}\,A`$ et $`X=U\Lambda U^*`$, on a $`f(X)=Uf(\Lambda)U^*=Up(\Lambda)U^*=p(X)`$ par
le théorème 1.5(b) ; donc $`f(A)=p(A)`$, $`f(A_c)=p(A_c)`$, et comme $`\varphi`$ est un homomorphisme unifère,
```math
\varphi(f(A))=p(\varphi(A))=\bigoplus_cp(A_c)=\bigoplus_cf(A_c),\qquad f(A_c)=U_cf(\Lambda_c)U_c^*
```
pour toute diagonalisation $`A_c=U_c\Lambda_cU_c^*`$ (théorème 1.3) : une diagonalisation par bloc.
(3) Les traces de $`\bigoplus_cM_{k_c}`$ sont les fonctionnelles $`\sum_cw_c\mathrm{Tr}\,A_c`$ (théorème 3.1).
(4) $`Z(\mathcal M)`$ est engendré par les projecteurs centraux minimaux $`P_c`$, et $`\varphi(P_c)`$ vaut $`1_{k_c}`$ dans le bloc
$`c`$ et $`0`$ ailleurs (théorème 2.6), donc $`\varphi(Z(\mathcal M))=\{\bigoplus_cz_c1_{k_c}:z\in\mathbb C^C\}`$.
(5) Proposition 2.8. (6) Corollaire 5.4 (blocs indépendants, chacun de loi de Haar sur $`U(k_c)`$) et théorème 5.3
(chaque bloc issu de la décomposition QR d'une matrice gaussienne). La dernière ligne n'est qu'un renvoi : les états
sont représentés par des densités (théorème 3.2), les canaux entre secteurs par des opérateurs de Kraus entre blocs
(proposition 6.4), et la dynamique fermée par $`e^{-itH}=\bigoplus_ce^{-itH_c}`$ (proposition 8.2).

*Dans la bibliothèque :* les produits et les fonctions comme `sqrt` et `apply_function` — [tutoriel, section 1.5](../../tutorial/fr/01_matrix_algebras.md#s-1-5) ; les traces, par exemple `Tr_blunt` — [tutoriel, section 3.1](../../tutorial/fr/03_traces_and_states.md#s-3-1) ; le centre, `TypeIAlgebra.central` — [tutoriel, section 2.5](../../tutorial/fr/02_von_neumann_algebras.md#s-2-5) ; l'espérance conditionnelle sur le centre, `channels.center_expectation` — [tutoriel, section 2.6](../../tutorial/fr/02_von_neumann_algebras.md#s-2-6) ; les unitaires de Haar, `random_unitary_operator` — [tutoriel, section 5.1](../../tutorial/fr/05_random_matrices.md#s-5-1).

*Calcul :* la ligne 2 du tableau sur $`\bigoplus_{c=1}^{64}M_{16}(\mathbb C)`$ — [tutoriel, section 0.2](../../tutorial/fr/00_overview.md#s-0-2) ; ce que coûte le calcul bloc par bloc — [tutoriel, section 0.1](../../tutorial/fr/00_overview.md#s-0-1).

<a id="s-0-2"></a>

## 0.2 Un premier calcul

*Calcul :* la racine carrée et la trace d'un élément positif aléatoire de $`\bigoplus_{c=1}^{64}M_{16}(\mathbb C)`$, calculées bloc par bloc et dans $`M_{1024}(\mathbb C)`$ — [tutoriel, section 0.2](../../tutorial/fr/00_overview.md#s-0-2).

<a id="s-0-3"></a>

## 0.3 Ce qui existe, et ce qui est nouveau

*Dans la bibliothèque :* les logiciels apparentés, et ce que la bibliothèque y ajoute — [tutoriel, section 0.3](../../tutorial/fr/00_overview.md#s-0-3).

<a id="s-0-4"></a>

## 0.4 Comment lire ce cours

| chapitres | contenu |
|---|---|
| 1–3 | les mathématiques : algèbres de matrices, théorème de structure, traces, états et entropie, avec des démonstrations complètes |
| 4–8 | d'où viennent ces algèbres et ce que l'on y calcule : symétries, matrices aléatoires, canaux, produits tensoriels, dynamique |
| 9–12 | applications : chaînes de spins, calcul quantique, optimisation, apprentissage automatique équivariant |
| 13 | fondements : algèbres de Jordan, et pourquoi la mécanique quantique utilise des *-algèbres |
| annexes | démonstrations complètes des résultats auxiliaires que les chapitres citent par leur numéro (par exemple le théorème A.24) |

Chacun des chapitres 1 à 13 énonce des résultats avec des démonstrations complètes. Les prérequis (algèbre
linéaire, analyse réelle, théorie de la mesure, probabilités élémentaires) sont énumérés dans le
[README](../../lecture_notes/fr/README.md#prerequisites) ; tout le reste est démontré dans les notes de cours, et les résultats
d'autres chapitres et des annexes sont cités par leur numéro. Les encadrés intitulés *Interprétation physique*
relèvent de l'interprétation physique et ne sont pas des théorèmes des notes de cours.

Chaque chapitre a un jumeau dans le [tutoriel](../../tutorial/fr/00_overview.md), qui présente le code de la bibliothèque, le confronte numériquement aux résultats démontrés ici et indique les outils standard auxquels le comparer. Les lignes *Calcul :* et *Dans la bibliothèque :* placées après un résultat renvoient à l'endroit du tutoriel où il est calculé ou utilisé.

**Suite :** [chapitre 1](01_matrix_algebras.md).
