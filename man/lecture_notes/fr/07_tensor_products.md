[English](../07_tensor_products.md) · [Español](../es/07_tensor_products.md) · **Français** · [Deutsch](../de/07_tensor_products.md) · [中文](../zh/07_tensor_products.md) · [日本語](../ja/07_tensor_products.md) · [Русский](../ru/07_tensor_products.md)

# 7. Produits tensoriels et intrication

Ce chapitre étudie le produit tensoriel de deux algèbres à secteurs (physiquement, l'algèbre de
deux systèmes pris ensemble ; voir les encadrés *Interprétation physique*). En présence de
secteurs, il y a deux produits : le produit simple, où chaque partie conserve sa propre charge, et
le produit *fusionné*, où seule la charge totale est conservée. Nous établissons la structure des
deux, définissons la trace partielle comme la duale d'un plongement et étudions l'intrication.

**Vous apprendrez** comment les secteurs se multiplient,
pourquoi la trace partielle est un canal et, une fois normalisée, une espérance conditionnelle
d'indice fini, quel est l'état réduit d'un état qui conserve la charge, comment l'intrication se
décompose en une partie de nombre et une partie configurationnelle, ainsi que la pureté et
l'entropie moyennes des états aléatoires. Les vérifications numériques de ce chapitre se trouvent
dans le [Tutoriel 7](../../tutorial/fr/07_tensor_products.md).

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

Pour des $`*`$-algèbres $`\mathcal M_1`$ sur $`H_1`$ et $`\mathcal M_2`$ sur $`H_2`$, le **produit tensoriel**
$`\mathcal M_1\otimes\mathcal M_2`$ est le sous-espace engendré par les opérateurs $`A\otimes B`$ sur
$`H_1\otimes H_2`$, avec $`A\in\mathcal M_1`$, $`B\in\mathcal M_2`$. Comme
$`(A\otimes B)(A'\otimes B')=AA'\otimes BB'`$, $`(A\otimes B)^*=A^*\otimes B^*`$ et $`1\otimes1=1`$, c'est
une $`*`$-algèbre, donc une algèbre de von Neumann (théorème 2.2).

<a id="r-7-1"></a>

**Théorème 7.1.** (a) $`H_1\otimes H_2=\bigoplus_{c,d}\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$, somme
orthogonale des images des projecteurs $`P_c\otimes Q_d`$, et
```math
\mathcal M_1\otimes\mathcal M_2=\bigoplus_{c,d}M_{k_c}\otimes M_{l_d}\cong\bigoplus_{c,d}M_{k_cl_d}(\mathbb C):
```
elle se compose de tous les opérateurs $`X=\sum_{c,d}(P_c\otimes Q_d)X(P_c\otimes Q_d)`$, avec des
blocs arbitraires sur les $`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$. (b) Ses projecteurs centraux minimaux
sont les $`C_1C_2`$ projecteurs $`P_c\otimes Q_d`$, et
$`\dim(\mathcal M_1\otimes\mathcal M_2)=\dim\mathcal M_1\dim\mathcal M_2`$.
(c) Pour des traces $`\mathrm{tr}_w`$ sur $`\mathcal M_1`$ et $`\mathrm{tr}_v`$ sur $`\mathcal M_2`$ (théorème 3.1),
$`\mathrm{tr}_w(A)\mathrm{tr}_v(B)=\mathrm{tr}_{w\otimes v}(A\otimes B)`$ avec
$`(w\otimes v)_{cd}=w_cv_d`$. En particulier, `Tr_blunt`, `Tr_norm` et `tau_vN` du produit sont les
produits des mêmes traces des facteurs.

*Démonstration.* (a) Les $`P_c\otimes Q_d`$ sont des projecteurs deux à deux orthogonaux de somme
$`1\otimes1`$, et l'image de $`P_c\otimes Q_d`$ est $`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$. Comme $`P_c`$ et
$`Q_d`$ sont centraux, $`(P_c\otimes Q_d)(A\otimes B)=AP_c\otimes BQ_d=(A\otimes B)(P_c\otimes Q_d)`$, donc
tout élément de $`\mathcal M_1\otimes\mathcal M_2`$ est diagonal par blocs. Réciproquement, avec les
unités matricielles $`E^c_{aa'}\in\mathcal M_1`$ et $`F^d_{bb'}\in\mathcal M_2`$ des blocs $`c`$ et $`d`$, on a
$`E^c_{aa'}\otimes F^d_{bb'}=|e_a\otimes e_b\rangle\langle e_{a'}\otimes e_{b'}|`$ dans
$`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$ ; ce sont les unités matricielles de $`M_{k_cl_d}(\mathbb C)`$, donc
tout opérateur porté par le bloc $`(c,d)`$ appartient à $`\mathcal M_1\otimes\mathcal M_2`$. (b) Le centre
d'une somme directe d'algèbres de matrices pleines est engendré par les projecteurs sur les termes
de la somme (chapitre 2, exercice 1), et $`\sum_{c,d}k_c^2l_d^2=\sum_ck_c^2\sum_dl_d^2`$. (c)
$`\mathrm{Tr}(X\otimes Y)=\sum_{a,b}X_{aa}Y_{bb}=\mathrm{Tr}\,X\mathrm{Tr}\,Y`$, donc
$`\mathrm{tr}_{w\otimes v}(A\otimes B)=\sum_{c,d}w_cv_d\mathrm{Tr}\,A_c\mathrm{Tr}\,B_d=\mathrm{tr}_w(A)\mathrm{tr}_v(B)`$.
Les poids $`1`$, $`1/k_c`$, $`1/(C_1k_c)`$ des trois traces (section 3.1) se multiplient pour donner $`1`$,
$`1/(k_cl_d)`$, $`1/(C_1C_2k_cl_d)`$, les poids des mêmes traces sur le produit. $`\square`$

*Dans la bibliothèque :* le produit est `tensor_product` et $`A\otimes B`$ est `kron` ; une vérification du théorème 7.1 sur $`(\mathbb C\oplus M_2)\otimes(M_2\oplus M_3)`$ — [tutoriel, section 7.1](../../tutorial/fr/07_tensor_products.md#s-7-1).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Si chacun de deux systèmes obéit à sa propre
> règle de supersélection — deux boîtes qui ne peuvent pas échanger de particules —, le système composé a un
> secteur pour chaque *couple* de charges, et aucune observable ne relie des couples différents. Un
> secteur composé compte avec le produit des poids (théorème 7.1(c)).

<a id="s-7-2"></a>

## 7.2 La trace partielle

**Définition.** Soit $`\rho\in\mathcal M_1\otimes\mathcal M_2`$. La **trace partielle**
$`\mathrm{Tr}_2\,\rho\in\mathcal M_1`$ est l'élément tel que
```math
\mathrm{Tr}\big((\mathrm{Tr}_2\,\rho)\,A\big)=\mathrm{Tr}\big(\rho\,(A\otimes1)\big)\qquad\text{for all }A\in\mathcal M_1 .
```
Elle existe et est unique d'après le théorème 3.2, appliqué à la forme linéaire
$`A\mapsto\mathrm{Tr}(\rho(A\otimes1))`$ et à la trace fidèle $`\mathrm{Tr}`$ de $`\mathcal M_1`$ : la
trace partielle est la duale, relativement aux traces brutes, du plongement $`A\mapsto A\otimes1`$.
$`\mathrm{Tr}_1`$ se définit de la même façon ; pour un état $`\rho`$, $`\mathrm{Tr}_2\,\rho`$ et
$`\mathrm{Tr}_1\,\rho`$ sont les **états réduits**. Pour les algèbres pleines $`\mathcal M_i=B(H_i)`$, c'est la
**trace partielle de l'espace de Hilbert** $`\mathrm{Tr}_{H_2}:B(H_1\otimes H_2)\to B(H_1)`$.

<a id="r-7-2"></a>

**Proposition 7.2.** (a) Pour $`\rho=\bigoplus_{c,d}\rho_{cd}`$,
```math
(\mathrm{Tr}_2\,\rho)_c=\sum_d\mathrm{tr}_{l_d}\,\rho_{cd},\qquad
(\mathrm{tr}_l\,X)_{aa'}=\sum_{b=1}^lX_{(a,b),(a',b)}\quad(X\text{ on }\mathbb C^k\otimes\mathbb C^l),
```
et $`\mathrm{Tr}_2\,\rho=\mathrm{Tr}_{H_2}\,\rho`$, où $`(\mathrm{Tr}_{H_2}\,X)_{ii'}=\sum_jX_{(i,j),(i',j)}`$.
(b) $`\mathrm{Tr}_2(A\otimes B)=(\mathrm{Tr}\,B)\,A`$.
(c) $`\mathrm{Tr}_2\big((A\otimes1)X(A'\otimes1)\big)=A\,(\mathrm{Tr}_2\,X)\,A'`$ pour $`A,A'\in\mathcal M_1`$.
(d) $`\mathrm{Tr}(\mathrm{Tr}_2\,X)=\mathrm{Tr}\,X`$, et $`\mathrm{Tr}_2`$ envoie les états sur des états.

*Démonstration.* (a) Avec $`A=\bigoplus_cA_c`$, on a $`A\otimes1=\bigoplus_{c,d}A_c\otimes1_{l_d}`$ et
```math
\mathrm{Tr}\big(\rho(A\otimes1)\big)=\sum_{c,d}\sum_{a,a',b}(\rho_{cd})_{(a,b),(a',b)}(A_c)_{a'a}
=\sum_c\mathrm{Tr}\Big(\Big(\sum_d\mathrm{tr}_{l_d}\,\rho_{cd}\Big)A_c\Big),
```
donc l'élément de droite vérifie la définition et est égal à $`\mathrm{Tr}_2\,\rho`$. Le même calcul
pour $`B(H_i)`$ donne la formule de $`\mathrm{Tr}_{H_2}`$ ; pour $`\rho`$ diagonale par blocs, les
coefficients $`\rho_{(i,j),(i',j)}`$ avec $`i,i'`$ dans des blocs différents sont nuls, donc
$`\mathrm{Tr}_{H_2}\,\rho`$ a les blocs de (a). (b) $`\mathrm{tr}_l(X\otimes Y)=(\mathrm{Tr}\,Y)X`$. (c) Pour
$`B\in\mathcal M_1`$,
$`\mathrm{Tr}\big(\mathrm{Tr}_2((A\otimes1)X(A'\otimes1))B\big)=\mathrm{Tr}\big(X(A'BA\otimes1)\big) =\mathrm{Tr}\big(A(\mathrm{Tr}_2\,X)A'B\big)`$ ; on conclut par unicité. (d) On prend $`A=1`$. Si $`\rho`$ est
un état, $`A\mapsto\mathrm{Tr}(\rho(A\otimes1))`$ est positive, puisque $`A^*A\otimes1=(A\otimes1)^*(A\otimes1)`$,
et vaut $`1`$ en $`A=1`$, donc sa densité $`\mathrm{Tr}_2\,\rho`$ est un état (théorème 3.2). $`\square`$

*Dans la bibliothèque :* $`\mathrm{Tr}_2`$ est `partial_trace` ; une vérification de la proposition 7.2(a), (b) par comparaison avec `numpy.einsum` — [tutoriel, section 7.2](../../tutorial/fr/07_tensor_products.md#s-7-2).

Les canaux font l'objet du chapitre 6 ; ici, seule la définition suivante est nécessaire : une
application linéaire $`\Phi`$ est **complètement positive** si $`\mathrm{id}_n\otimes\Phi`$ est positive
sur $`M_n(\mathbb C)\otimes(\cdot)`$ pour tout $`n`$, et c'est un **canal** si elle préserve de plus la
trace.

<a id="r-7-3"></a>

**Théorème 7.3.** Faisons parcourir à $`j`$ la base canonique de $`H_2`$ et posons
$`K_j=1\otimes\langle j|:H_1\otimes H_2\to H_1`$, $`K_j(x\otimes y)=\langle j,y\rangle x`$. Alors
```math
\mathrm{Tr}_{H_2}\,X=\sum_jK_jXK_j^*,\qquad\sum_jK_j^*K_j=1,\qquad\sum_jK_j^*AK_j=A\otimes1 .
```
Par conséquent, $`\mathrm{Tr}_2`$ est un canal de $`\mathcal M_1\otimes\mathcal M_2`$ dans $`\mathcal M_1`$, et sa
duale relativement aux traces brutes est $`A\mapsto A\otimes1`$.

*Démonstration.* $`K_j^*x=x\otimes e_j`$, donc $`\langle e_a,K_jXK_j^*e_{a'}\rangle=X_{(a,j),(a',j)}`$ ; la
somme sur $`j`$ donne la proposition 7.2(a). De plus,
$`\sum_jK_j^*K_j(x\otimes y)=\sum_j\langle j,y\rangle x\otimes e_j=x\otimes y`$
et $`\sum_jK_j^*AK_j=\sum_jA\otimes|j\rangle\langle j|=A\otimes1`$. Pour $`Y\ge0`$,
$`(\mathrm{id}_n\otimes\mathrm{Tr}_{H_2})(Y)=\sum_j(1_n\otimes K_j)Y(1_n\otimes K_j)^*\ge0`$, car
$`\langle x,TYT^*x\rangle=\langle T^*x,YT^*x\rangle\ge0`$ ; cela reste vrai sur la sous-algèbre
$`\mathcal M_1\otimes\mathcal M_2`$. Enfin,
$`\mathrm{Tr}\,\sum_jK_jXK_j^*A=\mathrm{Tr}\,X\sum_jK_j^*AK_j`$, ce qui, pour $`A=1`$, est la préservation de
la trace et, en général, la dualité. $`\square`$

*Dans la bibliothèque :* le canal du théorème 7.3 est `composite.partial_trace_channel`, avec sa duale — [tutoriel, section 7.2](../../tutorial/fr/07_tensor_products.md#s-7-2).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). L'état réduit $`\mathrm{Tr}_2\,\rho`$ est tout ce
> qu'apprend un observateur qui ne peut mesurer que le premier système — la définition dit qu'il reproduit toutes
> les valeurs moyennes des observables $`A\otimes1`$. Des états composés ayant les mêmes états réduits
> diffèrent par leurs corrélations.

**La trace partielle comme espérance conditionnelle.** Soit $`\mathcal M_2=M_m(\mathbb C)`$, de sorte
que les secteurs de $`\mathcal M_1\otimes M_m`$ sont ceux de $`\mathcal M_1`$ (théorème 7.1) ; soient
$`\iota(A)=A\otimes1_m`$ et $`\mathcal E=\frac1m\mathrm{Tr}_2:\mathcal M_1\otimes M_m\to\mathcal M_1`$, et
soient $`e_{ab}`$ les unités matricielles de $`M_m`$.

<a id="r-7-4"></a>

**Proposition 7.4.** (a) $`\iota\circ\mathcal E`$ est une espérance conditionnelle de
$`\mathcal M_1\otimes M_m`$ sur $`\mathcal M_1\otimes1`$ (section 2.6) et, pour toute trace fidèle
$`\mathrm{tr}_w`$ de $`\mathcal M_1\otimes M_m`$, c'est la seule qui préserve $`\mathrm{tr}_w`$. (b) Les
$`m^2`$ éléments $`u_{ab}=\sqrt m\,(1\otimes e_{ab})`$ vérifient
$`X=\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*X))`$ pour tout $`X`$, et
$`\sum_{a,b}u_{ab}u_{ab}^*=m^2\,1`$. (c) Toute famille $`(v_k)`$ telle que
$`X=\sum_kv_k\,\iota(\mathcal E(v_k^*X))`$ pour tout $`X`$ vérifie $`\sum_kv_kv_k^*=m^2\,1`$.

*Démonstration.* (a) $`\mathcal E`$ est positive (théorème 7.3), et $`\iota`$ aussi
($`T^*T\otimes1=(T\otimes1)^*(T\otimes1)`$) ; $`\mathcal E(A\otimes1)=A`$ d'après la proposition 7.2(b), et
la propriété de bimodule est la proposition 7.2(c). Notons encore $`\mathrm{tr}_w`$ la trace de mêmes
poids $`w_c`$ sur $`\mathcal M_1`$ ; alors $`\mathrm{tr}_w(Y\otimes1)=m\mathrm{tr}_w(Y)`$. La définition de
$`\mathrm{Tr}_2`$ appliquée à $`AP_c`$, multipliée par $`w_c`$ et sommée sur $`c`$, donne
$`\mathrm{tr}_w(X(A\otimes1))=\mathrm{tr}_w((\mathrm{Tr}_2\,X)A)`$, et $`A=1`$ donne
$`\mathrm{tr}_w(\iota\mathcal E(X))=\mathrm{tr}_w(X)`$. Si $`F=\iota\circ F_0`$ est une autre espérance
conditionnelle sur $`\mathcal M_1\otimes1`$ préservant $`\mathrm{tr}_w`$, alors, puisque
$`F(X(A\otimes1))=F(X)(A\otimes1)`$,
$`\mathrm{tr}_w(X(A\otimes1))=\mathrm{tr}_w(F(X)(A\otimes1))=m\mathrm{tr}_w(F_0(X)A)`$, donc
$`\mathrm{tr}_w((\mathrm{Tr}_2\,X-mF_0(X))A)=0`$ pour tout $`A\in\mathcal M_1`$, et $`F_0=\mathcal E`$ par fidélité
(prendre $`A=(\mathrm{Tr}_2\,X-mF_0(X))^*`$). (b) Écrivons $`X=\sum_{a,b}X_{ab}\otimes e_{ab}`$ avec
$`X_{ab}\in\mathcal M_1`$ (théorème 7.1(a)). D'après la proposition 7.2(b),
$`\mathcal E(u_{ab}^*X)=\frac{\sqrt m}m\sum_{b'}\mathrm{Tr}_2(X_{ab'}\otimes e_{bb'})=X_{ab}/\sqrt m`$, donc
$`u_{ab}\,\iota(\mathcal E(u_{ab}^*X))=X_{ab}\otimes e_{ab}`$ ; et
$`\sum_{a,b}u_{ab}u_{ab}^*=m\sum_{a,b}1\otimes e_{aa}=m^21`$. (c) $`\mathrm{Tr}_2(Y^*)=(\mathrm{Tr}_2\,Y)^*`$
d'après la proposition 7.2(a), donc l'adjoint du développement de $`Y=u_{ab}`$ est
$`u_{ab}^*=\sum_k\iota(\mathcal E(u_{ab}^*v_k))\,v_k^*`$, et, d'après (b),
$`\sum_{a,b}u_{ab}u_{ab}^*=\sum_k\big(\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*v_k))\big)v_k^*=\sum_kv_kv_k^*`$. $`\square`$

Une famille comme en (c) est une *quasi-base* de $`\mathcal E`$, et $`\sum_kv_kv_k^*`$, qui d'après (c)
n'en dépend pas et d'après (b) existe, est l'**indice** $`\mathrm{Ind}\,\mathcal E=m^2`$ de $`\mathcal E`$
(une définition ; le nom suit Y. Watatani, Mem. Amer. Math. Soc. 83, no. 424, 1990). La **matrice
d'inclusion** d'une inclusion de $`\bigoplus_cM_{k_c}`$ dans $`\bigoplus_sM_{n_s}`$ a pour coefficients
les $`\Lambda_{sc}`$, nombres de copies du bloc $`c`$ de la sous-algèbre dans le bloc $`s`$ de la plus
grande (le chapitre 9 l'utilise pour les tours d'algèbres). Pour un facteur,
$`M_k\otimes1\subset M_k\otimes M_m`$ : en ordonnant la base $`e_a\otimes e_b`$ de $`\mathbb C^{km}`$ d'abord
selon $`b`$, $`A\otimes1_m=\mathrm{diag}(A,\dots,A)`$ avec $`m`$ copies, donc $`\Lambda=(m)`$, et
$`m^2=\|\Lambda\|^2=\dim(M_k\otimes M_m)/\dim M_k`$. (Pour aller plus loin, sans que ce soit utilisé
ici : V. F. R. Jones, Invent. Math. 72, 1, 1983, sur l'indice des sous-facteurs, et F. M. Goodman,
P. de la Harpe et V. F. R. Jones, *Coxeter graphs and towers of algebras*, Springer, 1989, sur les
matrices d'inclusion.) Pour $`\mathcal M_1`$
quelconque, l'inclusion est la somme directe des $`C_1`$ inclusions
$`M_{k_c}\otimes1\subset M_{k_c}\otimes M_m`$, $`\Lambda=m1_{C_1}`$.

*Dans la bibliothèque :* `inclusions.from_tensor_product`, avec l'espérance conditionnelle, la quasi-base et l'indice de la proposition 7.4 — [tutoriel, section 7.2](../../tutorial/fr/07_tensor_products.md#s-7-2).

<a id="s-7-3"></a>

## 7.3 Fusion des charges

Supposons maintenant que le secteur $`c`$ de $`\mathcal M_1`$ porte une charge $`q_c`$, le secteur $`d`$ de
$`\mathcal M_2`$ une charge $`r_d`$, et qu'une **règle de fusion** $`f`$ donne la charge totale
$`f(q_c,r_d)`$ — l'addition pour des nombres de particules, l'addition modulo $`2`$ pour la parité
fermionique, ou toute autre application. Posons
```math
H_q=\bigoplus_{f(q_c,r_d)=q}\mathbb C^{k_c}\otimes\mathbb C^{l_d},\qquad
\Pi_q=\sum_{f(q_c,r_d)=q}P_c\otimes Q_d,\qquad n_q=\sum_{f(q_c,r_d)=q}k_cl_d .
```
Le **produit fusionné** $`\mathcal M_1\otimes_f\mathcal M_2=\{\Pi_q\}'`$ est formé des opérateurs sur
$`H_1\otimes H_2`$ qui conservent la charge totale.

<a id="r-7-5"></a>

**Proposition 7.5.** (a) $`\mathcal M_1\otimes_f\mathcal M_2=\bigoplus_qB(H_q)\cong\bigoplus_qM_{n_q}(\mathbb C)`$,
avec les projecteurs centraux minimaux $`\Pi_q`$. (b) $`\mathcal M_1\otimes\mathcal M_2\subseteq\mathcal M_1\otimes_f\mathcal M_2`$,
avec égalité si et seulement si chaque $`H_q`$ ne contient qu'un seul couple $`(c,d)`$. (c) Si les
charges sont des nombres réels et si $`f`$ est l'addition, alors $`\mathcal M_1\otimes_f\mathcal M_2=\{N\}'`$
avec $`N=N_1\otimes1+1\otimes N_2`$, $`N_1=\sum_cq_cP_c`$, $`N_2=\sum_dr_dQ_d`$.

*Démonstration.* (a) Chaque couple appartient à exactement un $`H_q`$, donc les $`\Pi_q`$ sont des
projecteurs deux à deux orthogonaux de somme $`1`$, et $`X`$ commute avec tous si et seulement si
$`X=\sum_q\Pi_qX\Pi_q`$, avec des blocs arbitraires dans $`B(H_q)`$ ; le centre est comme dans le
théorème 7.1(b). (b) Tout élément de $`\mathcal M_1\otimes\mathcal M_2`$ commute avec chaque
$`P_c\otimes Q_d`$, donc avec leurs sommes $`\Pi_q`$. Si chaque $`H_q`$ est un seul
$`\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$, les deux algèbres sont égales à
$`\bigoplus_{c,d}B(\mathbb C^{k_c}\otimes\mathbb C^{l_d})`$ ; si $`H_q`$ contient deux couples, $`B(H_q)`$
contient des opérateurs qui envoient le bloc d'un couple dans celui de l'autre, et qui ne sont pas
diagonaux par blocs. (c) $`N=\sum_{c,d}(q_c+r_d)P_c\otimes Q_d=\sum_qq\,\Pi_q`$, donc les $`\Pi_q`$ non
nuls sont les projecteurs spectraux de $`N`$. Ce sont des polynômes en $`N`$ (démonstration du
théorème 1.5(a)) et $`N`$ en est une combinaison linéaire, donc $`X`$ commute avec $`N`$ si et seulement
s'il commute avec chaque $`\Pi_q`$. $`\square`$

**Exemple : les bosons de cœur dur.** Sur $`L`$
sites, $`(\mathbb C^2)^{\otimes L}`$ avec le nombre d'occupation $`N=\sum_in^{(i)}`$, $`n=\mathrm{diag}(0,1)`$,
le sous-espace propre $`N=q`$ est engendré par les $`\binom Lq`$ vecteurs de la base produit ayant $`q`$
sites occupés, donc les opérateurs qui conservent le nombre, $`\{N\}'`$, forment
$`\bigoplus_qM_{\binom Lq}`$ (proposition 4.10(2)) : pour deux sites $`M_1\oplus M_2\oplus M_1`$ (charges
$`0,1,2`$), pour trois sites $`M_1\oplus M_3\oplus M_3\oplus M_1`$. D'après la proposition 7.5(c), le
produit fusionné est $`\{N\}'`$ pour le nombre d'occupation total de cinq sites, à l'ordre de la base
près, et en effet $`n_q=\sum_c\binom2c\binom3{q-c}=\binom5q`$ (identité de Vandermonde : une partie à
$`q`$ éléments de cinq sites est formée de $`c`$ des deux premiers et de $`q-c`$ des trois derniers).

*Dans la bibliothèque :* le produit fusionné est `fused_tensor_product` ; les bosons de cœur dur et la proposition 7.5 — [tutoriel, section 7.3](../../tutorial/fr/07_tensor_products.md#s-7-3).

L'algèbre fusionnée est de dimension $`\sum_q\binom5q^2=\binom{10}5=252`$ (encore Vandermonde), le
produit simple de dimension $`(1+4+1)(1+9+9+1)=120`$ (théorème 7.1(b)). D'après le théorème 7.1(a) et
la proposition 7.5(b), le produit simple est formé des éléments du produit fusionné qui sont
diagonaux par blocs pour les couples $`(c,d)`$ ; les $`132`$ dimensions restantes sont engendrées par
les opérateurs qui envoient un bloc de couple $`(c,d)`$ dans un autre bloc de couple $`(c',d')`$ du même
$`H_q`$, $`c'\ne c`$ (donc $`d'\ne d`$) : ils changent la charge de chaque partie en conservant la charge
totale, c'est-à-dire qu'ici ils déplacent des particules d'une partie à l'autre.

Un état de $`\mathcal M_1\otimes_f\mathcal M_2`$ peut présenter des cohérences entre les blocs de
couples différents, si bien que la proposition 7.2 ne s'applique pas. Soit $`E_1(X)=\sum_cP_cXP_c`$
sur $`B(H_1)`$ (l'espérance conditionnelle de la section 2.6) ; son image
est $`\mathcal M_1`$, et $`\mathrm{Tr}(E_1(X)A)=\mathrm{Tr}(XA)`$ pour $`A\in\mathcal M_1`$, puisque
$`A=\sum_cP_cAP_c`$.

<a id="r-7-6"></a>

**Théorème 7.6.** Soit $`\rho`$ la densité d'un état $`\omega`$ sur $`\mathcal M_1\otimes_f\mathcal M_2`$. La
densité de sa restriction $`A\mapsto\omega(A\otimes1)`$ à $`\mathcal M_1`$ est
```math
\rho_1=E_1\big(\mathrm{Tr}_{H_2}\,\rho\big).
```
Le bloc $`(c,c')`$ de $`\mathrm{Tr}_{H_2}\,\rho`$ est $`\sum_d\mathrm{tr}_{l_d}\,\rho_{(c,d),(c',d)}`$, la somme
portant sur les $`d`$ tels que $`f(q_c,r_d)=f(q_{c'},r_d)`$, où $`\rho_{(c,d),(c',d)}`$ est le bloc de
$`\rho`$ allant du couple $`(c',d)`$ au couple $`(c,d)`$. Par conséquent, $`E_1`$ agit trivialement sur tous
les $`\mathrm{Tr}_{H_2}\,\rho`$ si et seulement s'il n'existe pas de $`c\neq c'`$ et de $`d`$ tels que
$`f(q_c,r_d)=f(q_{c'},r_d)`$ ; c'est le cas par exemple si $`f`$ est l'addition et si les $`q_c`$ sont
deux à deux distincts.

*Démonstration.* $`A\otimes1\in\mathcal M_1\otimes_f\mathcal M_2`$ (proposition 7.5(b)) et
$`\omega(A\otimes1)=\mathrm{Tr}(\rho(A\otimes1))=\mathrm{Tr}\big((\mathrm{Tr}_{H_2}\,\rho)A\big) =\mathrm{Tr}\big(E_1(\mathrm{Tr}_{H_2}\,\rho)A\big)`$. Comme $`E_1(\mathrm{Tr}_{H_2}\,\rho)\in\mathcal M_1`$, c'est
la densité, par l'unicité dans le théorème 3.2. La formule par blocs est la proposition 7.2(a)
pour $`\mathrm{Tr}_{H_2}`$, avec $`i`$ dans le bloc $`c`$, $`i'`$ dans le bloc $`c'`$ et $`j`$ dans le bloc $`d`$ ; le
bloc de $`\rho`$ entre $`(c',d)`$ et $`(c,d)`$ est nul sauf si les deux couples appartiennent au même
$`H_q`$. S'il existe $`c\neq c'`$ et $`d`$ tels que $`f(q_c,r_d)=f(q_{c'},r_d)=q`$, prenons des vecteurs
unitaires $`x\in\mathbb C^{k_c}`$, $`x'\in\mathbb C^{k_{c'}}`$, $`y\in\mathbb C^{l_d}`$ et
$`\psi=(x+x')\otimes y/\sqrt2\in H_q`$. Alors $`|\psi\rangle\langle\psi|\in B(H_q)`$, et
$`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|=\frac12|x+x'\rangle\langle x+x'|`$ possède le bloc non nul
$`\frac12|x\rangle\langle x'|`$, que $`E_1`$ supprime. $`\square`$

D'après la formule par blocs du théorème 7.6, $`(\rho_1)_c=\sum_d\mathrm{tr}_{l_d}\,\rho_{(c,d),(c,d)}`$.

*Dans la bibliothèque :* `partial_trace` et `composite.partial_trace_channel` calculent $`\rho_1`$ pour toute règle de fusion ; une vérification du théorème 7.6 pour l'addition et pour la règle de parité — [tutoriel, section 7.3](../../tutorial/fr/07_tensor_products.md#s-7-3).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Dans la description BCS
> d'un supraconducteur, seule la parité du nombre d'électrons est conservée, de sorte que l'état global peut superposer les nombres de
> particules $`0`$ et $`2`$ d'une région. Un observateur dont les observables conservent le nombre de
> particules de la région ne peut pas voir cette cohérence : pour lui, l'état réduit est l'état
> rendu diagonal par blocs du théorème 7.6 (le contenu mathématique de cet encadré), un mélange
> de nombres de particules.

<a id="s-7-4"></a>

## 7.4 Intrication des états purs

<a id="r-7-7"></a>

**Lemme 7.7 (décomposition de Schmidt).** Tout vecteur unitaire $`\psi\in\mathbb C^a\otimes\mathbb C^b`$
peut s'écrire $`\psi=\sum_{i=1}^rs_i\,u_i\otimes v_i`$ avec $`s_i>0`$, $`\sum_is_i^2=1`$ et des familles
orthonormées $`u_i\in\mathbb C^a`$, $`v_i\in\mathbb C^b`$. Alors
$`\mathrm{Tr}_2|\psi\rangle\langle\psi|=\sum_is_i^2|u_i\rangle\langle u_i|`$
et $`\mathrm{Tr}_1|\psi\rangle\langle\psi|=\sum_is_i^2|v_i\rangle\langle v_i|`$ : les deux états réduits ont
les mêmes valeurs propres non nulles $`s_i^2`$.

*Démonstration.* Soit $`\Psi`$ la matrice $`a\times b`$ définie par $`\Psi_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,\psi\rangle`$.
D'après le théorème spectral, $`\Psi\Psi^*=\sum_{i=1}^a\lambda_i|u_i\rangle\langle u_i|`$ avec une base
orthonormée $`u_i`$ et $`\lambda_i=\|\Psi^*u_i\|^2\ge0`$, ordonnés de sorte que $`\lambda_i>0`$ exactement
pour $`i\le r`$. Posons $`s_i=\sqrt{\lambda_i}`$ et $`w_i=\Psi^*u_i/s_i`$ pour $`i\le r`$ ; alors
$`\langle w_i,w_j\rangle=\langle u_i,\Psi\Psi^*u_j\rangle/(s_is_j)=\delta_{ij}`$, et $`\Psi^*u_i=0`$ pour $`i>r`$.
Donc $`\Psi=\sum_i|u_i\rangle\langle u_i|\Psi=\sum_{i\le r}|u_i\rangle\langle\Psi^*u_i|=\sum_{i\le r}s_i|u_i\rangle\langle w_i|`$,
c'est-à-dire $`\Psi_{\alpha\beta}=\sum_is_i(u_i)_\alpha\overline{(w_i)_\beta}`$ : c'est l'assertion avec
$`v_i=\bar w_i`$, et $`\sum_is_i^2=\mathrm{Tr}\,\Psi\Psi^*=\|\psi\|^2=1`$. Enfin,
$`|\psi\rangle\langle\psi|=\sum_{i,j}s_is_j|u_i\rangle\langle u_j|\otimes|v_i\rangle\langle v_j|`$, et la
proposition 7.2(b) donne les états réduits. $`\square`$

L'**entropie d'intrication** de $`\psi`$ est
$`E(\psi)=S(\mathrm{Tr}_2|\psi\rangle\langle\psi|)=S(\mathrm{Tr}_1|\psi\rangle\langle\psi|)=-\sum_is_i^2\log s_i^2`$.
Elle est nulle si et seulement si $`r=1`$ (chaque terme $`-s_i^2\log s_i^2`$ est $`\ge0`$, et nul
seulement pour $`s_i^2=1`$), c'est-à-dire si et seulement si $`\psi=u\otimes v`$ (alors $`\Psi=uv^T`$ est de
rang $`1`$, et $`r=\mathrm{rank}\,\Psi\Psi^*=\mathrm{rank}\,\Psi`$), et $`E(\psi)\le\log\min(a,b)`$, avec
égalité si et seulement si l'état réduit sur le plus petit facteur est maximalement mélangé
(corollaire 3.6, appliqué aux deux états réduits).

En présence de charges, nous prenons un vecteur unitaire $`\psi`$ de charge totale définie,
c'est-à-dire appartenant à un secteur $`H_q`$ du produit fusionné (une superposition de différents
$`H_q`$ donne sur $`\mathcal M_1\otimes_f\mathcal M_2`$ le même état qu'un mélange, chapitre 3, exercice 3) ;
écrivons $`\psi=\sum_{(c,d)}\psi_{cd}`$ avec $`\psi_{cd}\in\mathbb C^{k_c}\otimes\mathbb C^{l_d}`$.

<a id="r-7-8"></a>

**Proposition 7.8 (intrication résolue en symétrie).** Supposons que $`f`$ soit l'addition, que les
$`q_c`$ soient deux à deux distincts et que les $`r_d`$ soient deux à deux distincts, de sorte que,
pour chaque $`c`$, au plus un $`d=d(c)`$ vérifie $`q_c+r_d=q`$. Posons $`\psi_c=\psi_{c\,d(c)}`$ (ou $`0`$),
$`p_c=\|\psi_c\|^2`$ et, pour $`p_c>0`$, $`\hat\psi_c=\psi_c/\sqrt{p_c}`$. Alors
$`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|`$ appartient à $`\mathcal M_1`$ (il commute avec $`N_1`$) et vaut
$`\bigoplus_cp_c\mathrm{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|`$ ; l'état réduit sur $`\mathcal M_2`$ a la
probabilité de secteur $`p_c`$ en $`d(c)`$ ; et
```math
E(\psi)=H(p)+\sum_cp_c\,E(\hat\psi_c),\qquad H(p)=-\sum_cp_c\log p_c .
```

*Démonstration.* D'après la formule par blocs du théorème 7.6, le bloc $`(c,c')`$ de
$`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|`$ est la somme des $`\mathrm{tr}_{l_d}|\psi_{cd}\rangle\langle\psi_{c'd}|`$
sur les $`d`$ tels que $`q_c+r_d=q=q_{c'}+r_d`$, ce qui impose $`c=c'`$ et $`d=d(c)`$. L'état réduit est
donc $`\bigoplus_c\mathrm{Tr}_2|\psi_c\rangle\langle\psi_c|`$, de traces de blocs $`p_c`$ ; en échangeant les
rôles des facteurs, on obtient l'énoncé pour $`\mathcal M_2`$. D'après le lemme 7.7, $`E(\psi)`$ est
l'entropie de cet état, et le théorème 3.3 avec $`\sigma_c=\mathrm{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|`$,
$`S(\sigma_c)=E(\hat\psi_c)`$, donne la formule. $`\square`$

$`H(p)`$ est l'**intrication de nombre** et $`\sum_cp_cE(\hat\psi_c)`$ l'**intrication
configurationnelle**.

*Calcul :* le lemme 7.7 et la proposition 7.8 pour un état aléatoire de deux particules sur cinq sites — [tutoriel, section 7.4](../../tutorial/fr/07_tensor_products.md#s-7-4).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Des parties soumises à la règle de
> supersélection du nombre de particules ne peuvent utiliser comme ressource que la partie configurationnelle : c'est
> l'*intrication de particules* de H. M. Wiseman et J. A. Vaccaro (Phys. Rev. Lett. 91, 097902,
> 2003). $`H(p)`$ est mesurée dans les atomes froids sous le nom d'*entropie de nombre* (A. Lukin et
> al., Science 364, 256, 2019).

<a id="s-7-5"></a>

## 7.5 États aléatoires et formule de Page

Soit $`g\in\mathbb C^N`$ à coordonnées indépendantes, dont les parties réelle et imaginaire sont
indépendantes de loi $`N(0,\tfrac12)`$ (un vecteur gaussien complexe standard, proposition A.16(a)), de
sorte que $`\mathbb E|g_i|^2=1`$ ; $`g\neq0`$ presque sûrement (théorème A.20). Un **état pur aléatoire de
Haar** est $`\psi=g/\|g\|`$. Sa loi est invariante par les unitaires : $`Ug`$ a la loi de $`g`$ pour tout
unitaire $`U`$ (proposition A.16(a)), et $`Ug/\|Ug\|=U\psi`$. D'après la proposition A.22, c'est l'unique
loi invariante par les unitaires sur la sphère unité et, comme le montre sa démonstration, c'est la
loi de la première colonne d'un unitaire de Haar (théorème 5.3).

<a id="r-7-9"></a>

**Lemme 7.9.** $`\|g\|`$ et $`\psi=g/\|g\|`$ sont indépendants.

*Démonstration.* (C'est une partie du théorème A.20(b).) Soient $`h\ge0`$ et $`\phi\ge0`$ des fonctions
boréliennes, respectivement sur $`(0,\infty)`$ et sur la sphère unité $`S^{2N-1}`$ de
$`\mathbb C^N\cong\mathbb R^{2N}`$. Comme $`g`$ a la densité $`\pi^{-N}e^{-|x|^2}`$ (proposition A.16(a)), les
coordonnées polaires (théorème A.18 avec $`n=2N`$ et la mesure $`\mu`$ qui y est définie) donnent
```math
\mathbb E\big[h(\|g\|)\phi(\psi)\big]=\pi^{-N}\!\int_{\mathbb R^{2N}}\!h(|x|)e^{-|x|^2}\phi\Big(\frac x{|x|}\Big)dx
=2N\pi^{-N}\int_0^\infty h(t)e^{-t^2}t^{2N-1}dt\int_{S^{2N-1}}\phi\,d\mu=:\alpha(h)\beta(\phi).
```
En prenant $`\phi=1`$, puis $`h=1`$, puis les deux, on obtient $`\mathbb Eh(\|g\|)=\alpha(h)\beta(1)`$,
$`\mathbb E\phi(\psi)=\alpha(1)\beta(\phi)`$ et $`\alpha(1)\beta(1)=1`$, d'où
$`\mathbb E[h(\|g\|)\phi(\psi)]=\mathbb Eh(\|g\|)\,\mathbb E\phi(\psi)`$. $`\square`$

<a id="r-7-10"></a>

**Lemme 7.10.** $`\mathbb E\|g\|^4=N(N+1)`$. Si $`G`$ est une matrice $`a\times b`$ à coefficients
indépendants distribués comme les $`g_i`$ (une matrice *de Ginibre*), alors
$`\mathbb E\mathrm{Tr}(GG^*)^2=ab(a+b)`$.

*Démonstration.* Avec $`g_i=x+iy`$ et $`x,y`$ indépendants de loi $`N(0,\frac12)`$, $`\mathbb Ex^4=3(\mathbb Ex^2)^2=\frac34`$
(théorème A.24 : trois appariements), donc
$`\mathbb E|g_i|^4=\mathbb Ex^4+2\mathbb Ex^2\mathbb Ey^2+\mathbb Ey^4=\tfrac34+\tfrac12+\tfrac34=2`$ et
$`\mathbb E\|g\|^4=\sum_i\mathbb E|g_i|^4+\sum_{i\ne j}\mathbb E|g_i|^2\mathbb E|g_j|^2=2N+N(N-1)`$. Ensuite,
$`\mathrm{Tr}(GG^*)^2=\sum_{i,i',j,j'}G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}`$. Les coefficients et leurs
conjugués sont des combinaisons linéaires complexes des coordonnées d'un même vecteur gaussien
centré (les parties réelles et imaginaires de tous les coefficients), avec
$`\mathbb E[G_{ij}G_{kl}]=\mathbb E[\bar G_{ij}\bar G_{kl}]=0`$ et $`\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}`$
(proposition A.16(a)). D'après la formule d'Isserlis–Wick (théorème A.24, forme complexe, avec les
trois appariements des quatre facteurs),
$`\mathbb E[G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}]=\delta_{ii'}\delta_{i'i}+0+\delta_{jj'}\delta_{jj'}=\delta_{ii'}+\delta_{jj'}`$.
En sommant sur $`i,i'\le a`$ et $`j,j'\le b`$, on obtient $`ab^2+a^2b`$. $`\square`$

<a id="r-7-11"></a>

**Théorème 7.11 (E. Lubkin, J. Math. Phys. 19, 1028, 1978).** Pour $`\psi\in\mathbb C^a\otimes\mathbb C^b`$
aléatoire de Haar et $`\rho_A=\mathrm{Tr}_2|\psi\rangle\langle\psi|`$,
```math
\mathbb E\mathrm{Tr}\,\rho_A^2=\frac{a+b}{ab+1}.
```

*Démonstration.* Soit $`G_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,g\rangle`$, une matrice de Ginibre ;
la matrice de $`\psi`$ est $`\Psi=G/\|g\|`$ et, d'après la proposition 7.2(a),
$`(\rho_A)_{\alpha\alpha'}=\sum_\beta\Psi_{\alpha\beta}\bar\Psi_{\alpha'\beta}`$, c'est-à-dire
$`\rho_A=\Psi\Psi^*=GG^*/\|g\|^2`$. Donc $`\mathrm{Tr}(GG^*)^2=\|g\|^4\mathrm{Tr}\,\rho_A^2`$, où
$`\mathrm{Tr}\,\rho_A^2\in[0,1]`$ est une fonction de $`\psi`$. D'après le lemme 7.9,
$`\mathbb E\mathrm{Tr}(GG^*)^2=\mathbb E\|g\|^4\,\mathbb E\mathrm{Tr}\,\rho_A^2`$, et d'après le lemme 7.10,
$`\mathbb E\mathrm{Tr}\,\rho_A^2=ab(a+b)/(ab(ab+1))`$. $`\square`$

Le chapitre 6 obtient la même moyenne à partir des moments d'ordre quatre d'un vecteur unitaire de
loi invariante par les unitaires (lemme 6.14), comme pureté de sortie d'un canal aléatoire
(théorème 6.15 et la remarque qui le suit).

<a id="r-7-12"></a>

**Corollaire 7.12.** $`\log\frac{ab+1}{a+b}\le\mathbb E\,S(\rho_A)\le\log\min(a,b)`$.

*Démonstration.* D'après le corollaire 3.6, $`S(\rho_A)\le\log a`$, et d'après le lemme 7.7,
$`S(\rho_A)=S(\mathrm{Tr}_1|\psi\rangle\langle\psi|)\le\log b`$. L'inégalité de Jensen finie pour la
fonction concave $`\log`$ (solution de l'exercice 2 du chapitre 3), avec les poids $`\lambda_i>0`$,
valeurs propres non nulles de $`\rho_A`$, donne $`\sum_i\lambda_i\log\lambda_i\le\log\sum_i\lambda_i^2`$,
c'est-à-dire $`S(\rho_A)\ge-\log X`$ avec $`X=\mathrm{Tr}\,\rho_A^2\in(0,1]`$. Soit
$`m=\mathbb EX=\frac{a+b}{ab+1}`$ (théorème 7.11). Le lemme 3.4 pour la fonction convexe $`-\log`$ sur
$`(0,\infty)`$ au point $`m`$ donne $`-\log X\ge-\log m-(X-m)/m`$ ; les deux membres sont bornés
($`X\ge1/a`$ par Cauchy–Schwarz), et en prenant les espérances on obtient
$`\mathbb E[-\log X]\ge-\log m`$. $`\square`$

Pour $`a`$ fixé et $`b\to\infty`$, les bornes se rejoignent : l'entropie moyenne d'une petite partie
d'un état pur aléatoire tend vers son maximum $`\log a`$. La moyenne exacte est donnée par la formule
de Page, démontrée dans l'annexe C.

<a id="r-7-13"></a>

**Théorème 7.13 (formule de Page).** Pour $`a\le b`$,
```math
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
```
De plus, $`\log a-\frac{a-1}{ab+1}-\frac{a-1}{2b}\le\mathbb E\,S(\rho_A)\le\log a-\frac{a-1}{2b}`$, donc
$`\big|\mathbb E\,S(\rho_A)-\log a+\frac a{2b}\big|\le\frac1{2b}`$ ; pour deux qubits, $`\mathbb E\,S(\rho_A)=\frac13`$.

*Démonstration.* C'est le théorème C.17 joint au corollaire C.18. La démonstration de l'annexe C part
de $`\rho_A=GG^*/\mathrm{Tr}\,GG^*`$ (démonstration du théorème 7.11), calcule la densité jointe des
valeurs propres de la matrice de Wishart $`GG^*`$ (théorème C.13) et de $`\rho_A`$ (théorème C.14), et
intègre $`x\log x`$ contre la densité de niveaux de Laguerre (lemmes C.15 et C.16). $`\square`$

*Calcul :* les théorèmes 7.11 et 7.13 et le corollaire 7.12 sur $`10^4`$ états purs aléatoires de $`\mathbb C^3\otimes\mathbb C^6`$ — [tutoriel, section 7.5](../../tutorial/fr/07_tensor_products.md#s-7-5).

La formule a été conjecturée par D. N. Page (Phys. Rev. Lett. 71, 1291, 1993) et démontrée par
S. K. Foong et S. Kanno (Phys. Rev. Lett. 72, 1148, 1994), J. Sánchez-Ruiz (Phys. Rev. E 52, 5653,
1995) et S. Sen (Phys. Rev. Lett. 77, 1, 1996).

Pour $`a>b`$, les rôles de $`a`$ et $`b`$ sont échangés. En effet, l'unitaire
$`e_\alpha\otimes e_\beta\mapsto e_\beta\otimes e_\alpha`$ de $`\mathbb C^a\otimes\mathbb C^b`$ sur $`\mathbb C^b\otimes\mathbb C^a`$
permute les coordonnées de $`g`$ ; il envoie donc un état pur aléatoire de Haar $`\psi`$ sur un état
pur aléatoire de Haar $`\psi'`$ de $`\mathbb C^b\otimes\mathbb C^a`$ ; sa matrice est $`\Psi^T`$, donc
$`\mathrm{Tr}_2|\psi'\rangle\langle\psi'|=\Psi^T\bar\Psi=\mathrm{Tr}_1|\psi\rangle\langle\psi|`$ (proposition 7.2(a)),
dont l'entropie est $`S(\rho_A)`$ d'après le lemme 7.7. Ainsi, avec $`a'=\min(a,b)`$ et $`b'=\max(a,b)`$,
le théorème 7.13 dit que l'entropie d'intrication moyenne reste en deçà de son maximum $`\log a'`$
d'au moins $`\frac{a'-1}{2b'}`$ et d'au plus $`\frac{a'-1}{2b'}+\frac{a'-1}{a'b'+1}`$, ce qui est moins d'un nat :
c'est $`0`$ pour $`b'=1`$ et, sinon, $`\frac{a'-1}{2b'}\le\frac12`$ et $`\frac{a'-1}{a'b'+1}<\frac1{b'}\le\frac12`$.
Ainsi, pour $`ab`$ fixé, en fonction de $`a`$, la moyenne est à moins d'un nat de
$`\log\min(a,b)=\min(\log a,\log(ab)-\log a)`$.

La loi de $`\rho_A`$ est la **mesure induite** (le nom est celui de K. Życzkowski et H.-J. Sommers,
J. Phys. A 34, 7111, 2001 ; le théorème C.14 donne la densité de ses valeurs propres).
`states.random_density_matrix(alg, rank=r)` forme $`G_cG_c^*`$ avec des blocs de Ginibre indépendants
$`G_c`$ de taille $`k_c\times r`$ et divise par $`\sum_c\mathrm{Tr}\,G_cG_c^*`$. Pour un seul bloc, c'est la loi
de $`\rho_A`$ ci-dessus avec $`b=r`$ (démonstration du théorème 7.11). Pour plusieurs blocs, c'est la loi
de la densité de la restriction à $`\mathcal M_1`$ de l'état vectoriel d'un état pur aléatoire de Haar
$`\psi=g/\|g\|`$ de $`H_1\otimes\mathbb C^r`$. En effet, soit $`G_{i\beta}=\langle e_i\otimes e_\beta,g\rangle`$, une
matrice de Ginibre $`D\times r`$ ($`D=\sum_ck_c`$), et soient $`G_c`$ ses lignes du bloc $`c`$ : des matrices
de Ginibre $`k_c\times r`$ indépendantes. Prenons $`\mathcal M_2=M_r`$ (un seul secteur) et une règle de
fusion constante ; il y a alors une seule charge totale, $`\Pi=1`$, et
$`\mathcal M_1\otimes_f\mathcal M_2=B(H_1\otimes\mathbb C^r)`$. D'après le théorème 7.6, la restriction a la
densité $`E_1(\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|)`$, et
$`\mathrm{Tr}_{H_2}|\psi\rangle\langle\psi|=GG^*/\|g\|^2`$ (proposition 7.2(a), comme dans la démonstration
du théorème 7.11), dont les blocs diagonaux sont les $`G_cG_c^*/\|g\|^2`$, avec
$`\|g\|^2=\sum_c\mathrm{Tr}\,G_cG_c^*`$.

*Dans la bibliothèque :* `states.random_density_matrix(alg, rank=r)` — [tutoriel, section 7.5](../../tutorial/fr/07_tensor_products.md#s-7-5).

<a id="r-7-14"></a>

**Proposition 7.14.** Pour $`\rho`$ tiré par `random_density_matrix(alg, rank=r)` sur
$`\bigoplus_cM_{k_c}`$, avec $`D=\sum_ck_c`$,
```math
\mathbb E\mathrm{Tr}\,\rho^2=\frac{\sum_ck_c(k_c+r)}{D(Dr+1)} .
```

*Démonstration.* $`\mathrm{Tr}\,\rho^2=\sum_c\mathrm{Tr}(G_cG_c^*)^2/\|g\|^4`$, où $`g\in\mathbb C^{Dr}`$ est la
collection de tous les coefficients, est une fonction de $`g/\|g\|`$ à valeurs dans $`[0,1]`$. Comme
dans le théorème 7.11, le lemme 7.9 donne
$`\sum_c\mathbb E\mathrm{Tr}(G_cG_c^*)^2=\mathbb E\|g\|^4\,\mathbb E\mathrm{Tr}\,\rho^2`$, et le lemme 7.10 donne
$`\mathbb E\mathrm{Tr}(G_cG_c^*)^2=k_cr(k_c+r)`$ et $`\mathbb E\|g\|^4=Dr(Dr+1)`$, donc
$`\mathbb E\mathrm{Tr}\,\rho^2=\sum_ck_cr(k_c+r)/(Dr(Dr+1))`$. $`\square`$

*Calcul :* la proposition 7.14 pour $`M_1\oplus M_2\oplus M_3`$, $`r=2`$ — [tutoriel, section 7.5](../../tutorial/fr/07_tensor_products.md#s-7-5).

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). L'entropie d'intrication moyenne
> d'un état pur aléatoire, tracée en fonction de la taille de la partie à dimension totale fixée
> (la courbe décrite après le théorème 7.13), s'appelle la *courbe de Page*. En physique, elle sert
> à décrire comment l'information quitte un trou noir dans le rayonnement de Hawking, et pourquoi
> les parties d'un système fermé paraissent thermiques ; ce sont des interprétations, non des
> résultats de ce cours.

<a id="s-7-6"></a>

## 7.6 Points de comparaison

*Dans la bibliothèque :* la comparaison avec d'autres logiciels — [tutoriel, section 7.6](../../tutorial/fr/07_tensor_products.md#s-7-6).

<a id="x-7"></a>

## Exercices

1. Calculer la duale de $`\mathrm{Tr}_2`$ relativement à $`\tau_{\mathrm{vN}}`$ sur
   $`\mathcal M_1\otimes\mathcal M_2`$ et sur $`\mathcal M_1`$.
2. (Sous-additivité.) Soit $`\rho`$ un état sur $`\mathcal M_1\otimes\mathcal M_2`$ dont les états réduits
   $`\rho_1,\rho_2`$ sont inversibles. Montrer que $`S(\rho)\le S(\rho_1)+S(\rho_2)`$, avec égalité si
   et seulement si $`\rho=\rho_1\otimes\rho_2`$.
3. Soit `site` l'algèbre $`\mathbb C\oplus\mathbb C`$ avec les charges $`0,1`$ (un site vide ou occupé),
   fusionnée avec elle-même par l'addition. Calculer $`E(\psi)`$ et ses deux parties pour
   $`\psi=(|01\rangle+|10\rangle)/\sqrt2`$, une particule partagée entre deux sites, et décrire la
   restriction de l'état vectoriel de $`(|00\rangle+|11\rangle)/\sqrt2`$ à l'algèbre fusionnée.

*Calcul :* les trois exercices, numériquement — [tutoriel, Exercices](../../tutorial/fr/07_tensor_products.md#x-7).

<details>
<summary>Solutions</summary>

1. Les poids sont $`1/(C_1k_c)`$ sur $`\mathcal M_1`$ et $`1/(C_1C_2k_cl_d)`$ sur le produit. Pour
   $`A\in\mathcal M_1`$, d'après la proposition 7.2(a),
   $`\tau_{\mathrm{vN}}((\mathrm{Tr}_2\,X)A)=\sum_c\frac1{C_1k_c}\sum_d\mathrm{Tr}(X_{cd}(A_c\otimes1))    =\sum_{c,d}\frac1{C_1C_2k_cl_d}\mathrm{Tr}\big(X_{cd}\,C_2l_d(A_c\otimes1)\big)`$, donc la duale est
   $`A\mapsto\bigoplus_{c,d}C_2l_d\,A_c\otimes1_{l_d}=A\otimes Z`$, avec l'élément central
   $`Z=\bigoplus_dC_2l_d1_{l_d}`$ de $`\mathcal M_2`$. Elle n'est unifère que si $`Z=1`$, c'est-à-dire
   $`\mathcal M_2=\mathbb C`$, et ne l'est pas dans l'exemple ($`Z=4\cdot1_2\oplus6\cdot1_3`$). Pour `Tr_norm`,
   le même calcul donne $`A\otimes\bigoplus_dl_d1_{l_d}`$, unifère seulement si tous les $`l_d=1`$. Les
   traces brutes font toujours de $`\mathrm{Tr}_2`$ la duale de $`A\mapsto A\otimes1`$ ; en général, c'est le
   cas de $`\mathrm{tr}_w`$ sur $`\mathcal M_1`$ avec les poids $`w_{cd}=w_c`$ sur le produit, par le même
   calcul.

2. $`\rho_1\otimes\rho_2`$ est un état (théorème 7.1(c)), inversible, de vecteurs propres
   $`u_i\otimes v_j`$ et de valeurs propres $`\lambda_i\mu_j`$, donc
   $`\log(\rho_1\otimes\rho_2)=\log\rho_1\otimes1+1\otimes\log\rho_2`$ (théorème 1.5(a)). Par définition
   de la trace partielle, $`\mathrm{Tr}\,\rho(\log\rho_1\otimes1)=\mathrm{Tr}\,\rho_1\log\rho_1`$, et de même
   pour le second facteur. Ainsi $`S(\rho\|\rho_1\otimes\rho_2)=S(\rho_1)+S(\rho_2)-S(\rho)`$, qui est
   $`\ge0`$, avec égalité si et seulement si $`\rho=\rho_1\otimes\rho_2`$ (corollaire 3.6). Cette
   différence est l'information mutuelle.

3. Les secteurs sont $`N=0,1,2`$, de tailles $`1,2,1`$. $`\psi`$ appartient à $`N=1`$, avec
   $`\psi_0=|01\rangle/\sqrt2`$, $`\psi_1=|10\rangle/\sqrt2`$, donc $`p=(\frac12,\frac12)`$ et, les blocs
   étant de dimension un, $`E(\hat\psi_c)=0`$ : d'après la proposition 7.8, $`E(\psi)=H(p)=\log2`$ est
   entièrement de l'intrication de nombre, et la partie configurationnelle est nulle.
   Le vecteur $`(|00\rangle+|11\rangle)/\sqrt2`$ superpose $`N=0`$ et $`N=2`$ ; sur l'algèbre fusionnée,
   diagonale par blocs, son état vectoriel est le mélange
   $`\frac12|00\rangle\langle00|+\frac12|11\rangle\langle11|`$ d'états produits (chapitre 3,
   exercice 3).

</details>

**Suite :** le [chapitre 8](08_dynamics.md) traite de la dynamique : les équations de Schrödinger et
de Lindblad dans une algèbre à secteurs, les sauts entre secteurs et la post-sélection.
