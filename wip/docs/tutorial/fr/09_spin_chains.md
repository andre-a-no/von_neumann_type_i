[English](../09_spin_chains.md) · [Español](../es/09_spin_chains.md) · **Français** · [Deutsch](../de/09_spin_chains.md) · [中文](../zh/09_spin_chains.md) · [日本語](../ja/09_spin_chains.md) · [Русский](../ru/09_spin_chains.md)

# 9. Chaînes de spins

Une chaîne de $L$ spins $\frac12$ vit sur $(\mathbb C^2)^{\otimes L}$, et l'algèbre linéaire dense en dimension
$2^L$ devient coûteuse au-delà de $L\approx14$ (une matrice dense a $4^L$ coefficients). Les symétries repoussent cette limite : l'aimantation découpe l'espace en $L+1$
secteurs (chapitre 4), les translations découpent à nouveau chaque secteur, et $SU(2)$ réduit encore les blocs.
Ce chapitre construit ces secteurs pour la chaîne XXZ, démontre ce qu'ils disent des états fondamentaux, confronte
la bibliothèque aux fermions libres et à l'ansatz de Bethe, et va au-delà des matrices denses (méthodes de Krylov)
et des chaînes propres (désordre).

**Vous apprendrez** pourquoi l'algèbre d'une chaîne est $\bigoplus_NM_{\binom LN}$ et pourquoi la tour de ces
algèbres a pour diagramme de Bratteli le triangle de Pascal et pour indice $4\cos^2\frac\pi{2L+4}$, comment on
dénombre les secteurs d'impulsion, la règle des signes de Marshall et ses conséquences pour l'impulsion et le spin
des états fondamentaux (Lieb–Mattis), la solution de Jordan–Wigner de la chaîne XX, les équations de Bethe à deux
magnons, pourquoi Lanczos converge de façon monotone, et comment mesurer les statistiques de niveaux et
l'intrication des fermions libres.

Dans tout ce chapitre, les sites sont $0,\dots,L-1$, la base de $(\mathbb C^2)^{\otimes L}$ est formée des
configurations $|b\rangle=|b_0\dots b_{L-1}\rangle$ avec $b_i=1$ pour un spin vers le haut (convention de la
bibliothèque : le bit $i$ d'une étiquette est le site $i$), $s^a_i$ ($a=x,y,z$) sont les matrices de spin
$\frac12\sigma_a$ sur le site $i$, $s^\pm_i=s^x_i\pm is^y_i$, $n_i=s^z_i+\frac12$, $S^a=\sum_is^a_i$, et
$N=S^z+\frac L2$ est le nombre de spins vers le haut. Une chaîne est **ouverte** (liaisons $(i,i+1)$, $i<L-1$) ou
**périodique** (avec aussi la liaison $(L-1,0)$ ; alors $L\ge3$ et les indices sont pris modulo $L$). $J>0$
(antiferromagnétique), $\hbar=1$, et $\log$ est le logarithme népérien.

## 9.1 La chaîne XXZ et ses secteurs $U(1)$

La **chaîne XXZ** d'anisotropie $\Delta\in\mathbb R$ et de champs $h_i\in\mathbb R$ est
$$
H=J\sum_{(i,j)}\Big[\tfrac12\big(s^+_is^-_j+s^-_is^+_j\big)+\Delta\,s^z_is^z_j\Big]+\sum_ih_is^z_i ,
$$
la somme portant sur les liaisons ; $\frac12(s^+_is^-_j+s^-_is^+_j)=s^x_is^x_j+s^y_is^y_j$ déplace un spin vers le
haut le long de la liaison. Pour $\Delta=1$, c'est la chaîne de Heisenberg $J\sum\mathbf s_i\cdot\mathbf s_j$ ; pour
$\Delta=0$, la chaîne XX.

**Proposition 9.1.** $[H,S^z]=0$. Par conséquent $H\in U(1)'\cong\bigoplus_{N=0}^LM_{\binom LN}(\mathbb C)$, le bloc
$N$ agissant sur l'espace engendré par les configurations à $N$ spins vers le haut.

*Démonstration.* Les termes diagonaux commutent avec la matrice diagonale $S^z$. Pour $i\neq j$, $s^+_is^-_j$ commute
avec $s^z_l$, $l\ne i,j$, et $[s^z_i+s^z_j,s^+_is^-_j]=[s^z_i,s^+_i]s^-_j+s^+_i[s^z_j,s^-_j]=s^+_is^-_j-s^+_is^-_j=0$
d'après $[s^z,s^\pm]=\pm s^\pm$ (section 4.5) ; pour $s^-_is^+_j$, on passe aux adjoints. Par conséquent $H$
commute avec $N=S^z+\frac L2$, qui a les valeurs propres entières $0,\dots,L$, de sous-espaces propres engendrés par
les configurations à $N$ spins vers le haut, de dimension $\binom LN$ ; le reste est la proposition 4.10(2) avec
$Q=N$. $\square$

`SpinChain(L).xxz` construit les blocs directement à partir des configurations de chaque secteur
(`chain.labels[N]`). Nous les comparons à la matrice complète $2^L\times2^L$ construite par produits de
Kronecker :

```python
import functools, itertools, math, time
import numpy as np
import torch
from scipy.optimize import brentq
from torch_vn_algebra import SpinChain, SU2Chain, DensityMatrix, inclusions, krylov

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
close = lambda a, b, tol=1e-10: bool((torch.as_tensor(a) - torch.as_tensor(b)).abs().max() < tol)
def op(L, ops):                             # product of one-site operators {site: 2 x 2}; site 0 = lowest bit
    return functools.reduce(torch.kron, [ops.get(j, torch.eye(2, dtype=f64)) for j in reversed(range(L))])

sp = torch.tensor([[0., 0.], [1., 0.]], dtype=f64)                    # s^+ : |0> = down -> |1> = up
sm, sz = sp.T.contiguous(), torch.diag(torch.tensor([-0.5, 0.5], dtype=f64))
L = 6
chain = SpinChain(L, boundary='periodic')
h = 0.3 * torch.randn(L, dtype=f64)                                  # random fields
H_full = sum(0.5 * (op(L, {i: sp, j: sm}) + op(L, {i: sm, j: sp})) + 0.7 * op(L, {i: sz, j: sz})
             for i, j in chain.bonds) + sum(h[i] * op(L, {i: sz}) for i in range(L))   # J = 1, Delta = 0.7
Sz = sum(op(L, {i: sz}) for i in range(L))
assert close(H_full @ Sz, Sz @ H_full)                               # Proposition 9.1
H = chain.xxz(J=1.0, Delta=0.7, h=h)                                 # an element of (+)_N M_binom(6, N)
for N, k in enumerate(chain.algebra.k_factors):
    lab = chain.labels[N]                                            # configurations of sector N
    assert close(H.matrix[0, N, :k, :k], H_full[lab][:, lab])
```

**La tour des chaînes.** Soit $\mathcal A_L=\bigoplus_NM_{\binom LN}$, agissant sur $L$ sites. Ajouter un site,
$x\mapsto x\otimes1_2$, envoie $\mathcal A_L$ dans $\mathcal A_{L+1}$, car $x\otimes1$ commute avec
$S^z\otimes1+1\otimes s^z$. On obtient une tour $\mathcal A_0\subset\mathcal A_1\subset\mathcal A_2\subset\cdots$.

**Proposition 9.2.** L'inclusion $\mathcal A_L\subset\mathcal A_{L+1}$ a pour matrice d'inclusion la matrice
$\Lambda_L$ de taille $(L+2)\times(L+1)$ telle que $\lambda_{N,N}=\lambda_{N,N-1}=1$, tous les autres coefficients étant
nuls ; les tailles des blocs obéissent à la règle de Pascal $\binom{L+1}N=\binom LN+\binom L{N-1}$, le diagramme de
Bratteli de la tour est le triangle de Pascal, et le graphe biparti de $\Lambda_L$ est un chemin à $2L+3$ sommets.

*Démonstration.* Les configurations à $N$ spins vers le haut sur $L+1$ sites sont celles qui ont $N$ spins vers le
haut sur les $L$ premiers sites et le nouveau spin vers le bas, et celles qui en ont $N-1$ avec le nouveau spin vers
le haut. Sur ces deux sous-espaces, $x\otimes1$ agit comme $x_N$ et $x_{N-1}$, donc le bloc $N$ de $\mathcal A_{L+1}$
contient une fois chacun les blocs $N$ et $N-1$ de $\mathcal A_L$ ; les dimensions s'additionnent selon la règle de
Pascal. Dans le graphe biparti, le secteur $N$ de $\mathcal A_{L+1}$ est relié aux secteurs $N-1$ et $N$ de
$\mathcal A_L$, ce qui donne le chemin $0'-0-1'-1-\cdots-L'-L-(L+1)'$ (les primes désignent $\mathcal A_{L+1}$). $\square$

Un chemin allant du sommet du diagramme de Bratteli au secteur $N$ du niveau $L$ est une configuration $b$ ; il y
en a $\binom LN$, la taille du bloc. Pour une inclusion connexe de matrice d'inclusion $\Lambda$, nous **définissons**
l'**indice** comme $\|\Lambda\|^2$ (norme d'opérateur) ; pour $M_k\otimes1\subset M_k\otimes M_m$, la matrice
d'inclusion est la matrice $1\times1$ $(m)$, donc l'indice est le $m^2$ de la section 7.2. `Inclusion.index` le
calcule. (Le nom est emprunté à la théorie de l'indice de V. Jones, *Invent. Math.* 72 (1983) 1 ; dans ce cours,
l'indice d'une inclusion d'algèbres de dimension finie est cette définition et rien de plus.)

**Lemme 9.3.** La matrice d'adjacence $A_n$ du chemin de sommets $1,\dots,n$ a pour valeurs propres
$2\cos\frac{\pi k}{n+1}$, $k=1,\dots,n$, de vecteurs propres $v_k(j)=\sin\frac{\pi jk}{n+1}$ ; par conséquent
$\|A_n\|=2\cos\frac\pi{n+1}$.

*Démonstration.* Posons $v_k(0)=v_k(n+1)=0$. Alors $(A_nv_k)(j)=v_k(j-1)+v_k(j+1)=2\cos\frac{\pi k}{n+1}\,v_k(j)$
d'après $\sin(a-b)+\sin(a+b)=2\sin a\cos b$, et $v_k\ne0$ car $v_k(1)=\sin\frac{\pi k}{n+1}\neq0$. Les $n$ valeurs
propres sont distinctes ($\cos$ est strictement décroissant sur $(0,\pi)$), ce sont donc toutes les valeurs
propres, et la norme de la matrice symétrique $A_n$ est le plus grand module, atteint en $k=1$ et $k=n$. $\square$

**Corollaire 9.4.** $[\mathcal A_{L+1}:\mathcal A_L]=\|\Lambda_L\|^2=4\cos^2\frac\pi{2L+4}$, qui croît vers $4$.

*Démonstration.* En ordonnant les sommets le long du chemin de la proposition 9.2, la matrice d'adjacence
$B=\begin{pmatrix}0&\Lambda_L\\\Lambda_L^T&0\end{pmatrix}$ devient $A_{2L+3}$ (une similitude par permutation). D'après le
lemme 1.1, $\|B\|^2=\|B^2\|$ et $\|\Lambda_L\Lambda_L^T\|=\|\Lambda_L^T\Lambda_L\|=\|\Lambda_L\|^2$, et $B^2=\Lambda_L\Lambda_L^T\oplus\Lambda_L^T\Lambda_L$, donc
$\|B\|=\|\Lambda_L\|$ ; le lemme 9.3 avec $n=2L+3$ donne la valeur. $\square$

```python
tower = inclusions.chain_tower(8)                                    # A_0 c A_1 c ... c A_8
print(inclusions.bratteli_diagram(tower).levels[:6])                 # Pascal's triangle
for L_, inc in enumerate(tower):
    sizes = [math.comb(L_ + 1, N) for N in range(L_ + 2)]
    assert inc.M.k_factors == sizes and inc.is_connected
    assert list(inc.matrix @ np.array(inc.N.k_factors)) == sizes      # unitality = Pascal's rule
    n = 2 * L_ + 3
    A_path = np.eye(n, k=1) + np.eye(n, k=-1)
    assert abs(np.linalg.eigvalsh(A_path).max() - 2 * math.cos(math.pi / (n + 1))) < 1e-12   # Lemma 9.3
    assert abs(inc.index() - 4 * math.cos(math.pi / (2 * L_ + 4)) ** 2) < 1e-12           # Corollary 9.4
```

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). $S^z$ est conservé dès que l'interaction est invariante par les rotations autour
> de l'axe $z$ ; $N$ compte les spins vers le haut, ou les particules après la transformation de Jordan–Wigner de
> la section 9.4. Dans la tour, chaque secteur de charge de la chaîne la plus longue est alimenté par deux
> secteurs de la plus courte.

## 9.2 Translations, impulsion et règle des signes de Marshall

Sur une chaîne périodique, soit $T$ la translation $T|b_0b_1\dots b_{L-1}\rangle=|b_{L-1}b_0\dots b_{L-2}\rangle$, qui
envoie le site $i$ sur le site $i+1$. $T$ est une matrice de permutation, donc unitaire, $T^L=1$ et
$Ts^a_iT^*=s^a_{i+1}$. Par conséquent $TS^zT^*=S^z$ et $THT^*=H$ pour un champ uniforme $h_i=h$ : $T$ envoie les
liaisons sur des liaisons. Ses valeurs propres sont donc de la forme $e^{ik}$, d'**impulsions** $k=2\pi m/L$, $m=0,\dots,L-1$.

**Proposition 9.5 (secteurs d'impulsion).** Fixons $N$. Pour une configuration $a$ à $N$ spins vers le haut, soit
$R_a$ la taille de son orbite $\{T^ja\}$. (a) $R_a$ divise $L$ ; pour $k=2\pi m/L$ avec $mR_a\equiv0\pmod L$, le
vecteur $|a,k\rangle=\sum_{j=0}^{R_a-1}e^{-ikj}T^j|a\rangle$ est non nul et $T|a,k\rangle=e^{ik}|a,k\rangle$, et ces $R_a$
vecteurs forment une base orthogonale de l'espace engendré par l'orbite. (b) La valeur propre $e^{2\pi im/L}$ de
$T$ sur le secteur $N$ a la multiplicité $d_{N,m}=\#\{\text{orbites }:mR_a\equiv0\pmod L\}$, et tout opérateur qui
commute avec $S^z$ et $T$ — par exemple la chaîne XXZ invariante par translation — est diagonal par blocs, dans
$\bigoplus_N\bigoplus_mM_{d_{N,m}}$.

*Démonstration.* (a) Le stabilisateur de $a$ dans le groupe cyclique $\{T^j\}\cong\mathbb Z_L$ est un sous-groupe, et
$R_a$ en est l'indice, un diviseur de $L$ ; les vecteurs $T^j|a\rangle$, $0\le j<R_a$, sont orthonormés et
$T^{R_a}|a\rangle=|a\rangle$. Si $e^{ikR_a}=1$, alors $T|a,k\rangle=e^{ik}\sum_{j=1}^{R_a}e^{-ikj}T^j|a\rangle=e^{ik}|a,k\rangle$,
puisque le terme $j=R_a$ est égal au terme $j=0$ ; le coefficient de $|a\rangle$ vaut $1$. Les
$m\in\{0,\dots,L-1\}$ admissibles sont les multiples de $L/R_a$, au nombre de $R_a$, avec des valeurs propres
distinctes de l'unitaire $T$ ; les vecteurs sont donc orthogonaux et engendrent l'espace de dimension $R_a$
engendré par l'orbite. (b) Le secteur $N$ est la somme orthogonale des espaces engendrés par ses orbites,
chacun invariant par $T$, et d'après (a) la valeur propre $e^{2\pi im/L}$ apparaît une fois dans l'espace engendré
par l'orbite de $a$ si $mR_a\equiv0\pmod L$, et pas du tout sinon. Un $X$
qui commute avec $S^z$ et $T$ commute avec leurs projecteurs spectraux $P_N$, $Q_m$ (théorème 1.5(d)), donc
$X=\sum_{N,m}P_NQ_mXP_NQ_m$. $\square$

`xxz_momentum` construit les blocs dans la base normalisée $|a,k\rangle$, avec un représentant $a$ par orbite. Nous
dénombrons les orbites indépendamment et comparons les spectres avec le secteur dense :

```python
L = 12
ring = SpinChain(L, boundary='periodic')
rot = lambda x: ((x << 1) | (x >> (L - 1))) & ((1 << L) - 1)        # T on labels: site i -> i + 1
for N in range(L + 1):
    orbits = {frozenset(itertools.accumulate(range(L - 1), lambda x, _: rot(x), initial=a))
              for a in ring.labels[N].tolist()}                       # {a, Ta, ..., T^{L-1} a}
    d = [sum((m * len(o)) % L == 0 for o in orbits) for m in range(L)]
    assert d == ring.momentum_dims(N)                                   # Proposition 9.5(b)
print(ring.momentum_dims(6))
Hk = ring.xxz_momentum(J=1.0, Delta=1.0, h=0.2, sector=6)               # 12 blocks instead of one of size 924
E_k = torch.sort(torch.cat([e[0] for e in Hk.eigenvalues()])).values
assert close(E_k, torch.linalg.eigvalsh(ring.xxz(1.0, 1.0, h=0.2, sector=6).matrix[0, 0]))
```

Le secteur $N=6$, de $924$ états, se découpe en blocs de tailles $80,75,78,76,78,75,80,75,78,76,78,75$.

**Lemme 9.6 (Perron–Frobenius, cas symétrique).** Soit $M$ une matrice réelle symétrique avec $M_{ij}\ge0$ pour
$i\ne j$, dont le graphe ($i\sim j$ si et seulement si $M_{ij}>0$) est connexe. Alors sa plus grande valeur propre
$\lambda$ est simple, avec un vecteur propre dont tous les coefficients sont $>0$.

*Démonstration.* $\lambda=\max_{\|x\|=1}\langle x,Mx\rangle$, et le maximum est atteint exactement sur les vecteurs
unitaires du sous-espace propre de $\lambda$ (développer dans une base orthonormée de vecteurs propres). Pour un tel
$x$, le vecteur $|x|$ des modules vérifie $\langle|x|,M|x|\rangle\ge\operatorname{Re}\langle x,Mx\rangle=\lambda$, car les
termes diagonaux coïncident et $M_{ij}|x_i||x_j|\ge\operatorname{Re}M_{ij}\bar x_ix_j$ pour $i\ne j$ ; donc $|x|$ est un
vecteur propre. Si $|x|_i=0$, alors $0=(\lambda-M_{ii})|x|_i=\sum_{j\ne i}M_{ij}|x|_j$ impose $|x|_j=0$ pour tous les
voisins $j$ de $i$, et par connexité $x=0$. Aucun vecteur propre non nul n'a donc de coefficient nul ; or un
sous-espace propre de dimension $\ge2$ contient un vecteur non nul avec $x_1=0$. Par conséquent $\lambda$ est simple,
et $|x|>0$ est un vecteur propre. $\square$

Nous démontrons la règle des signes pour toute disposition bipartie des couplages ; la chaîne en est un cas
particulier.

**Théorème 9.7 (règle des signes de Marshall).** Soit $A\sqcup B$ une partition des sites $\{0,\dots,L-1\}$ en deux
ensembles, et soit
$$
H=\sum_{i<j}\Big[\tfrac{J_{ij}}2\big(s^+_is^-_j+s^-_is^+_j\big)+K_{ij}\,s^z_is^z_j\Big]+\sum_ih_is^z_i
$$
avec $K_{ij},h_i$ réels et des $J_{ij}$ réels tels que $J_{ij}\ge0$ si $i,j$ sont dans des ensembles différents,
$J_{ij}\le0$ s'ils sont dans le même ensemble, et que le graphe $G$ sur les sites, avec une arête $i\sim j$ si et
seulement si $J_{ij}\ne0$, soit connexe. Soit $N_B(b)$ le nombre de spins vers le haut de $b$ sur $B$. Dans chaque
secteur $N$, la plus petite valeur propre de $H_N$ est simple, et son vecteur propre peut être choisi tel que
$(-1)^{N_B(b)}\psi(b)>0$ pour toute configuration $b$ à $N$ spins vers le haut.

En particulier, cela vaut pour la chaîne XXZ avec $J>0$ et $\Delta,h_i$ quelconques si la chaîne est ouverte, ou
périodique avec $L$ pair, avec pour $A$ les sites pairs et pour $B$ les sites impairs : alors $J_{ij}=J$,
$K_{ij}=J\Delta$ sur les liaisons et $0$ sinon, toute liaison relie un site pair et un site impair (pour l'anneau,
parce que $L$ est pair), et $G$ contient le chemin $0-1-\cdots-(L-1)$.

*Démonstration.* Soit $V|b\rangle=(-1)^{N_B(b)}|b\rangle$ ; $V$ est diagonal et $V=V^*=V^{-1}$. L'opérateur $s^+_is^-_j$
déplace un spin vers le haut de $j$ à $i$ ; il change $N_B$ de $\pm1$ si $i,j$ sont dans des ensembles différents
et le laisse inchangé sinon, donc $Vs^+_is^-_jV=\mp s^+_is^-_j$ selon le cas, et de même pour $s^-_is^+_j$. Deux
configurations distinctes $b,b'$ sont reliées par un tel terme si et seulement si $b'$ s'obtient à partir de $b$ en
échangeant les spins opposés de deux sites $i,j$, et alors $i,j$ sont déterminés par $b,b'$ (les sites où elles
diffèrent). Ainsi $M=-VH_NV$ a pour coefficients hors diagonale $J_{ij}/2$ (ensembles différents) ou $-J_{ij}/2$
(même ensemble), c'est-à-dire $|J_{ij}|/2\ge0$ dans les deux cas, strictement positifs si et seulement si $i\sim j$
dans $G$, et $0$ entre configurations qui diffèrent autrement. Le graphe de $M$ est connexe. En effet, notons
$(i\,j)$ l'échange des spins des sites $i,j$ ; pour $a,b,c$ distincts, $(a\,c)=(a\,b)(b\,c)(a\,b)$, donc, par
récurrence le long d'un chemin de $G$ de $i$ à $j$, tout $(i\,j)$ est une composée d'échanges le long d'arêtes de
$G$, et chacun de ces échanges ou bien fixe une configuration (spins égaux), ou bien suit une arête du graphe de
$M$. Si $b\ne b'$ ont $N$ spins vers le haut, il existe des sites $i$, vers le haut dans $b$ et vers le bas dans
$b'$, et $j$, vers le bas dans $b$ et vers le haut dans $b'$ ; $(i\,j)b$ diffère de $b'$ en deux sites de moins que
$b$, donc un nombre fini de telles étapes mène de $b$ à $b'$. D'après le lemme 9.6, la plus grande valeur propre de
$M$ — l'opposée de la plus petite de $H_N$ — est simple, avec un vecteur propre $\varphi>0$ ; on prend
$\psi=V\varphi$. $\square$

**Corollaire 9.8 (impulsion de l'état fondamental).** Sur une chaîne périodique avec $L$ pair, $J>0$ et un champ
uniforme, l'état fondamental du secteur $N$ a l'impulsion $k=0$ si $N$ est pair et $k=\pi$ si $N$ est impair.

*Démonstration.* $T$ commute avec $H_N$ et l'état fondamental $\psi=V\varphi$ est simple, donc $T\psi=e^{ik}\psi$. Or
$(T\psi)(b)=\psi(T^{-1}b)$, et $T^{-1}b$ porte sur le site $i$ le spin $b_{i+1}$ ; comme $L$ est pair, $i+1$ est pair
pour $i$ impair, donc $N_B(T^{-1}b)=N-N_B(b)$ et
$(VT\psi)(b)=(-1)^{N_B(b)}(-1)^{N-N_B(b)}\varphi(T^{-1}b)=(-1)^N\varphi(T^{-1}b)$. Par ailleurs $VT\psi=e^{ik}\varphi$.
$\varphi$ et $\varphi\circ T^{-1}$ sont tous deux positifs, donc $e^{ik}=(-1)^N$. $\square$

```python
for bc in ('open', 'periodic'):
    c = SpinChain(10, bc, complex_valued=False)
    for N in range(11):
        w, V = torch.linalg.eigh(c.xxz(1.0, 0.7, h=0.3 * torch.randn(10, dtype=f64), sector=N).matrix[0, 0])
        x = (-1.0) ** c.bits(N)[:, 1::2].sum(-1) * V[:, 0]                  # (-1)^{N_B(b)} psi(b)
        assert (x * x.sum().sign()).min() > 0 and (len(w) == 1 or w[1] - w[0] > 1e-8)   # Theorem 9.7
for N in range(1, L):
    lowest = [e[0, 0].item() for e in ring.xxz_momentum(1.0, 0.5, sector=N).eigenvalues()]
    m0 = int(np.argmin(lowest))
    assert m0 == (0 if N % 2 == 0 else L // 2)                              # Corollary 9.8: k = 0 or pi
    assert sorted(lowest)[1] - lowest[m0] > 1e-8
```

> **Interprétation physique** (interprétation physique, pas un théorème de ce cours). Dans un
> antiferromagnétique, le réseau fixe les signes de l'état fondamental : déplacer un spin à travers une liaison
> change le signe. Les physiciens s'en servent pour faire tourner le Monte-Carlo quantique sans problème du signe
> sur les antiferromagnétiques bipartis ; pour les chaînes frustrées, où un couplage antiferromagnétique relie
> deux sites du même ensemble, l'hypothèse du théorème 9.7 n'est pas satisfaite.

## 9.3 $SU(2)$ en $\Delta=1$

Pour $\Delta=1$ et $h=0$, $H=J\sum\mathbf s_i\cdot\mathbf s_j$ commute avec $S^x,S^y,S^z$ (proposition 4.15), donc
$H\cong\bigoplus_S1_{2S+1}\otimes H_S$ avec des blocs $H_S\in M_{m_S}$, $m_S=\binom L{L/2-S}-\binom L{L/2-S-1}$
(proposition 4.13), et le module `su2` construit $H_S$ sur les vecteurs de plus haut poids (proposition 4.15(b)).
Il en va de même pour tout $H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j$ avec des $J_{ij}$ réels : d'après la
proposition 4.15(a), chaque $\mathbf s_i\cdot\mathbf s_j$ appartient à $\{S_a\}'$. Notons $E_0(S)$ la plus petite valeur propre
de $H_S$, pour les spins $S$ qui apparaissent.

**Théorème 9.9 (Lieb–Mattis).** Soit $H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j$ avec $A,B,J_{ij}$ comme dans le
théorème 9.7 (donc $K_{ij}=J_{ij}$, $h_i=0$), $n_A=|A|\ge n_B=|B|$ et $S_*=\frac{n_A-n_B}2$.
(a) Pour $S_*\le S\le\frac L2$ (avec $S-\frac L2$ entier), l'état fondamental du secteur $S^z=S$ a le spin total $S$ ;
pour $0\le M\le S_*$ (avec $M-\frac L2$ entier), l'état fondamental du secteur $S^z=M$ a le spin total $S_*$.
(b) Par conséquent $E_0(S)$ est une valeur propre simple de $H_S$ pour $S\ge S_*$,
$E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)$, et $E_0(S)>E_0(S_*)$ pour $S<S_*$ : l'espace fondamental de $H$ est un unique
multiplet, de spin $S_*$.
Pour la chaîne de Heisenberg avec $J>0$, ouverte ou périodique avec $L$ pair, $A$ et $B$ sont les sites pairs et
impairs, $S_*=0$ pour $L$ pair et $S_*=\frac12$ pour $L$ impair (chaîne ouverte), et
$E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)$.

*Démonstration.* (a) Soient $M\ge0$ avec $M-\frac L2\in\mathbb Z$, $s=\max(M,S_*)$, et $\psi$ l'état fondamental du
secteur $S^z=M$, simple d'après le théorème 9.7 : $\psi=V\varphi$ avec $\varphi>0$. $\mathbf S^2$ commute avec $H$ et avec
$S^z$, donc il préserve l'espace fondamental de dimension un du secteur : $\psi$ a un spin total $S'$. Nous
construisons un vecteur $\chi$ du secteur avec $\mathbf S^2\chi=s(s+1)\chi$ et $\langle\psi,\chi\rangle\ne0$ ; comme des
vecteurs propres de $\mathbf S^2$ pour des valeurs propres distinctes sont orthogonaux, cela donne $S'=s$. Posons
$p=\frac L2-s$ ; alors $0\le p\le n_B$, car $s\ge S_*$ signifie $p\le\frac{n_A+n_B}2-\frac{n_A-n_B}2$. Choisissons des
$a_1,\dots,a_p\in A$ distincts et des $b_1,\dots,b_p\in B$ distincts, soit $R$ l'ensemble des $2s$ sites restants, et
soit
$$
\chi_0=\bigotimes_{l=1}^p\tfrac1{\sqrt2}\big(|{\uparrow}_{a_l}{\downarrow}_{b_l}\rangle-|{\downarrow}_{a_l}{\uparrow}_{b_l}\rangle\big)\otimes|{\uparrow\cdots\uparrow}\rangle_R ,
\qquad \chi=(S^-)^{s-M}\chi_0 .
$$
Chaque singulet est annulé par $s^\pm_{a_l}+s^\pm_{b_l}$ et les spins vers le haut de $R$ par $s^+_i$, donc
$S^+\chi_0=0$ ; avec $S^z\chi_0=s\chi_0$, cela donne $\mathbf S^2\chi_0=(S^-S^++(S^z)^2+S^z)\chi_0=s(s+1)\chi_0$
(section 4.5), et $\mathbf S^2\chi=s(s+1)\chi$ puisque $\mathbf S^2$ commute avec $S^-$. En écrivant
$S^-=\sum_l(s^-_{a_l}+s^-_{b_l})+S^-_R$, les termes commutent et les premiers annulent les singulets, donc $\chi$ est le
produit des singulets par $(S^-_R)^{s-M}|{\uparrow\cdots\uparrow}\rangle_R=(s-M)!\sum_D|D\rangle$, la somme portant sur les
configurations de $R$ ayant exactement $s-M$ spins vers le bas (il en existe, car $s-M\le2s=|R|$) ; ainsi $\chi\ne0$,
$S^z\chi=M\chi$, et les coefficients de $\chi$ sont réels. Son signe : un singulet contribue $+$ quand son spin vers
le haut est sur $A$ et $-$ quand il est sur $B$, et la partie sur $R$ est positive. Si $M=s$, les sites de $R$ sont
tous vers le haut ; si $M<s$, alors $s=S_*$, $p=n_B$ et $R\subseteq A$. Dans les deux cas, $N_B(b)$ est le nombre de
singulets dont le spin vers le haut est sur $B$, plus une constante, $|R\cap B|$ (si $M=s$) ou $0$, donc
$(-1)^{N_B(b)}\chi(b)$ a le même signe pour tous les $b$ tels que $\chi(b)\ne0$. Par conséquent
$\langle\psi,\chi\rangle=\langle\varphi,V\chi\rangle=\pm\sum_b\varphi(b)|\chi(b)|\ne0$.

(b) D'après la proposition 4.15(a), le corollaire 4.12 et le théorème 4.11, un multiplet de spin $S'$ rencontre le
secteur $S^z=M$ en un vecteur si $S'\ge|M|$ et pas du tout sinon, donc le spectre de $H$ sur le secteur $S^z=M\ge0$
est la réunion des spectres des $H_{S'}$, $S'\ge M$, avec multiplicités. D'après (a), sa plus petite valeur propre
est simple et appartient au spin $s=\max(M,S_*)$ : ainsi $E_0(s)$ est simple dans $H_s$ et $E_0(s)<E_0(S')$ pour tout
$S'\ge M$, $S'\ne s$. Pour $M=S\ge S_*$, cela donne la simplicité et $E_0(S)<E_0(S')$ pour $S'>S$ ; pour $M=S<S_*$,
cela donne $E_0(S_*)<E_0(S)$. Ainsi $E_0(S_*)$ est la plus petite valeur propre de $H$, atteinte seulement dans
$H_{S_*}$, où elle est simple. Pour la chaîne, les sites pairs sont plus nombreux que les impairs de $0$ ($L$ pair)
ou de $1$ ($L$ impair). $\square$

Le théorème 9.7 remonte à W. Marshall (*Proc. Roy. Soc. A* 232 (1955) 48) et le théorème 9.9 à E. Lieb et
D. Mattis (*J. Math. Phys.* 3 (1962) 749). Le **gap singulet–triplet** $E_0(1)-E_0(0)$ de l'anneau de Heisenberg
pair est positif d'après le théorème 9.9.

Nous comparons les tailles des blocs pour douze spins, vérifions le théorème 9.9 à l'aide des secteurs $U(1)$,
calculons le gap et vérifions les théorèmes 9.7 et 9.9 sur un graphe biparti qui n'est pas une chaîne ($n_A=5$,
$n_B=3$, $S_*=1$, couplages aléatoires) :

```python
sym = SU2Chain(12, boundary='periodic', complex_valued=False)
H_su2 = sym.heisenberg(J=1.0)                                       # (+)_S M_{m_S}, S = 0, ..., 6
u1_sizes = [math.comb(12, N) for N in range(13)]
print(sym.multiplicities, max(u1_sizes), sum(k ** 3 for k in u1_sizes), sum(m ** 3 for m in sym.multiplicities))
E0_S = sym.lowest_energies(H_su2)[0]                                # E_0(S), S = 0, ..., 6
E0_N = [torch.linalg.eigvalsh(ring.xxz(1.0, 1.0, sector=6 + S).matrix[0, 0])[0] for S in range(7)]
assert close(E0_S, torch.stack(E0_N))                               # Theorem 9.9: sector S^z = S has spin S
assert (E0_S[1:] > E0_S[:-1] + 1e-8).all()                          # E_0(S) strictly increasing
gap = {}
for n in (6, 8, 10, 12):
    s = SU2Chain(n, 'periodic', complex_valued=False)
    e = [torch.linalg.eigvalsh(s.heisenberg(spin=S).matrix[0, 0])[0].item() for S in (0, 1)]
    gap[n] = n * (e[1] - e[0])
print({n: round(g, 4) for n, g in gap.items()}, round(math.pi ** 2 / 2, 4))
assert 0 < gap[6] < gap[8] < gap[10] < gap[12] < math.pi ** 2 / 2
```

Les blocs de $SU(2)$ ont pour tailles $m_S=132,297,275,154,54,11,1$, contre $924$ pour le plus grand bloc de
$U(1)$, et le coût $\sum k^3$ de la diagonalisation tombe de $2.0\cdot10^9$ à $5.3\cdot10^7$. *Observation numérique :*
$L$ fois le gap croît lentement, $4.108,4.181,4.232,4.270$ pour $L=6,\dots,12$ ; nous ne démontrons rien sur son
comportement quand $L\to\infty$.

**Corrélations.** Les spins des sites forment des opérateurs vectoriels, et leurs composantes sphériques des
opérateurs tensoriels de rang un (lemme B.15(a),(b)). D'après le **théorème de Wigner–Eckart** (théorème B.16), leurs
éléments de matrice entre deux multiplets sont un coefficient de Clebsch–Gordan multiplié par un **élément de
matrice réduit** qui ne dépend pas des nombres quantiques magnétiques, et
$\mathbf s_i\cdot\mathbf s_j=\sum_q(-1)^qs^{(i)}_qs^{(j)}_{-q}$ (lemme B.15(a)). `spin_correlation` s'en sert pour construire
$\mathbf s_i\cdot\mathbf s_j$ sur chaque $H_S$ à partir des seuls vecteurs de plus haut poids ; la méthode `'direct'`
comprime à la place l'opérateur de liaison sur les vecteurs de plus haut poids, ce qui est correct d'après la
proposition 4.15(b) puisque $\mathbf s_i\cdot\mathbf s_j\in\{S_a\}'$. La correction de l'implémentation de la première
méthode n'est pas démontrée ici ; le code vérifie que les deux concordent entre elles et avec la définition,
évaluées sur l'état fondamental écrit dans le secteur $U(1)$ :

```python
gs = sym.ground_state(H_su2)
S0, coef = gs.spin.item(), gs.coefficients[0, :sym.multiplicity(gs.spin.item())]
pairs = [(0, r) for r in range(12)]
corr_we = sym.expectation(sym.spin_correlation(pairs), S0, coef)
corr_dir = sym.expectation(sym.spin_correlation(pairs, method='direct'), S0, coef)
psi = sym.embed(S0, coef)                                            # the ground state in the sector N = 6
lab, z = sym.u1.labels[6], sym.u1.bits(6).to(f64) - 0.5
def s_dot_s(i, j):                                                   # <psi| s_i . s_j |psi>, from the definition
    if i == j:
        return 0.75
    flip = ((lab >> i) & 1) != ((lab >> j) & 1)
    partner = sym.u1.index_of[lab[flip] ^ ((1 << i) | (1 << j))]
    return ((psi ** 2 * z[:, i] * z[:, j]).sum() + 0.5 * (psi[flip] * psi[partner]).sum()).item()
assert S0 == 0 and close(corr_we, corr_dir)
assert close(corr_we, torch.tensor([s_dot_s(0, r) for r in range(12)], dtype=f64))
assert ((-1.0) ** torch.arange(1, 12) * corr_we[1:] > 0).all()       # signs (-1)^r
print([round(x, 4) for x in corr_we.tolist()])
```

Les corrélations $\langle\mathbf s_0\cdot\mathbf s_r\rangle=0.75,-0.4489,0.1878,-0.1658,0.1209,-0.1224,0.1070,\dots$
alternent en signe, comme l'imposent les théorèmes 9.7 et 9.9 : le singulet est invariant par les rotations, qui
permutent les axes, et $\langle s^+_0s^+_r\rangle=0$ dans le secteur $S^z=0$, donc $\langle\mathbf s_0\cdot\mathbf s_r\rangle=3\langle s^x_0s^x_r\rangle=
\frac34\langle s^+_0s^-_r+s^-_0s^+_r\rangle=\frac34\sum_{b_0\ne b_r}\psi(b)\psi(b')$, où $b'$ s'obtient à partir de $b$ en
échangeant les spins des sites $0$ et $r$ ; comme $N_B(b')-N_B(b)$ vaut $\pm1$ pour $r$ impair et $0$ pour $r$ pair,
chaque terme a le signe $(-1)^r$. La décroissance est lente, comme $(-1)^r\sqrt{\log r}/r$ quand $L\to\infty$
(Affleck et al., *loc. cit.* ; énoncé sans démonstration).

## 9.4 Vérifications exactes : fermions libres et ansatz de Bethe

**Jordan–Wigner.** Sur une chaîne ouverte, posons $c_j=\big(\prod_{l<j}(-1)^{n_l}\big)s^-_j$.

**Théorème 9.10 (la chaîne XX est libre).** (a) $\{c_i,c_j^*\}=\delta_{ij}$, $\{c_i,c_j\}=0$ et $c_j^*c_j=n_j$.
(b) La chaîne XX ouverte s'écrit $H=\frac J2\sum_{j<L-1}(s^+_js^-_{j+1}+s^-_js^+_{j+1})=\sum_{i,j}h_{ij}c_i^*c_j$ avec
$h=\frac J2A_L$. (c) Soient $u_k$ des vecteurs propres réels orthonormés de $h$, de valeurs propres
$\varepsilon_k=J\cos\frac{\pi k}{L+1}$ (lemme 9.3), et $b_k=\sum_ju_k(j)c_j$. Pour $K=\{k_1<\dots<k_N\}$, les vecteurs
$b_K^*|0\rangle=b_{k_1}^*\cdots b_{k_N}^*|0\rangle$ ($|0\rangle$ : tous les spins vers le bas) forment une base orthonormée
de vecteurs propres du secteur $N$, d'énergies $\sum_{k\in K}\varepsilon_k$.

*Démonstration.* (a) Des opérateurs sur des sites différents commutent, et sur un même site $(-1)^n$ anticommute
avec $s^\pm$. Pour $i<j$, $c_i$ et $c_j^*$ (ou $c_j$) se rencontrent sur le site $i$ sous la forme $s^-_i$ et
$(-1)^{n_i}$ — une anticommutation — et commutent par ailleurs, donc ils anticommutent. Pour $i=j$, les chaînes de
signes sont de carré $1$ et $c_jc_j^*+c_j^*c_j=s^-_js^+_j+s^+_js^-_j=1$, $c_j^*c_j=s^+_js^-_j=n_j$. (b)
$c_j^*c_{j+1}=s^+_j(-1)^{n_j}s^-_{j+1}=s^+_js^-_{j+1}$, car $s^+_j\neq0$ seulement sur $n_j=0$ ; on passe aux adjoints.
(c) $\{b_k,b_l^*\}=\sum_ju_k(j)u_l(j)=\delta_{kl}$ et $\{b_k,b_l\}=0$, et $H=\sum_k\varepsilon_kb_k^*b_k$ par la
décomposition spectrale de $h$. D'après les CAR, $[b_l^*b_l,b_k^*]=\delta_{lk}b_k^*$, donc $[H,b_k^*]=\varepsilon_kb_k^*$ ;
avec $H|0\rangle=0$, cela donne $Hb_K^*|0\rangle=(\sum_K\varepsilon_k)b_K^*|0\rangle$. Comme $b_k|0\rangle=0$ (chaque $c_j$
abaisse un spin), faire passer les annihilateurs à droite à l'aide des CAR donne
$\langle b_{K'}^*0,b_K^*0\rangle=\delta_{KK'}$ par récurrence sur $N$. Chaque $c_j^*$ augmente $N$ de un, donc
$b_K^*|0\rangle$ appartient au secteur $|K|$ : $\binom LN$ vecteurs orthonormés dans un espace de cette dimension. $\square$

```python
L = 10
xx = SpinChain(L)
eps = [math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1)]
for N in range(L + 1):
    E = torch.linalg.eigvalsh(xx.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
    E_ff = sorted(sum(K) for K in itertools.combinations(eps, N))     # sums over N-element subsets
    assert close(E, torch.tensor(E_ff, dtype=f64))                    # all 1024 levels, Theorem 9.10
    assert abs(E[0].item() - sum(sorted(eps)[:N])) < 1e-12           # ground state: the N lowest modes
```

**Magnons.** Sur l'anneau avec $h=0$, l'état complètement polarisé $|0\rangle$ a l'énergie $E_F=J\Delta L/4$. Un spin
vers le haut en $x$ brise deux liaisons : $H|x\rangle=(E_F-J\Delta)|x\rangle+\frac J2(|x-1\rangle+|x+1\rangle)$ ($L\ge3$).

**Proposition 9.11 (un magnon).** Dans le secteur $N=1$ de l'anneau, $|k\rangle=\sum_xe^{-ikx}|x\rangle$,
$k=2\pi m/L$, vérifie $T|k\rangle=e^{ik}|k\rangle$ et $H|k\rangle=\big(E_F+J(\cos k-\Delta)\big)|k\rangle$.

*Démonstration.* $\sum_xe^{-ikx}(|x-1\rangle+|x+1\rangle)=\sum_y(e^{-ik(y+1)}+e^{-ik(y-1)})|y\rangle=2\cos k\,|k\rangle$, en
utilisant $e^{ikL}=1$ ; le même changement d'indice donne $T|k\rangle=e^{ik}|k\rangle$. $\square$

Deux magnons diffusent l'un sur l'autre. Notons $|x_1,x_2\rangle$, $x_1<x_2$, l'état aux spins vers le haut en
$x_1,x_2$ ; une paire de voisins brise deux liaisons, une paire séparée quatre.

**Proposition 9.12 (deux magnons ; équations de Bethe).** Soient $k_1,k_2\in\mathbb C$, $K=k_1+k_2$, et
$(A_{12},A_{21})\neq0$ tels que (i) $A_{12}(1+e^{iK}-2\Delta e^{ik_2})+A_{21}(1+e^{iK}-2\Delta e^{ik_1})=0$ et
(ii) $A_{21}=e^{ik_2L}A_{12}$, $A_{12}=e^{ik_1L}A_{21}$. Soient
$\psi(x_1,x_2)=A_{12}e^{i(k_1x_1+k_2x_2)}+A_{21}e^{i(k_2x_1+k_1x_2)}$ sur $\mathbb Z^2$ et
$\Psi=\sum_{0\le x_1<x_2<L}\psi(x_1,x_2)|x_1,x_2\rangle$. Si $\Psi\neq0$, alors
$H\Psi=\big(E_F+J(\cos k_1+\cos k_2-2\Delta)\big)\Psi$. Si les dénominateurs ne s'annulent pas, (i) et (ii)
s'écrivent
$$
e^{ik_1L}=-\frac{1+e^{iK}-2\Delta e^{ik_1}}{1+e^{iK}-2\Delta e^{ik_2}},\qquad
e^{ik_2L}=-\frac{1+e^{iK}-2\Delta e^{ik_2}}{1+e^{iK}-2\Delta e^{ik_1}} .
$$

*Démonstration.* Soit $E=E_F+J(\cos k_1+\cos k_2-2\Delta)$. (1) Pour tout $(x_1,x_2)\in\mathbb Z^2$, chaque onde plane
donne des valeurs $\psi(x_1\pm1,x_2)$ et $\psi(x_1,x_2\pm1)$ de somme $2(\cos k_1+\cos k_2)\psi(x_1,x_2)$. (2) Par
substitution, $\psi(x,x)+\psi(x+1,x+1)-2\Delta\psi(x,x+1)$ vaut $e^{iKx}$ fois le membre de gauche de (i), donc
s'annule pour tout $x$. (3) D'après (ii),
$\psi(x_2,x_1+L)=A_{12}e^{ik_2L}e^{i(k_2x_1+k_1x_2)}+A_{21}e^{ik_1L}e^{i(k_1x_1+k_2x_2)}=\psi(x_1,x_2)$, et en appliquant
ceci deux fois, $\psi(x_1+L,x_2+L)=\psi(x_1,x_2)$. Les représentants d'une configuration $\{y_1<y_2\}$ dans
$D=\{x_1<x_2<x_1+L\}$ sont $(y_1+nL,y_2+nL)$ et $(y_2+nL,y_1+(n+1)L)$, donc $\psi$ est une fonction des
configurations. Soit $(x_1,x_2)\in D$. Si $2\le x_2-x_1\le L-2$, les quatre sauts restent dans $D$ et le terme
diagonal vaut $E_F-2J\Delta$, donc $(H\Psi)(x_1,x_2)=(E_F-2J\Delta)\psi+\frac J2\cdot2(\cos k_1+\cos k_2)\psi=E\psi$
d'après (1). Si $x_2=x_1+1$, les sauts $(x_1+1,x_2)$, $(x_1,x_2-1)$ sont bloqués et le terme diagonal vaut
$E_F-J\Delta$ ; d'après (1) et (2),
$(H\Psi)(x_1,x_2)=(E_F-J\Delta)\psi+\frac J2\big(2(\cos k_1+\cos k_2)\psi-\psi(x_1,x_1)-\psi(x_1+1,x_1+1)\big)=E\psi$.
Si $x_2=x_1+L-1$, on utilise le représentant $(x_2,x_1+L)$, qui est du type précédent ($L\ge3$). $\square$

La proposition 9.12 dit seulement que les solutions donnent des vecteurs propres ; elle ne dit pas que tous les
vecteurs propres s'obtiennent ainsi. Écrivons $k_{1,2}=K/2\pm q$ et $c=\cos\frac K2$. Alors
$1+e^{iK}-2\Delta e^{ik_{1,2}}=2e^{iK/2}(c-\Delta e^{\pm iq})$, donc $A_{12}=c-\Delta e^{iq}$,
$A_{21}=-(c-\Delta e^{-iq})$ résolvent (i). Le produit des deux équations (ii) donne $e^{iKL}=1$, $K=2\pi n/L$, et
(ii) équivaut alors à $F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}=0$, une équation à une variable. Nous
cherchons ses racines $q\in(0,\pi)$, où $e^{-iqL/2}F(q)$ est réel ou imaginaire pur parce que $c-\Delta e^{\pm iq}$
sont complexes conjugués (états de diffusion), et sur les droites $q=i\kappa$, $\pi+i\kappa$, $\kappa>0$, où $F$ est
réel (états liés) ; pour $L=9$, $\Delta=0.6$, nous construisons $\Psi$ pour chaque racine et comparons avec le
secteur $N=2$ de $36$ états :

```python
L, Delta = 9, 0.6
ring9 = SpinChain(L, boundary='periodic')
H2 = ring9.xxz(J=1.0, Delta=Delta, sector=2).matrix[0, 0]
x1, x2 = np.triu_indices(L, 1)                                       # all pairs x1 < x2
rows = ring9.index_of[torch.tensor((1 << x1) | (1 << x2))]
def roots(g, grid):                                                  # sign changes of g on the grid, refined
    v = np.array([g(t) for t in grid])
    return [brentq(g, grid[i], grid[i + 1], xtol=1e-15) for i in np.nonzero(v[:-1] * v[1:] < 0)[0]]
E_bethe, vecs, n_bound = [], [], 0
for n in range(L):
    K, c = 2 * np.pi * n / L, math.cos(np.pi * n / L)
    F = lambda q: (-1) ** n * np.exp(1j * q * L) * (c - Delta * np.exp(-1j * q)) + c - Delta * np.exp(1j * q)
    G = lambda q: (lambda z: z.real + z.imag)(np.exp(-0.5j * q * L) * F(q))    # real, or imaginary
    bound = [s + 1j * t for s in (0, np.pi)                          # q = i kappa and pi + i kappa
             for t in roots(lambda t: F(s + 1j * t).real, np.linspace(1e-3, 20, 4000))]
    n_bound += len(bound)
    for q in roots(G, (np.arange(1, 2000) + 0.5 ** 0.5) * np.pi / 2000) + bound:
        psi = np.exp(0.5j * K * (x1 + x2)) * ((c - Delta * np.exp(1j * q)) * np.exp(1j * q * (x1 - x2))
                                              - (c - Delta * np.exp(-1j * q)) * np.exp(-1j * q * (x1 - x2)))
        v = torch.zeros(36, dtype=c128)
        v[rows] = torch.tensor(psi / np.linalg.norm(psi))
        E = 9 * Delta / 4 + 2 * (c * np.cos(q) - Delta)                 # E_F + J (cos k1 + cos k2 - 2 Delta)
        assert abs(E.imag) < 1e-12 and close(H2 @ v, E.real * v, 1e-9)   # Proposition 9.12
        E_bethe.append(E.real)
        vecs.append(v)
assert torch.linalg.matrix_rank(torch.stack(vecs)) == 36             # a complete eigenbasis
assert close(torch.sort(torch.tensor(E_bethe)).values, torch.linalg.eigvalsh(H2))
print(len(E_bethe), n_bound)
ring8, Phi = SpinChain(8, 'periodic'), torch.zeros(28, dtype=c128)     # even L: Phi = sum_x (-1)^x |x, x+1>
for x in range(8):
    Phi[ring8.index_of[(1 << x) | (1 << (x + 1) % 8)]] = (-1) ** x
assert close(ring8.xxz(1.0, Delta, sector=2).matrix[0, 0] @ Phi, (8 * Delta / 4 - Delta) * Phi)   # E_F - J Delta
```

Pour $L=9$, les $36$ états propres proviennent tous des équations de Bethe : $34$ états de diffusion et $2$ états
liés. La complétude est vérifiée ici, non démontrée, et elle est en défaut pour $L$ pair. Si $\Delta\ne0$,
$\Phi=\sum_x(-1)^x|x,x+1\rangle$ vérifie $T\Phi=-\Phi$ et $H\Phi=(E_F-J\Delta)\Phi$ (les sauts vers $|x-1,x+1\rangle$ et
$|x,x+2\rangle$ se compensent deux à deux). Un vecteur de Bethe vérifie $T\Psi=e^{-iK}\Psi$, puisque
$\psi(x_1-1,x_2-1)=e^{-iK}\psi(x_1,x_2)$ ; si $e^{iK}=-1$, alors $\cos k_1+\cos k_2=2\cos\frac K2\cos\frac{k_1-k_2}2=0$ et
son énergie vaut $E_F-2J\Delta$. Ainsi $\Phi$ est orthogonal à tous les vecteurs de Bethe (les dernières lignes du
code vérifient $H\Phi$ pour $L=8$). Avec $\Delta=0.6$ et $4\le L\le13$, la recherche ci-dessus trouve tous les états
pour $L$ impair, et tous les états sauf $\Phi$ pour $L$ pair. Pour $N>2$ magnons, l'ansatz comporte $N!$ termes et
$N$ équations de Bethe couplées (H. Bethe, *Z. Phys.* 71 (1931) 205 ; M. Karbach, G. Müller, *Computers in
Physics* 11 (1997) 36).

> **Interprétation physique.** Un magnon est une onde de spin : un spin retourné délocalisé sur l'anneau avec
> l'impulsion $k$. Deux magnons diffusent avec un déphasage ($q$ réel) ou forment un état lié ($q$ complexe).
> L'idée de Bethe fut que, pour $N$ magnons, la diffusion se factorise en phases à deux corps ; ses équations
> donnent par exemple l'énergie fondamentale par site $\frac14-\log2$ de la chaîne de Heisenberg infinie
> (L. Hulthén, *Ark. Mat. Astron. Fys.* 26A (1938) no 11 ; énoncé sans démonstration).

## 9.5 Lanczos et méthodes de Krylov

Le secteur $N=L/2$ est de dimension $12870$ pour $L=16$ et $184756$ pour $L=20$, mais $H$ est creuse ($O(L)$
coefficients non nuls par ligne, `xxz_sparse`), et les méthodes de Krylov n'ont besoin que de produits $Hv$. Soient
$H$ hermitienne de taille $n$ et $v\ne0$. L'**espace de Krylov** est
$\mathcal K_m=\operatorname{span}\{v,Hv,\dots,H^{m-1}v\}$.

**Proposition 9.13 (Lanczos).** Supposons $\dim\mathcal K_m=m$. Soit $q_1=v/\|v\|$ et, pour $j<m$, soit $q_{j+1}$
la composante normalisée de $Hq_j$ orthogonale à $\mathcal K_j$. Alors $q_1,\dots,q_m$ est une base orthonormée de
$\mathcal K_m$ ; avec $Q_m=(q_1\cdots q_m)$, $T_m=Q_m^*HQ_m$ est réelle, symétrique et tridiagonale, de diagonale
$\alpha_j=\langle q_j,Hq_j\rangle$ et de sous-diagonale $\beta_j=\langle q_{j+1},Hq_j\rangle>0$ ; et $HQ_m=Q_mT_m+r_me_m^T$,
où $r_m\perp\mathcal K_m$ s'annule si et seulement si $\mathcal K_{m+1}=\mathcal K_m$ (sinon $r_m=\beta_mq_{m+1}$,
$q_{m+1}=r_m/\|r_m\|$).

*Démonstration.* Supposons $\operatorname{span}\{q_1,\dots,q_j\}=\mathcal K_j$ (vrai pour $j=1$). Alors
$H\mathcal K_j\subseteq\mathcal K_j+\mathbb CHq_j$, puisque $Hq_i\in\mathcal K_{i+1}$ ; donc
$\mathcal K_{j+1}=\mathcal K_j+\mathbb CH^jv\subseteq\mathcal K_j+\mathbb CHq_j\subseteq\mathcal K_{j+1}$, et
$\mathcal K_{j+1}=\mathcal K_j+\mathbb CHq_j$. Pour $j<m$ la dimension croît, donc la composante de $Hq_j$ orthogonale à
$\mathcal K_j$ est non nulle et la récurrence se poursuit ; pour $j=m$ cette composante, $r_m$, s'annule si et
seulement si $\mathcal K_{m+1}=\mathcal K_m$. Pour $i>j+1$, $q_i\perp\mathcal K_{j+1}\ni Hq_j$, donc $(T_m)_{ij}=0$, et
$(T_m)_{ji}=0$ puisque $T_m$ est hermitienne. $\alpha_j$ est réel, et $\beta_j=\langle q_{j+1},Hq_j\rangle$ est la norme de
la composante, $>0$. Ainsi $Hq_j=\beta_{j-1}q_{j-1}+\alpha_jq_j+\beta_jq_{j+1}$ pour $j<m$ et
$Hq_m=\beta_{m-1}q_{m-1}+\alpha_mq_m+r_m$. $\square$

Les valeurs propres $\theta^{(m)}_1\le\dots\le\theta^{(m)}_m$ de $T_m$ sont les **valeurs de Ritz**.

**Théorème 9.14 (convergence monotone).** Soient $\lambda_1\le\dots\le\lambda_n$ les valeurs propres de $H$ et
$\dim\mathcal K_m=m$. Alors $\theta^{(m)}_j\ge\lambda_j$, et $\theta^{(m)}_j\ge\theta^{(m+1)}_j$ si
$\dim\mathcal K_{m+1}=m+1$ ($j\le m$). Si $v$ a une composante non nulle dans le sous-espace propre de $\lambda_1$ et
si $\mathcal K_{m+1}=\mathcal K_m$, alors $\theta^{(m)}_1=\lambda_1$.

*Démonstration.* D'après Courant–Fischer (algèbre linéaire standard),
$\theta^{(m)}_j=\min_{U}\max_{x\in U,\|x\|=1}\langle x,Hx\rangle$, le minimum portant sur les $U\subseteq\mathcal K_m$ de
dimension $j$, puisque $Q_m$ est une isométrie sur $\mathcal K_m$ et $\langle y,T_my\rangle=\langle Q_my,HQ_my\rangle$. La
même formule sur $U\subseteq\mathcal K_{m+1}$ donne $\theta^{(m+1)}_j$ et sur $U\subseteq\mathbb C^n$ donne $\lambda_j$ ;
minimiser sur davantage de sous-espaces donne une valeur plus petite. Si $\mathcal K_{m+1}=\mathcal K_m$, alors
$\mathcal K_m$ est $H$-invariant et contient $p(H)v$ pour tout polynôme $p$, en particulier $P_1v\ne0$ ($P_1$ le
projecteur spectral de $\lambda_1$, un polynôme en $H$ par interpolation de Lagrange, comme dans le théorème 1.5) ;
donc $\theta^{(m)}_1\le\langle x,Hx\rangle=\lambda_1$ en $x=P_1v/\|P_1v\|\in\mathcal K_m$. $\square$

La vitesse de convergence $\theta^{(m)}_1\to\lambda_1$ dépend du gap $\lambda_2-\lambda_1$ rapporté à l'étendue du spectre
(bornes de Kaniel–Paige–Saad via les polynômes de Tchebychev ; Y. Saad, *Numerical Methods for Large Eigenvalue
Problems*, SIAM 2011, chap. 6 ; énoncé sans démonstration). En virgule flottante, les $q_j$ perdent leur
orthogonalité ; `krylov.lanczos` réorthogonalise complètement. Pour la chaîne XX, l'énergie exacte de l'état
fondamental est connue (théorème 9.10) :

```python
L, N = 12, 6
xx = SpinChain(L, complex_valued=False)
Hs = xx.xxz_sparse(J=1.0, Delta=0.0, sector=N)                      # sparse: diagonal + hopping
alpha, beta, Q = krylov.lanczos(Hs.matvec, torch.randn(1, xx.sector_dim(N), dtype=f64), 40)
T = torch.diag(alpha[0]) + torch.diag(beta[0], 1) + torch.diag(beta[0], -1)
assert close(Q[0] @ Q[0].T, torch.eye(40, dtype=f64)) and close(Q[0] @ Hs.to_dense()[0] @ Q[0].T, T)
ritz = torch.tensor([torch.linalg.eigvalsh(T[:m, :m])[0].item() for m in range(1, 41)])   # T_m = leading block
E_xx = sum(sorted(math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1))[:N])
assert (ritz[1:] <= ritz[:-1] + 1e-12).all() and (ritz >= E_xx - 1e-12).all()             # Theorem 9.14
print([f"{(ritz[m - 1] - E_xx).item():.1e}" for m in (5, 10, 20, 40)])
```

L'erreur de la plus petite valeur de Ritz décroît de façon monotone, $0.42$, $3.1\cdot10^{-3}$ et $3.1\cdot10^{-6}$
après $5$, $10$ et $20$ pas, et atteint le niveau des erreurs d'arrondi après $40$, en dimension $924$.
`krylov.ground_state` redémarre Lanczos jusqu'à ce que le résidu $\|H\psi-E\psi\|$ soit inférieur à $10^{-10}$. Nous
le vérifions sur les énergies de l'anneau de Heisenberg tabulées à neuf chiffres dans `scripts/chains.py` (les
valeurs de référence de l'article compagnon), en dense jusqu'à $L=12$ et par Lanczos pour $L=16$ :

```python
RING = {4: -2.0, 6: -2.802775638, 8: -3.651093409, 10: -4.515446354, 12: -5.387390917,
        14: -6.263549534, 16: -7.142296361}
for n in (4, 6, 8, 10, 12):
    e0 = torch.linalg.eigvalsh(SpinChain(n, 'periodic').heisenberg(sector=n // 2).matrix[0, 0])[0].item()
    assert abs(e0 - RING[n]) < 5e-10
    print(n, round(e0 / n, 6))
t0 = time.time()
E16, psi16, res = krylov.ground_state(SpinChain(16, 'periodic', complex_valued=False).xxz_sparse(sector=8),
                                      generator=torch.Generator().manual_seed(0))
print(f"L = 16: E0 = {E16.item():.10f}, residual {res.item():.1e}, {time.time() - t0:.2f} s; "
      f"infinite chain: {0.25 - math.log(2):.6f}")
assert abs(E16.item() - RING[16]) < 5e-10 and res.item() < 1e-9
```

Pour $L=16$ (dimension $12870$), Lanczos prend nettement moins d'une seconde sur un seul thread CPU. L'énergie par
site, $-0.5,-0.467129,-0.456387,-0.451545,-0.448949$ pour $L=4,\dots,12$, s'approche par en dessous de la valeur
$\frac14-\log2=-0.443147$ de la chaîne infinie.

**Évolution temporelle.** Pour $e^{-i\tau H}v$, on projette sur $\mathcal K_m$ : $w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1$.

**Proposition 9.15.** Soit $\dim\mathcal K_m=m$. Pour tout $\tau$, $\|w_\tau\|=\|v\|$ et
$\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle$. Si $\mathcal K_m$ est $H$-invariant, $w_\tau=e^{-i\tau H}v$.

*Démonstration.* $Q_m$ est une isométrie et $e^{-i\tau T_m}$ est unitaire ($T_m$ est réelle symétrique), ce qui donne
la norme. Comme $Q_m^*HQ_m=T_m$ commute avec $e^{-i\tau T_m}$,
$\langle w_\tau,Hw_\tau\rangle=\|v\|^2\langle e^{-i\tau T_m}e_1,T_me^{-i\tau T_m}e_1\rangle=\|v\|^2\langle e_1,T_me_1\rangle=\langle v,Hv\rangle$. Si
$\mathcal K_m$ est invariant, $HQ_m=Q_mT_m$ (proposition 9.13, $r_m=0$), donc $H^jQ_m=Q_mT_m^j$ et
$e^{-i\tau H}Q_m=Q_me^{-i\tau T_m}$ ; on applique à $e_1$, avec $Q_me_1=v/\|v\|$. $\square$

Ainsi chaque pas de Krylov conserve la norme et l'énergie, qu'il soit précis ou non ; la précision exige $m$ grand
devant $\tau\|H\|$, et `krylov.evolve` subdivise l'intervalle de temps et estime l'erreur a posteriori
(M. Hochbruck, C. Lubich, *SIAM J. Numer. Anal.* 34 (1997) 1911 ; énoncé sans démonstration). Nous effectuons une
trempe (*quench*) à partir de l'état de Néel $|{\uparrow\downarrow\uparrow\downarrow}\cdots\rangle$ sur l'anneau de
Heisenberg et suivons l'aimantation alternée $m_s=\frac2L\sum_i(-1)^i\langle s^z_i\rangle$, en comparant avec le
propagateur dense pour $L=12$ :

```python
times, m_s = [0.0, 0.5, 1.0, 2.0, 4.0], {}
stag = lambda c, N, out: [((p[0].abs() ** 2) @ (c.bits(N).to(f64) - 0.5) @ torch.tensor(
    [(-1.0) ** i for i in range(c.L)], dtype=f64) * 2 / c.L).item() for p in out]
for n in (12, 16):
    c = SpinChain(n, 'periodic')
    Hn = c.xxz_sparse(J=1.0, Delta=1.0, sector=n // 2)
    psi0 = c.vector_in_sector('10' * (n // 2))[None]
    out = krylov.evolve(Hn, psi0, times)                               # (times, batch, dim)
    energy = torch.stack([(p.conj() * Hn.matvec(p)).sum(-1).real for p in out])
    assert close(out.norm(dim=-1), 1.0) and close(energy, energy[0])   # Proposition 9.15
    m_s[n] = stag(c, n // 2, out)
    print(n, [round(x, 4) for x in m_s[n]])
    if n == 12:
        w, U = torch.linalg.eigh(c.xxz(1.0, 1.0, sector=6).matrix[0, 0])
        for t, p in zip(times, out):
            assert close(p[0], U @ (torch.exp(-1j * t * w) * (U.mH @ psi0[0])))
print(f"{max(abs(a - b) for a, b in zip(m_s[12][:4], m_s[16][:4])):.1e}")   # L = 12 vs 16 up to t = 2
```

$m_s$ décroît de $1$ à $-0.183$ en $t=2$, en passant par $0.770$ et $0.279$, pour les deux tailles à
$1.7\cdot10^{-5}$ près ; en $t=4$, elles diffèrent ($0.026$ contre $0.035$) : la taille finie du plus petit anneau
devient visible.

## 9.6 Désordre : statistiques de niveaux et intrication

La **chaîne de Heisenberg en champ aléatoire** a des champs indépendants $h_i$ uniformes sur $[-W,W]$. Elle conserve
$S^z$ (proposition 9.1) et génériquement rien d'autre : les champs brisent les translations, les réflexions, le
retournement des spins et $SU(2)$. Les statistiques de niveaux se calculent donc par secteur $N$ (théorème 5.15 et
section 5.4 : mélanger les secteurs simule des niveaux indépendants). Pour un désordre faible, les niveaux se
repoussent comme dans le GOE, $\langle r\rangle\approx0.5307$ (section 5.3) ; pour un désordre fort, les états
propres se localisent (**localisation à N corps**) et les niveaux deviennent indépendants,
$\langle r\rangle=2\log2-1\approx0.3863$ (proposition 5.14). Nous calculons $\langle r\rangle$ dans la moitié centrale du
secteur $N=L/2$, avec une erreur standard sur les réalisations indépendantes (chacune donne un rapport moyen) :

```python
gen = torch.Generator().manual_seed(1)
r_mean = {}
for n, R in ((8, 200), (10, 100)):
    c = SpinChain(n, 'periodic', complex_valued=False)
    k = c.sector_dim(n // 2)
    for W in (0.5, 1, 2, 4, 8, 16):
        E = torch.linalg.eigvalsh(c.random_field_heisenberg(W, R, sector=n // 2, generator=gen).matrix[:, 0])
        r = SpinChain.level_spacing_ratio(E[:, k // 4: 3 * k // 4])     # one value per realisation
        r_mean[n, W] = (r.mean().item(), (r.std() / math.sqrt(R)).item())
    print(n, '  '.join(f"W={W}: {m:.3f}+-{s:.3f}" for (n_, W), (m, s) in r_mean.items() if n_ == n))
assert abs(r_mean[10, 0.5][0] - 0.5307) < 0.015 and abs(r_mean[10, 16][0] - 0.3863) < 0.015
assert r_mean[10, 1][0] > r_mean[10, 2][0] > r_mean[10, 4][0] > r_mean[10, 8][0]
```

Pour $L=10$, le rapport moyen passe de $0.530$ ($W=0.5$) à $0.384$ ($W=16$), en passant par $0.484$ ($W=2$) et
$0.410$ ($W=4$), avec des erreurs standard de $0.003$–$0.004$. Pour $L=8$, la valeur en $W=0.5$ ($0.507$) est
inférieure à celle en $W=1$ ($0.525$) : en $W=0$ la chaîne est intégrable et possède toutes les symétries
ci-dessus, et une valeur GOE exige un désordre grand devant l'espacement des niveaux. Le crossover est large et
dérive avec $L$ ; une analyse de taille finie le situe près de $W\approx3.7$ (D. Luitz, N. Laflorencie, F. Alet,
*Phys. Rev. B* 91 (2015) 081103), et la survie d'une phase localisée quand $L\to\infty$ est débattue (J. Šuntajs,
J. Bonča, T. Prosen, L. Vidmar, *Phys. Rev. E* 102 (2020) 062144).

**Intrication des fermions libres.** Pour les états propres de la chaîne XX ouverte, l'entropie d'intrication des
$n$ premiers sites est déterminée par une matrice $n\times n$.

**Théorème 9.16 (formule de Peschel).** Soient $\Phi=b_{k_1}^*\cdots b_{k_N}^*|0\rangle$ comme dans le théorème 9.10,
$C_{ij}=\langle\Phi,c_i^*c_j\Phi\rangle$ et $C_A=(C_{ij})_{i,j<n}$. Alors $0\le C_A\le1$, et l'entropie d'intrication des
sites $A=\{0,\dots,n-1\}$ est, $\nu$ parcourant les valeurs propres de $C_A$ comptées avec multiplicité,
$$
S_A=\sum_\nu \eta(\nu),\qquad \eta(\nu)=-\nu\log\nu-(1-\nu)\log(1-\nu)\quad(0\log0=0).
$$

*Démonstration.* Pour un vecteur réel $x$, posons $b(x)=\sum_jx_jc_j$ ; alors $b(x)^*b(y)^*=-b(y)^*b(x)^*$ et
$\{b(x),b(y)^*\}=\langle x,y\rangle$. Soient $U\subseteq\mathbb R^L$ l'espace engendré par les $u_{k_a}$, et $P_A$,
$P_B=1-P_A$ les projecteurs de coordonnées sur $\{0,\dots,n-1\}$ et sur son complémentaire. Diagonaliser la forme
$\langle P_Ax,P_Ay\rangle$ sur $U$ donne une base orthonormée $\phi_1,\dots,\phi_N$ de $U$ avec
$\langle P_A\phi_a,P_A\phi_b\rangle=\nu_a\delta_{ab}$, d'où $\langle P_B\phi_a,P_B\phi_b\rangle=(1-\nu_a)\delta_{ab}$. Comme
$(x_1,\dots,x_N)\mapsto b(x_1)^*\cdots b(x_N)^*|0\rangle$ est multilinéaire et alternée,
$\Phi=\pm\prod_a\big(\sqrt{\nu_a}f_a^*+\sqrt{1-\nu_a}\,g_a^*\big)|0\rangle$ avec $f_a=b(P_A\phi_a/\sqrt{\nu_a})$,
$g_a=b(P_B\phi_a/\sqrt{1-\nu_a})$ (en omettant les termes de coefficient nul). En développant,
$\Phi=\sum_I\pm\prod_{a\in I}\sqrt{\nu_a}\prod_{a\notin I}\sqrt{1-\nu_a}\,F_I^*G_{I^c}^*|0\rangle$, la somme portant sur les
parties $I$, avec $F_I^*$, $G_{I^c}^*$ les produits ordonnés. Pour $j\ge n$, la chaîne de signes de $c_j^*$ contient
tout $A$ : $c_j^*=(-1)^{N_A}\otimes\tilde c_j^*$ avec des opérateurs de Jordan–Wigner $\tilde c_j$ des sites $\ge n$,
tandis que pour $i<n$, $c_i^*=\hat c_i^*\otimes1$. Comme $(-1)^{N_A}$ vaut $1$ sur $|0_A\rangle$,
$F_I^*G_{I^c}^*|0\rangle=\xi_I\otimes\chi_{I^c}$ avec $\xi_I=\hat F_I^*|0_A\rangle$, $\chi_{I^c}=\tilde G_{I^c}^*|0_B\rangle$, qui
sont des familles orthonormées d'après le calcul avec les CAR du théorème 9.10. C'est une décomposition de Schmidt
(lemme 7.7) de poids $p_I=\prod_{a\in I}\nu_a\prod_{a\notin I}(1-\nu_a)$, la loi de variables de Bernoulli
indépendantes, donc $S_A=-\sum_Ip_I\log p_I=\sum_a\eta(\nu_a)$. Enfin $C_{ij}=\sum_a\phi_a(i)\phi_a(j)$ (développer
$c_i^*c_j$ sur les $b(\phi_a)$), donc $C_A=\sum_a(P_A\phi_a)(P_A\phi_a)^T$ a les valeurs propres $\nu_a$ sur les vecteurs
orthogonaux $P_A\phi_a$ et $0$ ailleurs, et $\eta(0)=0$. $\square$

`entanglement_entropy` calcule l'état réduit secteur par secteur (il commute avec l'aimantation du bloc,
proposition 7.8) et son entropie. L'état fondamental du secteur $N$ est $\Phi$ avec les $N$ modes les plus bas,
simple parce que les $\varepsilon_k$ sont distincts :

```python
L, N = 12, 6
xx = SpinChain(L)
w, V = torch.linalg.eigh(xx.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
rho = DensityMatrix(xx.sector_algebra(N), matrix=torch.outer(V[:, 0], V[:, 0].conj())[None, None])
eps_k, u = np.linalg.eigh(np.eye(L, k=1) / 2 + np.eye(L, k=-1) / 2)  # h = A_L / 2, ascending
C = u[:, :N] @ u[:, :N].T                                            # correlation matrix of the N lowest modes
bern = lambda nu: -(nu * np.log(nu) + (1 - nu) * np.log(1 - nu))
S_peschel = [bern(np.clip(np.linalg.eigvalsh(C[:n, :n]), 1e-15, 1 - 1e-15)).sum() for n in range(1, L)]
S_lib = [xx.entanglement_entropy(rho, n, sector=N).item() for n in range(1, L)]
assert np.allclose(S_lib, S_peschel, atol=1e-9)                       # Theorem 9.16
print([round(s, 4) for s in S_lib[:6]])
```

Les entropies $0.6931,0.5199,0.7456,0.6283,0.7719,0.6560$ pour $n=1,\dots,6$ oscillent avec la parité de $n$ et
croissent comme $\frac16\log L$ à $n/L$ fixé, la loi de la théorie conforme des champs pour une chaîne ouverte
(I. Peschel, *J. Phys. A* 36 (2003) L205 ; P. Calabrese, J. Cardy, *J. Stat. Mech.* (2004) P06002 ; énoncé sans
démonstration).

## 9.7 Points de comparaison

QuSpin (Python) construit des hamiltoniens de chaînes de spins dans des secteurs d'aimantation, d'impulsion, de
parité et d'inversion de spin, et utilise les solveurs creux de valeurs propres et les exponentielles de matrices
de SciPy ; des codes dédiés de diagonalisation exacte (XDiag, SPINPACK, HΦ) sont optimisés pour les grands
problèmes creux, et, avec les symétries du réseau et de grosses machines, de tels codes ont atteint environ
40–50 spins ; ITensor, TeNPy et TensorKit.jl (avec des arbres de fusion $SU(2)$) traitent les longues chaînes par
états produits de matrices. Ce qu'apporte la bibliothèque est modeste et spécifique : la chaîne en tant que
`TypeIAlgebra` $\bigoplus_NM_{\binom LN}$, de sorte que les états, canaux, entropies et dynamiques des chapitres 3–8
s'appliquent secteur par secteur ; la tour d'inclusions avec son indice ; des blocs $SU(2)$ issus des vecteurs de
plus haut poids, avec des corrélations par Wigner–Eckart ; des lots de réalisations du désordre en un seul appel ;
et des états réduits qui ne quittent jamais les secteurs. Elle n'a pas de secteurs de réflexion ni d'inversion de
spin, n'a de secteurs d'impulsion que pour les chaînes XXZ uniformes, des méthodes denses jusqu'à $L\approx16$, des
méthodes de Krylov jusqu'à $L\approx24$ (environ $26$ sur un gros GPU), et pas de réseaux de tenseurs.

## Exercices

1. (Secteurs d'impulsion à la main) Pour $L=4$ et $N=2$, trouver les orbites des six configurations sous $T$ et les
   dimensions $d_{2,m}$, $m=0,\dots,3$.
2. (Règles de somme) Soit $\psi$ l'état fondamental de l'anneau de Heisenberg avec $L$ pair. Montrer que
   $E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$ et $\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0$.
3. (Trop de symétrie) Montrer que, dans le secteur $N=L/2$, le vecteur dont tous les coefficients valent $1$ est un
   multiple de $(S^-)^{L/2}|{\uparrow\cdots\uparrow}\rangle$ et un vecteur propre de l'anneau de Heisenberg d'énergie
   $JL/4$, et que tout espace de Krylov $\mathcal K_m(H,v)$ est contenu dans chaque sous-espace $H$-invariant
   contenant $v$. Que renvoie `krylov.ground_state` lorsqu'on part de ce vecteur ?

<details>
<summary>Solutions</summary>

1. $0011\to1001\to1100\to0110\to0011$ est une orbite avec $R=4$, et $0101\to1010\to0101$ une orbite avec $R=2$.
   D'après la proposition 9.5, la première contribue à tout $m$, la seconde aux $m$ tels que $2m\equiv0\pmod4$,
   c'est-à-dire $m=0,2$ : $d_{2,m}=2,1,2,1$.

```python
assert SpinChain(4, 'periodic').momentum_dims(2) == [2, 1, 2, 1]
```

2. L'état fondamental est un singulet simple (théorème 9.9) et $T$ commute avec $H$, donc $T\psi=e^{ik}\psi$. Comme
   $\mathbf s_i\cdot\mathbf s_{i+1}=T^i(\mathbf s_0\cdot\mathbf s_1)T^{*i}$, $\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle$ est le même
   pour les $L$ liaisons, et $E_0=J\sum_i\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$. Comme
   $\mathbf S^2\psi=0$,
   $0=\langle\mathbf S^2\rangle=\sum_{i,j}\langle\mathbf s_i\cdot\mathbf s_j\rangle=L\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle$, là encore
   par invariance par translation.

```python
assert abs(12 * corr_we[1].item() - gs.energy.item()) < 1e-10 and abs(corr_we.sum().item()) < 1e-10
```

3. $S^-$ abaisse l'un des spins vers le haut ; en l'appliquant $L/2$ fois à l'état complètement polarisé, on
   atteint chaque configuration à $L/2$ spins vers le haut selon $(L/2)!$ ordres, tous avec le coefficient $1$.
   $|{\uparrow\cdots\uparrow}\rangle$ a ses $L$ liaisons parallèles, d'énergie $JL/4$, et $S^-$ commute avec $H$, donc
   le vecteur est un vecteur propre de même énergie. Si $W$ est invariant et $v\in W$, alors $H^jv\in W$ pour tout
   $j$. Ici $\mathcal K_m=\mathbb Cv$ pour tout $m$, et `ground_state` renvoie $E=JL/4$, le sommet du spectre (chaque
   $\mathbf s_i\cdot\mathbf s_j\le\frac14$) ; un vecteur de départ aléatoire évite cela avec probabilité un.

```python
c10 = SpinChain(10, 'periodic', complex_valued=False)
E_sym, _, _ = krylov.ground_state(c10.xxz_sparse(1.0, 1.0, sector=5), v0=torch.ones(1, 252, dtype=f64))
assert abs(E_sym.item() - 2.5) < 1e-10 and E_sym.item() > RING[10] + 7     # far from E_0
```

</details>

**Suite :** le [chapitre 10](10_quantum_computing.md) traite du calcul quantique avec des secteurs : sous-espaces
sans décohérence et sous-systèmes sans bruit, correction d'erreurs dans les algèbres d'opérateurs, et mélangeurs
qui respectent un poids de Hamming conservé — l'algèbre $\bigoplus_wM_{\binom nw}$ de ce chapitre dans un nouveau
rôle.
