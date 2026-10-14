[English](../05_random_matrices.md) · [Español](../es/05_random_matrices.md) · **Français** · [Deutsch](../de/05_random_matrices.md) · [中文](../zh/05_random_matrices.md) · [日本語](../ja/05_random_matrices.md) · [Русский](../ru/05_random_matrices.md)

# Tutoriel 5. Matrices aléatoires dans une algèbre

Ce tutoriel accompagne le [chapitre 5](../../lecture_notes/fr/05_random_matrices.md) des notes de cours, qui construit la mesure de Haar
sur le groupe unitaire de $\mathcal M=\bigoplus_cM_{k_c}$ ([corollaire 5.4](../../lecture_notes/fr/05_random_matrices.md#r-5-4)), calcule ses moments d'ordre un et
deux ([proposition 5.7](../../lecture_notes/fr/05_random_matrices.md#r-5-7), [théorème 5.10](../../lecture_notes/fr/05_random_matrices.md#r-5-10)) et montre pourquoi les statistiques de niveaux doivent être
calculées secteur par secteur ([théorème 5.15](../../lecture_notes/fr/05_random_matrices.md#r-5-15)). Nous tirons des unitaires de Haar et des ensembles gaussiens avec
la bibliothèque et vérifions ces résultats par Monte-Carlo : une moyenne doit coïncider avec la valeur exacte à $5$
erreurs types près. Les résultats qui ne reposent sur aucun théorème sont signalés comme observations numériques.

Dans tout ce chapitre, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ avec $D=\sum_ck_c$, les
$P_c$ sont les projecteurs centraux minimaux, et $E$ est l'espérance conditionnelle sur le centre
de la proposition 2.8. $U(k)$ désigne le groupe des matrices unitaires $k\times k$ et
$U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}$ le groupe unitaire de $\mathcal M$. Les variables
aléatoires sont définies sur un espace probabilisé d'espérance $\mathbb E$ ; l'espérance d'une
matrice aléatoire se prend coefficient par coefficient. Une variable *gaussienne complexe
standard* a des parties réelle et imaginaire indépendantes de loi $\sim N(0,\tfrac12)$, de sorte
que $\mathbb E|z|^2=1$ (c'est ce que tire `torch.randn` pour un dtype complexe, d'après sa
documentation ; le premier bloc de code vérifie les moments d'ordre deux) ; une *matrice de
Ginibre* a des coefficients indépendants gaussiens complexes standard.

<a id="s-5-1"></a>

## 5.1 Mesure de Haar

*Notes de cours :* [section 5.1](../../lecture_notes/fr/05_random_matrices.md#s-5-1) ; [théorème 5.3](../../lecture_notes/fr/05_random_matrices.md#r-5-3) (le facteur unitaire de la décomposition QR d'une matrice de Ginibre, avec les phases fixées par $\Lambda=\mathrm{diag}(R'_{ii}/|R'_{ii}|)$, est distribué selon la loi de Haar), [corollaire 5.4](../../lecture_notes/fr/05_random_matrices.md#r-5-4) (les blocs d'un unitaire de Haar de $U(\mathcal M)$ sont indépendants et distribués selon la loi de Haar) et [corollaire 5.5](../../lecture_notes/fr/05_random_matrices.md#r-5-5) ($\mathbb EU_{ij}=0$, $\mathbb E|U_{ij}|^2=1/k$). La mesure de Haar de $SU(k)$ est la loi de $e^{-i\theta/k}U$, $\theta=\arg\det U$.

Sans la correction $\Lambda$ du [théorème 5.3](../../lecture_notes/fr/05_random_matrices.md#r-5-3), la loi de $Q'$ dépend de la convention
de signe de la routine QR. Pour la routine qu'appelle `torch.linalg.qr` (LAPACK), le code
ci-dessous trouve $\mathbb E\mathrm{Re}\,\mathrm{Tr}\,Q'$ éloignée de $0$ de nombreuses
erreurs types, alors qu'un unitaire de Haar vérifie $\mathbb EU=0$ (corollaire 5.5) : ce $Q'$
n'est donc pas distribué selon la loi de Haar. C'est un énoncé sur cette routine, vérifié par le
code, et non un théorème.

La bibliothèque tire des unitaires de Haar de $U(\mathcal M)$ avec `random_unitary_operator`, au
moyen d'une seule décomposition QR par lots par taille de bloc, pour tous les secteurs et tous les
échantillons à la fois (sur GPU, par Gram–Schmidt, qui donne le même $Q$ grâce à l'unicité du
lemme 5.2). La mesure `'haar'` correspond au théorème 5.3 dans chaque bloc ; pour une algèbre
complexe, `'haar_su'` multiplie par $e^{-i\theta/k}$ comme ci-dessus ; `'coe'` est l'ensemble
$S=W^TW$ de l'exercice 3 (les autres mesures de la bibliothèque ne sont pas utilisées dans ce
cours). Nous vérifions les corollaires 5.4 et 5.5, et que l'échantillonneur de la bibliothèque est
exactement la QR corrigée du théorème 5.3 :

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinChain, channels

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def mc(samples, exact, name):
    """Monte Carlo mean of real samples with its standard error; must agree within 5 standard errors."""
    mean, err = samples.mean().item(), samples.std().item() / math.sqrt(samples.numel())
    print(f"{name}: {mean:.4f} +- {err:.4f}   (theory {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12
    return mean

alg = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')        # M_1 (+) M_2 (+) M_3, D = 6
B = 100_000
U = alg.random_unitary_operator(B)                                   # Haar on U(1) x U(2) x U(3)
assert close(U.adjoint() @ U, alg.identity())                        # unitary elements of the algebra
u = [U.matrix[:, c, :k, :k] for c, k in enumerate(alg.k_factors)]   # the blocks U_c
for c, k in enumerate(alg.k_factors):                                # Corollary 5.5
    mc(u[c][:, 0, 0].real, 0.0, f"E Re (U_{c})_11")
    mc(u[c][:, 0, 0].abs() ** 2, 1 / k, f"E |(U_{c})_11|^2")
mc(u[1][:, 0, 0].abs() ** 2 * u[2][:, 0, 0].abs() ** 2, 1 / 6, "E |(U_1)_11|^2 |(U_2)_11|^2")   # independent blocks

z = torch.randn(B, dtype=dt)                                         # torch.randn draws standard complex Gaussians:
mc(z.real ** 2, 0.5, "E (Re z)^2")                                   # Re z, Im z ~ N(0, 1/2), uncorrelated
mc(z.imag ** 2, 0.5, "E (Im z)^2")
mc(z.real * z.imag, 0.0, "E Re z Im z")

SU = alg.random_unitary_operator(4, measure='haar_su').matrix
assert torch.allclose(torch.linalg.det(SU[:, 2, :3, :3]), torch.ones(4, dtype=dt))

torch.manual_seed(1)
Z = torch.randn(B, 2, 2, dtype=dt)                                   # Ginibre matrices
torch.manual_seed(1)
V = alg.random_unitary(2, batch_size=B)                              # the library, same random numbers
Q, R = torch.linalg.qr(Z)                                            # LAPACK: diag(R) real, either sign
d = R.diagonal(dim1=-2, dim2=-1)
assert close(V, Q * (d / d.abs()).unsqueeze(-2))                     # Theorem 5.3: V = Q Lambda
raw = Q.diagonal(dim1=-2, dim2=-1).sum(-1).real                      # Re Tr Q without the correction
print(f"E Re Tr Q without the phase correction: {raw.mean():.4f} +- {raw.std() / math.sqrt(B):.4f}")
assert abs(raw.mean()) > 20 * raw.std() / math.sqrt(B)               # far from 0: not Haar
```

*Observation numérique.* Sans la correction, le code trouve
$\mathbb E\mathrm{Re}\,\mathrm{Tr}\,Q'\approx-0.84$ pour $k=2$, alors qu'un unitaire de Haar
vérifie $\mathbb EU=0$.

<a id="s-5-2"></a>

## 5.2 Moments et calcul de Weingarten

*Notes de cours :* [section 5.2](../../lecture_notes/fr/05_random_matrices.md#s-5-2) ; [proposition 5.7](../../lecture_notes/fr/05_random_matrices.md#r-5-7) ($\mathbb E[UXU^*]=E(X)$), [lemme 5.8](../../lecture_notes/fr/05_random_matrices.md#r-5-8) (le commutant de tous les $V\otimes V$ est $\mathrm{span}\{1,F\}$, où $F$ est l'échange), [théorème 5.10](../../lecture_notes/fr/05_random_matrices.md#r-5-10) (la formule de Weingarten d'ordre deux) et [corollaire 5.11](../../lecture_notes/fr/05_random_matrices.md#r-5-11).

Nous vérifions la proposition 5.7 par rapport à `channels.center_expectation`, le lemme 5.8 par
algèbre linéaire (le commutant de trois $V\otimes V$ tirés au hasard, calculé comme un noyau,
contient $F$ et est de dimension $2$, donc, d'après le lemme 5.8, c'est déjà
$\mathrm{span}\{1,F\}$), et le théorème 5.10 et le corollaire 5.11 par Monte-Carlo avec le
bloc $U(3)$ des échantillons précédents :

```python
X = alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])      # a fixed element
T = (U @ X @ U.adjoint()).matrix                                    # B samples of U X U*
mean, err = T.mean(0), T.std(0) / math.sqrt(B)
EX = channels.center_expectation(alg)(X).matrix[0]
random = err > 1e-12                                                 # entries that are not constant
z = ((mean - EX).abs()[random] / err[random]).max().item()
print(f"twirl over U(M) vs E(X): largest deviation {z:.2f} standard errors")
assert ((mean - EX).abs() <= 5 * err + 1e-12).all()                 # Proposition 5.7

def commutant_dim(ops, tol=1e-9):
    """Dimension of {X : XA = AX for all A in ops}, as the null space of X -> (AX - XA)_A."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    ops = [A.contiguous() for A in ops]                              # torch.kron needs contiguous input
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return int((torch.linalg.svdvals(L) < tol).sum())

k = 3
F = torch.zeros(k * k, k * k, dtype=dt)
for i in range(k):
    for j in range(k):
        F[j * k + i, i * k + j] = 1                                  # F (e_i (x) e_j) = e_j (x) e_i
VV = [torch.kron(W, W) for W in alg.random_unitary(k, batch_size=3).contiguous()]
assert all(close(F @ W, W @ F) for W in VV)
assert commutant_dim(VV) == 2                                        # Lemma 5.8: span{1, F}

u3 = u[2]                                                            # Haar on U(3)
mc(u3[:, 0, 0].abs() ** 4, 2 / (k * (k + 1)), "E|U_11|^4")
mc(u3[:, 0, 0].abs() ** 2 * u3[:, 0, 1].abs() ** 2, 1 / (k * (k + 1)), "E|U_11|^2 |U_12|^2")
mc((u3[:, 0, 0] * u3[:, 1, 1] * (u3[:, 0, 1] * u3[:, 1, 0]).conj()).real,
   -1 / ((k - 1) * k * (k + 1)), "E U_11 U_22 conj(U_12 U_21)")

n_s = 20_000                                                         # Theorem 5.10 for a random X
X2 = torch.randn(k * k, k * k, dtype=dt)
trX, trFX = X2.diagonal().sum(), (F @ X2).diagonal().sum()
alpha, beta = (trX - trFX / k) / (k ** 2 - 1), (trFX - trX / k) / (k ** 2 - 1)
UU = torch.einsum('bij,bkl->bikjl', u3[:n_s], u3[:n_s]).reshape(n_s, k * k, k * k)   # U (x) U
T2 = UU @ X2 @ UU.mH
mean, err = T2.mean(0), T2.std(0) / math.sqrt(n_s)
assert ((mean - alpha * torch.eye(k * k, dtype=dt) - beta * F).abs() <= 5 * err).all()
```

<a id="s-5-3"></a>

## 5.3 Statistiques de niveaux

*Notes de cours :* [section 5.3](../../lecture_notes/fr/05_random_matrices.md#s-5-3) : le GOE $H=(G+G^T)/2$, le GUE $H=(G+G^*)/\sqrt2$, et le rapport d'espacements $r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})$ des espacements $s_i=\lambda_{i+1}-\lambda_i$ ; [proposition 5.13](../../lecture_notes/fr/05_random_matrices.md#r-5-13) (la surmise de Wigner pour les matrices $2\times2$) et [proposition 5.14](../../lecture_notes/fr/05_random_matrices.md#r-5-14) (le rapport moyen vaut $2\ln2-1$ pour des niveaux indépendants, $4-2\sqrt3$ et $\frac{2\sqrt3}\pi-\frac12$ pour le GOE et le GUE $3\times3$, où $\mathbb Er=I_1/I_0$ avec $I_j=\int_0^1x^jf_\beta(x)\,dx$).

La fonction `SpinChain.level_spacing_ratio` (dans `chains`) renvoie le rapport moyen de chaque
ligne de niveaux. Nous vérifions les propositions 5.13 et 5.14 et relevons la valeur du GOE pour
$n=64$ :

```python
ratio = SpinChain.level_spacing_ratio
n2 = 200_000
alg_r = TypeIAlgebra([2], [2], complex_valued=False, precision='double')     # real symmetric 2 x 2
alg_c = TypeIAlgebra([2], [2], precision='double')
G = torch.randn(n2, 1, 2, 2, dtype=torch.float64)
Gc = torch.randn(n2, 1, 2, 2, dtype=dt)
goe2 = alg_r.operator((G + G.mT) / 2, is_self_adjoint=True)
gue2 = alg_c.operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
cdf = {1: lambda t: 1 - math.exp(-math.pi * t ** 2 / 4),
       2: lambda t: math.erf(2 * t / math.sqrt(math.pi)) - 4 * t / math.pi * math.exp(-4 * t ** 2 / math.pi)}
for beta, H, mean_s in ((1, goe2, math.sqrt(math.pi)), (2, gue2, 4 / math.sqrt(math.pi))):
    lam = H.eigenvalues()[0]                                         # (n2, 2), ascending
    S = (lam[:, 1] - lam[:, 0]) / mean_s
    mc(S, 1.0, f"beta={beta}: E S")
    for t in (0.25, 1.0, 2.0):
        mc((S <= t).double(), cdf[beta](t), f"beta={beta}: P(S <= {t})")

levels = torch.rand(4000, 50, dtype=torch.float64)                   # independent levels
mc(ratio(levels), 2 * math.log(2) - 1, "<r>, independent uniform levels")      # Proposition 5.14(a)

n3 = 100_000                                                         # Proposition 5.14(b): 3 x 3
G = torch.randn(n3, 1, 3, 3, dtype=torch.float64)
Gc = torch.randn(n3, 1, 3, 3, dtype=dt)
goe3 = TypeIAlgebra([3], [3], complex_valued=False, precision='double').operator((G + G.mT) / 2, is_self_adjoint=True)
gue3 = TypeIAlgebra([3], [3], precision='double').operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
h = 1e-4
xs = (torch.arange(int(1 / h), dtype=torch.float64) + 0.5) * h      # midpoint rule on [0, 1]
for beta, H, exact in ((1, goe3, 4 - 2 * math.sqrt(3)), (2, gue3, 2 * math.sqrt(3) / math.pi - 0.5)):
    f = (xs + xs ** 2) ** beta / (1 + xs + xs ** 2) ** (1 + 1.5 * beta)
    assert abs((xs * f).sum() / f.sum() - exact) < 1e-7              # E r = I_1 / I_0
    r3 = ratio(H.eigenvalues()[0])                                   # one ratio per matrix
    mc(r3, exact, f"beta={beta}: <r>, 3 x 3")
    mc((r3 <= 0.5).double(), (f * (xs <= 0.5)).sum().item() / f.sum().item(), f"beta={beta}: P(r <= 1/2), 3 x 3")

n = 64                                                               # numerical observation, no theorem
alg64 = TypeIAlgebra([n], [n], complex_valued=False, precision='double')
G = torch.randn(800, 1, n, n, dtype=torch.float64)
lam = alg64.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()[0]
r_goe = ratio(lam[:, n // 4: 3 * n // 4])                            # 800 samples of the mean ratio
print(f"<r>, GOE(64), middle half of the spectrum: {r_goe.mean():.4f} +- {r_goe.std() / math.sqrt(800):.4f}")

def same_mean(a, b, name):
    """Two independent Monte Carlo samples of the same law: means agree within 5 combined standard errors."""
    ea, eb = a.std().item() / math.sqrt(a.numel()), b.std().item() / math.sqrt(b.numel())
    print(f"{name}: {a.mean():.4f} vs {b.mean():.4f} +- {math.hypot(ea, eb):.4f}")
    assert abs(a.mean().item() - b.mean().item()) <= 5 * math.hypot(ea, eb)
```

*Observation numérique.* Le code trouve $\langle r\rangle\approx0.53$ pour la moitié centrale du
spectre du GOE avec $n=64$, en dessous de la valeur exacte $3\times3$, $4-2\sqrt3\approx0.536$, et
bien au-dessus de la valeur $2\ln2-1\approx0.386$ des niveaux indépendants. Ce nombre est un
résultat du calcul, non un théorème ; plus bas, il sert de référence pour des blocs GOE de même
taille $64$, dont les rapports d'espacements ont exactement la même loi.

<a id="s-5-4"></a>

## 5.4 Pourquoi il ne faut pas mélanger les secteurs

*Notes de cours :* [section 5.4](../../lecture_notes/fr/05_random_matrices.md#s-5-4) ; d'après le [théorème 5.15](../../lecture_notes/fr/05_random_matrices.md#r-5-15), les niveaux de blocs indépendants différents ne se repoussent pas : pour deux blocs GOE, le plus petit espacement $\delta$ du spectre fusionné vérifie $P(\delta<\varepsilon)$ au moins d'ordre $\varepsilon$, alors qu'à l'intérieur d'un bloc GOE cette probabilité est au plus d'ordre $\varepsilon^2$.

Avec davantage de secteurs, l'effet est plus marqué. Nous ne démontrons pas de théorème limite pour
de nombreux secteurs ; l'expérience qui suit est une observation numérique. Dans la bibliothèque,
les blocs sont séparés dès le départ : `eigenvalues()` renvoie un spectre par secteur, et les
mélanger exige un `torch.cat` explicite. Nous plaçons $C$ blocs GOE indépendants de taille $64$
dans une même algèbre et comparons le rapport moyen par secteur à celui du spectre fusionné. Chaque
secteur est une matrice GOE de taille $64$, donc ses rapports ont exactement la loi de `r_goe`
ci-dessus, et le code le vérifie ; le comportement du spectre fusionné est seulement observé :

```python
k, merged_r = 64, {}
for C in (1, 2, 4, 8, 16):
    alg_C = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double')
    G = torch.randn(800 // C, C, k, k, dtype=torch.float64)
    lam = alg_C.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()    # one spectrum per sector
    per_sector = torch.stack([ratio(l[:, k // 4: 3 * k // 4]) for l in lam])
    merged = torch.cat(lam, dim=-1).sort(dim=-1).values                       # the spectrum on C^D
    r_m = ratio(merged[:, k * C // 4: 3 * k * C // 4])
    merged_r[C] = r_m.mean().item()
    print(f"C = {C:2d}: <r> per sector {per_sector.mean():.4f}, merged {merged_r[C]:.4f}"
          f" +- {r_m.std() / math.sqrt(len(r_m)):.4f}")
    same_mean(per_sector.flatten(), r_goe, f"C = {C:2d}: sectors vs GOE(64)")   # every sector is GOE(64)
assert merged_r[1] > merged_r[2] > merged_r[4] > merged_r[16]       # numerical observation
assert abs(merged_r[16] - (2 * math.log(2) - 1)) < 0.01              # numerical observation
```

*Observation numérique.* Chaque secteur a $\langle r\rangle\approx0.53$, la valeur GOE(64) de la
section 5.3 (cet accord est le test d'un fait : les secteurs sont des matrices GOE(64)). Le spectre
fusionné tombe à $0.42$ pour deux secteurs et à environ $0.39$, à moins de $0.01$ de la valeur
$2\ln2-1$ des niveaux indépendants (proposition 5.14(a)), pour $8$ et $16$ secteurs. Que le
spectre fusionné s'approche de niveaux indépendants quand $C$ croît, c'est ce qu'observe le code ;
ce n'est pas démontré dans ce cours.

Il en va de même dans un modèle physique. La chaîne de Heisenberg en champs aléatoires
(chapitre 9) conserve l'aimantation, de sorte que son hamiltonien appartient à
$\bigoplus_{N=0}^LM_{\binom LN}$ (proposition 9.1). Ses blocs ne sont *pas* indépendants — ils
dépendent tous des mêmes champs aléatoires —, si bien que le théorème 5.15 ne s'applique pas tel
quel ; l'expérience ci-dessous montre le même effet, à titre d'observation numérique.

```python
chain = SpinChain(10, complex_valued=False)                         # sectors N = 0, ..., 10
lam = chain.random_field_heisenberg(W=0.5, batch_size=20).eigenvalues()
k5 = chain.sector_dim(5)                                            # 252
r_sector = ratio(lam[5][:, k5 // 4: 3 * k5 // 4]).mean().item()
spectrum = torch.cat(lam, dim=-1).sort(dim=-1).values               # all 1024 levels mixed
r_all = ratio(spectrum[:, 256:768]).mean().item()
print(f"Heisenberg chain, L = 10: <r> in sector N = 5: {r_sector:.4f}, all sectors mixed: {r_all:.4f}")
assert r_sector > 0.49 and r_all < 0.42                             # numerical observation
```

*Observation numérique.* Dans le secteur $N=5$, le code trouve $\langle r\rangle\approx0.53$,
proche de la valeur GOE(64) ci-dessus, et le spectre où tous les secteurs sont mélangés donne
$0.39$, proche de $2\ln2-1$.

<a id="s-5-5"></a>

## 5.5 Points de comparaison

Unitaires de Haar : `scipy.stats.unitary_group` et `ortho_group`, les unitaires aléatoires de Haar
de QuTiP, et toute QR de `numpy`/`torch` suivie de la correction des phases. Calcul de
Weingarten : les paquets symboliques IntU (Mathematica ; Z. Puchała, J. A. Miszczak), RTNI
(Mathematica et Python ; M. Fukuda, R. König, I. Nechita, J. Phys. A 52 (2019), 425303) et le
paquet Python Haarpy calculent symboliquement les moments de tout ordre ; le théorème 5.10 en est
le cas d'ordre deux, et le théorème 5.12 (démontré dans l'annexe B) la formule générale qu'ils
implémentent. Les statistiques de niveaux se calculent d'ordinaire en quelques lignes de `numpy`,
et les codes de diagonalisation exacte comme QuSpin fournissent les spectres secteur par secteur.
Ce qu'apporte la bibliothèque, c'est le groupe unitaire d'une algèbre en tant qu'objet unique : des
blocs de Haar indépendants de toutes tailles tirés par des appels vectorisés, le twirl identifié à
l'espérance conditionnelle `channels.center_expectation`, et des spectres renvoyés par secteur, de
sorte que mélanger les secteurs n'est jamais le comportement par défaut.

<a id="x-5"></a>

## Exercices

Les parties numériques des [exercices](../../lecture_notes/fr/05_random_matrices.md#x-5) des notes de cours :
1. Pour $U$ distribué selon la loi de Haar sur $U(\mathcal M)$, $\mathcal M=M_1\oplus M_2\oplus M_3$, vérifier que
   $\mathbb E|\mathrm{Tr}\,U|^2=C=3$ et $\mathbb E|\mathrm{Tr}\,U|^4=17$.
2. Calculer numériquement la matrice de Weingarten d'ordre $3$ pour $k=3$ ([théorème 5.12](../../lecture_notes/fr/05_random_matrices.md#r-5-12)) et vérifier que
   $\mathbb E|U_{11}|^6=\frac1{10}$.
3. L'**ensemble orthogonal circulaire** (`measure='coe'`) est $S=W^TW$ avec $W$ distribué selon la loi de Haar sur
   $U(k)$. Vérifier que le COE de la bibliothèque est $W^TW$, que $S$ est unitaire et symétrique, que
   $\mathbb E\mathrm{Tr}\,S^m=0$ pour $m=1,2,3$ et que la densité moyenne des phases propres est uniforme.
   Calculer $\langle r\rangle$ pour les phases propres du COE avec $k=32$ et le comparer à la valeur du GOE de la
   section 5.3 (une observation numérique).

<details>
<summary>Solutions</summary>

*Exercice 1.*

```python
tr = U.trace                                                        # Tr U = sum_c Tr U_c, shape (B,)
mc(tr.abs() ** 2, 3.0, "E|Tr U|^2")
mc(tr.abs() ** 4, 17.0, "E|Tr U|^4")
```

*Exercice 2.*

```python
import itertools

def n_cycles(p):                                                    # number of cycles of a permutation
    seen, count = set(), 0
    for i in range(len(p)):
        count += i not in seen
        while i not in seen:
            seen.add(i)
            i = p[i]
    return count

perms = list(itertools.permutations(range(3)))                      # S_3; s.index(j) = s^{-1}(j)
Gram = torch.tensor([[3.0 ** n_cycles([s.index(p[i]) for i in range(3)]) for p in perms] for s in perms],
                    dtype=torch.float64)                            # k^{#cycles(s^{-1} p)}, k = 3
Wg3 = torch.linalg.inv(Gram)                                        # the Weingarten matrix, n = k = 3
assert abs(Wg3.sum().item() - 6 / (3 * 4 * 5)) < 1e-12
mc(u3[:, 0, 0].abs() ** 6, 0.1, "E|U_11|^6")
```

*Exercice 3.*

   *Observation numérique.* Le code trouve $\langle r\rangle\approx0.53$ pour le COE avec $k=32$,
   à moins de $0.015$ de la valeur GOE(64) de la section 5.3. Cet accord est un résultat du calcul ;
   aucun théorème de ce cours ne relie les statistiques d'espacements du COE à celles du GOE (pour
   les ensembles circulaires, voir F. J. Dyson, J. Math. Phys. 3 (1962), 140–156).

```python
alg_s = TypeIAlgebra([32], [32], precision='double')
torch.manual_seed(2)
S = alg_s.random_unitary_operator(400, measure='coe').matrix[:, 0]
torch.manual_seed(2)
W = alg_s.random_unitary(32, batch_size=400)                         # same random numbers
assert close(S, W.mT @ W)                                            # the library's COE is W^T W
assert close(S, S.mT) and close(S @ S.mH, torch.eye(32, dtype=dt))
for m in (1, 2, 3):                                                  # E Tr S^m = 0
    trm = torch.linalg.matrix_power(S, m).diagonal(dim1=-2, dim2=-1).sum(-1)
    mc(trm.real, 0.0, f"E Re Tr S^{m}")
    mc(trm.imag, 0.0, f"E Im Tr S^{m}")
phases = torch.linalg.eigvals(S).angle().sort(dim=-1).values         # in (-pi, pi]
mc(((phases > 0) & (phases <= 1)).sum(-1).double(), 32 / (2 * math.pi), "E #eigenphases in (0, 1]")
s = torch.diff(torch.cat([phases, phases[:, :1] + 2 * math.pi], dim=-1), dim=-1)   # around the circle
r = torch.minimum(s, s.roll(-1, -1)) / torch.maximum(s, s.roll(-1, -1))
r_coe = r.mean(-1)
print(f"<r>, COE(32) eigenphases: {r_coe.mean():.4f}; GOE(64), Section 5.3: {r_goe.mean():.4f}")
assert abs(r_coe.mean() - r_goe.mean()) < 0.015                      # numerical observation
```

</details>

**Suite :** [Tutoriel 6](06_quantum_channels.md) : les canaux quantiques dans la bibliothèque.
