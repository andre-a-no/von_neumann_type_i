[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · **Français** · [Deutsch](../de/11_optimisation.md) · [中文](../zh/11_optimisation.md) · [日本語](../ja/11_optimisation.md) · [Русский](../ru/11_optimisation.md)

# Tutoriel 11. Optimisation et inégalités de traces

Ce tutoriel accompagne le [chapitre 11](../../lecture_notes/fr/11_optimisation.md) des notes de cours, qui démontre les extrêmes de trois
fonctionnelles de traces et explique pourquoi l'échantillonnage d'unitaires ne les trouve pas. Ici, les paramétrages
d'ensembles de contraintes de la bibliothèque (`optimize.PositiveParam`, `UnitaryParam`, `SelfAdjointParam`) et son
optimiseur multi-départs (`optimize.extremize`) sont confrontés à ces théorèmes. Avec $\Delta$ le contraste de Michelson
et $\ell(\delta)=(1-\delta)/(1+\delta)$, les trois fonctionnelles sont :

| | fonctionnelle | variables |
|---|---|---|
| Exp. 1 | $z=\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)$ | $X,Y\ge0$, $U\in\mathcal M$ unitaire |
| Exp. 2 | $z=\mathrm{Tr}\lvert XY\rvert-\mathrm{Tr}(\lvert X\rvert\,\lvert Y\rvert)$ | $X=UX_0$, $Y=VY_0$, $X_0,Y_0\ge0$, $U,V$ unitaires |
| Exp. 3 | $z=\mathrm{Tr}(Y\lvert X\rvert Y)-\mathrm{Tr}\lvert YXY\rvert$ | $X=X^*$, $Y\ge0$ |

avec $|X|$, $|Y|$ normalisés par $\lambda_{\max}=1$ et leurs contrastes comme paramètres. Les **enveloppes** sont
$\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}$ et $\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}$, toutes les
autres variables étant libres.

<a id="s-11-1"></a>

## 11.1 Ensembles de contraintes et leurs paramétrages

*Notes de cours :* [section 11.1](../../lecture_notes/fr/11_optimisation.md#s-11-1) ; [lemme 11.1](../../lecture_notes/fr/11_optimisation.md#r-11-1), [proposition 11.2](../../lecture_notes/fr/11_optimisation.md#r-11-2).

`optimize.PositiveParam`, avec un contraste prescrit $\delta$, réalise le paramétrage de la
[section 11.1](../../lecture_notes/fr/11_optimisation.md#s-11-1) des notes de cours : pour chaque départ, il tire deux emplacements distincts $s_\pm$ (`spectrum.pin_max`,
`spectrum.pin_min` ; chaque couple ordonné a la probabilité $1/(D(D-1))$) et un unitaire $U_0$ selon la loi de Haar
(`unitary.base`, tiré par `random_unitary_operator`, [lemme 11.10](../../lecture_notes/fr/11_optimisation.md#r-11-10)), et envoie les paramètres $\theta$
(`spectrum.theta`) et $G$ (`unitary.G`) sur
$$
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell(\delta),\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
$$
`UnitaryParam` est l'application $G\mapsto U_0e^{G-G^*}$ seule, et `SelfAdjointParam` ajoute des signes,
$X=U\mathrm{diag}(\sigma\lambda)U^*$. D'après la [proposition 11.2](../../lecture_notes/fr/11_optimisation.md#r-11-2), l'image est l'ensemble des opérateurs de
la couche de contraste $\mathcal P_\delta$ ayant $1$ et $\ell(\delta)$ aux emplacements tirés. Un contraste libre dans $[a,b]$
est un paramètre de plus ; la valeur par défaut de la bibliothèque est $[a,b]=[0,0.999]$, et le code ci-dessous la conserve.

Le code vérifie la contrainte, ainsi que la proposition 11.2 de manière constructive : les paramètres de la
démonstration reproduisent une cible.

```python
import math
import scipy.linalg
import torch
from torch_vn_algebra import TypeIAlgebra, optimize as opt

torch.manual_seed(0)
alg = TypeIAlgebra([3, 3], [3, 3], precision='double')          # M_3 (+) M_3, D = 6
ell = lambda delta: (1 - delta) / (1 + delta)                     # lambda_min for contrast delta
d = torch.tensor([0.0, 0.25, 0.5, 0.9], dtype=torch.float64)
x = opt.PositiveParam(alg, 4, delta=d).operator()                 # four starts, one contrast each
assert torch.allclose(x.michelson_contrast, d, atol=1e-12)        # the constraint holds exactly
assert torch.allclose(x.lambda_max, torch.ones(4, dtype=torch.float64), atol=1e-12)

P = opt.PositiveParam(alg, 1, delta=torch.tensor([0.5], dtype=torch.float64))
lo = ell(0.5)
mu = lo + (1 - lo) * torch.rand(1, 2, 3, dtype=torch.float64)    # target spectrum in [lo, 1] ...
(cp, ip), (cm, im) = P.spectrum.pin_max[0].tolist(), P.spectrum.pin_min[0].tolist()
mu[0, cp, ip], mu[0, cm, im] = 1.0, lo                            # ... with 1 and lo in the drawn slots
W = alg.random_unitary_operator(1).matrix                         # target eigenvectors
target = W @ torch.diag_embed(mu.to(W.dtype)) @ W.mH
with torch.no_grad():                                             # the parameters of the proof
    P.spectrum.theta.copy_(torch.asin((2 * (mu - lo) / (1 - lo) - 1).clamp(-1, 1)))
    V = P.unitary.base.mH @ W                                     # U0^* W = exp(S), G = S / 2
    S = torch.stack([torch.from_numpy(scipy.linalg.logm(V[0, c].numpy())) for c in range(2)])
    P.unitary.G.copy_(S[None] / 2)
err = (P.operator().matrix - target).abs().max().item()
print(f"slots of 1 and lo: {(cp, ip)}, {(cm, im)};  reconstruction error {err:.1e}")
assert err < 1e-10
```

<a id="s-11-2"></a>

## 11.2 L'inégalité de traces de von Neumann

*Notes de cours :* [section 11.2](../../lecture_notes/fr/11_optimisation.md#s-11-2) ; [théorème 11.4](../../lecture_notes/fr/11_optimisation.md#r-11-4), [corollaire 11.5](../../lecture_notes/fr/11_optimisation.md#r-11-5), [proposition 11.6](../../lecture_notes/fr/11_optimisation.md#r-11-6), [corollaire 11.7](../../lecture_notes/fr/11_optimisation.md#r-11-7).

D'après [notes de cours, théorème 11.4](../../lecture_notes/fr/11_optimisation.md#r-11-4) (inégalité de traces de von Neumann),
$|\mathrm{Tr}(AB)|\le\sum_is_i(A)s_i(B)$, et d'après le [corollaire 11.5](../../lecture_notes/fr/11_optimisation.md#r-11-5), $\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i$ pour $P,Q$ autoadjointes de valeurs
propres $p,q$, et $\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i$ pour $P,Q\ge0$. Pour Exp. 2 à spectres fixés, la
[proposition 11.6](../../lecture_notes/fr/11_optimisation.md#r-11-6) donne $\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}(y^\downarrow_{c,i}-y^\uparrow_{c,i})$, atteint pour $X_0=\mathrm{diag}(x^\downarrow)$, $V=J$ la permutation qui renverse l'ordre et $Y_0=\mathrm{diag}(y^\uparrow)$
(resp. $\mathrm{diag}(y^\downarrow)$), et le [corollaire 11.7](../../lecture_notes/fr/11_optimisation.md#r-11-7) donne l'enveloppe $\sup z=-\inf z=\sum_c\lfloor k_c/2\rfloor(1-\ell_X)(1-\ell_Y)$.

Ici $1-\ell_X=2\delta/(1+\delta)$ pour $\Delta(|X|)=\delta$ ; pour $k_c$ impair, la valeur propre du milieu est sans
importance. Le code vérifie le théorème 11.4, le corollaire 11.5 et les configurations extrémales de la
proposition 11.6 pour $k=3$, $C=2$.

```python
M5 = TypeIAlgebra([5], [5], precision='double')
A = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))         # Ginibre, no structure
B = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))
sA, sB = torch.linalg.svdvals(A.matrix[:, 0]), torch.linalg.svdvals(B.matrix[:, 0])   # decreasing
assert ((A @ B).trace.abs() <= (sA * sB).sum(-1) + 1e-12).all()             # Theorem 11.4
Ph, Qh = (M5.operator_from_eigenvalues(lambda k: torch.randn(2000, k), batch_size=2000, force_self_adjoint=True)
          for _ in range(2))
p, q = Ph.eigenvalues()[0], Qh.eigenvalues()[0]                             # increasing order
t = (Ph @ Qh).trace.real
assert (t <= (p * q).sum(-1) + 1e-12).all() and (t >= (p * q.flip(-1)).sum(-1) - 1e-12).all()   # Cor. 11.5(a)
Pp, Qp = Ph @ Ph, Qh @ Qh                                                   # positive, spectra p^2, q^2
nuc = (Pp @ Qp).trace_norm()
assert ((Pp @ Qp).trace.real <= nuc + 1e-12).all()                          # Cor. 11.5(b)
assert (nuc <= ((p ** 2).sort(-1)[0] * (q ** 2).sort(-1)[0]).sum(-1) + 1e-10).all()

alg3 = TypeIAlgebra([3, 3], [3, 3], precision='double')
x_dn = torch.tensor([1.0, 0.6, ell(0.5)], dtype=torch.float64)              # Delta = 0.5, middle value 0.6
y_dn = torch.tensor([1.0, 0.7, ell(0.8)], dtype=torch.float64)              # Delta = 0.8, middle value 0.7
J = torch.eye(3, dtype=torch.float64).flip(0)                               # reversal permutation
blocks = lambda m: alg3.from_blocks([m, m])                                 # the same block in both sectors

def z_exp2(X0, V, Y0):
    return ((X0 @ V @ Y0).trace_norm() - (X0 @ Y0).trace.real).item()

X0 = blocks(torch.diag(x_dn))
z_max = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn.flip(0))))
z_min = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn)))
exact = 2 * (x_dn * (y_dn - y_dn.flip(0))).sum().item()                     # Proposition 11.6
envelope = 2 * (1 - ell(0.5)) * (1 - ell(0.8))                              # Corollary 11.7, floor(3/2) = 1
print(f"Exp. 2, k = 3, C = 2: z_max = {z_max:.6f}, z_min = {z_min:.6f}, formula {envelope:.6f}")
assert abs(z_max - exact) < 1e-12 and abs(z_min + exact) < 1e-12 and abs(exact - envelope) < 1e-12
```

<a id="s-11-3"></a>

## 11.3 Une inégalité de traces avec cas d'égalité

*Notes de cours :* [section 11.3](../../lecture_notes/fr/11_optimisation.md#s-11-3) ; [lemme 11.8](../../lecture_notes/fr/11_optimisation.md#r-11-8), [proposition 11.9](../../lecture_notes/fr/11_optimisation.md#r-11-9).

D'après [notes de cours, proposition 11.9](../../lecture_notes/fr/11_optimisation.md#r-11-9), $z=\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\ge0$ dans Exp. 3, avec égalité si et seulement si
$YX_+Y^2X_-Y=0$, en particulier si $XY=YX$. Pour $X=\mathrm{diag}(1,-1,0)$ et $Y$ de lignes $(1,0,a)$, $(0,1,0)$, $(a,0,1)$, $|a|<1$, l'égalité a lieu bien
que $XY\ne YX$ pour $a\ne0$. Le code vérifie les deux (avec $a=0.5$) ; `X.abs()` est $|X|$ et `trace_norm()` est
$\|\cdot\|_1$.

```python
def z_exp3(X, Y):
    return (Y @ X.abs() @ Y).trace.real - (Y @ X @ Y).trace_norm()

Xs = alg3.operator_from_eigenvalues(lambda k: torch.randn(500, k), batch_size=500, force_self_adjoint=True)
Yp = alg3.operator_from_eigenvalues(lambda k: torch.rand(500, k), batch_size=500, force_positive=True)
assert (z_exp3(Xs, Yp) >= -1e-12).all()                                     # Proposition 11.9

M3 = TypeIAlgebra([3], [3], precision='double')
Xe = M3.operator(torch.diag(torch.tensor([1.0, -1.0, 0.0]))[None, None])
Ye = M3.operator(torch.tensor([[1.0, 0, 0.5], [0, 1, 0], [0.5, 0, 1]])[None, None])
print(f"equality without commutation: z = {z_exp3(Xe, Ye).item():.1e}, "
      f"||XY - YX|| = {(Xe @ Ye - Ye @ Xe).operator_norm().item():.2f}")
assert abs(z_exp3(Xe, Ye).item()) < 1e-12 and (Xe @ Ye - Ye @ Xe).operator_norm().item() > 0.1
```

<a id="s-11-4"></a>

## 11.4 Échantillonnage contre optimisation

*Notes de cours :* [section 11.4](../../lecture_notes/fr/11_optimisation.md#s-11-4) ; [lemme 11.10](../../lecture_notes/fr/11_optimisation.md#r-11-10), [proposition 11.11](../../lecture_notes/fr/11_optimisation.md#r-11-11), [corollaire 11.12](../../lecture_notes/fr/11_optimisation.md#r-11-12), [lemme 11.13](../../lecture_notes/fr/11_optimisation.md#r-11-13).

Quelle est la taille typique de $|\mathrm{Tr}(AU)|$, comparée à son maximum $\|A\|_1$ ? La bibliothèque tire $U$
comme le facteur $Q$ de la décomposition QR $Z=QR$, $R_{ii}>0$, d'une **matrice de Ginibre** $Z$ (coefficients
indépendants de densité $\pi^{-1}e^{-|z|^2}$ ; recette de Mezzadri), indépendamment dans chaque bloc.

D'après [notes de cours, proposition 11.11](../../lecture_notes/fr/11_optimisation.md#r-11-11), $\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k$ pour $U\in U(k)$ de loi de Haar, et d'après le
[corollaire 11.12](../../lecture_notes/fr/11_optimisation.md#r-11-12), pour $X,Y\ge0$ à spectres dans $[a,b]$, $a>0$, la moyenne quadratique de $|\mathrm{Tr}(XUY)|$ vaut
au plus $\min(1/\sqrt k,b^2/(a^2k))$ fois son maximum sur $U$, et $\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le b^4/(a^4k^2)$.

Ainsi, $N$ échantillons de Haar trouvent une violation de $|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)$ avec une
probabilité d'au plus $Nb^4/(a^4k^2)$ (une réunion de $N$ événements), bien qu'il en existe une pour tout couple qui
ne commute pas (proposition 3.9). La borne n'est informative que pour $k>\sqrt N\,b^2/a^2$ (pour $b/a=6$ et $N=2000$ :
$k>1609$) ; numériquement, les échantillons manquent la violation dès des $k$ bien plus petits (ci-dessous). Le code
vérifie la proposition 11.11 par Monte-Carlo (la moyenne empirique est à moins de cinq erreurs standard de la valeur
exacte, pour ce calcul à graine fixée) et affiche, pour des spectres uniformes dans $[0.2,1.2]$ ($b/a=6$), le rapport du
corollaire 11.12, $\sup z$ et le meilleur $z$ parmi $2000$ échantillons, tous deux rapportés à $\mathrm{Tr}(XY)$.

```python
def random_pair(alg_k):
    pos = lambda k: 0.2 + torch.rand(1, k)                                  # spectrum in [a, b] = [0.2, 1.2]
    return tuple(alg_k.operator_from_eigenvalues(pos, force_positive=True) for _ in range(2))

M6 = TypeIAlgebra([6], [6], precision='double')
Xf, Yf = random_pair(M6)
N = 20_000
t2 = (Yf @ Xf @ M6.random_unitary_operator(N)).trace.abs() ** 2
m2, err, exact2 = t2.mean().item(), t2.std().item() / math.sqrt(N), (Yf @ Xf).frobenius_norm().item() ** 2 / 6
print(f"E|Tr(AU)|^2 = {m2:.4f} +- {err:.4f}   (exact ||A||_2^2 / k = {exact2:.4f})")
assert abs(m2 - exact2) < 5 * err                                           # Proposition 11.11

print(" k   rms/max   sup z/TrXY   best of 2000 samples")
for k in (2, 4, 8, 16, 32):
    Mk = TypeIAlgebra([k], [k], precision='double')
    Xk, Yk = random_pair(Mk)
    trXY, nuc = (Xk @ Yk).trace.real.item(), (Yk @ Xk).trace_norm().item()
    ratio = (Yk @ Xk).frobenius_norm().item() / (math.sqrt(k) * nuc)
    best = (Xk @ Mk.random_unitary_operator(2000) @ Yk).trace.abs().max().item()
    print(f"{k:2d}   {ratio:.4f}    {(nuc - trXY) / trXY:.2e}     {(best - trXY) / trXY:+.3f}")
    assert ratio <= min(1 / math.sqrt(k), 36 / k) + 1e-12                   # Corollary 11.12
    assert nuc > trXY and best < nuc and (k < 4 or best < trXY)            # from k = 4 on: no violation found
```

*Observation numérique* (le tableau affiché ; affirmée seulement pour ce calcul à graine fixée) : le rapport reste
sous la borne démontrée $\min(1/\sqrt k,b^2/(a^2k))$ et est à peu près divisé par deux chaque fois que $k$ double,
c'est-à-dire qu'il décroît comme $1/k$, et à partir de $k=4$ le meilleur des $2000$ échantillons n'atteint même pas
$z=0$. Une recherche par gradient, en revanche, s'approche des deux extrêmes dans le calcul ci-dessous (là encore une
observation numérique). Pour le couple `Xf, Yf`, les deux extrêmes en $U$ sont connus :

D'après [notes de cours, lemme 11.13](../../lecture_notes/fr/11_optimisation.md#r-11-13), $\min_U|\mathrm{Tr}(AU)|=0$ pour toute $A\in M_k(\mathbb C)$, $k\ge2$.

Ainsi $\inf_Uz=-\mathrm{Tr}(XY)$ et $\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)$ (proposition 3.9). Un seul appel à
`extremize` approche les deux, avec $12$ départs qui maximisent $z$ et $12$ qui maximisent $-z$ (à une précision
relative de $10^{-3}$ dans ce calcul, affirmée ci-dessous). D'après la proposition 11.11, la moyenne quadratique de
$|\mathrm{Tr}(XUY)|$ vaut $\|YX\|_2/\sqrt k$ ; avec les bornes de la démonstration du corollaire 11.12 (spectres
dans $[a,b]$), $\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4$ et $\mathrm{Tr}(XY)\ge ka^2$, donc l'inégalité de Markov
donne $\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)$ pour $\varepsilon>0$ :
pour $k$ grand, $z$ est proche de son infimum avec grande probabilité.
*Observation numérique* (affirmée pour ce calcul à graine fixée) : l'échantillonnage arrive à moins de $1\,\%$ de
l'infimum, mais n'atteint même pas $z=0$, alors que le supremum est positif :

```python
trXY = (Xf @ Yf).trace.real
sup_exact, inf_exact = ((Yf @ Xf).trace_norm() - trXY).item(), -trXY.item()
U6 = opt.UnitaryParam(M6, 24)                                               # 24 starts
pm = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(12)   # 12 maximise z, 12 maximise -z
best = opt.extremize(lambda: pm * ((Xf @ U6.operator() @ Yf).trace.abs() - trXY), [U6],
                     maximize=True, rounds=6, steps=30).reshape(2, 12).max(1)[0]
best_up, best_down = best[0].item(), -best[1].item()
z_mc = (Xf @ M6.random_unitary_operator(N) @ Yf).trace.abs() - trXY
print(f"sup z: exact {sup_exact:.6e}, ascent {best_up:.6e}, best sample {z_mc.max().item():+.4f}")
print(f"inf z: exact {inf_exact:.6f}, descent {best_down:.6f}, worst sample {z_mc.min().item():+.6f}")
assert best_up <= sup_exact + 1e-12 and best_down >= inf_exact - 1e-12    # attained values: one-sided
assert sup_exact - best_up < 1e-3 * sup_exact and best_down - inf_exact < 1e-3 * abs(inf_exact)
assert z_mc.min().item() - inf_exact < 1e-2 * abs(inf_exact)               # sampling finds the typical side ...
assert z_mc.max().item() < 0                                               # ... but misses the violation
```

<a id="s-11-5"></a>

## 11.5 Recherche multi-départs des enveloppes

*Notes de cours :* [section 11.5](../../lecture_notes/fr/11_optimisation.md#s-11-5) ; [proposition 11.14](../../lecture_notes/fr/11_optimisation.md#r-11-14), [corollaire 11.7](../../lecture_notes/fr/11_optimisation.md#r-11-7), [proposition 11.9](../../lecture_notes/fr/11_optimisation.md#r-11-9).

`optimize.extremize(objective, params, maximize)` exécute Adam en `rounds` tours de `steps` pas (taux
d'apprentissage `lr * decay**round`) sur la perte $\mp\sum_b\text{objective}_b$, recentre les cartes unitaires après
chaque tour et renvoie, **pour chaque départ**, la plus grande (la plus petite, pour `maximize=False`) valeur de
l'objectif parmi celles calculées aux itérés du gradient et une fois au point final. Chaque élément du lot est un
départ indépendant avec ses propres emplacements, signes, $U_0$ et contraste : Adam agit coordonnée par coordonnée,
donc si la $b$-ième composante de l'objectif ne dépend que des $b$-ièmes paramètres (comme dans tous les objectifs
ci-dessous), les départs n'interagissent pas. Les signes de `SelfAdjointParam` sont discrets : entre les tours,
`discrete_step` essaie huit changements de signe individuels par départ (chacun en un emplacement uniformément
aléatoire), l'un après l'autre, et conserve chaque changement qui améliore l'objectif.

En arithmétique exacte, chaque itéré est admissible (la proposition 11.2 et les remarques qui la suivent), donc chaque
valeur renvoyée est la valeur de $z$ en un point admissible ; en virgule flottante, c'est cette valeur aux erreurs
d'arrondi près (les assertions ci-dessous tolèrent $10^{-12}$). Par conséquent, une maximisation donne une **borne
inférieure** du supremum, une minimisation une borne supérieure de l'infimum, et rien de plus — une recherche locale
peut s'arrêter en un extremum local. L'autre sens exige une démonstration, comme le corollaire 11.7 et la
proposition 11.9, ou une relaxation convexe (section 11.6). Une enveloppe calculée est confrontée (i) à des formes
closes ; (ii) à des configurations construites à la main, telles que des sommes directes de solutions plus petites,
en dessous desquelles elle ne doit pas tomber ; (iii) à des échantillons structurés (spectres à deux niveaux, bases
de permutation), dont les valeurs sont d'autres valeurs atteintes.

D'après [notes de cours, proposition 11.14](../../lecture_notes/fr/11_optimisation.md#r-11-14), si tous les $k_c\ge2$, l'enveloppe inférieure d'Exp. 1 est $\underline z(\delta)=-\big(D-2\delta/(1+\delta)\big)$.

Le code vérifie la proposition 11.14 sur $M_2\oplus M_2$, le corollaire 11.7 sur $M_3\oplus M_3$ (les deux contrastes
prescrits, $k=3$ impair) et la proposition 11.9 ; suprema et infima partagent un même appel, comme ci-dessus. Les
assertions unilatérales (aucune valeur renvoyée au-delà d'un extrême démontré) sont des conséquences de ces
théorèmes ; les assertions selon lesquelles la recherche arrive à la tolérance indiquée de l'extrême sont des
observations numériques de ce calcul à graine fixée. Le mot-clé `is_self_adjoint` choisit l'algorithme de la norme
trace (valeurs singulières, ou valeurs propres pour $YXY$) sans test à chaque pas. La dernière vérification est de
type (ii) : le $z$ d'Exp. 3 est additif sur les sommes directes, et le contraste d'une somme directe est le plus grand
contraste de ses termes (exercice 3), donc, pour les enveloppes d'Exp. 3 avec $\Delta(Y)$ libre dans $[0,0.999]$
comme dans le code,
$$
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
$$
en particulier, la somme de deux valeurs atteintes sur $M_2$ aux contrastes $\delta$ et $\delta'\le\delta$ est une
valeur atteinte sur $M_2\oplus M_2$ au contraste $\delta$.

```python
torch.manual_seed(1)
S = 6                                                                       # starts per grid point
grid = torch.tensor([0.0, 0.5, 0.9], dtype=torch.float64)
alg22 = TypeIAlgebra([2, 2], [2, 2], precision='double')                    # D = 4

# Exp. 1, lower envelope along Delta(X); Delta(Y) free in [0, 0.999]
X1, Y1, U1 = (opt.PositiveParam(alg22, 3 * S, delta=grid.repeat_interleave(S)),
              opt.PositiveParam(alg22, 3 * S), opt.UnitaryParam(alg22, 3 * S))
def z1():
    x, y = X1.operator(), Y1.operator()
    return (x @ U1.operator() @ y).trace.abs() - (x @ y).trace.real
inf1 = opt.extremize(z1, [X1, Y1, U1], maximize=False, rounds=4, steps=30).reshape(3, S).min(1)[0]
exact1 = -(4 - 2 * grid / (1 + grid))                                       # Proposition 11.14
print("Exp. 1  inf:", [f"{v:.4f}" for v in inf1.tolist()], " exact:", [f"{v:.4f}" for v in exact1.tolist()])
assert (inf1 >= exact1 - 1e-12).all() and (inf1 - exact1).max() < 4e-3

# Exp. 2 on M_3 (+) M_3, Delta(X0) = 0.5, Delta(Y0) = 0.8: sup and inf in one call
sign = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(S)
X2, Y2, V2 = (opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.5, dtype=torch.float64)),
              opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.8, dtype=torch.float64)),
              opt.UnitaryParam(alg3, 2 * S))
def z2():
    x0, y0 = X2.operator(), Y2.operator()
    M = alg3.operator((x0 @ V2.operator() @ y0).matrix, is_self_adjoint=False)
    return sign * (M.trace_norm() - (x0 @ y0).trace.real)
best2 = opt.extremize(z2, [X2, Y2, V2], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]
sup2, inf2 = best2[0].item(), -best2[1].item()
print(f"Exp. 2  sup {sup2:.6f}, inf {inf2:.6f}, exact +-{envelope:.6f}")
assert sup2 <= envelope + 1e-12 and inf2 >= -envelope - 1e-12               # Corollary 11.7: one-sided ...
assert envelope - sup2 < 1e-3 * envelope and inf2 + envelope < 1e-3 * envelope   # ... and reached

# Exp. 3: maximise sign * z at the contrasts `deltas` of |X|; Delta(Y) free in [0, 0.999]
def exp3(alg_, deltas, sign):
    X3, Y3 = opt.SelfAdjointParam(alg_, len(deltas), delta=deltas), opt.PositiveParam(alg_, len(deltas))
    def z3():
        y = Y3.operator()
        M = alg_.operator((y @ X3.operator() @ y).matrix, is_self_adjoint=True)
        return sign * ((y @ X3.abs_operator() @ y).trace.real - M.trace_norm())
    return opt.extremize(z3, [X3, Y3], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]

best3 = exp3(alg22, torch.full((2 * S,), 0.5, dtype=torch.float64), sign)
sup_M22, inf3 = best3[0].item(), -best3[1].item()
s_M2 = exp3(TypeIAlgebra([2], [2], precision='double'),
            torch.tensor([0.0, 0.5], dtype=torch.float64).repeat_interleave(S), 1.0)
direct_sum = (s_M2[1] + s_M2.max()).item()                                 # delta = 0.5, best delta' <= 0.5
conj = [(1 - dl) * (2 * 0.999 / 1.999) ** 2 for dl in (0.0, 0.5)]          # unproven closed form, d = 0.999
print(f"Exp. 3  inf {inf3:.1e};  sup for M_2 at delta = 0, 0.5: {s_M2[0]:.5f}, {s_M2[1]:.5f}"
      f"  (closed form {conj[0]:.5f}, {conj[1]:.5f})")
print(f"        sup for M_2 (+) M_2 at 0.5: search {sup_M22:.5f}, direct sum {direct_sum:.5f}")
assert inf3 >= -1e-12 and inf3 < 1e-6                                       # Proposition 11.9
assert sup_M22 >= direct_sum - 5e-3 * direct_sum                           # (ii), up to 0.5 %
```

*Observations numériques* (les valeurs affichées par ce calcul). La recherche sur $M_2\oplus M_2$ reste en dessous de
la somme directe d'environ $0.2\,\%$ : elle n'a pas entièrement convergé, et la meilleure borne inférieure de
$\overline z_{M_2\oplus M_2}(0.5)$ est le plus grand des deux nombres — détecter cela est le but de la
vérification (ii). Les valeurs pour $M_2$ coïncident à $10^{-4}$ près avec
$(1-\delta)\big(2d/(1+d)\big)^2=(1-\delta)(1-\ell(d))^2$, où $d=0.999$ est le plus grand contraste admis pour $Y$.
Cette forme close n'est **pas** un théorème de ce cours : l'article qui accompagne la bibliothèque l'a observée
numériquement (dans des algèbres réelles) sans la démontrer, et le code l'affiche sans rien affirmer à son sujet.
Pour des blocs plus grands, nous ne connaissons pas de forme close des enveloppes supérieures d'Exp. 1 et d'Exp. 3 ;
les courbes calculées sont alors des bornes inférieures (des valeurs atteintes), rien de plus.

<a id="s-11-6"></a>

## 11.6 Points de comparaison

**Manopt**, **pymanopt** et **geoopt** fournissent des variétés matricielles (variétés de Stiefel, groupes orthogonaux
et unitaires) avec des rétractions telles que QR, polaire ou de Cayley, et des solveurs riemanniens (gradients
conjugués, régions de confiance, Adam riemannien), une alternative à Adam dans une carte recentrée (nous ne faisons
ici aucune comparaison de performances). **cvxpy** transmet des problèmes convexes en variables opérateurs à des
solveurs SDP, qui renvoient un optimum accompagné d'un certificat dual, tous deux à la tolérance du solveur près ;
une relaxation — un problème convexe sur un ensemble admissible plus grand — majore le supremum, ce qui est le sens
inaccessible à une recherche locale. **NCTSSOS**, **ncpol2sdpa** et **NCSOStools** bornent les valeurs propres
et les traces de polynômes non commutatifs en toutes dimensions à la fois (section 0.3). La bibliothèque apporte des
paramétrages exacts des contraintes spectrales et de contraste de $\bigoplus_cM_{k_c}$ (les extrêmes peuvent se
trouver dans n'importe quel secteur), des unitaires qui respectent les blocs, et un seul calcul par lots pour de
nombreux départs et valeurs de contrainte, sur CPU ou GPU. Elle produit des candidats extrémaux et des bornes
inférieures pour les suprema ; démonstrations ou relaxations doivent fournir le reste.

<a id="x-11"></a>

## Exercices

[Exercices 1 et 2](../../lecture_notes/fr/11_optimisation.md#x-11) des notes de cours, numériquement : la transformation de Cayley $C(S)=(1-S)(1+S)^{-1}$ d'un
$S$ antiautoadjoint est unitaire, n'a pas de valeur propre $-1$, et $S=(1-C)(1+C)^{-1}$ (exercice 1) ; pour $O$ distribué selon la
loi de Haar sur $O(k)$, tiré par `random_unitary_operator` d'une algèbre réelle (`complex_valued=False`),
$\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k$ (exercice 2).

<details>
<summary>Solutions</summary>

```python
Sk = torch.randn(4, 4, dtype=torch.complex128)
Sk = Sk - Sk.mH                                                   # skew-adjoint
I4 = torch.eye(4, dtype=torch.complex128)
Cay = (I4 - Sk) @ torch.linalg.inv(I4 + Sk)
assert torch.allclose(Cay.mH @ Cay, I4, atol=1e-12)
assert torch.allclose((I4 - Cay) @ torch.linalg.inv(I4 + Cay), Sk, atol=1e-10)
assert torch.linalg.eigvals(Cay).add(1).abs().min() > 1e-6        # -1 is never an eigenvalue
```

```python
Mr = TypeIAlgebra([5], [5], complex_valued=False, precision='double')
Ar = Mr.operator(torch.randn(1, 1, 5, 5, dtype=torch.float64))
tr2 = (Ar @ Mr.random_unitary_operator(20_000)).trace ** 2
assert abs(tr2.mean().item() - Ar.frobenius_norm().item() ** 2 / 5) < 5 * tr2.std().item() / math.sqrt(20_000)
```

</details>

**Suite :** [Tutoriel 12](12_equivariant_learning.md) : les couches équivariantes comme éléments d'un commutant, dans la bibliothèque.
