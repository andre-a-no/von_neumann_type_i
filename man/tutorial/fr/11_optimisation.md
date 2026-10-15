[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · **Français** · [Deutsch](../de/11_optimisation.md) · [中文](../zh/11_optimisation.md) · [日本語](../ja/11_optimisation.md) · [Русский](../ru/11_optimisation.md)

# Tutoriel 11. Optimisation et inégalités de traces

Ce tutoriel accompagne le [chapitre 11](../../lecture_notes/fr/11_optimisation.md) des notes de cours, qui démontre les extrêmes de trois
fonctionnelles de traces et explique pourquoi l'échantillonnage d'unitaires ne les trouve pas. Ici, les paramétrages
d'ensembles de contraintes de la bibliothèque (`optimize.PositiveParam`, `UnitaryParam`, `SelfAdjointParam`) et son
optimiseur multi-départs (`optimize.extremize`) sont confrontés à ces théorèmes. Avec $`\Delta`$ le contraste de Michelson
et $`\ell(\delta)=(1-\delta)/(1+\delta)`$, les trois fonctionnelles sont :

| | fonctionnelle | variables |
|---|---|---|
| Exp. 1 | $`z=\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)`$ | $`X,Y\ge0`$, $`U\in\mathcal M`$ unitaire |
| Exp. 2 | $`z=\mathrm{Tr}\lvert XY\rvert-\mathrm{Tr}(\lvert X\rvert\,\lvert Y\rvert)`$ | $`X=UX_0`$, $`Y=VY_0`$, $`X_0,Y_0\ge0`$, $`U,V`$ unitaires |
| Exp. 3 | $`z=\mathrm{Tr}(Y\lvert X\rvert Y)-\mathrm{Tr}\lvert YXY\rvert`$ | $`X=X^*`$, $`Y\ge0`$ |

avec $`|X|`$, $`|Y|`$ normalisés par $`\lambda_{\max}=1`$ et leurs contrastes comme paramètres. Les **enveloppes** sont
$`\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}`$ et $`\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}`$, toutes les
autres variables étant libres.

<a id="s-11-1"></a>

## 11.1 Ensembles de contraintes et leurs paramétrages

*Notes de cours :* [section 11.1](../../lecture_notes/fr/11_optimisation.md#s-11-1) ; [lemme 11.1](../../lecture_notes/fr/11_optimisation.md#r-11-1), [proposition 11.2](../../lecture_notes/fr/11_optimisation.md#r-11-2).

`optimize.PositiveParam`, avec un contraste prescrit $`\delta`$, réalise le paramétrage de la
[section 11.1](../../lecture_notes/fr/11_optimisation.md#s-11-1) des notes de cours : pour chaque départ, il tire deux emplacements distincts $`s_\pm`$ (`spectrum.pin_max`,
`spectrum.pin_min` ; chaque couple ordonné a la probabilité $`1/(D(D-1))`$) et un unitaire $`U_0`$ selon la loi de Haar
(`unitary.base`, tiré par `random_unitary_operator`, [lemme 11.10](../../lecture_notes/fr/11_optimisation.md#r-11-10)), et envoie les paramètres $`\theta`$
(`spectrum.theta`) et $`G`$ (`unitary.G`) sur
```math
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell(\delta),\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
```
`UnitaryParam` est l'application $`G\mapsto U_0e^{G-G^*}`$ seule, et `SelfAdjointParam` ajoute des signes,
$`X=U\mathrm{diag}(\sigma\lambda)U^*`$. D'après la [proposition 11.2](../../lecture_notes/fr/11_optimisation.md#r-11-2), l'image est l'ensemble des opérateurs de
la couche de contraste $`\mathcal P_\delta`$ ayant $`1`$ et $`\ell(\delta)`$ aux emplacements tirés. Un contraste libre dans $`[a,b]`$
est un paramètre de plus ; la valeur par défaut de la bibliothèque est $`[a,b]=[0,0.999]`$, et le code ci-dessous la conserve.

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
$`|\mathrm{Tr}(AB)|\le\sum_is_i(A)s_i(B)`$, et d'après le [corollaire 11.5](../../lecture_notes/fr/11_optimisation.md#r-11-5), $`\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i`$ pour $`P,Q`$ autoadjointes de valeurs
propres $`p,q`$, et $`\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i`$ pour $`P,Q\ge0`$. Pour Exp. 2 à spectres fixés, la
[proposition 11.6](../../lecture_notes/fr/11_optimisation.md#r-11-6) donne $`\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}(y^\downarrow_{c,i}-y^\uparrow_{c,i})`$, atteint pour $`X_0=\mathrm{diag}(x^\downarrow)`$, $`V=J`$ la permutation qui renverse l'ordre et $`Y_0=\mathrm{diag}(y^\uparrow)`$
(resp. $`\mathrm{diag}(y^\downarrow)`$), et le [corollaire 11.7](../../lecture_notes/fr/11_optimisation.md#r-11-7) donne l'enveloppe $`\sup z=-\inf z=\sum_c\lfloor k_c/2\rfloor(1-\ell_X)(1-\ell_Y)`$.

Ici $`1-\ell_X=2\delta/(1+\delta)`$ pour $`\Delta(|X|)=\delta`$ ; pour $`k_c`$ impair, la valeur propre du milieu est sans
importance. Le code vérifie le théorème 11.4, le corollaire 11.5 et les configurations extrémales de la
proposition 11.6 pour $`k=3`$, $`C=2`$.

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

D'après [notes de cours, proposition 11.9](../../lecture_notes/fr/11_optimisation.md#r-11-9), $`z=\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\ge0`$ dans Exp. 3, avec égalité si et seulement si
$`YX_+Y^2X_-Y=0`$, en particulier si $`XY=YX`$. Pour $`X=\mathrm{diag}(1,-1,0)`$ et $`Y`$ de lignes $`(1,0,a)`$, $`(0,1,0)`$, $`(a,0,1)`$, $`-1<a<1`$, l'égalité a lieu bien
que $`XY\ne YX`$ pour $`a\ne0`$. Le code vérifie les deux (avec $`a=0.5`$) ; `X.abs()` est $`|X|`$ et `trace_norm()` est
$`\|\cdot\|_1`$.

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

Quelle est la taille typique de $`|\mathrm{Tr}(AU)|`$, comparée à son maximum $`\|A\|_1`$ ? La bibliothèque tire $`U`$
comme le facteur $`Q`$ de la décomposition QR $`Z=QR`$, $`R_{ii}>0`$, d'une **matrice de Ginibre** $`Z`$ (coefficients
indépendants de densité $`\pi^{-1}e^{-|z|^2}`$ ; recette de Mezzadri), indépendamment dans chaque bloc.

D'après [notes de cours, proposition 11.11](../../lecture_notes/fr/11_optimisation.md#r-11-11), $`\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k`$ pour $`U\in U(k)`$ de loi de Haar, et d'après le
[corollaire 11.12](../../lecture_notes/fr/11_optimisation.md#r-11-12), pour $`X,Y\ge0`$ à spectres dans $`[a,b]`$, $`a>0`$, la moyenne quadratique de $`|\mathrm{Tr}(XUY)|`$ vaut
au plus $`\min(1/\sqrt k,b^2/(a^2k))`$ fois son maximum sur $`U`$, et $`\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le b^4/(a^4k^2)`$.

Ainsi, $`N`$ échantillons de Haar trouvent une violation de $`|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)`$ avec une
probabilité d'au plus $`Nb^4/(a^4k^2)`$ (une réunion de $`N`$ événements), bien qu'il en existe une pour tout couple qui
ne commute pas (proposition 3.9). La borne n'est informative que pour $`k>\sqrt N\,b^2/a^2`$ (pour $`b/a=6`$ et $`N=2000`$ :
$`k>1609`$) ; numériquement, les échantillons manquent la violation dès des $`k`$ bien plus petits (ci-dessous). Le code
vérifie la proposition 11.11 par Monte-Carlo (la moyenne empirique est à moins de cinq erreurs standard de la valeur
exacte, pour ce calcul à graine fixée) et affiche, pour des spectres uniformes dans $`[0.2,1.2]`$ ($`b/a=6`$), le rapport du
corollaire 11.12, $`\sup z`$ et le meilleur $`z`$ parmi $`2000`$ échantillons, tous deux rapportés à $`\mathrm{Tr}(XY)`$.

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
sous la borne démontrée $`\min(1/\sqrt k,b^2/(a^2k))`$ et est à peu près divisé par deux chaque fois que $`k`$ double,
c'est-à-dire qu'il décroît comme $`1/k`$, et à partir de $`k=4`$ le meilleur des $`2000`$ échantillons n'atteint même pas
$`z=0`$. Une recherche par gradient, en revanche, s'approche des deux extrêmes dans le calcul ci-dessous (là encore une
observation numérique). Pour le couple `Xf, Yf`, les deux extrêmes en $`U`$ sont connus :

D'après [notes de cours, lemme 11.13](../../lecture_notes/fr/11_optimisation.md#r-11-13), $`\min_U|\mathrm{Tr}(AU)|=0`$ pour toute $`A\in M_k(\mathbb C)`$, $`k\ge2`$.

Ainsi $`\inf_Uz=-\mathrm{Tr}(XY)`$ et $`\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)`$ (proposition 3.9). Un seul appel à
`extremize` approche les deux, avec $`12`$ départs qui maximisent $`z`$ et $`12`$ qui maximisent $`-z`$ (à une précision
relative de $`10^{-3}`$ dans ce calcul, affirmée ci-dessous). D'après la proposition 11.11, la moyenne quadratique de
$`|\mathrm{Tr}(XUY)|`$ vaut $`\|YX\|_2/\sqrt k`$ ; avec les bornes de la démonstration du corollaire 11.12 (spectres
dans $`[a,b]`$), $`\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4`$ et $`\mathrm{Tr}(XY)\ge ka^2`$, donc l'inégalité de Markov
donne $`\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)`$ pour $`\varepsilon>0`$ :
pour $`k`$ grand, $`z`$ est proche de son infimum avec grande probabilité.
*Observation numérique* (affirmée pour ce calcul à graine fixée) : l'échantillonnage arrive à moins de $`1\,\%`$ de
l'infimum, mais n'atteint même pas $`z=0`$, alors que le supremum est positif :

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

*Notes de cours :* [section 11.5](../../lecture_notes/fr/11_optimisation.md#s-11-5) ; [proposition 11.14](../../lecture_notes/fr/11_optimisation.md#r-11-14), [corollaire 11.7](../../lecture_notes/fr/11_optimisation.md#r-11-7), [proposition 11.9](../../lecture_notes/fr/11_optimisation.md#r-11-9), [proposition 11.15](../../lecture_notes/fr/11_optimisation.md#r-11-15), [corollaire 11.16](../../lecture_notes/fr/11_optimisation.md#r-11-16).

`optimize.extremize(objective, params, maximize)` exécute Adam en `rounds` tours de `steps` pas (taux
d'apprentissage `lr * decay**round`) sur la perte $`\mp\sum_b\text{objective}_b`$, recentre les cartes unitaires après
chaque tour et renvoie, **pour chaque départ**, la plus grande (la plus petite, pour `maximize=False`) valeur de
l'objectif parmi celles calculées aux itérés du gradient et une fois au point final. Chaque élément du lot est un
départ indépendant avec ses propres emplacements, signes, $`U_0`$ et contraste : Adam agit coordonnée par coordonnée,
donc si la $`b`$-ième composante de l'objectif ne dépend que des $`b`$-ièmes paramètres (comme dans tous les objectifs
ci-dessous), les départs n'interagissent pas. Les signes de `SelfAdjointParam` sont discrets : entre les tours,
`discrete_step` essaie huit changements de signe individuels par départ (chacun en un emplacement uniformément
aléatoire), l'un après l'autre, et conserve chaque changement qui améliore l'objectif.

En arithmétique exacte, chaque itéré est admissible (la proposition 11.2 et les remarques qui la suivent), donc chaque
valeur renvoyée est la valeur de $`z`$ en un point admissible ; en virgule flottante, c'est cette valeur aux erreurs
d'arrondi près (les assertions ci-dessous tolèrent $`10^{-12}`$). Par conséquent, une maximisation donne une **borne
inférieure** du supremum, une minimisation une borne supérieure de l'infimum, et rien de plus — une recherche locale
peut s'arrêter en un extremum local. L'autre sens exige une démonstration, comme le corollaire 11.7 et la
proposition 11.9, ou une relaxation convexe (section 11.6). Une enveloppe calculée est confrontée (i) à des formes
closes ; (ii) à des configurations construites à la main, telles que des sommes directes de solutions plus petites,
en dessous desquelles elle ne doit pas tomber ; (iii) à des échantillons structurés (spectres à deux niveaux, bases
de permutation), dont les valeurs sont d'autres valeurs atteintes.

D'après [notes de cours, proposition 11.14](../../lecture_notes/fr/11_optimisation.md#r-11-14), si tous les $`k_c\ge2`$, l'enveloppe inférieure d'Exp. 1 est $`\underline z(\delta)=-\big(D-2\delta/(1+\delta)\big)`$.

Le code vérifie la proposition 11.14 sur $`M_2\oplus M_2`$, le corollaire 11.7 sur $`M_3\oplus M_3`$ (les deux contrastes
prescrits, $`k=3`$ impair) et la proposition 11.9 ; suprema et infima partagent un même appel, comme ci-dessus. Les
assertions unilatérales (aucune valeur renvoyée au-delà d'un extrême démontré) sont des conséquences de ces
théorèmes ; les assertions selon lesquelles la recherche arrive à la tolérance indiquée de l'extrême sont des
observations numériques de ce calcul à graine fixée. Le mot-clé `is_self_adjoint` choisit l'algorithme de la norme
trace (valeurs singulières, ou valeurs propres pour $`YXY`$) sans test à chaque pas. La dernière vérification est de
type (ii) : le $`z`$ d'Exp. 3 est additif sur les sommes directes, et le contraste d'une somme directe est le plus grand
contraste de ses termes (exercice 3), donc, pour les enveloppes d'Exp. 3 avec $`\Delta(Y)`$ libre dans $`[0,0.999]`$
comme dans le code,
```math
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
```
en particulier, la somme de deux valeurs atteintes sur $`M_2`$ aux contrastes $`\delta`$ et $`\delta'\le\delta`$ est une
valeur atteinte sur $`M_2\oplus M_2`$ au contraste $`\delta`$.

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
closed = [(1 - dl) * (2 * 0.999 / 1.999) ** 2 for dl in (0.0, 0.5)]        # Corollary 11.16, d = 0.999
print(f"Exp. 3  inf {inf3:.1e};  sup for M_2 at delta = 0, 0.5: {s_M2[0]:.5f}, {s_M2[1]:.5f}"
      f"  (closed form {closed[0]:.5f}, {closed[1]:.5f})")
print(f"        sup for M_2 (+) M_2 at 0.5: search {sup_M22:.5f}, direct sum {direct_sum:.5f}")
assert inf3 >= -1e-12 and inf3 < 1e-6                                       # Proposition 11.9
assert sup_M22 >= direct_sum - 5e-3 * direct_sum                           # (ii), up to 0.5 %
assert (s_M2 <= torch.tensor(closed, dtype=torch.float64) + 1e-12).all()   # Corollary 11.16: one-sided
assert sup_M22 <= (2 - 0.5) * (2 * 0.999 / 1.999) ** 2 + 1e-12
```

*Observations numériques* (les valeurs affichées par ce calcul). La recherche sur $`M_2\oplus M_2`$ reste en dessous de
la somme directe d'environ $`0.2\,\%`$ : elle n'a pas entièrement convergé, et la meilleure borne inférieure de
$`\overline z_{M_2\oplus M_2}(0.5)`$ est le plus grand des deux nombres — détecter cela est le but de la
vérification (ii). Les valeurs pour $`M_2`$ coïncident à $`10^{-4}`$ près avec
$`(1-\delta)\big(2d/(1+d)\big)^2=(1-\delta)(1-\ell(d))^2`$, où $`d=0.999`$ est le plus grand contraste admis pour $`Y`$.
Cette forme close, d'abord observée numériquement dans l'article qui accompagne la bibliothèque, est l'enveloppe
supérieure du [corollaire 11.16](../../lecture_notes/fr/11_optimisation.md#r-11-16) (plus bas) ; le code n'affirme que le sens
démontré (aucune valeur au-dessus). D'après le même corollaire, la somme directe est optimale :
$`\overline z_{M_2\oplus M_2}(0.5)=1.5\,(2d/(1+d))^2`$.
Pour des blocs de taille quelconque, les enveloppes supérieures d'Exp. 1 et d'Exp. 3 sont données par le
[théorème 11.20](../../lecture_notes/fr/11_optimisation.md#r-11-20) ([section 11.7](#s-11-7)).

D'après [notes de cours, proposition 11.15](../../lecture_notes/fr/11_optimisation.md#r-11-15), sur $`M_2`$, pour des spectres
$`x_1\ge x_2`$ de $`X`$ (de $`|X|`$ dans Exp. 3) et $`y_1\ge y_2`$ de $`Y`$, $`\sup z`$ vaut
$`(x_1-x_2)^2(y_1-y_2)^2/(4(x_1+x_2)(y_1+y_2))`$ dans Exp. 1 et $`\frac{2x_1x_2}{x_1+x_2}(y_1-y_2)^2`$ dans Exp. 3. Avec
$`e_1`$ un vecteur propre de $`X`$ pour sa valeur propre de module $`x_1`$ et $`f_1`$ un vecteur propre de $`Y`$ pour $`y_1`$,
le supremum est atteint quand $`\gamma=|\langle e_1,f_1\rangle|^2`$ vaut $`\gamma_1=\frac14+\frac S{2q}`$ (Exp. 1, avec
$`|\mathrm{Tr}(XUY)|=\|YX\|_1`$) ou $`\gamma_3=\frac Sq`$ (Exp. 3, $`X=\mathrm{diag}(x_1,-x_2)`$), où $`S=x_1y_2+x_2y_1`$
et $`q=(x_1+x_2)(y_1+y_2)`$. D'après le
[corollaire 11.16](../../lecture_notes/fr/11_optimisation.md#r-11-16), pour $`C`$ blocs de taille $`2`$ avec $`\Delta(|X|)=\delta`$ et
$`\Delta(Y)=\delta'`$, Exp. 1 a $`\sup z=C\delta^2\delta'^2/((1+\delta)(1+\delta'))`$ et
$`\inf z=-(2C-1+\ell(\delta)\ell(\delta'))`$, et Exp. 3 a $`\sup z=(C-\delta)\big(2\delta'/(1+\delta')\big)^2`$. Les deux
suprema croissent avec $`\delta'`$ : avec $`\Delta(Y)`$ libre dans $`[0,d]`$, les enveloppes supérieures sont les valeurs en
$`\delta'=d`$, et quand $`d\to1`$ elles tendent vers $`C\delta^2/(2(1+\delta))`$ et $`C-\delta`$, les valeurs pour un $`Y`$
singulier ($`\Delta(Y)=1`$).

Le code construit les maximiseurs réels de la démonstration pour $`C=1,2`$ (`from_blocks`), vérifie les contrastes et les
valeurs, vérifie les bornes de la proposition 11.15 sur des paires complexes aléatoires, atteint les limites $`d\to1`$ avec des
projecteurs de rang un dans chaque bloc, et lance une recherche multi-départs pour Exp. 1 sur $`M_2`$ le long de $`\Delta(X)`$
avec $`\Delta(Y)`$ libre dans $`[0,0.999]`$, qui doit rester sous la forme close (un théorème) et, dans ce calcul, s'en approche
à $`10^{-3}`$ près.

```python
def rot(gamma):                                                             # real basis f with |<e_1, f_1>|^2 = gamma
    t = math.acos(math.sqrt(gamma))
    return torch.tensor([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]], dtype=torch.float64)

def gammas(x1, x2, y1, y2):                                                 # gamma_1 (Exp. 1), gamma_3 (Exp. 3)
    S, q = x1 * y2 + x2 * y1, (x1 + x2) * (y1 + y2)
    return 0.25 + S / (2 * q), S / q

env1 = lambda C, dl, dp: C * dl ** 2 * dp ** 2 / ((1 + dl) * (1 + dp))     # Corollary 11.16(a)
env3 = lambda C, dl, dp: (C - dl) * (2 * dp / (1 + dp)) ** 2               # Corollary 11.16(b)

def maximisers(C, dl, dp):
    """Exp. 1: (X, U, Y) and Exp. 3: (X, Y) attaining Corollary 11.16 on C blocks of size 2."""
    a2 = TypeIAlgebra([2] * C, [2] * C, complex_valued=False, precision='double')
    lx, ly = ell(dl), ell(dp)
    R = rot(gammas(1.0, lx, 1.0, ly)[0])
    X0 = torch.diag(torch.tensor([1.0, lx], dtype=torch.float64))
    Y0 = R @ torch.diag(torch.tensor([1.0, ly], dtype=torch.float64)) @ R.T
    W, _, Vh = torch.linalg.svd(Y0 @ X0)                                    # Y X = W Sigma V^T, U = V W^T
    exp1 = [a2.from_blocks([m] * C) for m in (X0, Vh.T @ W.T, Y0)]
    xs = [lx] + [1.0] * (C - 1)                                             # |X|: {1, ell} in block 1, {1, 1} else
    Rs = [rot(gammas(1.0, x2, 1.0, ly)[1]) for x2 in xs]
    X3 = a2.from_blocks([torch.diag(torch.tensor([1.0, -x2], dtype=torch.float64)) for x2 in xs])
    Y3 = a2.from_blocks([Rc @ torch.diag(torch.tensor([1.0, ly], dtype=torch.float64)) @ Rc.T for Rc in Rs])
    return exp1, (X3, Y3)

print(" C  delta  delta'   Exp. 1 z / formula      Exp. 3 z / formula")
for C in (1, 2):
    for dl, dp in ((0.5, 0.8), (0.9, 0.3), (0.4, 1.0)):
        (X, U, Y), (X3, Y3) = maximisers(C, dl, dp)
        assert abs(X.michelson_contrast.item() - dl) < 1e-12 and abs(Y.michelson_contrast.item() - dp) < 1e-12
        assert abs(X3.abs().michelson_contrast.item() - dl) < 1e-12 and abs(Y3.michelson_contrast.item() - dp) < 1e-12
        z1v = ((X @ U @ Y).trace.abs() - (X @ Y).trace.real).item()
        z3v = ((Y3 @ X3.abs() @ Y3).trace.real - (Y3 @ X3 @ Y3).trace_norm()).item()
        print(f" {C}   {dl:.1f}    {dp:.1f}    {z1v:.6f} / {env1(C, dl, dp):.6f}    {z3v:.6f} / {env3(C, dl, dp):.6f}")
        assert abs(z1v - env1(C, dl, dp)) < 1e-12 and abs(z3v - env3(C, dl, dp)) < 1e-12
print("limits d -> 1 at delta = 0.4:", env1(1, 0.4, 1.0), 0.4 ** 2 / (2 * 1.4), env3(1, 0.4, 1.0), 1 - 0.4)

# Proposition 11.15 on random complex pairs in M_2: no pair exceeds the bounds
M2c = TypeIAlgebra([2], [2], precision='double')
x = torch.rand(4000, 2, dtype=torch.float64).sort(-1, descending=True)[0]
y = torch.rand(4000, 2, dtype=torch.float64).sort(-1, descending=True)[0]
V, Wc = M2c.random_unitary_operator(4000).matrix, M2c.random_unitary_operator(4000).matrix
Xr = M2c.operator(V @ torch.diag_embed(x.to(V.dtype)).unsqueeze(1) @ V.mH, is_positive=True)
Xi = M2c.operator(V @ torch.diag_embed((x * torch.tensor([1.0, -1.0])).to(V.dtype)).unsqueeze(1) @ V.mH,
                  is_self_adjoint=True)
Yr = M2c.operator(Wc @ torch.diag_embed(y.to(V.dtype)).unsqueeze(1) @ Wc.mH, is_positive=True)
(x1, x2), (y1, y2) = x.unbind(-1), y.unbind(-1)
sup1 = (Yr @ Xr).trace_norm() - (Xr @ Yr).trace.real                        # sup over U (Theorem 1.9)
z3r = (Yr @ Xi.abs() @ Yr).trace.real - M2c.operator((Yr @ Xi @ Yr).matrix, is_self_adjoint=True).trace_norm()
assert (sup1 <= (x1 - x2) ** 2 * (y1 - y2) ** 2 / (4 * (x1 + x2) * (y1 + y2)) + 1e-12).all()
assert (z3r <= 2 * x1 * x2 / (x1 + x2) * (y1 - y2) ** 2 + 1e-12).all()

# multistart search for Exp. 1 on M_2 along Delta(X), Delta(Y) free in [0, 0.999]
torch.manual_seed(2)
M2r = TypeIAlgebra([2], [2], complex_valued=False, precision='double')
gx = torch.tensor([0.3, 0.6, 0.9], dtype=torch.float64)
Xp, Yp1, Up = opt.PositiveParam(M2r, 3 * S, delta=gx.repeat_interleave(S)), opt.PositiveParam(M2r, 3 * S), \
    opt.UnitaryParam(M2r, 3 * S)
def z1m():
    xm, ym = Xp.operator(), Yp1.operator()
    return (xm @ Up.operator() @ ym).trace.abs() - (xm @ ym).trace.real
sup1m = opt.extremize(z1m, [Xp, Yp1, Up], maximize=True, rounds=4, steps=30).reshape(3, S).max(1)[0]
exact1m = env1(1, gx, 0.999)                                                # Corollary 11.16(c)
print("Exp. 1 sup on M_2:", [f"{v:.5f}" for v in sup1m.tolist()], " closed form:", [f"{v:.5f}" for v in exact1m.tolist()])
assert (sup1m <= exact1m + 1e-12).all() and ((exact1m - sup1m) / exact1m).max() < 1e-3
```

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
paramétrages exacts des contraintes spectrales et de contraste de $`\bigoplus_cM_{k_c}`$ (les extrêmes peuvent se
trouver dans n'importe quel secteur), des unitaires qui respectent les blocs, et un seul calcul par lots pour de
nombreux départs et valeurs de contrainte, sur CPU ou GPU. Elle produit des candidats extrémaux et des bornes
inférieures pour les suprema ; démonstrations ou relaxations doivent fournir le reste.

<a id="s-11-7"></a>

## 11.7 Enveloppes pour des blocs de taille quelconque

*Notes de cours :* [section 11.7](../../lecture_notes/fr/11_optimisation.md#s-11-7) ; [lemme 11.17](../../lecture_notes/fr/11_optimisation.md#r-11-17), [lemme 11.18](../../lecture_notes/fr/11_optimisation.md#r-11-18), [proposition 11.19](../../lecture_notes/fr/11_optimisation.md#r-11-19), [théorème 11.20](../../lecture_notes/fr/11_optimisation.md#r-11-20).

D'après le [théorème 11.20 des notes de cours](../../lecture_notes/fr/11_optimisation.md#r-11-20), avec
$`N=\sum_c\lfloor k_c/2\rfloor`$, $`\Delta(X)=\delta`$ (dans Exp. 3 : $`\Delta(|X|)=\delta`$) et $`\Delta(Y)=\delta'`$, les suprema
d'Exp. 1 et d'Exp. 3 sont
```math
\sup z_1=N\,\frac{\delta^2\delta'^2}{(1+\delta)(1+\delta')},\qquad
\sup z_3=(N-\epsilon\,\delta)\Big(\frac{2\delta'}{1+\delta'}\Big)^2,
```
où $`\epsilon=1`$ si tous les $`k_c`$ sont pairs et $`\epsilon=0`$ sinon.
L'étape clé est la [proposition 11.19](../../lecture_notes/fr/11_optimisation.md#r-11-19) : pour $`X=X_+\oplus(-1)`$ sur
$`\mathbb C^m\oplus\mathbb C^m`$ et $`Y=b+(1-b)Q`$, $`Q`$ le projecteur sur l'image de $`\binom C\Sigma`$, l'opérateur $`YXY`$ a
les valeurs propres d'une matrice $`\mathcal L`$ affine en $`A=C^2`$, et $`z\le(1-b)^2\sum_i2x_i/(1+x_i)`$. Le premier bloc le
vérifie sur $`M_6`$ ($`m=3`$, $`X_+`$ de valeurs propres $`0.3,1,1`$ dans une base aléatoire) : les spectres de $`YXY`$
(`eigenvalues`) et de $`\mathcal L`$ coïncident pour des angles aléatoires, aucun angle ne dépasse la borne, et les angles
du cas d'égalité l'atteignent.

```python
# Proposition 11.19 on M_6, m = 3: X = X_+ (+) (-1), Y = b + (1 - b) Q, Q the projection onto the range of (C; Sigma)
torch.manual_seed(4)
m, b = 3, 0.2
M6 = TypeIAlgebra([2 * m], [2 * m], precision='double')
W = TypeIAlgebra([m], [m], precision='double').random_unitary_operator(1).matrix[0, 0]
x = torch.tensor([0.3, 1.0, 1.0], dtype=torch.float64)                    # X_+ = 1 - (1 - a) e e^*, a = 0.3
Xp = W @ torch.diag(x).to(W.dtype) @ W.mH

def model(theta, Xp):
    v = torch.cat([torch.diag(torch.cos(theta)), torch.diag(torch.sin(theta))]).to(Xp.dtype)
    Xm = torch.block_diag(Xp, -torch.eye(m, dtype=Xp.dtype))
    Ym = b * torch.eye(2 * m, dtype=Xp.dtype) + (1 - b) * v @ v.mH
    return M6.operator(Xm[None, None], is_self_adjoint=True), M6.operator(Ym[None, None], is_positive=True)

def lin(theta, Xp):                                                        # the matrix L of the proof
    w, V = torch.linalg.eigh(Xp)
    f = lambda g: V @ torch.diag(g(w)).to(V.dtype) @ V.mH
    A = torch.diag(torch.cos(theta) ** 2).to(Xp.dtype)
    K12, X12 = f(lambda t: (1 + t).sqrt()), f(torch.sqrt)
    top = (1 - b ** 2) * K12 @ A @ K12 + b ** 2 * Xp - torch.eye(m, dtype=Xp.dtype)
    return torch.cat([torch.cat([top, b * X12], 1), torch.cat([b * X12, torch.zeros_like(Xp)], 1)])

bound = (1 - b) ** 2 * (2 * x / (1 + x)).sum().item()                       # Proposition 11.19
worst = 0.0
for _ in range(200):
    theta = torch.rand(m, dtype=torch.float64) * math.pi / 2
    X, Y = model(theta, Xp)
    YXY = M6.operator((Y @ X @ Y).matrix, is_self_adjoint=True)
    assert torch.allclose(YXY.eigenvalues()[0][0], torch.linalg.eigvalsh(lin(theta, Xp)), atol=1e-12)
    worst = max(worst, ((Y @ X.abs() @ Y).trace.real - YXY.trace_norm()).item())
th_star = torch.acos(((1 + b * x) / ((1 + x) * (1 + b))).sqrt())            # equality case, X_+ diagonal
X, Y = model(th_star, torch.diag(x).to(W.dtype))
z_star = ((Y @ X.abs() @ Y).trace.real - (Y @ X @ Y).trace_norm()).item()
print(f"Proposition 11.19: largest z at 200 random angles {worst:.6f} <= {bound:.6f}; at the maximiser {z_star:.12f}")
assert worst <= bound + 1e-12 and abs(z_star - bound) < 1e-12
```

Le deuxième bloc construit avec `from_blocks` les maximiseurs de la démonstration du théorème 11.20 pour plusieurs tailles
de blocs : des blocs $`2\times2`$ comme dans la proposition 11.15 et des blocs $`1\times1`$ ; dans Exp. 3, la valeur propre
$`\ell(\delta)`$ de $`|X|`$ est sur un bloc $`2\times2`$ si tous les $`k_c`$ sont pairs, sur un bloc $`1\times1`$ sinon ; dans
Exp. 1, l'unitaire est l'adjoint du facteur polaire de $`YX`$ (`polar`). Les contrastes et les valeurs doivent coïncider avec
le théorème à $`10^{-12}`$ près.

```python
def sup_exact(exp, ks, dl, dp):                                             # Theorem 11.20
    N = sum(k // 2 for k in ks)
    if exp == 1:
        return N * dl ** 2 * dp ** 2 / ((1 + dl) * (1 + dp))
    return (N - dl * all(k % 2 == 0 for k in ks)) * (2 * dp / (1 + dp)) ** 2

def maximisers_any(ks, dl, dp):
    """Exp. 1 (X, U, Y) and Exp. 3 (X, Y) attaining Theorem 11.20: 2 x 2 maximisers of Proposition 11.15 and 1 x 1
    blocks; in Exp. 3 the eigenvalue ell(delta) of |X| on a 2 x 2 block if all k_c are even, else on a 1 x 1 block."""
    lx, ly = ell(dl), ell(dp)
    even = all(k % 2 == 0 for k in ks)
    two = lambda d: torch.diag(torch.tensor(d, dtype=torch.float64))
    R1 = rot(gammas(1.0, lx, 1.0, ly)[0])
    blocks = {name: [] for name in ('X1', 'Y1', 'X3', 'Y3')}
    placed = False
    for k in ks:
        b1 = {name: [] for name in blocks}
        for j in range(k // 2):
            x3 = lx if even and not placed else 1.0
            placed = placed or x3 == lx
            R3 = rot(gammas(1.0, x3, 1.0, ly)[1])
            b1['X1'].append(two([1.0, lx]))
            b1['Y1'].append(R1 @ two([1.0, ly]) @ R1.T)
            b1['X3'].append(two([1.0, -x3]))
            b1['Y3'].append(R3 @ two([1.0, ly]) @ R3.T)
        if k % 2:
            b1['X1'].append(two([1.0]))
            b1['Y1'].append(two([1.0]))
            b1['X3'].append(two([1.0 if placed else lx]))
            b1['Y3'].append(two([1.0]))
            placed = True
        for name in blocks:
            blocks[name].append(torch.block_diag(*b1[name]))
    a = TypeIAlgebra(ks, ks, complex_valued=False, precision='double')
    X1, Y1, X3, Y3 = (a.from_blocks(blocks[name]) for name in ('X1', 'Y1', 'X3', 'Y3'))
    U1 = (Y1 @ X1).polar()[0].adjoint()                                      # U = W^*, Y X = W |Y X|
    return (X1, U1, Y1), (X3, Y3)

print(" sizes      delta  delta'   Exp. 1 z / Theorem 11.20    Exp. 3 z / Theorem 11.20")
for ks in ([4, 4], [3, 2], [2, 2, 2], [16], [3]):
    for dl, dp in ((0.5, 0.8), (0.9, 0.999)):
        (X, U, Y), (X3, Y3) = maximisers_any(ks, dl, dp)
        assert abs(X.michelson_contrast.item() - dl) < 1e-12 and abs(Y.michelson_contrast.item() - dp) < 1e-12
        assert abs(X3.abs().michelson_contrast.item() - dl) < 1e-12 and abs(Y3.michelson_contrast.item() - dp) < 1e-12
        z1v = ((X @ U @ Y).trace.abs() - (X @ Y).trace.real).item()
        z3v = ((Y3 @ X3.abs() @ Y3).trace.real - (Y3 @ X3 @ Y3).trace_norm()).item()
        print(f" {str(ks):10s} {dl:.1f}    {dp:.3f}   {z1v:.6f} / {sup_exact(1, ks, dl, dp):.6f}"
              f"       {z3v:.6f} / {sup_exact(3, ks, dl, dp):.6f}")
        assert abs(z1v - sup_exact(1, ks, dl, dp)) < 1e-12 and abs(z3v - sup_exact(3, ks, dl, dp)) < 1e-12
```

Des points de départ aléatoires de Haar s'approchent rarement de ces maximiseurs quand il y a beaucoup de blocs. Avec
`init='diagonal'` pour $`X`$, `init='pairs'` pour $`Y`$ (la base propre relative est alors une somme directe de rotations
$`2\times2`$), `spectrum_init` (spectres à deux niveaux ; dans Exp. 3, $`|X|`$ proche de $`1`$ et une valeur propre positive et
une négative de $`X`$ par paire) et `init='rotations'` pour $`U`$ (des rotations des mêmes paires, de déterminant $`+1`$ : dans une algèbre réelle, une carte $`U_0e^{G-G^*}`$ reste
dans la composante connexe de $`U_0`$, et le maximiseur d'Exp. 1 est de déterminant $`+1`$ dans chaque secteur), la moitié
des départs commence sur cette structure, l'autre moitié en des points aléatoires de Haar (`init_fraction` ; le masque
`X.structured` fait partir $`Y`$ et $`U`$ structurés aux mêmes départs). Les recherches ci-dessous restent sous le
théorème 11.20 (affirmé : un théorème) et, dans cette exécution à graine fixée, s'en approchent à $`10^{-3}`$ près (une
observation numérique).

```python
# multistart search along Delta(X), resp. Delta(|X|), with Delta(Y) free in [0, 0.999]; half of the starts structured
def search_any(exp, ks, deltas, starts=8, seed=0):
    torch.manual_seed(seed)
    a = TypeIAlgebra(ks, ks, complex_valued=False, precision='double')
    d = torch.tensor(deltas, dtype=torch.float64).repeat_interleave(starts)
    cls = opt.SelfAdjointParam if exp == 3 else opt.PositiveParam
    Xs = cls(a, len(d), delta=d, init='diagonal', init_fraction=0.5, spectrum_init='high' if exp == 3 else 'pairs')
    Ys = opt.PositiveParam(a, len(d), init='pairs', init_fraction=Xs.structured, spectrum_init='pairs')
    Us = opt.UnitaryParam(a, len(d), init='rotations', init_fraction=Xs.structured)
    def z():
        x, y = Xs.operator(), Ys.operator()
        if exp == 1:
            return (x @ Us.operator() @ y).trace.abs() - (x @ y).trace.real
        M = a.operator((y @ x @ y).matrix, is_self_adjoint=True)
        return (y @ Xs.abs_operator() @ y).trace.real - M.trace_norm()
    params = [Xs, Ys, Us] if exp == 1 else [Xs, Ys]
    vals = opt.extremize(z, params, maximize=True, rounds=6, steps=40).reshape(len(deltas), starts)
    return vals.nan_to_num(nan=-math.inf).max(1)[0]

for exp, ks in ((3, [4, 4]), (3, [3]), (1, [2] * 6)):
    deltas = [0.3, 0.8]
    found = search_any(exp, ks, deltas)
    exact = torch.tensor([sup_exact(exp, ks, dl, 0.999) for dl in deltas], dtype=torch.float64)
    print(f"Exp. {exp} on {ks}: search {[round(v, 6) for v in found.tolist()]}, Theorem 11.20 "
          f"{[round(v, 6) for v in exact.tolist()]}")
    assert (found <= exact + 1e-12).all() and ((exact - found) / exact).max() < 1e-3
```

<a id="x-11"></a>

## Exercices

[Exercices 1 et 2](../../lecture_notes/fr/11_optimisation.md#x-11) des notes de cours, numériquement : la transformation de Cayley $`C(S)=(1-S)(1+S)^{-1}`$ d'un
$`S`$ antiautoadjoint est unitaire, n'a pas de valeur propre $`-1`$, et $`S=(1-C)(1+C)^{-1}`$ (exercice 1) ; pour $`O`$ distribué selon la
loi de Haar sur $`O(k)`$, tiré par `random_unitary_operator` d'une algèbre réelle (`complex_valued=False`),
$`\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k`$ (exercice 2).

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
