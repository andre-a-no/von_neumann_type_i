[English](../08_dynamics.md) · [Español](../es/08_dynamics.md) · **Français** · [Deutsch](../de/08_dynamics.md) · [中文](../zh/08_dynamics.md) · [日本語](../ja/08_dynamics.md) · [Русский](../ru/08_dynamics.md)

# Tutoriel 8. Dynamique

Ce tutoriel accompagne le [chapitre 8](../../lecture_notes/fr/08_dynamics.md) des notes de cours, qui démontre qu'un hamiltonien
de $`\mathcal M`$ conserve les poids des secteurs ([théorème 8.3](../../lecture_notes/fr/08_dynamics.md#r-8-3)), que les générateurs de forme GKSL produisent
des semi-groupes dynamiques quantiques ([théorème 8.6](../../lecture_notes/fr/08_dynamics.md#r-8-6)), et traite des sauts entre secteurs, de la post-sélection et
des états stationnaires. Nous utilisons les solveurs exacts et de Runge–Kutta de `dynamics` pour les équations de
Schrödinger, de von Neumann et de Lindblad, et comparons les résultats à SciPy et à des générateurs GKSL sur $`M_D`$.

<a id="s-8-1"></a>

## 8.1 Systèmes fermés

*Notes de cours :* [section 8.1](../../lecture_notes/fr/08_dynamics.md#s-8-1) ; [proposition 8.2](../../lecture_notes/fr/08_dynamics.md#r-8-2) (le propagateur $`U(t)=e^{-itH}`$), [théorème 8.3](../../lecture_notes/fr/08_dynamics.md#r-8-3) (conservation de la charge).

D'après la proposition 8.2(1), l'évolution se scinde en $`C`$ problèmes indépendants de tailles $`k_c`$,
un par bloc. La bibliothèque calcule $`U(t)`$ avec une décomposition spectrale par taille de bloc
(`dynamics.propagator`), et `dynamics.von_neumann` fait évoluer une densité jusqu'à de nombreux
instants à la fois dans la base propre, où $`\rho(t)_{ij}=e^{-it(\lambda_i-\lambda_j)}\rho_{ij}`$.
Nous comparons avec `scipy.linalg.expm` :

```python
import math
import numpy as np
import scipy.integrate
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, InterSectorChannel, channels, dynamics, states

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
alg = TypeIAlgebra([2, 3, 4], [2, 3, 4], precision='double')     # M_2 (+) M_3 (+) M_4, D = 9
H = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
t = 0.7
U = dynamics.propagator(H, t)                                      # exp(-itH), block by block
for c, k in enumerate(alg.k_factors):                              # independent reference
    ref = sla.expm(-1j * t * H.matrix[0, c, :k, :k].numpy())
    assert np.abs(U.matrix[0, c, :k, :k].numpy() - ref).max() < 1e-12
assert (U.adjoint() @ U - alg.identity()).matrix.abs().max() < 1e-12

rho0 = states.random_density_matrix(alg)
times = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.von_neumann(rho0, H, times)                        # exact, all times at once
for s, rho_s in zip(times, traj):
    U_s = dynamics.propagator(H, s)
    assert (rho_s.matrix - (U_s @ rho0 @ U_s.adjoint()).matrix).abs().max() < 1e-12
    assert (rho_s.sector_probabilities() - rho0.sector_probabilities()).abs().max() < 1e-12   # Theorem 8.3
A = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
A_t = U.adjoint() @ A @ U                                          # Heisenberg picture
assert abs(rho0.expectation(A_t) - dynamics.von_neumann(rho0, H, [t])[0].expectation(A)) < 1e-12
```

Un vecteur unitaire $`\psi=\bigoplus_c\psi_c`$ définit sur $`\mathcal M`$ l'état $`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$, qui
évolue par $`U(t)`$ ([notes de cours, section 8.1](../../lecture_notes/fr/08_dynamics.md#s-8-1)), de sorte que `dynamics.schrodinger` (vecteurs dans la
disposition par blocs `(batch, C, k_max)`) et `von_neumann` doivent coïncider sur lui. Pour un $`H(t)`$ dépendant du temps, la
bibliothèque utilise la méthode de Runge–Kutta classique à pas fixe (`dynamics.rk4`,
`schrodinger_rk4`) : un pas multiplie $`\psi`$ par un polynôme en des valeurs de $`H(t)\in\mathcal M`$, il
ne déplace donc jamais d'amplitude entre secteurs, mais il n'est pas unitaire : pour $`H`$ constant, un
pas de taille $`h`$ multiplie un vecteur propre de $`H`$ de valeur propre $`\lambda`$ par $`T(-ih\lambda)`$,
$`T(z)=\sum_{k=0}^4z^k/k!`$, et $`|T(ix)|^2=1-\frac{x^6}{72}+\frac{x^8}{576}<1`$ pour $`0<|x|<2\sqrt2`$
(développer). Ainsi RK4 ne conserve pas exactement les poids des secteurs : pour $`H`$ constant, une
composante de secteur ayant un coefficient non nul sur un vecteur propre avec $`0<|h\lambda|<2\sqrt2`$
perd du poids à chaque pas. L'ampleur de leur variation au cours d'un calcul est une observation
numérique (ci-dessous). (Pour un seul grand secteur, `krylov.evolve` n'a besoin que de produits
matrice–vecteur ; voir le [Tutoriel 9, section 9.5](../../tutorial/fr/09_spin_chains.md#s-9-5).)

```python
rows = torch.tensor([[1, 1, 0, 0], [1, 1, 1, 0], [1, 1, 1, 1]], dtype=c128)   # entries of C^2 (+) C^3 (+) C^4
psi0 = torch.randn(1, alg.C, alg.k_max, dtype=c128) * rows
psi0 = psi0 / psi0.norm()
psi = dynamics.schrodinger(psi0, H, times)                         # shape (T, batch, C, k_max)
p_vec = (psi.abs() ** 2).sum(-1)                                   # ||P_c psi(t)||^2
assert (p_vec - p_vec[0]).abs().max() < 1e-12
mix = DensityMatrix(alg, matrix=psi0[..., :, None] * psi0[..., None, :].conj())   # the state of psi0 on M
outer = psi[..., :, None] * psi[..., None, :].conj()
assert all((o - r.matrix).abs().max() < 1e-12 for o, r in zip(outer, dynamics.von_neumann(mix, H, times)))

V = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
H_of_t = lambda s: H.matrix + math.cos(3 * s) * V.matrix          # driven, in M at every time
psi_rk = dynamics.schrodinger_rk4(psi0, H_of_t, times, substeps=50, algebra=alg)
drift = ((psi_rk.abs() ** 2).sum(-1) - p_vec[0]).abs().max().item()
print(f"RK4, time-dependent H: largest change of a sector weight {drift:.1e}")
assert drift < 1e-8 and (dynamics.schrodinger_rk4(psi0, H, times, substeps=50) - psi).abs().max() < 1e-7

```

Observation numérique : dans ce calcul, les poids des secteurs varient de $`1.1\cdot10^{-9}`$ avec RK4
(la valeur affichée) et de moins que la tolérance $`10^{-12}`$ des `assert` avec les solveurs exacts,
pour lesquels le théorème 8.3 donne une conservation exacte.

<a id="s-8-2"></a>

## 8.2 Semi-groupes dynamiques quantiques

*Notes de cours :* [section 8.2](../../lecture_notes/fr/08_dynamics.md#s-8-2) ; [théorème 8.6](../../lecture_notes/fr/08_dynamics.md#r-8-6) (générateurs GKSL et développement de Dyson), [théorème 8.7](../../lecture_notes/fr/08_dynamics.md#r-8-7).

La bibliothèque propose quatre routines de Lindblad : `lindblad_rhs` renvoie le second membre sous
forme de fonction pour `rk4` ; `lindblad_evolve` intègre avec RK4 ; `lindblad_superoperator`
renvoie, pour chaque secteur, la matrice $`k_c^2\times k_c^2`$ de $`\mathcal L`$ agissant sur
$`\mathrm{vec}\,\rho_c`$ rangé ligne par ligne ; et `lindblad_channel` calcule $`e^{t\mathcal L}`$
exactement, secteur par secteur, sous forme de `Channel`. Nous vérifions le théorème 8.6 pour des
sauts arbitraires sur $`M_3`$ : la matrice de Choi de $`e^{t\mathcal L}`$ est positive (théorème 6.2),
l'application est TP, la loi de semi-groupe est satisfaite, RK4 concorde avec l'exponentielle, et la
partie due aux sauts $`e^{t\mathcal L}-S_0(t)=\sum_{n\ge1}S_n(t)`$ est elle-même CP :

```python
M3 = TypeIAlgebra([3], [3], precision='double')
H3 = M3.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
jumps = [M3.operator(torch.randn(1, 1, 3, 3, dtype=c128) / 2) for _ in range(2)]     # arbitrary L_j
T = lambda s: dynamics.lindblad_channel(M3, H3, jumps, t=s)        # exp(sL), exact
S = lambda s: T(s).superoperator()[0][0]                           # its 9 x 9 matrix on row-major vec
choi = lambda X: X.reshape(3, 3, 3, 3).permute(2, 0, 3, 1).reshape(9, 9)   # J_{(i,a),(j,b)} = S_{(a,b),(i,j)}
assert (choi(S(0.8)) - T(0.8).choi()[0][0]).abs().max() < 1e-12
assert torch.linalg.eigvalsh(choi(S(0.8))).min() > -1e-12         # completely positive
assert T(0.8).trace_preservation_error().max() < 1e-12             # trace preserving
assert (S(0.3) @ S(0.5) - S(0.8)).abs().max() < 1e-12              # semigroup

rho3 = states.random_density_matrix(M3)
rk = dynamics.lindblad_evolve(rho3, H3, jumps, [0.0, 0.8], substeps=40)[-1]
assert (rk.matrix - T(0.8)(rho3).matrix).abs().max() < 1e-7        # RK4 against exp(tL)
assert abs(dynamics.lindblad_rhs(H3, jumps)(0.0, rho3.matrix)[0, 0].trace()) < 1e-12   # Tr L(rho) = 0

H_eff = H3.matrix - 0.5j * sum(L.adjoint().matrix @ L.matrix for L in jumps)
S0 = channels.Channel(M3, torch.linalg.matrix_exp(-0.8j * H_eff)[:, None]).superoperator()[0][0]
jump_part = torch.linalg.eigvalsh(choi(S(0.8) - S0))
print(f"smallest eigenvalue of the Choi matrix of the jump part: {jump_part.min().item():.1e}")
assert jump_part.min() > -1e-12                                    # sum_{n >= 1} S_n(t) is CP
```

<a id="s-8-3"></a>

## 8.3 Équations de Lindblad dans une algèbre à secteurs

*Notes de cours :* [section 8.3](../../lecture_notes/fr/08_dynamics.md#s-8-3) : le générateur $`\mathcal L`$ avec des sauts à l'intérieur des secteurs et entre eux, qu'implémente `dynamics.lindblad_rhs` ; [proposition 8.8](../../lecture_notes/fr/08_dynamics.md#r-8-8), [proposition 8.9](../../lecture_notes/fr/08_dynamics.md#r-8-9) (poids des secteurs).

**Exemple : perte de particules** ([notes de cours, section 8.3](../../lecture_notes/fr/08_dynamics.md#s-8-3)). Deux modes bosoniques avec au plus deux
particules, $`\mathcal M=M_1\oplus M_2\oplus M_3`$ (secteurs $`N=0,1,2`$), le dimère de Bose–Hubbard
$`H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)`$, et chaque mode perd des particules au taux $`\gamma`$ : un
`InterSectorChannel` d'opérateurs de Kraus $`a_1,a_2`$, restreints à $`N\to N-1`$. D'après la
[proposition 8.9(2)](../../lecture_notes/fr/08_dynamics.md#r-8-9), les poids des secteurs obéissent à $`p'=Qp`$ avec $`Q=R-\mathrm{diag}(\sum_dR_{dc})_c`$ pour le taux $`1`$,
$`R`$ = `transition_matrix()`, et en partant de $`N=2`$ ils suivent la loi binomiale $`p_2=q^2`$, $`p_1=2q(1-q)`$, $`p_0=(1-q)^2`$,
$`q=e^{-\gamma t}`$, quels que soient les effets de l'effet tunnel entre modes (hopping) et de l'interaction. Le code vérifie
$`\Phi^*(1)=0\oplus1\oplus2`$, c'est-à-dire $`a_1^*a_1+a_2^*a_2=N`$ sur les blocs utilisés, et compare $`p(t)`$ aux deux formules.

```python
fock = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double', charges=[0, 1, 2])   # N = 0, 1, 2
r2 = math.sqrt(2)
a1 = {(0, 1): [[1.0, 0.0]], (1, 2): [[r2, 0.0, 0.0], [0.0, 1.0, 0.0]]}           # a_1 : N -> N - 1
a2 = {(0, 1): [[0.0, 1.0]], (1, 2): [[0.0, 1.0, 0.0], [0.0, 0.0, r2]]}           # a_2 : N -> N - 1
blocks = lambda a: {dc: torch.tensor(K, dtype=f64) for dc, K in a.items()}
loss = InterSectorChannel.from_blocks(fock, fock, {dc: [blocks(a1)[dc], blocks(a2)[dc]] for dc in a1})
J, U_int = 1.0, 0.8                                                # hopping and interaction
H_bh = fock.from_blocks([torch.zeros(1, 1, dtype=f64), torch.tensor([[0.0, J], [J, 0.0]], dtype=f64),
                         torch.tensor([[U_int, r2 * J, 0.0], [r2 * J, 0.0, r2 * J], [0.0, r2 * J, U_int]], dtype=f64)])

assert (loss.adjoint()(fock.identity()).matrix - fock.central([0.0, 1.0, 2.0]).matrix).abs().max() < 1e-12
R = loss.transition_matrix()[0]                                    # r_dc (columns: from sector c)
Q = R - torch.diag(R.sum(0))                                       # Proposition 8.9(2), for rate 1
print(Q)

gamma = 0.5
psi2 = torch.randn(3, dtype=c128)
psi2 = psi2 / psi2.norm()
m0 = torch.zeros(1, 3, 3, 3, dtype=c128)
m0[0, 2] = torch.outer(psi2, psi2.conj())                          # a pure state with N = 2
ts = [0.0, 0.5, 1.0, 2.0, 4.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m0), H_bh, [loss], ts, rates=[gamma], substeps=100)
for s, rho_s in zip(ts, traj):
    q = math.exp(-gamma * s)
    p = rho_s.sector_probabilities()[0]
    assert (p - torch.tensor([(1 - q) ** 2, 2 * q * (1 - q), q ** 2], dtype=f64)).abs().max() < 1e-8
    assert (p - torch.from_numpy(sla.expm(gamma * s * Q.numpy())[:, 2])).abs().max() < 1e-8
neg = min(w.min().item() for rho_s in traj for w in rho_s.eigenvalues())
print(f"smallest eigenvalue of rho(t) under RK4: {neg:.1e}")
assert neg > -1e-7                                                 # positive up to the RK4 error

sm = torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=c128)           # |g><e|; H = 0, rate 3, one step h = 1
decay = lambda s, X: 3.0 * (sm @ X @ sm.mH - 0.5 * (sm.mH @ sm @ X + X @ sm.mH @ sm))
step = dynamics.rk4(decay, torch.diag(torch.tensor([0.0, 1.0], dtype=c128)), [0.0, 1.0], substeps=1)[-1]
assert (step - torch.diag(torch.tensor([-3 / 8, 11 / 8], dtype=c128))).abs().max() < 1e-14   # not positive
```

RK4 ne préserve pas nécessairement la positivité : pour $`\rho'=\mathcal L\rho`$, un pas de taille $`h`$
est $`\rho\mapsto T(h\mathcal L)\rho`$ avec le polynôme $`T`$ de la section 8.1 (insérer les quatre
étapes), qui préserve la trace (puisque $`\mathrm{Tr}\,\mathcal L(X)=0`$ pour tout $`X`$) mais n'est
pas nécessairement positif. Par exemple, sur $`M_2`$ de base $`|g\rangle,|e\rangle`$, avec $`H=0`$ et un saut
$`|g\rangle\langle e|`$ de taux $`\gamma`$, $`\mathcal L`$ envoie $`\mathrm{diag}(a,b)`$ sur
$`\gamma b\mathrm{diag}(1,-1)`$, donc un pas avec $`h\gamma=3`$ envoie $`\mathrm{diag}(0,1)`$ sur
$`\mathrm{diag}(1-T(-3),T(-3))=\mathrm{diag}(-\frac38,\frac{11}8)`$. Observation numérique :
dans ce calcul, la plus petite valeur propre de $`\rho(t)`$ vaut $`\approx-2\cdot10^{-8}`$ (la valeur
affichée, avec des pas allant jusqu'à $`0.02`$). `lindblad_channel` est exact mais ne traite que les
sauts à l'intérieur des secteurs. Enfin, nous vérifions la proposition 8.8(2) pour une application
aléatoire entre les secteurs et la proposition 8.8(3) pour la perte de particules, en les comparant à
des générateurs GKSL sur $`M_6`$ construits avec numpy :

```python
Phi = channels.random_inter_sector_channel(fock, fock, kraus_rank=2)
rho = states.random_density_matrix(fock)
on_M = dynamics.lindblad_rhs(None, [Phi], algebra=fock)(0.0, rho.matrix)[0].numpy()
ks, K = fock.k_factors, Phi.kraus[0].numpy()                       # K[l, d, c] = K_l^{dc}, zero padded
embed = lambda Y: sla.block_diag(*[Y[c, :k, :k] for c, k in enumerate(ks)])    # (+)_c Y_c in M_6
Ls = [np.block([[K[l, d, c, :kd, :kc] for c, kc in enumerate(ks)] for d, kd in enumerate(ks)])
      for l in range(2)]                                           # L_l = sum_{d,c} J_l^{dc}
X = embed(rho.matrix[0].numpy())
full = sum(L @ X @ L.conj().T - 0.5 * (L.conj().T @ L @ X + X @ L.conj().T @ L) for L in Ls)   # GKSL on M_6
E_full = full * embed(np.ones((3, 3, 3)))                          # E(X) = sum_d P_d X P_d
assert np.abs(E_full - embed(on_M)).max() < 1e-12                  # Proposition 8.8(2)
print(f"coherences between sectors created on M_6: {np.abs(full - E_full).max():.2f}")

zero = lambda d, c: np.zeros((ks[d], ks[c]))
L_loss = [np.block([[blocks(a)[(d, c)].numpy() if (d, c) in a else zero(d, c) for c in range(3)] for d in range(3)])
          for a in (a1, a2)]                                       # a_1, a_2 on M_6: one block per block row and column
full_loss = sum(L @ X @ L.T - 0.5 * (L.T @ L @ X + X @ L.T @ L) for L in L_loss)
loss_on_M = dynamics.lindblad_rhs(None, [loss], algebra=fock)(0.0, rho.matrix)[0].numpy()
assert np.abs(full_loss - embed(loss_on_M)).max() < 1e-12         # Proposition 8.8(3): no compression needed
```

<a id="s-8-4"></a>

## 8.4 Post-sélection et évolution sans saut

*Notes de cours :* [section 8.4](../../lecture_notes/fr/08_dynamics.md#s-8-4) ; [proposition 8.10](../../lecture_notes/fr/08_dynamics.md#r-8-10) (la probabilité sans saut $`w_0`$), [proposition 8.11](../../lecture_notes/fr/08_dynamics.md#r-8-11) (post-sélection sur « aucune particule perdue »).

La fonction `states.lueders_update(rho, P)` de la bibliothèque renvoie
$`P\rho P/\mathrm{Tr}(P\rho P)`$ et $`\mathrm{Tr}(\rho P)`$. Supposons maintenant que les
deux modes fuient à des taux différents $`\gamma_1>\gamma_2`$ ; alors $`A=\gamma_1n_1+\gamma_2n_2`$
n'est pas central (dans le secteur $`N=2`$, c'est $`\mathrm{diag}(2\gamma_1,\gamma_1+\gamma_2,2\gamma_2)`$
dans la base $`|20\rangle,|11\rangle,|02\rangle`$), la proposition 8.9(2) ne s'applique pas, et l'état
post-sélectionné sur $`N=2`$ est l'évolution normalisée par $`H_{\mathrm{eff}}`$, dont la partie
anti-hermitienne $`-\frac i2A`$ n'est pas un multiple de $`1`$ dans ce secteur. Nous vérifions la proposition 8.11 pour $`n=0`$ et, en calculant
l'intégrale de $`S_1`$ par la méthode de Simpson, pour $`n=1`$ :

```python
g1, g2 = 0.6, 0.2                                                  # mode 1 leaks faster
loss1 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a1).items()})
loss2 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a2).items()})
m1 = torch.zeros(1, 3, 3, 3, dtype=c128)
m1[0, 2, 0, 0] = 1.0                                               # |20>: both particles in mode 1
ts = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m1), H_bh, [loss1, loss2], ts, rates=[g1, g2],
                                substeps=100)

n_1 = [np.diag(x) for x in ([0.0], [1.0, 0.0], [2.0, 1.0, 0.0])]   # n_1 in each sector; n_2 = N - n_1
H_eff = [H_bh.matrix[0, N, :N + 1, :N + 1].numpy() - 0.5j * (g1 * n + g2 * (N * np.eye(N + 1) - n))
         for N, n in enumerate(n_1)]

def S0(s, N, X):                                  # e^{-isH_eff} X e^{isH_eff^*} in sector N, s a time or an array
    E = torch.linalg.matrix_exp(torch.from_numpy(-1j * np.asarray(s)[..., None, None] * H_eff[N])).numpy()
    return E @ X @ E.conj().swapaxes(-1, -2)

rho_20 = m1[0, 2].numpy()
for s, rho_s in zip(ts, traj):
    w0 = np.trace(S0(s, 2, rho_20)).real                           # probability of no loss
    post, prob = states.lueders_update(rho_s, fock.central([0.0, 0.0, 1.0]))   # post-select N = 2
    assert abs(prob.item() - w0) < 1e-8                            # Proposition 8.11, n = 0
    assert np.abs(post.matrix[0, 2].numpy() - S0(s, 2, rho_20) / w0).max() < 1e-7
closed = dynamics.von_neumann(DensityMatrix(fock, matrix=m1), H_bh, [2.0])[0]   # no losses at all
dist = 0.5 * np.abs(np.linalg.eigvalsh(S0(2.0, 2, rho_20) / w0 - closed.matrix[0, 2].numpy())).sum()
print(f"t = 2: P(no loss) = {w0:.3f} (without hopping: {math.exp(-2 * g1 * 2):.3f}); "
      f"trace distance post-selected/closed system: {dist:.3f}")

K1, K2 = blocks(a1)[(1, 2)].numpy(), blocks(a2)[(1, 2)].numpy()   # S_1(2) rho = int_0^2 S_0(2-s) Phi S_0(s) rho ds
jump = lambda X: g1 * K1 @ X @ K1.T + g2 * K2 @ X @ K2.T
s_grid = np.linspace(0.0, 2.0, 201)
S1 = scipy.integrate.simpson(S0(2.0 - s_grid, 1, jump(S0(s_grid, 2, rho_20))), x=s_grid, axis=0)
assert np.abs(S1 - traj[-1].matrix[0, 1, :2, :2].numpy()).max() < 1e-8     # Proposition 8.11, n = 1
```

Observation numérique : à $`t=2`$, le poids $`p_2(2)=w_0(2)`$ du secteur $`N=2`$ (« aucune perte ») vaut
$`0.230`$ (affiché ci-dessus). Sans effet tunnel ($`J=0`$), il vaudrait exactement
$`e^{-2\gamma_1t}=e^{-2.4}=0.091`$, puisqu'alors $`|20\rangle`$ est un vecteur propre de
$`H_{\mathrm{eff}}`$ de valeur propre $`U-i\gamma_1`$ et $`w_0(t)=|e^{-it(U-i\gamma_1)}|^2`$ (proposition
8.10) ; ici, les poids des secteurs dépendent donc bel et bien du hamiltonien. L'état post-sélectionné
est à distance de trace $`0.135`$ (affichée ci-dessus) de l'état du système fermé : conditionner par
$`N=2`$ modifie l'état, au-delà de l'évolution unitaire.

<a id="s-8-5"></a>

## 8.5 États stationnaires et relaxation

*Notes de cours :* [section 8.5](../../lecture_notes/fr/08_dynamics.md#s-8-5) ; [corollaire 8.13](../../lecture_notes/fr/08_dynamics.md#r-8-13), [proposition 8.14](../../lecture_notes/fr/08_dynamics.md#r-8-14) (émission spontanée).

D'après la [proposition 8.14](../../lecture_notes/fr/08_dynamics.md#r-8-14), l'émission spontanée sur $`M_2`$ ($`H=\omega|e\rangle\langle e|`$, un saut
$`|g\rangle\langle e|`$ de taux $`\gamma`$) a $`e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}}`$, le spectre
$`\{0,-\gamma,-\gamma/2\pm i\omega\}`$ et le trou spectral $`\gamma/2`$. Nous vérifions le spectre avec `lindblad_superoperator`, le
canal et la distance de trace avec `lindblad_channel`, ainsi que le
[corollaire 8.13](../../lecture_notes/fr/08_dynamics.md#r-8-13) pour des sauts aléatoires à l'intérieur des secteurs de $`M_2\oplus M_3`$ :

```python
q2 = TypeIAlgebra([2], [2], precision='double')                    # basis |g>, |e>
omega, gamma = 1.3, 0.4
H_at = q2.from_blocks([torch.tensor([[0.0, 0.0], [0.0, omega]], dtype=f64)])
sigma_minus = q2.from_blocks([torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=f64)])   # |g><e|
gen = dynamics.lindblad_superoperator(q2, H_at, [sigma_minus], rates=[gamma])[0][0]
spec = sorted(torch.linalg.eigvals(gen).tolist(), key=lambda z: (z.real, z.imag))
print([f"{z.real:+.2f}{z.imag:+.2f}i" for z in spec])
assert np.allclose(spec, [-gamma, -gamma / 2 - 1j * omega, -gamma / 2 + 1j * omega, 0.0], atol=1e-12)

psi_at = torch.tensor([0.6, 0.8], dtype=c128)                      # rho_ee(0) = 0.64, rho_ge(0) = 0.48
rho_at = DensityMatrix(q2, matrix=torch.outer(psi_at, psi_at.conj())[None, None])
ground = DensityMatrix(q2, matrix=torch.tensor([[[[1.0, 0.0], [0.0, 0.0]]]], dtype=c128))
for s in (0.5, 2.0, 10.0):
    ch = dynamics.lindblad_channel(q2, H_at, [sigma_minus], t=s, rates=[gamma])
    UA = channels.unitary_channel(dynamics.propagator(H_at, s)).compose(
        channels.amplitude_damping_channel(q2, 1 - math.exp(-gamma * s)))
    assert (ch.superoperator()[0] - UA.superoperator()[0]).abs().max() < 1e-12
    r = ch(rho_at).matrix[0, 0]
    assert abs(r[1, 1].real - 0.64 * math.exp(-gamma * s)) < 1e-12
    assert abs(r[0, 1] - 0.48 * np.exp((1j * omega - gamma / 2) * s)) < 1e-12
    exact = math.sqrt(0.64 ** 2 * math.exp(-2 * gamma * s) + 0.48 ** 2 * math.exp(-gamma * s))
    assert abs(states.trace_distance(ch(rho_at), ground).item() - exact) < 1e-12

two = TypeIAlgebra([2, 3], [2, 3], precision='double')
H_two = two.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
L_two = [two.operator(torch.randn(1, 2, 3, 3, dtype=c128)) for _ in range(2)]
kernel = [int((torch.linalg.eigvals(G_c[0]).abs() < 1e-10).sum())
          for G_c in dynamics.lindblad_superoperator(two, H_two, L_two)]
print("stationary states per sector:", kernel)
assert all(dim >= 1 for dim in kernel)                             # Corollary 8.13
```

Observation numérique : pour ces sauts aléatoires, le code trouve exactement un état stationnaire par
secteur. D'après le théorème D.6, appliqué secteur par secteur ([notes de cours, section 8.5](../../lecture_notes/fr/08_dynamics.md#s-8-5)), c'est le cas si chaque
secteur a un état stationnaire fidèle et si les seules matrices qui commutent avec $`H_c`$ et avec tous les $`L_{j,c}`$, $`L_{j,c}^*`$
sont les multiples de $`1`$ ; le code vérifie ces deux hypothèses (en virgule flottante, pour les taux $`g_j=1`$ utilisés
ci-dessus), puis teste la conclusion.

```python
for c, k in enumerate(two.k_factors):                              # Theorem D.6 in every sector
    G_c = dynamics.lindblad_superoperator(two, H_two, L_two)[c][0].numpy()   # on row-major vec(rho_c)
    w_c, V_c = np.linalg.eig(G_c)
    sigma_c = V_c[:, np.argmin(abs(w_c))].reshape(k, k)
    sigma_c = sigma_c / np.trace(sigma_c)
    assert np.abs(G_c @ sigma_c.reshape(-1)).max() < 1e-12            # stationary
    assert np.abs(sigma_c - sigma_c.conj().T).max() < 1e-10 and np.linalg.eigvalsh(sigma_c).min() > 1e-3   # faithful
    ops = [H_two.matrix[0, c, :k, :k].numpy()] + [L.matrix[0, c, :k, :k].numpy() for L in L_two]
    ops += [A.conj().T for A in ops[1:]]
    comm = np.concatenate([np.kron(A, np.eye(k)) - np.kron(np.eye(k), A.T) for A in ops])   # X -> [A, X]
    assert int((np.linalg.svd(comm, compute_uv=False) < 1e-9).sum()) == 1   # commutant = C 1
assert kernel == [1, 1]                                            # the conclusion of Theorem D.6
```

<a id="s-8-6"></a>

## 8.6 Points de comparaison

QuTiP (`sesolve`, `mesolve` pour l'équation de Lindblad avec des solveurs d'EDO adaptatifs,
`mcsolve` pour les trajectoires quantiques, `steadystate`, `liouvillian`), dynamiqs (JAX ; solveurs
de Schrödinger et de Lindblad par lots et différentiables, sur GPU), QuantumOptics.jl et SciPy
(`scipy.linalg.expm`, et `scipy.sparse.linalg.expm_multiply` pour $`e^{t\mathcal L}\rho`$ avec un
générateur vectorisé creux) résolvent les mêmes équations sur un seul espace de Hilbert, et doivent
donner les mêmes nombres pour les problèmes ci-dessus. Sur $`\mathcal M`$, on leur passe des matrices
$`D\times D`$ diagonales par blocs ; le générateur vectorisé vit alors sur un espace de dimension
$`D^2`$ (creux, mais avec les secteurs couplés dans une seule résolution), à moins que l'utilisateur
ne découpe le problème en secteurs à la main. La bibliothèque intègre la structure par blocs :
propagateurs exacts et canaux de Lindblad secteur par secteur (superopérateurs $`k_c^2\times k_c^2`$),
sauts entre secteurs sous forme d'`InterSectorChannel` avec la compression de la proposition 8.8,
lots de hamiltoniens et d'opérateurs de saut (pour les moyennes sur le désordre) et différentiation
automatique (ces deux derniers points aussi dans dynamiqs). Elle n'a pas de contrôle adaptatif du
pas (RK4 est à pas fixe et seulement approximativement positif), pas de solveur de trajectoires, et
`lindblad_channel` n'accepte pas les sauts entre secteurs ; pour ceux-ci, QuTiP ou dynamiqs sur
$`M_D`$ font référence. Outre RK4, `lindblad_evolve` a deux intégrateurs en temps dont les pas sont des canaux pour toute taille de
pas (section 8.7) : le schéma de décomposition d'ordre deux `method='kraus'` et `method='exact'`, qui traite aussi les sauts
entre secteurs, au moyen de la matrice de $`\mathcal L`$ sur $`\mathcal M`$ (pour de petits secteurs).

<a id="s-8-7"></a>

## 8.7 Pas de temps qui sont des canaux

*Notes de cours :* [section 8.7](../../lecture_notes/fr/08_dynamics.md#s-8-7) ; [proposition 8.15](../../lecture_notes/fr/08_dynamics.md#r-8-15) (RK4 n'est pas positif), [théorème 8.16](../../lecture_notes/fr/08_dynamics.md#r-8-16) (un schéma de décomposition dont les pas sont des canaux).

`lindblad_evolve` a trois intégrateurs en temps, choisis par `method` : `'rk4'` (la valeur par défaut, inchangée) est le
pas $`R_h=T_4(h\mathcal L)`$ ; `'kraus'` est le schéma de décomposition $`\tilde\Psi_h`$ du théorème 8.16, construit à partir
de $`E_h=e^{-ihH_{\mathrm{eff}}/2}`$, de la partie due aux sauts $`\Phi_{\mathrm{tot}}`$ (sauts à l'intérieur des secteurs et
entre eux) et de la correction $`M_h^{-1/2}`$ ; `'exact'` est $`e^{h\mathcal L}`$, calculé à partir de la matrice de
$`\mathcal L`$ sur $`\bigoplus_cM_{k_c}`$ (complétée par des zéros jusqu'à $`C k_{\max}^2`$ coordonnées ; pour un générateur
constant, une exponentielle par longueur d'intervalle, et `substeps` n'est pas utilisé). Nous vérifions d'abord la
proposition 8.15 sur le système à deux niveaux avec $`\omega=\gamma=1`$ et $`h=1`$ : un pas de RK4 est $`T_4(h\mathcal L)`$
avec la matrice de `lindblad_superoperator`, $`\tau-|a|^2<0`$, et l'état pur
$`\cos\frac\theta2|g\rangle+\sin\frac\theta2|e\rangle`$, $`\theta=0.2746`$, est envoyé sur une matrice ayant une valeur propre
négative ; les deux autres méthodes donnent des états.

```python
T4 = lambda z: 1 + z + z ** 2 / 2 + z ** 3 / 6 + z ** 4 / 24
om, gm, h = 1.0, 1.0, 1.0                                          # omega = gamma: Proposition 8.15(3)
qb = TypeIAlgebra([2], [2], precision='double')                    # basis |g>, |e>
H_qb = qb.from_blocks([torch.tensor([[0.0, 0.0], [0.0, om]], dtype=f64)])
L_qb = [qb.from_blocks([torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=f64)])]   # |g><e|
step = lambda rho, method: dynamics.lindblad_evolve(rho, H_qb, L_qb, [0.0, h], rates=[gm], substeps=1,
                                                    method=method)[-1]
gen = dynamics.lindblad_superoperator(qb, H_qb, L_qb, rates=[gm])[0][0].numpy()
R_h = sum(np.linalg.matrix_power(h * gen, n) / math.factorial(n) for n in range(5))   # T_4(hL)
X = states.random_density_matrix(qb)
assert np.abs(step(X, 'rk4').matrix[0, 0].numpy().reshape(-1) - R_h @ X.matrix[0, 0].numpy().reshape(-1)).max() < 1e-14

tau, a = T4(-h * gm), T4(h * (1j * om - gm / 2))
u = gm ** 2 + 4 * om ** 2
c = [-gm * u * (4 * om ** 2 - 3 * gm ** 2) / 384, u ** 2 * (4 * om ** 2 - 5 * gm ** 2) / 4608, gm * u ** 3 / 9216,
     -u ** 4 / 147456]
assert abs(tau - abs(a) ** 2 - h ** 5 * sum(ci * h ** i for i, ci in enumerate(c))) < 1e-15   # Proposition 8.15(3)
print(f"tau - |a|^2 = {tau - abs(a) ** 2:.4f}")

th = 0.2746
psi_th = torch.tensor([math.cos(th / 2), math.sin(th / 2)], dtype=c128)
pure = DensityMatrix(qb, matrix=torch.outer(psi_th, psi_th.conj())[None, None])
low = {m: torch.linalg.eigvalsh(step(pure, m).matrix[0, 0]).min().item() for m in dynamics.LINDBLAD_METHODS}
print({m: f"{v:.2e}" for m, v in low.items()})
assert low['rk4'] < -8e-5 and low['kraus'] > 0 and low['exact'] > 0
x_star = max(r.real for r in np.roots([1, -4, 12, -24]) if abs(r.imag) < 1e-12)
print(f"x_* = {x_star:.4f}")
```

Observation numérique : ici $`\tau-|a|^2=-0.0091`$ et la plus petite valeur propre après un pas de RK4 vaut
$`-8.5\cdot10^{-5}`$ (les valeurs affichées) ; $`\theta=0.2746`$ est, à quatre chiffres près, l'angle qui la rend minimale. La
proposition 8.15 porte sur un qubit ; dans une algèbre à secteurs, le même défaut apparaît dans les matrices de Choi des
composantes de secteur $`R_{dc}`$, qui sont toutes positives si et seulement si le pas est CP (proposition 6.4). La fonction
`step_map` ci-dessous reconstruit un pas à partir de `lindblad_evolve` sur des matrices hermitiennes (la sortie est rendue
autoadjointe, donc le pas est évalué sur $`E_{ij}+E_{ji}`$ et $`i(E_{ij}-E_{ji})`$ puis étendu par linéarité). Pour la perte
de particules de la section 8.3, nous vérifions le théorème 8.16(1) et la proposition 8.8(1) : les pas de `'kraus'` et de
`'exact'` sont CP et TP pour des tailles de pas de $`0.1`$ à $`1000`$, tandis que le pas de RK4 n'est déjà plus CP pour
$`h=0.1`$.

```python
def step_map(alg, H, jumps, rates, h, method):
    """{(d, c): matrix of the sector component R_dc of one step, on row-major vec}"""
    ks_, out = alg.k_factors, {}
    for c, kc in enumerate(ks_):
        E = torch.zeros(kc * kc, alg.C, alg.k_max, alg.k_max, dtype=c128)
        for m in range(kc * kc):
            E[m, c, m // kc, m % kc] = 1.0
        Y = [dynamics.lindblad_evolve(alg.operator(Z), H, jumps, [0.0, h], rates=rates, substeps=1,
                                      method=method)[-1].matrix for Z in (E + E.mH, 1j * (E - E.mH))]
        RE = 0.5 * (Y[0] - 1j * Y[1])                                  # R(E_ij) by linearity
        for d, kd in enumerate(ks_):
            out[d, c] = RE[:, d, :kd, :kd].reshape(kc * kc, kd * kd).T.numpy()
    return out

def choi_min(alg, S):                                              # smallest eigenvalue of all J(R_dc)
    ks_ = alg.k_factors
    return min(np.linalg.eigvalsh(S[d, c].reshape(ks_[d], ks_[d], ks_[c], ks_[c]).transpose(2, 0, 3, 1)
                                  .reshape(ks_[c] * ks_[d], -1)).min() for d, c in S)

def tp_error(alg, S):                                              # Tr sum_d R_dc(X) = Tr X
    ks_ = alg.k_factors
    return max(np.abs(sum(np.trace(S[d, c].reshape(kd, kd, -1), axis1=0, axis2=1) for d, kd in enumerate(ks_))
                      - np.eye(kc).reshape(-1)).max() for c, kc in enumerate(ks_))

for h in (0.1, 3.0, 1000.0):
    for m in dynamics.LINDBLAD_METHODS:
        S = step_map(fock, H_bh, [loss], [0.5], h, m)
        print(f"h = {h:6}, {m:5}: smallest Choi eigenvalue {choi_min(fock, S):+.1e}, trace error {tp_error(fock, S):.0e}")
        if m != 'rk4':
            assert choi_min(fock, S) > -1e-12 and tp_error(fock, S) < 1e-12        # a channel
assert choi_min(fock, step_map(fock, H_bh, [loss], [0.5], 0.1, 'rk4')) < -5e-5
```

Observation numérique : pour $`h=0.1`$, la plus petite valeur propre d'une matrice de Choi du pas de RK4 vaut
$`-9.7\cdot10^{-5}`$ ; pour $`h=3`$, où certaines valeurs propres $`z`$ de $`h\mathcal L`$ vérifient $`|T_4(z)|>1`$ (RK4 est
instable), elle vaut $`-9.6\cdot10^{2}`$ (valeurs affichées ci-dessus) ; les deux autres méthodes restent des canaux aux
arrondis près. Pour $`h=0.1`$, toutes les valeurs propres vérifient $`|T_4(z)|\le1`$ : la perte de positivité n'est pas une
instabilité. La renormalisation par $`M_h^{-1/2}`$ est calculée à partir des valeurs propres de $`M_h`$, dont le
conditionnement peut croître exponentiellement avec $`h`$ lorsque les taux de décroissance d'un secteur diffèrent ; lorsqu'il
dépasse $`\varepsilon^{-1/4}`$ ($`\varepsilon`$ la précision machine), la bibliothèque effectue le pas comme deux pas de taille
$`h/2`$, récursivement, de sorte que le pas calculé préserve la trace aux arrondis près.

Enfin, la précision. `'exact'` coïncide avec $`e^{t\hat{\mathcal L}}`$ du générateur $`\hat{\mathcal L}`$ sur $`M_6`$ avec les
sauts $`\sqrt\gamma a_1,\sqrt\gamma a_2`$ (proposition 8.8(1) et (3), calculé avec SciPy à partir des matrices `L_loss` de la
section 8.3), et la distance de trace des deux autres méthodes à celui-ci décroît comme $`h^2`$ (théorème 8.16(3)) et comme
$`h^4`$, respectivement :

```python
rho_f = states.random_density_matrix(fock)
ref = dynamics.lindblad_evolve(rho_f, H_bh, [loss], [0.0, 2.0], rates=[0.5], method='exact')[-1]
H6, I6 = embed(H_bh.matrix[0].numpy()), np.eye(6)
G6 = -1j * (np.kron(H6, I6) - np.kron(I6, H6.T))
for L6 in L_loss:
    L6 = math.sqrt(0.5) * L6
    G6 = G6 + np.kron(L6, L6.conj()) - 0.5 * (np.kron(L6.T @ L6, I6) + np.kron(I6, (L6.T @ L6).T))
exact6 = sla.expm(2.0 * G6) @ embed(rho_f.matrix[0].numpy()).reshape(-1)
assert np.abs(exact6 - embed(ref.matrix[0].numpy()).reshape(-1)).max() < 1e-12

for m, subs in (('kraus', (4, 8, 16, 32)), ('rk4', (4, 8, 16, 32))):
    err = [states.trace_distance(dynamics.lindblad_evolve(rho_f, H_bh, [loss], [0.0, 2.0], rates=[0.5], substeps=n,
                                                          method=m)[-1], ref).item() for n in subs]
    ratios = [e0 / e1 for e0, e1 in zip(err[:-1], err[1:])]
    print(m, [f"{e:.1e}" for e in err], "ratios", [f"{r:.1f}" for r in ratios])
    assert all((3.5 < r < 4.5) if m == 'kraus' else (12 < r < 20) for r in ratios)
```

Observation numérique : diviser $`h`$ par deux divise l'erreur de `'kraus'` par environ $`4`$ et celle de RK4 par environ
$`16`$ (les rapports affichés) : le schéma de décomposition achète la positivité pour toute taille de pas au prix d'un ordre
plus bas. Pour un $`H`$ dépendant du temps ou des sauts dépendant du temps, `'kraus'` et `'exact'` figent le générateur au
milieu de chaque pas, de sorte que chaque pas reste un canal d'après le théorème 8.16(1) ; l'ordre de ce schéma à point
milieu figé n'est pas démontré dans les notes de cours.

<a id="x-8"></a>

## Exercices

[Exercices 2 et 4](../../lecture_notes/fr/08_dynamics.md#x-8) des notes de cours : les valeurs propres
$`-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}`$, $`\delta=(\gamma_1-\gamma_2)/4`$, de $`H_{\mathrm{eff}}`$ dans le secteur $`N=1`$ de la
section 8.4, avec un unique bloc de Jordan de taille $`2`$ en $`J=|\delta|`$ ; et le spectre $`\{0,0,-2\gamma\pm i\omega\}`$ du
déphasage pur sur $`M_2`$ ($`L=\sigma_z`$). Vérifier les deux numériquement.

<details>
<summary>Solutions</summary>

```python
lam = np.linalg.eigvals(H_eff[1])
delta = (g1 - g2) / 4
root = np.sqrt(J ** 2 - delta ** 2 + 0j)
assert np.allclose(sorted(lam, key=np.real), sorted(-0.25j * (g1 + g2) + np.array([root, -root]), key=np.real))
H_ep = np.array([[-0.5j * g1, delta], [delta, -0.5j * g2]])         # J = delta: the exceptional point
Y = H_ep + 0.25j * (g1 + g2) * np.eye(2)                           # H_ep minus its only eigenvalue
assert np.abs(Y @ Y).max() < 1e-15 and np.linalg.matrix_rank(Y) == 1   # one Jordan block of size 2
```

```python
deph = dynamics.lindblad_superoperator(q2, H_at, [q2.from_blocks([torch.diag(torch.tensor([1.0, -1.0], dtype=f64))])],
                                       rates=[gamma])[0][0]
ev = torch.linalg.eigvals(deph)
assert int((ev.abs() < 1e-12).sum()) == 2
assert np.allclose(sorted(ev.tolist(), key=lambda z: (z.real, z.imag))[:2],
                   [-2 * gamma - 1j * omega, -2 * gamma + 1j * omega], atol=1e-12)
```

</details>

**Suite :** [Tutoriel 9. Chaînes de spins](09_spin_chains.md).
