[English](../10_quantum_computing.md) · [Español](../es/10_quantum_computing.md) · **Français** · [Deutsch](../de/10_quantum_computing.md) · [中文](../zh/10_quantum_computing.md) · [日本語](../ja/10_quantum_computing.md) · [Русский](../ru/10_quantum_computing.md)

# Tutoriel 10. Calcul quantique

Ce tutoriel accompagne le [chapitre 10](../../lecture_notes/fr/10_quantum_computing.md) des notes de cours, qui décrit le bruit d'un ordinateur
quantique par la $`*`$-algèbre de ses opérateurs de Kraus. Nous calculons des algèbres de points fixes et des sous-systèmes
sans bruit ([théorème 10.2](../../lecture_notes/fr/10_quantum_computing.md#r-10-2)), vérifions le théorème de Knill–Laflamme ([théorème 10.5](../../lecture_notes/fr/10_quantum_computing.md#r-10-5)) et la correction
d'une étiquette de secteur ([proposition 10.6](../../lecture_notes/fr/10_quantum_computing.md#r-10-6)) sur le code de répétition, et simulons un ansatz à mélangeur XY
avec amortissement d'amplitude et post-sélection dans l'algèbre $`\bigoplus_wM_{\binom nw}`$ des secteurs de poids de Hamming
([propositions 10.8](../../lecture_notes/fr/10_quantum_computing.md#r-10-8)–[10.10](../../lecture_notes/fr/10_quantum_computing.md#r-10-10)).

<a id="s-10-1"></a>

## 10.1 Algèbres de bruit, sous-systèmes sans bruit et points fixes

*Notes de cours :* [section 10.1](../../lecture_notes/fr/10_quantum_computing.md#s-10-1) ; [théorème 10.2](../../lecture_notes/fr/10_quantum_computing.md#r-10-2), [corollaire 10.3](../../lecture_notes/fr/10_quantum_computing.md#r-10-3), [proposition 10.4](../../lecture_notes/fr/10_quantum_computing.md#r-10-4).

D'après [notes de cours, corollaire 10.3](../../lecture_notes/fr/10_quantum_computing.md#r-10-3), le déphasage collectif $`\Phi(\rho)=\tfrac12(U\rho U^*+U^*\rho U)`$,
$`U=e^{-i\theta Z_{\mathrm{tot}}}`$, sur $`n`$ qubits a pour algèbre d'interaction $`\mathcal A=\mathrm{span}\{P_0,\dots,P_n\}`$
($`P_w`$ le projecteur sur le sous-espace $`H_w`$ engendré par les vecteurs de base de poids de Hamming $`w`$), et
$`\mathrm{Fix}(\Phi)=\mathcal A'=\bigoplus_wB(H_w)`$ ([théorème 10.2](../../lecture_notes/fr/10_quantum_computing.md#r-10-2)).

Les commutants sont des noyaux de $`T\mapsto(AT-TA)_A`$ (chapitre 2), les points fixes celui de
$`\sum_iE_i\otimes\bar E_i-1`$ (la matrice de $`\Phi-\mathrm{id}`$ sur les matrices vectorisées ligne par
ligne). Pour $`n=4`$, les trois espaces du corollaire sont de dimensions $`5`$, $`70=1+16+36+16+1`$ et
$`70`$ :

```python
import itertools
import math
import torch

torch.manual_seed(0)
dt = torch.complex128
I2, X = torch.eye(2, dtype=dt), torch.tensor([[0, 1], [1, 0]], dtype=dt)
Y, Z = torch.tensor([[0, -1j], [1j, 0]], dtype=dt), torch.tensor([[1, 0], [0, -1]], dtype=dt)

def on(op, j, n):                             # op on qubit j of n (j = 0 is the leftmost factor)
    out = torch.ones(1, 1, dtype=dt)
    for i in range(n):
        out = torch.kron(out, op if i == j else I2)
    return out

def null_space(L, tol=1e-9):                  # orthonormal columns spanning {x : L x = 0}
    w, V = torch.linalg.eigh(L.mH @ L)
    return V[:, w < tol]

def commutant_basis(ops):
    d, ops = ops[0].shape[0], [A.resolve_conj().contiguous() for A in ops]
    eye = torch.eye(d, dtype=dt)
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return [v.reshape(d, d) for v in null_space(L).T]

def fixed_points(kraus):                      # vec(E X E^*) = (E (x) conj E) vec(X)
    d = kraus[0].shape[0]
    S = sum(torch.kron(E.resolve_conj(), E.conj().resolve_conj()) for E in kraus)
    return [v.reshape(d, d) for v in null_space(S - torch.eye(d * d, dtype=dt)).T]

n, theta = 4, 0.3
U = torch.linalg.matrix_exp(-1j * theta * sum(on(Z, j, n) for j in range(n)))
kraus = [U / math.sqrt(2), U.mH / math.sqrt(2)]                  # unital: U U^* = U^* U = 1
comm = commutant_basis(kraus + [E.mH for E in kraus])             # A'
A_int, fix = commutant_basis(comm), fixed_points(kraus)           # A = A'' (Theorem 2.2), Fix(Phi)
print(len(A_int), len(comm), len(fix))
assert (len(A_int), len(comm), len(fix)) == (5, 70, 70)
span = torch.linalg.svdvals(torch.stack([F.reshape(-1) for F in comm + fix]))
assert int((span > 1e-9).sum()) == 70                             # Fix(Phi) = A' (Theorem 10.2)
```

$`\bigoplus_wM_{\binom nw}`$ est l'algèbre d'une `SpinChain` (chapitre 9), dont les secteurs sont les
nombres de spins vers le haut, c'est-à-dire les poids de Hamming. Là, les opérateurs de Kraus sont
centraux, $`U=\sum_we^{-i\theta(n-2w)}P_w`$, et le canal est l'identité : la bibliothèque stocke la
partie sans bruit.

```python
from torch_vn_algebra import SpinChain, TypeIAlgebra, DensityMatrix, channels, states, dynamics

M = SpinChain(n).algebra                                          # (+)_w M_binom(4, w)
print(M.k_factors)
assert sum(k * k for k in M.k_factors) == len(fix)
phases = [torch.exp(torch.tensor(-1j * theta * (n - 2 * w), dtype=dt)) for w in range(n + 1)]
Phi_M = channels.kraus_channel([M.central([z / math.sqrt(2) for z in phases]),
                                M.central([z.conj() / math.sqrt(2) for z in phases])])
rho = states.random_density_matrix(M, batch_size=3)
assert (Phi_M(rho).matrix - rho.matrix).abs().max() < 1e-12      # no state of M is disturbed
```

Les rotations collectives $`e^{-i\theta\,\mathbf n\cdot\mathbf S}`$, $`S_\alpha=\tfrac12\sum_j\sigma_{\alpha,j}`$, appartiennent à la
$`*`$-algèbre $`\mathcal B`$ engendrée par $`1,S_x,S_y,S_z`$. Pour $`n=4`$, la [section 10.1](../../lecture_notes/fr/10_quantum_computing.md#s-10-1) des notes de cours montre que
$`\dim\mathcal B=35`$, $`\dim\mathcal B'=14`$ et $`\dim K=2`$ pour $`K=\ker S_x\cap\ker S_y\cap\ker S_z`$ ; le code trouve les mêmes
nombres et vérifie la [proposition 10.4](../../lecture_notes/fr/10_quantum_computing.md#r-10-4) : un canal dont les opérateurs de Kraus sont des rotations collectives
aléatoires fixe un état à support dans $`K`$, et déplace un état générique.

```python
S = [sum(on(P, j, n) for j in range(n)) / 2 for P in (X, Y, Z)]
comm_S = commutant_basis(S)
print(len(comm_S), len(commutant_basis(comm_S)))
assert len(comm_S) == 14 and len(commutant_basis(comm_S)) == 35
K = null_space(torch.cat(S))                                      # columns: a basis of K
assert K.shape[1] == 2
axes = torch.nn.functional.normalize(torch.randn(3, 3, dtype=torch.float64), dim=1)
rotations = [torch.linalg.matrix_exp(-1j * t * sum(c * Sa for c, Sa in zip(ax, S))) / math.sqrt(3)
             for t, ax in zip(torch.rand(3, dtype=torch.float64) * math.pi, axes)]
apply = lambda kr, r: sum(E @ r @ E.mH for E in kr)
for G in (K @ torch.randn(2, 2, dtype=dt), torch.randn(16, 16, dtype=dt)):  # supported in K / generic
    r = G @ G.mH / (G @ G.mH).trace()
    assert ((apply(rotations, r) - r).abs().max() < 1e-12) == (G.shape[1] == 2)
```

<a id="s-10-2"></a>

## 10.2 Correction d'erreurs : le théorème de Knill–Laflamme

*Notes de cours :* [section 10.2](../../lecture_notes/fr/10_quantum_computing.md#s-10-2) ; [théorème 10.5](../../lecture_notes/fr/10_quantum_computing.md#r-10-5).

L'exemple du bit-flip de la [section 10.2](../../lecture_notes/fr/10_quantum_computing.md#s-10-2) des notes de cours : pour $`\mathcal C=\mathrm{span}\{|000\rangle,|111\rangle\}`$ et
« au plus un bit inversé », la condition de Knill–Laflamme du [théorème 10.5](../../lecture_notes/fr/10_quantum_computing.md#r-10-5) est vérifiée avec
$`\Lambda=\mathrm{diag}(p_0,\dots,p_3)`$, et la récupération $`\mathcal R(\sigma)=\sum_{k=0}^3PX_k\sigma X_kP`$ restaure tout
état du code. Si chaque qubit est inversé indépendamment avec probabilité $`p`$, la condition est en défaut, et la
fidélité de $`|0_L\rangle`$ après récupération vaut $`1-3p^2+2p^3=0.972`$ pour $`p=0.1`$. La fonction auxiliaire `kl_matrix` renvoie
$`\Lambda`$, ou `None` si la condition est en défaut. Dans la bibliothèque, $`M_8`$ est une algèbre à un seul bloc, et
`R @ Phi` est la composée $`\mathcal R\circ\Phi`$ :

```python
M8 = TypeIAlgebra([8], [8], precision='double')
op8 = lambda T: M8.operator(T.reshape(1, 1, 8, 8))
state8 = lambda v: DensityMatrix(M8, matrix=torch.outer(v, v.conj()).reshape(1, 1, 8, 8))
ket = lambda bits: torch.eye(8, dtype=dt)[int(bits, 2)]
def kl_matrix(kraus, P, tol=1e-12):
    """(lambda_ij) with P E_i^* E_j P = lambda_ij P, or None if the Knill-Laflamme condition fails."""
    lam = torch.tensor([[(P @ a.mH @ b @ P).trace() / P.trace() for b in kraus] for a in kraus])
    ok = all((P @ a.mH @ b @ P - lam[i, j] * P).abs().max() < tol
             for i, a in enumerate(kraus) for j, b in enumerate(kraus))
    return lam if ok else None
def independent(op, prob):                    # op on each qubit independently with probability prob
    return [math.sqrt(prob ** sum(s) * (1 - prob) ** (3 - sum(s)))
            * torch.linalg.multi_dot([on(op, j, 3) if s[j] else torch.eye(8, dtype=dt) for j in range(3)])
            for s in itertools.product([0, 1], repeat=3)]

P = torch.outer(ket('000'), ket('000')) + torch.outer(ket('111'), ket('111'))
flips = [torch.eye(8, dtype=dt)] + [on(X, j, 3) for j in range(3)]
E_one = [math.sqrt(pk) * Xk for pk, Xk in zip((0.7, 0.1, 0.1, 0.1), flips)]
assert torch.allclose(kl_matrix(E_one, P), torch.diag(torch.tensor([0.7, 0.1, 0.1, 0.1], dtype=dt)))
Phi_one = channels.kraus_channel([op8(E) for E in E_one])
R = channels.kraus_channel([op8(P @ Xk) for Xk in flips])        # syndrome measurement + correction
assert R.is_trace_preserving()
a, b = torch.nn.functional.normalize(torch.randn(2, dtype=dt), dim=0)
rho_L = state8(a * ket('000') + b * ket('111'))
assert abs((R @ Phi_one)(rho_L).fidelity(rho_L).item() - 1) < 1e-12       # Theorem 10.5
p = 0.1                                                           # independent flips: KL fails
E_ind = independent(X, p)
assert kl_matrix(E_ind, P) is None
zero_L = state8(ket('000'))
F = (R @ channels.kraus_channel([op8(E) for E in E_ind]))(zero_L).fidelity(zero_L).item()
print(f"fidelity after recovery {F:.4f}, exact {1 - 3 * p**2 + 2 * p**3:.4f}")
assert abs(F - (1 - 3 * p**2 + 2 * p**3)) < 1e-12
```

<a id="s-10-3"></a>

## 10.3 Correction d'erreurs par algèbres d'opérateurs

*Notes de cours :* [section 10.3](../../lecture_notes/fr/10_quantum_computing.md#s-10-3) ; [proposition 10.6](../../lecture_notes/fr/10_quantum_computing.md#r-10-6).

L'exemple de la [section 10.3](../../lecture_notes/fr/10_quantum_computing.md#s-10-3) des notes de cours : le code de répétition soumis à des inversions de phase
indépendantes (probabilité $`q`$ par qubit), suivies d'au plus une inversion de bit. Le qubit n'est pas protégé (la
condition de Knill–Laflamme est en défaut), mais le bit classique l'est : la condition 2 de la
[proposition 10.6](../../lecture_notes/fr/10_quantum_computing.md#r-10-6) est vérifiée pour $`P_0=|000\rangle\langle000|`$, $`P_1=|111\rangle\langle111|`$, et le vote majoritaire
$`N_0=\sum_{|b|\le1}|b\rangle\langle b|`$, $`N_1=1-N_0`$ lit le bit. Après la récupération de la section 10.2, la fidélité de
$`|{+_L}\rangle`$ vaut $`(1+(1-2q)^3)/2=0.8645`$ pour $`q=0.05`$.

```python
q = 0.05
Phi2 = Phi_one @ channels.kraus_channel([op8(E) for E in independent(Z, q)])   # dephasing, then a flip
E2 = [Phi2.kraus[0, i, 0] for i in range(Phi2.kraus_rank)]        # its 32 Kraus operators
P0, P1 = torch.outer(ket('000'), ket('000')), torch.outer(ket('111'), ket('111'))
assert kl_matrix(E2, P) is None                                   # the qubit is not protected
assert max((P0 @ e.mH @ f @ P1).abs().max().item() for e in E2 for f in E2) == 0      # condition 2
N0 = torch.diag(torch.tensor([float(bin(i).count('1') <= 1) for i in range(8)], dtype=dt))
for N, Pd in ((N0, P0), (torch.eye(8, dtype=dt) - N0, P1)):       # P Phi^*(N_d) P = P_d
    assert (P @ sum(e.mH @ N @ e for e in E2) @ P - Pd).abs().max() < 1e-12
votes = states.born_probabilities(Phi2(rho_L), [op8(N0), op8(torch.eye(8, dtype=dt) - N0)])[0]
assert torch.allclose(votes, torch.stack([abs(a) ** 2, abs(b) ** 2]))      # the classical bit survives
plus_L = state8((ket('000') + ket('111')) / math.sqrt(2))
F_plus = (R @ Phi2)(plus_L).fidelity(plus_L).item()
print(f"fidelity of |+_L> {F_plus:.4f}, exact {(1 + (1 - 2 * q)**3) / 2:.4f}")
assert abs(F_plus - (1 + (1 - 2 * q) ** 3) / 2) < 1e-12
```

<a id="s-10-4"></a>

## 10.4 Mélangeurs sous contraintes

*Notes de cours :* [section 10.4](../../lecture_notes/fr/10_quantum_computing.md#s-10-4) ; [proposition 10.8](../../lecture_notes/fr/10_quantum_computing.md#r-10-8), [proposition 10.9](../../lecture_notes/fr/10_quantum_computing.md#r-10-9), [proposition 10.10](../../lecture_notes/fr/10_quantum_computing.md#r-10-10).

L'ansatz à opérateurs alternés quantique alterne $`e^{-i\gamma H_C}`$, où $`H_C`$ est un hamiltonien de coût diagonal,
et $`e^{-i\beta H_M}`$, où $`H_M=\tfrac14\sum_{(i,j)}(X_iX_j+Y_iY_j)`$ est le mélangeur XY sur un anneau. D'après
[notes de cours, proposition 10.8](../../lecture_notes/fr/10_quantum_computing.md#r-10-8), tous ces opérateurs appartiennent à $`\mathcal M=\{P_0,\dots,P_n\}'\cong\bigoplus_wM_{\binom nw}`$,
donc l'ansatz ne quitte jamais le secteur de poids de Hamming de son état initial.

Dans la bibliothèque, $`\mathcal M`$ est `SpinChain(n).algebra` (bit $`1`$ = spin vers le haut, donc
$`S^z_i=-\tfrac12Z_i`$), et $`H_M`$ est son hamiltonien XXZ
$`J\sum_{(i,j)}\big(\tfrac12(S^+_iS^-_j+S^-_iS^+_j)+\Delta S^z_iS^z_j\big)`$ avec $`J=1`$, $`\Delta=0`$, puisque
$`S^+_iS^-_j+S^-_iS^+_j=\tfrac12(X_iX_j+Y_iY_j)`$ ; nous le vérifions sur la matrice dense $`64\times64`$.
Pour $`H_C`$, nous prenons la bissection du graphe formé des triangles $`\{0,1,2\}`$, $`\{3,4,5\}`$ reliés
par les arêtes $`(2,3)`$, $`(5,0)`$ : $`H_C=\sum_{(i,j)\in E}\tfrac12(1-Z_iZ_j)`$, $`Z_iZ_j=4S^z_iS^z_j`$, compte
les arêtes coupées. La meilleure bissection coupe $`2`$ arêtes (les deux triangles ; tout autre ensemble de trois sommets
coupe les deux triangles, et un triangle coupé a deux arêtes coupées) ; nous partons de $`|101010\rangle`$, qui en
coupe $`6`$. Enfin, tout l'ansatz est recalculé avec des matrices denses $`64\times64`$, indépendamment de
la bibliothèque.

```python
n, w = 6, 3
chain = SpinChain(n, 'periodic')
M = chain.algebra                                                 # (+)_w M_binom(6, w), padded to 20
H_M = chain.xxz(J=1.0, Delta=0.0)
H_dense = sum(on(X, i, n) @ on(X, j, n) + on(Y, i, n) @ on(Y, j, n) for i, j in chain.bonds) / 4
spec = torch.sort(torch.cat([e[0] for e in H_M.eigenvalues()]))[0]
assert torch.allclose(spec, torch.linalg.eigvalsh(H_dense))       # same spectrum, sector by sector
edges = chain.bonds + [(0, 2), (3, 5)]                            # two triangles joined by (2,3), (5,0)
H_C = M.zero()
for i, j in edges:
    H_C = H_C + (M.identity() - chain.sz(i) @ chain.sz(j) * 4.0) * 0.5
layers = [dynamics.propagator(H_M, beta) @ dynamics.propagator(H_C, gamma)
          for gamma, beta in ((0.6, 1.0), (1.4, 1.4))]
rho0 = ideal = chain.basis_state('101010')
for U_l in layers:
    ideal = channels.unitary_channel(U_l)(ideal)
assert abs(ideal.sector_probabilities()[0, w].item() - 1) < 1e-12           # Proposition 10.8

HC_dense = sum(torch.eye(64, dtype=dt) - on(Z, i, n) @ on(Z, j, n) for i, j in edges) / 2
psi = torch.eye(64, dtype=dt)[int('101010', 2)]                   # the same circuit on 6 full qubits
for gamma, beta in ((0.6, 1.0), (1.4, 1.4)):
    psi = torch.linalg.matrix_exp(-1j * beta * H_dense) @ torch.linalg.matrix_exp(-1j * gamma * HC_dense) @ psi
weight = torch.tensor([bin(b).count('1') for b in range(64)])
assert abs(psi[weight == w].norm().item() - 1) < 1e-12           # Proposition 10.8, dense
assert abs((psi.conj() @ HC_dense @ psi).real.item() - states.expectation(ideal, H_C).real.item()) < 1e-10
```

L'**amortissement d'amplitude** du qubit $`i`$ ($`|1\rangle\to|0\rangle`$ avec probabilité $`g`$) a pour
opérateurs de Kraus $`K_0^{(i)}=|0\rangle\langle0|_i+\sqrt{1-g}\,|1\rangle\langle1|_i`$ et
$`K_1^{(i)}=\sqrt g\,|0\rangle\langle1|_i`$, qui abaisse le poids : c'est une application entre secteurs,
`chain.site_amplitude_damping(i, g)` (un `InterSectorChannel`).

D'après [notes de cours, proposition 10.9](../../lecture_notes/fr/10_quantum_computing.md#r-10-9), si chaque couche du circuit est suivie de l'amortissement
de tous les qubits, le poids $`w`$ de l'état initial est encore obtenu avec probabilité $`(1-g)^{wp}`$ après $`p`$ couches,
et l'état post-sélectionné est l'état sans bruit.

Un `InterSectorChannel` ne conserve que les blocs diagonaux de sa sortie ; ici il n'y en a pas
d'autres (la démonstration de la proposition 10.9 s'applique mot pour mot à l'amortissement d'un seul qubit,
avec $`s\in\{0,1\}`$), donc rien n'est perdu. Avec $`g=0.05`$ et $`p=2`$, le taux d'acceptation vaut $`0.95^6=0.7351`$ :

```python
g = 0.05
damping = [chain.site_amplitude_damping(i, g) for i in range(n)]
noisy = rho0
for U_l in layers:
    noisy = channels.unitary_channel(U_l)(noisy)
    for D in damping:                                             # one qubit after the other
        noisy = D(noisy)
print(noisy.sector_probabilities()[0])                            # weight leaks to w = 2, 1, 0
P_w = M.central([1.0 if v == w else 0.0 for v in range(n + 1)])
post, accept = noisy.condition_on(P_w)                            # measure the weight, keep w = 3
cut = lambda r: states.expectation(r, H_C).real.item()
print(f"P(weight 3) {accept.item():.4f}, fidelity {noisy.fidelity(ideal).item():.4f}, "
      f"post-selected {post.fidelity(ideal).item():.4f}; "
      f"cut: ideal {cut(ideal):.4f}, noisy {cut(noisy):.4f}, post-selected {cut(post):.4f}")
assert abs(accept.item() - (1 - g) ** (w * len(layers))) < 1e-12              # Proposition 10.9
assert abs(post.fidelity(ideal).item() - 1) < 1e-10 and abs(cut(post) - cut(ideal)) < 1e-10
```

Observation numérique (les valeurs affichées) : la partie non admissible abaisse la coupe moyenne de
$`4.3789`$ à $`4.3123`$, une amélioration factice qui tromperait un optimiseur des angles ; après
post-sélection, la valeur est exacte (proposition 10.9).

D'après [notes de cours, proposition 10.10](../../lecture_notes/fr/10_quantum_computing.md#r-10-10), sous l'équation de Lindblad avec amortissement d'amplitude
(taux $`\gamma_1`$) et déphasage $`L_k=S^z_k`$ (taux $`\gamma_2`$), un état partant de $`H_w`$ conserve le poids $`w`$ avec
probabilité $`e^{-\gamma_1wt}`$, et l'état post-sélectionné ne dépend pas de $`\gamma_1`$ ; pour $`\gamma_2=0`$, la
post-sélection est parfaite. L'article
d'accompagnement (`examples/paper/ex9_qaoa_xy_mixer.py`) fait évoluer $`|111000\rangle`$ sous $`H_M`$
pendant $`T=2`$ ; nous reproduisons sa ligne $`\gamma_1=\gamma_2=0.1`$ : $`P(w=3)=e^{-0.6}=0.5488`$
(proposition 10.10) et, comme observation numérique du code ci-dessous, une fidélité de $`0.4805`$ et une
fidélité post-sélectionnée de $`0.8755`$ ; d'après la proposition 10.10, cette dernière est la même que
pour $`\gamma_1=0`$.

```python
T = 2.0
jumps = [chain.lowering(i) for i in range(n)] + [chain.sz(i) for i in range(n)]
rho0 = chain.basis_state('111000')                                # the initial state of the paper
ideal_T = dynamics.von_neumann(rho0, H_M, [T])[0]
post_T = {}
for g1, g2 in ((0.1, 0.0), (0.0, 0.1), (0.1, 0.1)):          # the last one is printed
    rho_T = dynamics.lindblad_evolve(rho0, H_M, jumps, [0.0, T], rates=[g1] * n + [g2] * n, substeps=80)[-1]
    post_T[g1, g2], accept_T = rho_T.condition_on(P_w)
    assert abs(accept_T.item() - math.exp(-w * g1 * T)) < 1e-8                 # Proposition 10.10
print(f"P(weight 3) {accept_T.item():.4f}, fidelity {rho_T.fidelity(ideal_T).item():.4f}, "
      f"post-selected {post_T[0.1, 0.1].fidelity(ideal_T).item():.4f}")
assert abs(post_T[0.1, 0.0].fidelity(ideal_T).item() - 1) < 1e-6               # gamma_2 = 0: perfect
assert post_T[0.1, 0.1].trace_distance(post_T[0.0, 0.1]).item() < 1e-6       # independent of gamma_1
```

<a id="s-10-5"></a>

## 10.5 Points de comparaison

Qiskit (`qiskit.quantum_info`, modèles de bruit d'Aer, `QAOAAnsatz` avec un mélangeur personnalisé),
Cirq, PennyLane (`qml.qaoa.xy_mixer`, `default.mixed`) et QuTiP/`qutip-qip` simulent tout cela sur $`n`$
qubits complets et sont les bons outils pour les circuits ; stim échantillonne des circuits de
Clifford avec un bruit de type Pauli sur des milliers de qubits (l'amortissement d'amplitude n'est
pas un canal de Pauli : les canaux de Pauli sont unitaux, l'amortissement d'amplitude ne l'est pas). Nous ne connaissons dans ces logiciels aucune fonction qui calcule l'algèbre
d'interaction d'un canal, son commutant ou toute son algèbre de points fixes (`steadystate` de QuTiP
renvoie un seul état stationnaire d'un lindbladien) ; pour les décompositions numériques par blocs
des $`*`$-algèbres de matrices, voir Murota et al., *Japan J. Indust. Appl. Math.* 27 (2010). La
bibliothèque ajoute le stockage dans $`\bigoplus_wM_{\binom nw}`$, les applications entre secteurs et la
post-sélection comme projecteur central, par lots. Un secteur requiert $`\binom nw^2`$ coefficients au
lieu de $`4^n`$ ($`1.7\cdot10^8`$ contre $`4.3\cdot10^9`$ pour $`n=16`$, $`w=8`$), mais l'algèbre complète est
complétée par des zéros jusqu'à $`(n+1)\binom n{\lfloor n/2\rfloor}^2`$ coefficients ($`2.8\cdot10^9`$), une
économie modeste. Elle n'est pas un simulateur de circuits et n'a pas de formalisme des
stabilisateurs.

<a id="x-10"></a>

## Exercices

[Exercices 1–3](../../lecture_notes/fr/10_quantum_computing.md#x-10) des notes de cours, numériquement : pour l'amortissement d'amplitude d'un qubit,
$`\mathcal A'`$, $`\mathrm{Fix}(\Phi)`$ et $`\mathrm{Fix}(\Phi^*)`$ sont de dimension un (exercice 1) ; pour le canal du
[corollaire 10.3](../../lecture_notes/fr/10_quantum_computing.md#r-10-3) avec $`\theta=\pi/2`$ et $`n=4`$, $`\dim\mathrm{Fix}(\Phi)=128`$ (exercice 2) ; la récupération
$`\mathcal R`$ du code de bit-flip corrige l'erreur cohérente $`e^{-i\phi X_1}`$ (exercice 3).

<details>
<summary>Solutions</summary>

```python
E_ad = [torch.tensor([[1, 0], [0, math.sqrt(0.7)]], dtype=dt), torch.tensor([[0, math.sqrt(0.3)], [0, 0]], dtype=dt)]
fix_ad = fixed_points(E_ad)
assert len(fix_ad) == 1 and len(commutant_basis(E_ad + [E.mH for E in E_ad])) == 1
assert torch.allclose(fix_ad[0] / fix_ad[0][0, 0], torch.tensor([[1, 0], [0, 0]], dtype=dt))
assert len(fixed_points([E.mH for E in E_ad])) == 1                     # Fix(Phi^*) = C 1
```

```python
U2 = torch.linalg.matrix_exp(-1j * (math.pi / 2) * sum(on(Z, j, 4) for j in range(4)))
assert len(fixed_points([U2 / math.sqrt(2), U2.mH / math.sqrt(2)])) == 128
assert sum(math.comb(4, v) ** 2 for v in range(5)) == math.comb(8, 4) == 70
```

```python
for phi in (0.2, 1.0, 2.5):
    G1 = torch.linalg.matrix_exp(-1j * phi * on(X, 0, 3))
    out = (R @ channels.kraus_channel([op8(G1)]))(rho_L)
    assert (out.matrix - rho_L.matrix).abs().max() < 1e-12       # corrected exactly
```

</details>

**Suite :** [Tutoriel 11](11_optimisation.md) : optimisation sur les groupes unitaires d'une algèbre, et inégalités de traces.
