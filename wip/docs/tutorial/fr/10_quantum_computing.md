[English](../10_quantum_computing.md) · [Español](../es/10_quantum_computing.md) · **Français** · [Deutsch](../de/10_quantum_computing.md) · [中文](../zh/10_quantum_computing.md) · [日本語](../ja/10_quantum_computing.md) · [Русский](../ru/10_quantum_computing.md)

# 10. Calcul quantique

Le bruit sur un ordinateur quantique est un canal, et un canal engendre une algèbre : la
$*$-algèbre de ses opérateurs de Kraus. Sa structure (chapitre 2) indique quelle information survit
intacte au bruit ; une condition sur les produits d'opérateurs de Kraus, faisant intervenir une
autre algèbre, indique quelle information une récupération peut restaurer. La dernière partie
traite des algorithmes qui conservent une charge : leurs états vivent dans l'algèbre des secteurs
de charge, et le bruit dominant les fait passer d'un secteur à l'autre.

**Vous apprendrez** pourquoi un canal unifère et son dual ont les mêmes points fixes (le commutant
du bruit), le théorème de Knill–Laflamme avec une démonstration complète, comment ne protéger que
l'étiquette d'un secteur, et à quoi ressemblent le mélangeur XY de l'ansatz à opérateurs alternés
quantique, l'amortissement d'amplitude et la post-sélection dans l'algèbre $\bigoplus_wM_{\binom nw}$
des secteurs de poids de Hamming.

Dans tout ce chapitre, un **canal** sur $M_D(\mathbb C)$ est une application $\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*$
dont les **opérateurs de Kraus** vérifient $\sum_iE_i^*E_i=1$ (préservation de la trace) ; il est
**unifère** si de plus $\sum_iE_iE_i^*=1$. Son **dual** $\Phi^*(T)=\sum_iE_i^*TE_i$ vérifie
$\operatorname{Tr}(\Phi(\rho)T)=\operatorname{Tr}(\rho\,\Phi^*(T))$ (section 6.2). (D'après le théorème
6.2, toute application complètement positive est de cette forme ; ici les canaux sont *donnés* par
des opérateurs de Kraus.) L'**algèbre d'interaction** $\mathcal A$ de $\Phi$ est la $*$-algèbre
engendrée par $1$ et les $E_i$, l'espace engendré par tous les produits des $E_i$ et des $E_i^*$ ;
une matrice commute avec $\mathcal A$ si et seulement si elle commute avec tous les $E_i$ et $E_i^*$.
Un état $\rho$ est **à support** dans un sous-espace $K$ de projecteur $P_K$ si $\rho=P_K\rho P_K$. Sur
$n$ qubits, $D=2^n$ ; $X,Y,Z$ sont les matrices de Pauli, $Z=\operatorname{diag}(1,-1)$, $X_j$ est $X$ sur le qubit $j$ (numérotés $1,\dots,n$ ;
la section 10.4 suit la bibliothèque, sites $0,\dots,n-1$), $|b\rangle$, $b\in\{0,1\}^n$, la base
de calcul, $|b|$ le **poids de Hamming**, et $P_w$ le projecteur sur $H_w=\operatorname{span}\{|b\rangle:|b|=w\}$.

## 10.1 Algèbres de bruit, sous-systèmes sans bruit et points fixes

D'après le théorème 2.6, il existe un unitaire $V$ tel que $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$
et $V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$ : chaque opérateur de Kraus n'agit que sur le
premier facteur tensoriel de chaque terme de la somme.

**Proposition 10.1 (sous-systèmes sans bruit ; sous-espaces sans décohérence).** (a) Écrivons
$V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}$, et soit $\rho=V(\sigma\otimes\tau)V^*$, placé dans le terme
$c$, pour des états $\sigma$ sur $\mathbb C^{k_c}$ et $\tau$ sur $\mathbb C^{m_c}$. Alors
$\Phi(\rho)=V(\Phi_c(\sigma)\otimes\tau)V^*$ avec le canal $\Phi_c(\sigma)=\sum_iE_{i,c}\sigma E_{i,c}^*$ sur
$M_{k_c}$ : l'état $\tau$ n'est pas touché. De plus, $\Phi^*(T)=T$ pour tout $T\in\mathcal A'$.
(b) Si $K$ est un sous-espace tel que $E_iP_K=e_iP_K$, $e_i\in\mathbb C$, pour tout $i$, alors
$\Phi(\rho)=\rho$ pour tout état $\rho$ à support dans $K$.

*Démonstration.* (a) $\Phi(\rho)=V\sum_i(E_{i,c}\sigma E_{i,c}^*\otimes\tau)V^*$, et $\Phi_c$ préserve la
trace car le terme $c$ de $V^*(\sum_iE_i^*E_i)V=1$ est $(\sum_iE_{i,c}^*E_{i,c})\otimes1_{m_c}$. Pour
$T\in\mathcal A'$, $\Phi^*(T)=T\sum_iE_i^*E_i=T$. (b) $P_KE_i^*=(E_iP_K)^*=\bar e_iP_K$, donc
$E_i\rho E_i^*=E_iP_K\rho P_KE_i^*=|e_i|^2\rho$, et $\sum_i|e_i|^2P_K=P_K\sum_iE_i^*E_iP_K=P_K$. $\square$

Un terme avec $k_c=1$ relève du cas (b) : c'est un **sous-espace sans décohérence** (Zanardi–Rasetti
1997, Lidar–Chuang–Whaley 1998). Pour $k_c>1$, le facteur $\mathbb C^{m_c}$ est un **sous-système
sans bruit** (Knill–Laflamme–Viola 2000) : le bruit agit, mais seulement sur l'autre facteur. La
proposition 10.1 vaut pour tout canal ; pour les canaux unifères, les observables préservées sont
*exactement* celles de $\mathcal A'$ :

**Théorème 10.2 (points fixes des canaux unifères).** Soit $\Phi$ unifère. Alors
$$
\operatorname{Fix}(\Phi)=\{T\in M_D(\mathbb C):\Phi(T)=T\}=\mathcal A'=\operatorname{Fix}(\Phi^*).
$$
En particulier, $\operatorname{Fix}(\Phi)$ est une algèbre de von Neumann, $*$-isomorphe à $\bigoplus_cM_{m_c}$.

*Démonstration.* La première égalité est le théorème 6.10 pour l'algèbre $M_D(\mathbb C)$ (un seul
bloc) : $T$ est un point fixe si et seulement s'il commute avec tous les $E_i$ et $E_i^*$,
c'est-à-dire si et seulement si $T\in\mathcal A'$. Le dual $\Phi^*$ est encore un canal unifère : ses
opérateurs de Kraus $E_i^*$ vérifient $\sum_iE_iE_i^*=1$ (préservation de la trace) et
$\sum_iE_i^*E_i=1$ (unifère), et ils engendrent la même $*$-algèbre $\mathcal A$. Le théorème 6.10
appliqué à $\Phi^*$ donne donc $\operatorname{Fix}(\Phi^*)=\mathcal A'$. La dernière assertion est le
théorème 2.6. $\square$

Pour tout canal, $\mathcal A'\subseteq\operatorname{Fix}(\Phi^*)$ (proposition 10.1(a)), mais sans
l'hypothèse unifère le théorème est en défaut : pour l'amortissement d'amplitude,
$\operatorname{Fix}(\Phi)=\mathbb C|0\rangle\langle0|$, alors que $\mathcal A'=\operatorname{Fix}(\Phi^*)=\mathbb C1$
(section 6.4 et exercice 1). (Kribs 2003, Arias–Gheondea–Gudder 2002 ; en dimension infinie, il faut
davantage d'hypothèses.)

**Corollaire 10.3 (déphasage collectif).** Soient $U=e^{-i\theta Z_{\mathrm{tot}}}$,
$Z_{\mathrm{tot}}=\sum_jZ_j$, avec $e^{2i\theta d}\neq1$ pour $d=1,\dots,n$, et
$\Phi(\rho)=\tfrac12(U\rho U^*+U^*\rho U)$. Alors $\mathcal A=\operatorname{span}\{P_0,\dots,P_n\}$, chaque
$H_w$ est un sous-espace sans décohérence, et $\operatorname{Fix}(\Phi)=\bigoplus_wB(H_w)\cong\bigoplus_wM_{\binom nw}$,
de dimension $\sum_w\binom nw^2$.

*Démonstration.* $Z_{\mathrm{tot}}|b\rangle=(n-2|b|)|b\rangle$, donc $U=\sum_we^{-i\theta(n-2w)}P_w$, de
valeurs propres distinctes (leurs quotients sont les $e^{2i\theta(w-w')}$, $0<|w-w'|\le n$). Les
polynômes de Lagrange de la démonstration du théorème 1.5(a) donnent $P_w=p_w(U)\in\mathcal A$ ;
réciproquement, $U,U^*$ appartiennent à l'algèbre $\operatorname{span}\{P_w\}$. Donc
$\mathcal A=\operatorname{span}\{P_w\}$ ($k_c=1$, $m_c=\binom nw$), et la proposition 10.1(b) s'applique à
chaque $H_w$. Si $T$ commute avec tous les $P_w$, alors $P_vTP_{v'}=TP_vP_{v'}=0$ pour $v\neq v'$,
c'est-à-dire que $T$ est diagonal par blocs, et réciproquement ; le théorème 10.2 donne le reste.
$\square$

Les commutants sont des noyaux de $T\mapsto(AT-TA)_A$ (chapitre 2), les points fixes celui de
$\sum_iE_i\otimes\bar E_i-1$ (la matrice de $\Phi-\mathrm{id}$ sur les matrices vectorisées ligne par
ligne). Pour $n=4$, les trois espaces du corollaire sont de dimensions $5$, $70=1+16+36+16+1$ et
$70$ :

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

$\bigoplus_wM_{\binom nw}$ est l'algèbre d'une `SpinChain` (chapitre 9), dont les secteurs sont les
nombres de spins vers le haut, c'est-à-dire les poids de Hamming. Là, les opérateurs de Kraus sont
centraux, $U=\sum_we^{-i\theta(n-2w)}P_w$, et le canal est l'identité : la bibliothèque stocke la
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

> **Interprétation physique.** Un champ fluctuant, le même sur tous les qubits, décale la phase de
> $|b\rangle$ proportionnellement à $n-2|b|$. Les états d'un même poids de Hamming acquièrent une
> phase commune, de sorte que $|01\rangle,|10\rangle$ codent un qubit insensible à ce bruit (ions
> piégés : Kielpinski et al., *Science* 291, 2001), tandis que les cohérences entre poids
> s'amortissent : le bruit impose une règle de supersélection (section 3.3).

**Décohérence collective.** Si le champ fluctue aussi en direction, les opérateurs de Kraus sont
des rotations $e^{-i\theta\,\mathbf n\cdot\mathbf S}$, $S_\alpha=\tfrac12\sum_j\sigma_{\alpha,j}$, $\sigma\in\{X,Y,Z\}$,
qui appartiennent à la $*$-algèbre $\mathcal B$ engendrée par $1,S_x,S_y,S_z$ (théorème 1.5(a)).

**Proposition 10.4 (sous-espace singulet).** Soit $K=\ker S_x\cap\ker S_y\cap\ker S_z$. Tout canal
dont les opérateurs de Kraus appartiennent à $\mathcal B$ fixe tout état à support dans $K$.

*Démonstration.* $\mathcal B$ est engendrée linéairement par les mots en $S_x,S_y,S_z$ (et $1$) ; les
mots de longueur $\ge1$ annulent $K$, donc chaque élément de $\mathcal B$ agit sur $K$ comme un
scalaire. On applique la proposition 10.1(b). $\square$

Pour $n=4$, le code trouve $\dim K=2$, $\dim\mathcal B'=14$, $\dim\mathcal B=35$, comme pour
$\mathcal B\cong M_5\oplus M_3\oplus M_1$ avec les multiplicités $1,3,2$ (spin total $2,1,0$), ce que
le chapitre 4 (en préparation) établit : quatre qubits portent un qubit sans décohérence ($k=1$,
$m=2$) et un qutrit sans bruit ($k=3$, $m=3$).

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

## 10.2 Correction d'erreurs : le théorème de Knill–Laflamme

Un **code** est un sous-espace $\mathcal C\subseteq\mathbb C^D$ de projecteur $P$. Il **corrige**
$\Phi$ s'il existe un canal $\mathcal R$ (la **récupération**) tel que $\mathcal R(\Phi(\rho))=\rho$ pour
tout état $\rho$ à support dans $\mathcal C$.

**Théorème 10.5 (Knill–Laflamme).** $\mathcal C$ corrige $\Phi$ si et seulement s'il existe des
nombres $\lambda_{ij}$ tels que
$$
PE_i^*E_jP=\lambda_{ij}P\qquad\text{for all }i,j.
$$

*Démonstration.* *Condition nécessaire.* Soient $R_m$ les opérateurs de Kraus de $\mathcal R$ ;
alors $\mathcal R\circ\Phi$ a pour opérateurs de Kraus les $R_mE_i$. Pour un vecteur unitaire
$\psi\in\mathcal C$, $\sum_{m,i}R_mE_i|\psi\rangle\langle\psi|E_i^*R_m^*=|\psi\rangle\langle\psi|$, et pour
$\phi\perp\psi$ la valeur moyenne en $\phi$ donne $\sum_{m,i}|\langle\phi,R_mE_i\psi\rangle|^2=0$, donc
$R_mE_i\psi\in\mathbb C\psi$. Ainsi $R_mE_iP=c_{mi}P$ : une application linéaire $T$ telle que
$Tx\in\mathbb Cx$ pour tout $x$ est un multiple de $1$ (si $Tx=ax$, $Ty=by$, $T(x+y)=c(x+y)$, avec $x,y$
indépendants, alors $a=c=b$). Comme $\sum_mR_m^*R_m=1$,
$PE_i^*E_jP=\sum_m(R_mE_iP)^*(R_mE_jP)=\sum_m\bar c_{mi}c_{mj}P$.

*Condition suffisante.* $\Lambda=(\lambda_{ij})$ est hermitienne, car
$\lambda_{ji}P=(PE_i^*E_jP)^*=\bar\lambda_{ij}P$, et positive, car
$\sum_{ij}\bar c_ic_j\lambda_{ij}P=PF^*FP\ge0$ pour $F=\sum_jc_jE_j$. Écrivons $\Lambda=UDU^*$, avec
$U$ unitaire, $D=\operatorname{diag}(d_k)$, $d_k\ge0$ (théorème 1.3), et $F_k=\sum_iu_{ik}E_i$. Alors
$\sum_kF_k\rho F_k^*=\sum_{ij}(UU^*)_{ij}E_i\rho E_j^*=\Phi(\rho)$ et $PF_k^*F_lP=(U^*\Lambda U)_{kl}P=d_k\delta_{kl}P$.
Si $d_k=0$, alors $\|F_kPx\|^2=0$ pour tout $x$, donc $F_kP=0$. Si $d_k>0$, posons
$V_k=F_kP/\sqrt{d_k}$ ; alors $V_k^*V_l=\delta_{kl}P$, donc $V_kP=V_k$, les $Q_k=V_kV_k^*$ sont des
projecteurs deux à deux orthogonaux et $Q=\sum_kQ_k$ est un projecteur. La récupération
$$
\mathcal R(\sigma)=\sum_{k:d_k>0}V_k^*\sigma V_k+(1-Q)\sigma(1-Q)
$$
est un canal, puisque $\sum_kV_kV_k^*+(1-Q)^2=1$. Pour $\rho=P\rho P$, $\Phi(\rho)=\sum_{d_k>0}d_kV_k\rho V_k^*$ ;
comme $V_l^*V_k=\delta_{kl}P$ et $(1-Q)V_k=(1-Q)Q_kV_k=0$, $\mathcal R(\Phi(\rho))=\sum_kd_k\rho$, et
$\sum_kd_k=\operatorname{Tr}\Lambda=1$ car $\sum_i\lambda_{ii}P=P\sum_iE_i^*E_iP=P$. $\square$

(Knill–Laflamme 1997 ; Nielsen–Chuang, théorème 10.1.) $\mathcal R$ mesure le **syndrome** — lequel
des sous-espaces orthogonaux $Q_k\mathbb C^D$ contient l'état — et ramène l'état dans le code par une
rotation.

**Exemple : le code de bit-flip.** $\mathcal C=\operatorname{span}\{|000\rangle,|111\rangle\}$ et « au plus
un bit inversé », $E_0=\sqrt{p_0}\,1$, $E_j=\sqrt{p_j}\,X_j$. Pour $i\neq j$, $X_iX_j$ (avec $X_0=1$)
inverse un ou deux bits et envoie $\mathcal C$ dans $\mathcal C^\perp$ ; donc
$\Lambda=\operatorname{diag}(p_0,\dots,p_3)$, $V_k=X_kP$, $Q=1$ et $\mathcal R(\sigma)=\sum_{k=0}^3PX_k\sigma X_kP$.
Si chaque qubit est inversé indépendamment avec probabilité $p$, la condition est en défaut, et
$\mathcal R$ envoie $|0_L\rangle=|000\rangle$ sur $|1_L\rangle$ dès que $\ge2$ inversions se produisent et le
laisse en place sinon : la fidélité de $|0_L\rangle$ après récupération est la probabilité d'au plus
une inversion, $(1-p)^3+3p(1-p)^2=1-3p^2+2p^3=0.972$ pour $p=0.1$. Dans la bibliothèque, $M_8$ est une
algèbre à un seul bloc, et `R @ Phi` est la composée $\mathcal R\circ\Phi$ :

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

## 10.3 Correction d'erreurs par algèbres d'opérateurs

Le point de vue de Heisenberg permet de protéger moins que le code entier. Soit $\mathfrak A$ une
algèbre de von Neumann sur $\mathcal C=P\mathbb C^D$ (une $*$-algèbre d'opérateurs sur $\mathcal C$
contenant $P$), les observables à protéger. $\mathfrak A$ est **corrigeable** pour $\Phi$ s'il existe
un canal $\mathcal R$ tel que
$$
P\,(\mathcal R\circ\Phi)^*(A)\,P=A\qquad\text{for all }A\in\mathfrak A,
$$
c'est-à-dire $\operatorname{Tr}(\mathcal R(\Phi(\rho))A)=\operatorname{Tr}(\rho A)$ pour tous les états
$\rho$ à support dans $\mathcal C$ et tout $A\in\mathfrak A$. Pour $\mathfrak A=B(\mathcal C)$, c'est la
correction de la section 10.2 : si $\sigma=\mathcal R(\Phi(\rho))$ vérifie $P\sigma P=\rho$, alors
$\operatorname{Tr}(1-P)\sigma(1-P)=0$, donc $\sigma^{1/2}(1-P)=0$ et $\sigma=\rho$. L'extrême opposé est une
algèbre commutative : seule une étiquette classique est protégée.

**Proposition 10.6 (correction de l'étiquette d'un secteur).** Soient $P=\sum_{c=1}^CP_c$, avec des
projecteurs $P_c\neq0$ deux à deux orthogonaux, et $\mathfrak A=\operatorname{span}\{P_1,\dots,P_C\}$. Les
assertions suivantes sont équivalentes :
1. il existe $N_1,\dots,N_C\ge0$, $\sum_cN_c=1$, tels que $\operatorname{Tr}(N_c\Phi(\rho))=1$ pour tout $c$
   et tout état $\rho$ à support dans $P_c\mathbb C^D$ ;
2. $P_cE_i^*E_jP_d=0$ pour tous $i,j$ et $c\neq d$, c'est-à-dire $PE_i^*E_jP\in\mathfrak A'$, le
   commutant de $\mathfrak A$ dans $B(\mathcal C)$ ;
3. $\mathfrak A$ est corrigeable pour $\Phi$.

*Démonstration.* $T\in B(\mathcal C)$ commute avec tous les $P_c$ si et seulement si $P_cTP_d=0$ pour
$c\neq d$ (comme dans le corollaire 10.3).

(2)$\Rightarrow$(1). Les sous-espaces $S_c=\operatorname{span}\{E_iP_cx\}$ sont deux à deux orthogonaux,
car $\langle E_iP_cx,E_jP_dy\rangle=\langle x,P_cE_i^*E_jP_dy\rangle=0$. Soient $N_c$ le projecteur sur
$S_c$ pour $c<C$ et $N_C=1-\sum_{c<C}N_c$, qui est aussi un projecteur. Pour $\rho$ à support dans
$P_c\mathbb C^D$, $\Phi(\rho)=\sum_iE_iP_c\rho P_cE_i^*$ a son image dans $S_c$, et $N_c$ est l'identité
sur $S_c$ (pour $c=C$ parce que $S_C\perp S_{c'}$, $c'<C$), donc
$\operatorname{Tr}(N_c\Phi(\rho))=\operatorname{Tr}\Phi(\rho)=1$.
(1)$\Rightarrow$(2). Pour un vecteur unitaire $x\in P_c\mathbb C^D$ et $d\neq c$, les nombres
$\operatorname{Tr}(N_{d'}\Phi(|x\rangle\langle x|))\ge0$ ont pour somme $1$ et celui pour $d'=c$ vaut $1$,
donc $0=\operatorname{Tr}(N_d\Phi(|x\rangle\langle x|))=\sum_i\|N_d^{1/2}E_ix\|^2$ ; ainsi $N_dE_iP_c=0$
et, en sommant sur $d$, $N_cE_iP_c=E_iP_c$. Par conséquent $P_cE_i^*E_jP_d=P_cE_i^*N_cE_jP_d=0$ pour
$c\neq d$.
(1),(2)$\Rightarrow$(3). Soient $\omega_c=P_c/\operatorname{Tr}P_c$ et $\mathcal R(\sigma)=\sum_c\operatorname{Tr}(N_c\sigma)\,\omega_c$
(« mesurer, puis préparer »), un canal ayant pour opérateurs de Kraus les
$\sqrt\mu\,|u\rangle\langle f|N_c^{1/2}$ pour les couples propres $(\mu,u)$ de $\omega_c$ et les vecteurs
$f$ d'une base orthonormée. Comme $\operatorname{Tr}(\mathcal R(\sigma)P_d)=\operatorname{Tr}(N_d\sigma)$,
$\mathcal R^*(P_d)=N_d$. En utilisant $N_dE_iP_c=\delta_{cd}E_iP_c$ (étape précédente), puis 2,
$$
P(\mathcal R\circ\Phi)^*(P_d)P=P\,\Phi^*(N_d)\,P=\sum_i\sum_{c,c'}P_{c'}E_i^*N_dE_iP_c=\sum_i\sum_{c'}P_{c'}E_i^*E_iP_d=P_d\Big(\sum_iE_i^*E_i\Big)P_d=P_d .
$$
(3)$\Rightarrow$(1). Soient $N_c=\mathcal R^*(P_c)$ pour $c<C$ et $N_C=\mathcal R^*(1-\sum_{c<C}P_c)$ ; ils
sont positifs et de somme $\mathcal R^*(1)=1$. Pour $\rho$ à support dans $P_c\mathbb C^D$ et $d<C$,
$\operatorname{Tr}(N_d\Phi(\rho))=\operatorname{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(P_d)P)=\operatorname{Tr}(\rho P_d)=\delta_{cd}$,
donc aussi $\operatorname{Tr}(N_C\Phi(\rho))=\delta_{cC}$. $\square$

**Théorème 10.7 (Bény–Kempf–Kribs 2007 ; énoncé sans démonstration).** Une algèbre de von Neumann
$\mathfrak A$ sur $\mathcal C$ est corrigeable pour $\Phi$ si et seulement si $PE_i^*E_jP\in\mathfrak A'$
pour tous $i,j$, où $\mathfrak A'$ est le commutant de $\mathfrak A$ dans $B(\mathcal C)$.

Pour $\mathfrak A=B(\mathcal C)$, $\mathfrak A'=\mathbb CP$ et l'on retrouve le théorème 10.5 ; la
proposition 10.6 est le cas d'une $\mathfrak A$ commutative. Pour
$\mathfrak A=M_k\otimes1_m$ sur $\mathcal C=\mathbb C^k\otimes\mathbb C^m$, on obtient les codes
d'opérateurs (codes de sous-systèmes) de Kribs–Laflamme–Poulin (*Phys. Rev. Lett.* 94, 180501, 2005),
qui contiennent les sous-systèmes sans bruit de la section 10.1 avec $\mathcal R=\mathrm{id}$. Pour la
démonstration, voir C. Bény, A. Kempf, D. W. Kribs, *Phys. Rev. Lett.* 98, 100502 (2007) et
*Phys. Rev. A* 76, 042303 (2007).

**Exemple.** Le code de répétition soumis à des inversions de phase indépendantes (probabilité $q$
par qubit), suivies d'au plus une inversion de bit. Les opérateurs de Kraus sont des multiples des
$X_kZ^s$, $s\in\{0,1\}^3$, et chaque $E_i^*E_j$ inverse au plus deux bits, donc la condition 2 est
vérifiée pour $P_0=|000\rangle\langle000|$, $P_1=|111\rangle\langle111|$ : le bit classique survit, lu par
vote majoritaire, $N_0=\sum_{|b|\le1}|b\rangle\langle b|$, $N_1=1-N_0$. Le qubit, lui, ne survit pas :
$PZ_1P$ n'est pas un multiple de $P$, et après la récupération $\mathcal R$ de la section 10.2 l'état
$|{+_L}\rangle=(|000\rangle+|111\rangle)/\sqrt2$ subit l'inversion de phase logique $Z_1Z_2Z_3$ exactement
lorsqu'un nombre impair d'inversions de phase s'est produit ; sa fidélité vaut
$(1+(1-2q)^3)/2=0.8645$ pour $q=0.05$.

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

> **Interprétation physique.** Protéger une algèbre commutative, c'est protéger de l'information
> classique (un bit, l'étiquette d'un secteur de supersélection), ce qui est bien moins coûteux que
> de protéger un qubit : un bruit qui commute avec les projecteurs des secteurs est inoffensif. Les
> « codes hybrides » correspondent à $\bigoplus_cM_{k_c}$.

## 10.4 Mélangeurs sous contraintes

L'ansatz à opérateurs alternés quantique (Hadfield et al., *Algorithms* 12, 34, 2019) prépare
$|\beta,\gamma\rangle=e^{-i\beta_pH_M}e^{-i\gamma_pH_C}\cdots e^{-i\beta_1H_M}e^{-i\gamma_1H_C}|b_0\rangle$, où $H_C$ est
un hamiltonien de coût diagonal et $H_M$ un **mélangeur**. Pour une contrainte $|b|=w$ (bissection de
graphe, portefeuilles de $w$ actifs), on prend $|b_0|=w$ et le **mélangeur XY** sur un anneau,
$H_M=\tfrac14\sum_{(i,j)}(X_iX_j+Y_iY_j)$.

**Proposition 10.8.** $H_M$, $H_C$, $e^{-i\beta H_M}$ et $e^{-i\gamma H_C}$ appartiennent à
$\mathcal M=\{P_0,\dots,P_n\}'=\bigoplus_wB(H_w)\cong\bigoplus_wM_{\binom nw}$, donc $|\beta,\gamma\rangle\in H_w$ :
l'ansatz ne quitte jamais le sous-espace admissible.

*Démonstration.* Sur les qubits $i,j$, $X\otimes X+Y\otimes Y=2(|01\rangle\langle10|+|10\rangle\langle01|)$
(utiliser $Y|0\rangle=i|1\rangle$, $Y|1\rangle=-i|0\rangle$) : $X_iX_j+Y_iY_j$ échange les bits $i,j$ lorsqu'ils
diffèrent et annule $|b\rangle$ sinon, donc il envoie chaque $H_w$ dans lui-même. Un $H_C$ diagonal
commute avec les $P_w$ diagonaux, et les exponentielles aussi (théorème 1.5(d)). $\square$

Dans la bibliothèque, $\mathcal M$ est `SpinChain(n).algebra` (bit $1$ = spin vers le haut, donc
$S^z_i=-\tfrac12Z_i$), et $H_M$ est son hamiltonien XXZ
$J\sum_{(i,j)}\big(\tfrac12(S^+_iS^-_j+S^-_iS^+_j)+\Delta S^z_iS^z_j\big)$ avec $J=1$, $\Delta=0$, puisque
$S^+_iS^-_j+S^-_iS^+_j=\tfrac12(X_iX_j+Y_iY_j)$ ; nous le vérifions sur la matrice dense $64\times64$.
Pour $H_C$, nous prenons la bissection du graphe formé des triangles $\{0,1,2\}$, $\{3,4,5\}$ reliés
par les arêtes $(2,3)$, $(5,0)$ : $H_C=\sum_{(i,j)\in E}\tfrac12(1-Z_iZ_j)$, $Z_iZ_j=4S^z_iS^z_j$, compte
les arêtes coupées. La meilleure bissection coupe $2$ arêtes ; nous partons de $|101010\rangle$, qui en
coupe $6$. Enfin, tout l'ansatz est recalculé avec des matrices denses $64\times64$, indépendamment de
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

L'**amortissement d'amplitude** du qubit $i$ ($|1\rangle\to|0\rangle$ avec probabilité $g$) a pour
opérateurs de Kraus $K_0^{(i)}=|0\rangle\langle0|_i+\sqrt{1-g}\,|1\rangle\langle1|_i$ et
$K_1^{(i)}=\sqrt g\,|0\rangle\langle1|_i$, qui abaisse le poids : c'est une application entre secteurs,
`chain.site_amplitude_damping(i, g)` (un `InterSectorChannel`).

**Proposition 10.9 (fuite et post-sélection).** Soient $\rho_0$ à support dans $H_w$, $U_l$ des
unitaires de $\mathcal M$, $\mathcal D$ l'amortissement de tous les qubits, et
$\rho_l=\mathcal D(U_l\rho_{l-1}U_l^*)$. Alors chaque $\rho_l$ est diagonal par blocs, de poids $\le w$,
et
$$
P_w\rho_pP_w=(1-g)^{wp}\;U\rho_0U^*,\qquad U=U_p\cdots U_1 :
$$
le poids $w$ est obtenu avec probabilité $(1-g)^{wp}$, et l'état post-sélectionné est l'état sans
bruit.

*Démonstration.* Les opérateurs de Kraus de $\mathcal D$ sont les $K_s=\prod_iK^{(i)}_{s_i}$,
$s\in\{0,1\}^n$ ; $K_s$ envoie $H_v$ dans $H_{v-|s|}$, c'est-à-dire $P_vK_s=K_sP_{v+|s|}$ ($P_v=0$ pour
$v>n$). Si $\sigma$ est diagonal par blocs de poids $\le w$, alors
$P_vK_s\sigma K_s^*P_{v'}=K_sP_{v+|s|}\sigma P_{v'+|s|}K_s^*$ s'annule pour $v\neq v'$ et pour $v>w-|s|$ :
$\mathcal D(\sigma)$ est diagonal par blocs de poids $\le w$, et son bloc de poids $w$ provient du seul
$s=0$. Comme $K_0|b\rangle=(1-g)^{|b|/2}|b\rangle$, $P_w\mathcal D(\sigma)P_w=(1-g)^wP_w\sigma P_w$. Les $U_l$
commutent avec les $P_v$ ; on conclut par récurrence sur $l$. $\square$

Un `InterSectorChannel` ne conserve que les blocs diagonaux de sa sortie ; ici il n'y en a pas
d'autres, donc rien n'est perdu. Avec $g=0.05$ et $p=2$, le taux d'acceptation vaut $0.95^6=0.7351$ :

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

La partie non admissible abaisse la coupe moyenne de $4.3789$ à $4.3123$, une amélioration factice
qui tromperait un optimiseur des angles ; après post-sélection, la valeur est exacte.

**Temps continu.** Il en va de même pour l'équation de Lindblad (théorème 8.6, [chapitre 8](08_dynamics.md))
$$
\dot\rho=-i[H,\rho]+\gamma_1\sum_i\Big(S^-_i\rho S^+_i-\tfrac12\{S^+_iS^-_i,\rho\}\Big)+\gamma_2\sum_k\Big(L_k\rho L_k^*-\tfrac12\{L_k^*L_k,\rho\}\Big)
$$
avec $H\in\mathcal M$ et des sauts $L_k\in\mathcal M$ à l'intérieur des secteurs, par exemple le
déphasage pur $L_k=S^z_k$. Ici $\sum_iS^+_iS^-_i=\sum_i|1\rangle\langle1|_i=N$, l'**opérateur nombre**,
$NP_v=vP_v$ : il est central.

**Proposition 10.10.** Si $\rho(0)$ est à support dans $H_w$, alors $P_w\rho(t)P_w=e^{-\gamma_1wt}\sigma(t)$,
où $\sigma$ est solution de la même équation avec $\gamma_1=0$ et $\sigma(0)=\rho(0)$. Le taux
d'acceptation vaut donc $e^{-\gamma_1wt}$, et l'état post-sélectionné $\sigma(t)$ ne dépend pas de
$\gamma_1$.

*Démonstration.* Soit $B_{vv'}=P_v\rho P_{v'}$. Comme $H,L_k,N$ commutent avec les $P_v$ et que
$P_vS^-_i=S^-_iP_{v+1}$,
$\dot B_{vv'}=-i[H,B_{vv'}]+\gamma_1\sum_iS^-_iB_{v+1,v'+1}S^+_i-\tfrac{\gamma_1}2(v+v')B_{vv'}+\gamma_2\mathcal L_2(B_{vv'})$,
où $\mathcal L_2$ est la partie de déphasage. Les blocs avec $\max(v,v')>w$ vérifient un système
linéaire fermé de valeurs initiales nulles, donc ils sont nuls (unicité dans le lemme 8.1). Alors
$\dot B_{ww}=-i[H,B_{ww}]-\gamma_1wB_{ww}+\gamma_2\mathcal L_2(B_{ww})$, de sorte que
$e^{\gamma_1wt}B_{ww}$ est solution de l'équation de $\sigma$ avec la même valeur initiale, et est égal
à $\sigma$ (lemme 8.1). Cette équation est encore de forme GKSL, donc $\operatorname{Tr}\sigma(t)=1$
(théorème 8.6). $\square$

Comparer avec la proposition 8.11, qui admet un $\sum_jL_j^*L_j$ non central mais aucun saut à
l'intérieur des secteurs : ici le terme de perte $\gamma_1N$ est central, si bien que le facteur sans
saut est le scalaire $e^{-\gamma_1wt}$ et que le déphasage à l'intérieur du secteur est permis.

Pour $\gamma_2=0$, $\sigma(t)=e^{-iHt}\rho(0)e^{iHt}$ : la post-sélection est parfaite. L'article
d'accompagnement (`examples/paper/ex9_qaoa_xy_mixer.py`) fait évoluer $|111000\rangle$ sous $H_M$
pendant $T=2$ ; nous reproduisons sa ligne $\gamma_1=\gamma_2=0.1$ : $P(w=3)=e^{-0.6}=0.5488$, fidélité
$0.4805$, post-sélectionnée $0.8755$ (comme pour $\gamma_1=0$).

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

> **Interprétation physique.** La relaxation d'énergie ($T_1$), l'une des principales erreurs des
> qubits supraconducteurs, est un saut vers un secteur inférieur. La mesure du poids de Hamming —
> un projecteur central, inoffensif à l'intérieur d'un secteur (section 3.3) — la détecte, et la
> post-sélection l'élimine au prix du taux d'acceptation. Le déphasage commute avec les projecteurs
> des secteurs et n'est pas détecté.

## 10.5 Points de comparaison

Qiskit (`qiskit.quantum_info`, modèles de bruit d'Aer, `QAOAAnsatz` avec un mélangeur personnalisé),
Cirq, PennyLane (`qml.qaoa.xy_mixer`, `default.mixed`) et QuTiP/`qutip-qip` simulent tout cela sur $n$
qubits complets et sont les bons outils pour les circuits ; stim échantillonne des circuits de
Clifford avec un bruit de type Pauli sur des milliers de qubits (l'amortissement d'amplitude n'est
pas un canal de Pauli). Nous ne connaissons dans ces logiciels aucune fonction qui calcule l'algèbre
d'interaction d'un canal, son commutant ou toute son algèbre de points fixes (`steadystate` de QuTiP
renvoie un seul état stationnaire d'un lindbladien) ; pour les décompositions numériques par blocs
des $*$-algèbres de matrices, voir Murota et al., *Japan J. Indust. Appl. Math.* 27 (2010). La
bibliothèque ajoute le stockage dans $\bigoplus_wM_{\binom nw}$, les applications entre secteurs et la
post-sélection comme projecteur central, par lots. Un secteur requiert $\binom nw^2$ coefficients au
lieu de $4^n$ ($1.7\cdot10^8$ contre $4.3\cdot10^9$ pour $n=16$, $w=8$), mais l'algèbre complète est
complétée par des zéros jusqu'à $(n+1)\binom n{\lfloor n/2\rfloor}^2$ coefficients ($2.8\cdot10^9$), une
économie modeste. Elle n'est pas un simulateur de circuits et n'a pas de formalisme des
stabilisateurs.

## Exercices

1. (Le théorème 10.2 exige l'hypothèse unifère ; comparer avec la section 6.4.) Pour
   l'amortissement d'amplitude d'un qubit, $E_0=|0\rangle\langle0|+\sqrt{1-g}\,|1\rangle\langle1|$,
   $E_1=\sqrt g\,|0\rangle\langle1|$, $0<g<1$, montrer que
   $\mathcal A'=\mathbb C1$, $\operatorname{Fix}(\Phi)=\mathbb C|0\rangle\langle0|$ et $\operatorname{Fix}(\Phi^*)=\mathbb C1$.
2. Montrer que $\sum_w\binom nw^2=\binom{2n}n$. Pour le canal du corollaire 10.3 avec $\theta=\pi/2$
   (qui viole son hypothèse), que valent $\mathcal A$ et
   $\operatorname{Fix}(\Phi)$ ? Vérifier $\dim\operatorname{Fix}(\Phi)$ pour $n=4$.
3. (Linéarité de la correction d'erreurs.) Soit $PE_a^*E_bP=\lambda_{ab}P$ pour un nombre fini
   d'opérateurs $E_a$ (pas nécessairement les opérateurs de Kraus d'un canal), et construisons
   $\mathcal R$ à partir de $\Lambda$ comme dans la démonstration du théorème 10.5. Montrer que
   $\mathcal R$ corrige tout canal dont les opérateurs de Kraus appartiennent à
   $\operatorname{span}\{E_a\}$. En déduire que le code de bit-flip corrige l'erreur cohérente
   $e^{-i\phi X_1}$.

<details>
<summary>Solutions</summary>

1. $E_1$, $E_1^*$ et leurs produits donnent toutes les unités matricielles, donc $\mathcal A=M_2$ et
   $\mathcal A'=\mathbb C1$. Pour $T=(x_{kl})$, $\Phi(T)=\begin{pmatrix}x_{00}+gx_{11}&\sqrt{1-g}\,x_{01}\\\sqrt{1-g}\,x_{10}&(1-g)x_{11}\end{pmatrix}=T$
   si et seulement si $x_{11}=x_{01}=x_{10}=0$, et
   $\Phi^*(T)=\begin{pmatrix}x_{00}&\sqrt{1-g}\,x_{01}\\\sqrt{1-g}\,x_{10}&(1-g)x_{11}+gx_{00}\end{pmatrix}$,
   qui est égal à $T$ si et seulement si $T=x_{00}1$. ($\Phi$ n'est pas unifère, et $\Phi^*$ est
   unifère mais ne préserve pas la trace.)

```python
E_ad = [torch.tensor([[1, 0], [0, math.sqrt(0.7)]], dtype=dt), torch.tensor([[0, math.sqrt(0.3)], [0, 0]], dtype=dt)]
fix_ad = fixed_points(E_ad)
assert len(fix_ad) == 1 and len(commutant_basis(E_ad + [E.mH for E in E_ad])) == 1
assert torch.allclose(fix_ad[0] / fix_ad[0][0, 0], torch.tensor([[1, 0], [0, 0]], dtype=dt))
assert len(fixed_points([E.mH for E in E_ad])) == 1                     # Fix(Phi^*) = C 1
```

2. Dénombrer les parties à $n$ éléments de $\{1,\dots,2n\}$ selon le nombre $w$ de leurs éléments
   situés dans la première moitié : $\binom{2n}n=\sum_w\binom nw\binom n{n-w}=\sum_w\binom nw^2$. Pour
   $\theta=\pi/2$, $U=(-i)^n\sum_w(-1)^wP_w$ ne voit que la parité de $w$, donc
   $\mathcal A=\operatorname{span}\{P_{\mathrm{even}},P_{\mathrm{odd}}\}$ avec $P_{\mathrm{even}}=\sum_{w\text{ even}}P_w$,
   et d'après le théorème 10.2, $\operatorname{Fix}(\Phi)$ est formé des opérateurs diagonaux par
   blocs pour la parité : dimension $2\cdot(2^{n-1})^2=128$ pour $n=4$.

```python
U2 = torch.linalg.matrix_exp(-1j * (math.pi / 2) * sum(on(Z, j, 4) for j in range(4)))
assert len(fixed_points([U2 / math.sqrt(2), U2.mH / math.sqrt(2)])) == 128
assert sum(math.comb(4, v) ** 2 for v in range(5)) == math.comb(8, 4) == 70
```

3. Avec $F_k=\sum_au_{ak}E_a$, $\operatorname{span}\{E_a\}=\operatorname{span}\{F_k\}$, donc chaque opérateur de
   Kraus du canal s'écrit $G_m=\sum_kh_{mk}F_k$, et $G_mP=\sum_{d_k>0}h_{mk}\sqrt{d_k}V_k$. Comme
   $V_l^*V_k=\delta_{lk}P$ et $(1-Q)V_k=0$, $\mathcal R(G_m\rho G_m^*)=\sum_l|h_{ml}|^2d_l\,\rho$ pour
   $\rho=P\rho P$, donc $\mathcal R(\Phi(\rho))=c\rho$ ; en prenant les traces, $c=1$. Pour
   $E_a\in\{1,X_1,X_2,X_3\}$, $\Lambda=1$, $V_k=X_kP$, et $\mathcal R$ est la récupération de la
   section 10.2 ; $e^{-i\phi X_1}=\cos\phi\,1-i\sin\phi\,X_1$ appartient à l'espace engendré.

```python
for phi in (0.2, 1.0, 2.5):
    G1 = torch.linalg.matrix_exp(-1j * phi * on(X, 0, 3))
    out = (R @ channels.kraus_channel([op8(G1)]))(rho_L)
    assert (out.matrix - rho_L.matrix).abs().max() < 1e-12       # corrected exactly
```

</details>

**Suite :** chapitre 11 (en préparation) : optimisation sur les groupes unitaires d'une algèbre, et
inégalités de traces.
