[English](../06_quantum_channels.md) · [Español](../es/06_quantum_channels.md) · **Français** · [Deutsch](../de/06_quantum_channels.md) · [中文](../zh/06_quantum_channels.md) · [日本語](../ja/06_quantum_channels.md) · [Русский](../ru/06_quantum_channels.md)

# Tutoriel 6. Canaux quantiques

Ce tutoriel accompagne le [chapitre 6](../../lecture_notes/fr/06_quantum_channels.md) des notes de cours, qui démontre les descriptions des
applications complètement positives entre sommes directes d'algèbres de matrices — matrices de Choi,
opérateurs de Kraus, dilatation de Stinespring — et étudie les duales, les points fixes, les matrices de
transition et les canaux aléatoires. Nous voyons ici comment la bibliothèque stocke les canaux et
vérifions numériquement ces résultats.

Dans tout ce chapitre, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ et
$\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)$, de projecteurs centraux minimaux $P_c$ et
$Q_d$ ; $A=\bigoplus_cA_c$ ; $M_{k'\times k}$ est l'espace des matrices complexes $k'\times k$ ;
les $E_{ij}=|i\rangle\langle j|$ sont les unités matricielles ; $\mathrm{Tr}$ est la trace brute
$\sum_c\mathrm{Tr}\,A_c$, et $\mathrm{tr}_w=\sum_cw_c\mathrm{Tr}\,A_c$ une trace fidèle
(théorème 3.1). Les matrices par blocs agissent sur $\mathbb C^n\otimes\mathbb C^k$, l'indice de
bloc venant en premier.

<a id="s-6-1"></a>

## 6.1 Applications complètement positives et théorème de Choi

*Notes de cours :* [section 6.1](../../lecture_notes/fr/06_quantum_channels.md#s-6-1) : [proposition 6.1](../../lecture_notes/fr/06_quantum_channels.md#r-6-1) (la transposition), [théorème 6.2](../../lecture_notes/fr/06_quantum_channels.md#r-6-2) (Choi), [corollaire 6.3](../../lecture_notes/fr/06_quantum_channels.md#r-6-3) (rang de Kraus) et [proposition 6.4](../../lecture_notes/fr/06_quantum_channels.md#r-6-4) (applications entre sommes directes).

La bibliothèque comporte deux classes. `Channel` stocke les applications **qui préservent les
secteurs** ($C'=C$, $k'=k$, $\Phi_{dc}=0$ pour $d\ne c$), dont les opérateurs de Kraus
$K_l=\bigoplus_cK^{cc}_l$ appartiennent à $\mathcal M$, sous la forme d'un tenseur de forme
`(batch, r, C, k_max, k_max)` ; `InterSectorChannel` stocke l'application générale de la
[proposition 6.4](../../lecture_notes/fr/06_quantum_channels.md#r-6-4), de forme `(batch, r, C', C, k'_max, k_max)`. `Channel.choi()` renvoie
$J(\Phi_{cc})$ pour chaque secteur et `Channel.superoperator()` la matrice $S_c$ de $\Phi_{cc}$ sur
les vecteurs rangés ligne par ligne, $\mathrm{vec}(X)_{(a,b)}=X_{ab}$. Toutes deux contiennent
les nombres $\Phi(E_{ij})_{ab}=J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$, dans un ordre différent.
`Channel.from_superoperator` lève `ValueError` si une matrice de Choi n'est pas hermitienne,
$\|J-J^*\|>2\cdot10^{-10}s$ (en double précision), où $s=\max(1,\|J\|)$ avec la norme spectrale ;
sinon, elle remplace $J$ par sa partie hermitienne $(J+J^*)/2$, la diagonalise comme dans la
démonstration du [théorème 6.2](../../lecture_notes/fr/06_quantum_channels.md#r-6-2), lève `ValueError` si une valeur propre est inférieure à $-10^{-10}s$,
et met à $0$ les valeurs propres négatives restantes. Comme $J^*=J$ si et seulement si
$\Phi(X^*)=\Phi(X)^*$ pour tout $X$ ([notes de cours, section 6.1](../../lecture_notes/fr/06_quantum_channels.md#s-6-1)), aux tolérances d'arrondi près,
elle refuse exactement les applications qui ne préservent pas les adjoints et, parmi les autres,
exactement celles qui ne sont pas CP ; par exemple, $X\mapsto iX$ et la transposition sont refusées.

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, states, channels
from torch_vn_algebra.channels import Channel, InterSectorChannel

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def unit(i, j, k):
    X = torch.zeros(k, k, dtype=dt)
    X[i, j] = 1
    return X

def choi(f, k):
    """Choi matrix sum_ij E_ij (x) f(E_ij) of a linear map f on k x k matrices."""
    return sum(torch.kron(unit(i, j, k), f(unit(i, j, k)).contiguous()) for i in range(k) for j in range(k))

def to_superoperator(J, k):
    """S[(a,b),(i,j)] = J[(i,a),(j,b)], with a leading batch axis."""
    return J.reshape(k, k, k, k).permute(1, 3, 0, 2).reshape(1, k * k, k * k)

ev = torch.linalg.eigvalsh(choi(lambda X: X.T, 2))
print(ev)                                                 # the Choi matrix of the transpose is the swap
assert abs(ev[0] + 1) < 1e-12                             # Proposition 6.1

alg = TypeIAlgebra([2, 3], [2, 3], precision='double')    # M_2 (+) M_3
Phi = channels.random_channel(alg, kraus_rank=2)          # Kraus tensor of shape (1, 2, 2, 3, 3)

def block(Phi, c):
    """The block Phi_cc of a sector-preserving map, as a function on k_c x k_c matrices."""
    ks = Phi.algebra.k_factors
    def f(X):
        blocks = [torch.zeros(k, k, dtype=dt) for k in ks]
        blocks[c] = X
        return Phi(Phi.algebra.from_blocks(blocks)).matrix[0, c, :ks[c], :ks[c]]
    return f

J, S = Phi.choi(), Phi.superoperator()                    # one matrix per sector, (1, k_c^2, k_c^2)
rho = states.random_density_matrix(alg)
out = Phi(rho)
for c, k in enumerate(alg.k_factors):
    assert close(J[c][0], choi(block(Phi, c), k))         # the convention of Theorem 6.2
    w = torch.linalg.eigvalsh(J[c][0])
    assert w.min() > -1e-12 and int((w > 1e-10).sum()) == 2      # J >= 0, Kraus rank 2 (Corollary 6.3)
    assert close(S[c][0], to_superoperator(J[c][0], k)[0])
    assert close(S[c][0] @ rho.matrix[0, c, :k, :k].reshape(-1), out.matrix[0, c, :k, :k].reshape(-1))
assert Phi.is_trace_preserving()
assert close(Channel.from_superoperator(alg, S)(rho), out)       # back from the superoperator
try:
    Channel.from_superoperator(alg, [to_superoperator(choi(lambda X: X.T, k), k) for k in (2, 3)])
    raise AssertionError("the transpose was accepted")
except ValueError as e:
    print(e)                                              # superoperator is not completely positive
try:
    Channel.from_superoperator(alg, [1j * torch.eye(k * k, dtype=dt)[None] for k in (2, 3)])
    raise AssertionError("X -> iX was accepted")
except ValueError as e:
    print(e)                                              # superoperator does not preserve Hermiticity
```

<a id="s-6-2"></a>

## 6.2 Dilatation de Stinespring et liberté dans le choix des opérateurs de Kraus

*Notes de cours :* [section 6.2](../../lecture_notes/fr/06_quantum_channels.md#s-6-2) : [théorème 6.5](../../lecture_notes/fr/06_quantum_channels.md#r-6-5) (Stinespring) et [théorème 6.7](../../lecture_notes/fr/06_quantum_channels.md#r-6-7) (liberté unitaire). Pour le bloc $c=1$ du canal de la section 6.1, le code vérifie que $V$ est une isométrie, que $\Phi(\rho)=\mathrm{Tr}_E(V\rho V^*)$ et $\Phi^*(B)=V^*(1\otimes B)V$, et que les opérateurs de Kraus $\sum_lu_{jl}K_l$ donnent le même canal.

Dans le code, l'indice d'environnement de $V$ est le premier indice du tenseur de Kraus :

```python
K = Phi.kraus[0, :, 1, :3, :3]                            # the two Kraus operators of block c = 1
V = K.reshape(2 * 3, 3)                                   # V = sum_l e_l (x) K_l : C^3 -> C^2 (x) C^3
assert close(V.conj().T @ V, torch.eye(3, dtype=dt))      # TP <=> V is an isometry
r1 = rho.matrix[0, 1, :3, :3]
TrE = torch.einsum('lalb->ab', (V @ r1 @ V.conj().T).reshape(2, 3, 2, 3))
assert close(TrE, out.matrix[0, 1, :3, :3])               # Phi(rho) = Tr_E(V rho V^*)
B1 = torch.randn(3, 3, dtype=dt)
Bop = alg.from_blocks([torch.zeros(2, 2, dtype=dt), B1])
assert close(V.conj().T @ torch.kron(torch.eye(2, dtype=dt), B1) @ V,
             Phi.adjoint()(Bop).matrix[0, 1, :3, :3])    # Phi^*(B) = V^*(1 (x) B)V

u = alg.random_unitary(3)                                 # Theorem 6.7 with n = 3: pad with a zero operator
Kpad = torch.cat([Phi.kraus, torch.zeros_like(Phi.kraus[:, :1])], dim=1)
PhiL = Channel(alg, torch.einsum('jl,blcxy->bjcxy', u, Kpad))
assert close(PhiL(rho), out) and all(close(a, b) for a, b in zip(PhiL.choi(), J))
```

<a id="s-6-3"></a>

## 6.3 Applications duales relativement à une trace

*Notes de cours :* [section 6.3](../../lecture_notes/fr/06_quantum_channels.md#s-6-3) : l'application duale $\Phi^\sharp$ relativement à $(w,v)$, [proposition 6.8](../../lecture_notes/fr/06_quantum_channels.md#r-6-8) et [corollaire 6.9](../../lecture_notes/fr/06_quantum_channels.md#r-6-9).

Dans la bibliothèque, `Channel.adjoint()` est $\Phi^*$, et `InterSectorChannel.adjoint(trace=...)` est
$\Phi^\sharp$ pour la même trace nommée sur les deux algèbres (`'blunt'`, `'norm'`, `'tau_vN'` ; les
poids $1/(Ck_c)$ et $1/(C'k'_d)$ de `tau_vN` tiennent aussi compte du nombre de secteurs). Vérifions
la [proposition 6.8](../../lecture_notes/fr/06_quantum_channels.md#r-6-8) pour une application aléatoire $M_2\oplus M_3\to M_1\oplus M_2\oplus M_2$ et les
trois traces, et le [corollaire 6.9(a)](../../lecture_notes/fr/06_quantum_channels.md#r-6-9) pour la perte de particule $\Phi(x)=(x_0+\mathrm{Tr}\,x_1)\oplus0$
sur $M_1\oplus M_2$, qui est TP pour les traces brutes mais pas pour `Tr_norm`
([notes de cours, section 6.3](../../lecture_notes/fr/06_quantum_channels.md#s-6-3)) :

```python
alg_out = TypeIAlgebra([1, 2, 2], [1, 2, 2], precision='double')
Psi = channels.random_inter_sector_channel(alg, alg_out, kraus_rank=2)     # Kraus tensor (1, 2, 3, 2, 2, 3)
x = alg.operator(torch.randn(1, 2, 3, 3, dtype=dt))
y = alg_out.operator(torch.randn(1, 3, 2, 2, dtype=dt))
weights = {'blunt': lambda k, C: 1.0, 'norm': lambda k, C: 1 / k, 'tau_vN': lambda k, C: 1 / (C * k)}
trace = {'blunt': lambda X: X.Tr_blunt(), 'norm': lambda X: X.Tr_norm(), 'tau_vN': lambda X: X.tau_vN()}

def z(a, t):                                              # the central element sum_c w_c P_c of trace t
    return a.central(torch.tensor([weights[t](k, a.C) for k in a.k_factors], dtype=torch.float64))

for t in trace:
    dual = Psi.adjoint(trace=t)
    assert (trace[t](Psi(x) @ y) - trace[t](x @ dual(y))).abs().max() < 1e-10   # the definition
    assert close(dual(y), z(alg, t).inverse() @ Psi.adjoint()(z(alg_out, t) @ y))   # Proposition 6.8
dev = [(Psi.adjoint(trace=t)(alg_out.identity()) - alg.identity()).matrix.abs().max().item() for t in trace]
print(dev)                         # TP only for the blunt traces: the other duals are not unital
assert Psi.is_trace_preserving() and dev[0] < 1e-12 and min(dev[1:]) > 0.1

alg01 = TypeIAlgebra([1, 2], [1, 2], precision='double')
loss = InterSectorChannel.from_blocks(alg01, alg01, {(0, 0): [torch.ones(1, 1)],
                                                     (0, 1): [torch.tensor([[1.0, 0.0]]), torch.tensor([[0.0, 1.0]])]})
assert loss.is_trace_preserving()
assert close(loss.adjoint(trace='norm')(alg01.identity()), alg01.central([1.0, 2.0]))
s = states.random_density_matrix(alg01)
assert abs(loss(s).Tr_norm() - (s @ alg01.central([1.0, 2.0])).Tr_norm()) < 1e-12   # Corollary 6.9(a)
print(s.Tr_norm().real.item(), loss(s).Tr_norm().real.item())      # Tr_norm is not conserved
```

<a id="s-6-4"></a>

## 6.4 Décohérence

*Notes de cours :* [section 6.4](../../lecture_notes/fr/06_quantum_channels.md#s-6-4) : [théorème 6.10](../../lecture_notes/fr/06_quantum_channels.md#r-6-10) (points fixes des canaux unifères), [théorème 6.11](../../lecture_notes/fr/06_quantum_channels.md#r-6-11) (les canaux unifères augmentent l'entropie) et [proposition 6.12](../../lecture_notes/fr/06_quantum_channels.md#r-6-12) (la matrice de transition).

Les canaux standard de la bibliothèque agissent à l'intérieur de chaque bloc ($p,\gamma\in[0,1]$) :

| canal | $\Phi(\rho)_c$ | opérateurs de Kraus dans le bloc $c$ |
|---|---|---|
| `dephasing_channel` | $(1-p)\rho_c+p\mathrm{diag}\,\rho_c$ | $\sqrt{1-p}\,1$, $\sqrt p\,E_{ii}$ |
| `depolarizing_channel` | $(1-p)\rho_c+p\mathrm{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1-p}\,1$, $\sqrt{p/k_c}\,E_{ij}$ |
| `center_expectation` | $\mathrm{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1/k_c}\,E_{ij}$ |
| `amplitude_damping_channel` | les niveaux $i>0$ se désexcitent vers $0$ | $E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}$, $\sqrt\gamma\,E_{0i}$ |
| `random_mixed_unitary_channel` | $\sum_jp_jU_j\rho U_j^*$ | $\sqrt{p_j}\,U_j$ |

Pour chaque canal du tableau sur $M_2\oplus M_3$, le code calcule la dimension de ses points fixes (le
noyau de $S_c-1$) et celle du commutant de ses opérateurs de Kraus : elles coïncident pour les canaux
unifères ([théorème 6.10](../../lecture_notes/fr/06_quantum_channels.md#r-6-10)), et pour l'amortissement d'amplitude les dimensions sont égales mais les
espaces différents. Il vérifie ensuite le [théorème 6.11](../../lecture_notes/fr/06_quantum_channels.md#r-6-11) sur $100$ états aléatoires, et la
[proposition 6.12](../../lecture_notes/fr/06_quantum_channels.md#r-6-12) pour l'application `Psi` de la section 6.3 et pour la perte de particule.

```python
def fixed_dims(Phi):
    """dim {X in M_{k_c} : Phi_cc(X) = X} for every sector: the null space of S_c - 1."""
    return [int((torch.linalg.svdvals(S_c[0] - torch.eye(k * k, dtype=dt)) < 1e-9).sum())
            for S_c, k in zip(Phi.superoperator(), Phi.algebra.k_factors)]

def commutant_dims(Phi):
    """dim of the commutant of {K_l, K_l^*} in every block (row-major vec, as in Chapter 2)."""
    out = []
    for c, k in enumerate(Phi.algebra.k_factors):
        I, Ks = torch.eye(k, dtype=dt), Phi.kraus[0, :, c, :k, :k]
        ops = [A.contiguous() for A in torch.cat([Ks, Ks.mH])]
        L = torch.cat([torch.kron(A, I) - torch.kron(I, A.T.contiguous()) for A in ops])
        out.append(k * k - int((torch.linalg.svdvals(L) > 1e-9).sum()))
    return out

zoo = {'dephasing': channels.dephasing_channel(alg, 0.3),
       'depolarizing': channels.depolarizing_channel(alg, 0.3),
       'centre': channels.center_expectation(alg),
       'mixed unitary': channels.random_mixed_unitary_channel(alg, 2),
       'amplitude damping': channels.amplitude_damping_channel(alg, 0.4)}
for name, ch in zoo.items():
    print(f"{name:18s} unital {ch.is_unital()!s:5s} fixed {fixed_dims(ch)} commutant {commutant_dims(ch)}")
    assert ch.is_trace_preserving()
    if ch.is_unital():
        assert fixed_dims(ch) == commutant_dims(ch)                         # Theorem 6.10
assert fixed_dims(zoo['dephasing']) == [2, 3] and fixed_dims(zoo['depolarizing']) == [1, 1]
E00 = alg.from_blocks([unit(0, 0, 2), unit(0, 0, 3)])
E01 = alg.from_blocks([unit(0, 1, 2), unit(0, 1, 3)])
assert close(zoo['amplitude damping'](E00), E00) and not close(E00 @ E01, E01 @ E00)   # equal dims, different spaces
assert fixed_dims(zoo['amplitude damping']) == commutant_dims(zoo['amplitude damping']) == [1, 1]

rhos = states.random_density_matrix(alg, batch_size=100)
for name in ('dephasing', 'depolarizing', 'centre', 'mixed unitary'):                 # Theorem 6.11
    assert (zoo[name](rhos).von_neumann_entropy() >= rhos.von_neumann_entropy() - 1e-12).all()
    assert (zoo[name](rhos).purity() <= rhos.purity() + 1e-12).all()
mm = states.maximally_mixed_state(alg)
S_after = channels.amplitude_damping_channel(alg, 1.0)(mm).von_neumann_entropy().item()
print(mm.von_neumann_entropy().item(), S_after)          # not unital: the entropy drops
assert abs(S_after + 0.4 * math.log(0.4) + 0.6 * math.log(0.6)) < 1e-12

T = Psi.transition_matrix()[0]                            # (C', C) = (3, 2)
p, q = rho.sector_probabilities()[0], Psi(rho).sector_probabilities()[0]
assert torch.allclose(T.sum(0), torch.ones(2, dtype=T.dtype))                       # Proposition 6.12
assert torch.allclose(q, Psi.transition_matrix(rho)[0] @ p)
classical = channels.center_expectation(alg_out) @ Psi @ channels.center_expectation(alg)
Tp = T @ p
assert close(classical(rho), alg_out.central([Tp[d] / k for d, k in enumerate(alg_out.k_factors)]))
assert torch.allclose(loss.transition_matrix()[0], torch.tensor([[1.0, 1.0], [0.0, 0.0]], dtype=T.dtype))
```

Sous l'amortissement d'amplitude complet, l'entropie de l'état maximalement mélangé de
$M_2\oplus M_3$ tombe de $\log5\approx1.609$ à $H(\tfrac25,\tfrac35)\approx0.673$ : un canal non
unifère peut purifier.

Sur $M_2$, le déphasage $D_p$ et l'amortissement d'amplitude $A_\gamma$ commutent, et leur composée
multiplie les cohérences par $(1-p)\sqrt{1-\gamma}$ ; avec $1-p=e^{-t/T_\varphi}$ et $1-\gamma=e^{-t/T_1}$, cela
donne $1/T_2=1/(2T_1)+1/T_\varphi$ ([notes de cours, section 6.4](../../lecture_notes/fr/06_quantum_channels.md#s-6-4), *Cohérences d'un qubit*). Le code
vérifie ces identités et les lois de composition $D_{p(t)}D_{p(s)}=D_{p(t+s)}$, $A_{\gamma(t)}A_{\gamma(s)}=A_{\gamma(t+s)}$ :

```python
q2 = TypeIAlgebra([2], [2], precision='double')
p_, g_ = 0.3, 0.4
Dp, Ag = channels.dephasing_channel(q2, p_), channels.amplitude_damping_channel(q2, g_)
r2 = states.random_density_matrix(q2)
assert close(Dp(Ag(r2)), Ag(Dp(r2)))                                         # D_p A_gamma = A_gamma D_p
coh = Dp(Ag(r2)).matrix[0, 0, 0, 1] / r2.matrix[0, 0, 0, 1]
assert abs(coh - (1 - p_) * math.sqrt(1 - g_)) < 1e-12                       # coherence factor
T1, Tphi, t1, t2 = 2.0, 3.0, 0.7, 1.1
gam = lambda t: 1 - math.exp(-t / T1)
assert close(channels.amplitude_damping_channel(q2, gam(t1))(channels.amplitude_damping_channel(q2, gam(t2))(r2)),
             channels.amplitude_damping_channel(q2, gam(t1 + t2))(r2))       # A_gamma(t) A_gamma(s) = A_gamma(t+s)
pdeph = lambda t: 1 - math.exp(-t / Tphi)
assert close(channels.dephasing_channel(q2, pdeph(t1))(channels.dephasing_channel(q2, pdeph(t2))(r2)),
             channels.dephasing_channel(q2, pdeph(t1 + t2))(r2))           # D_p(t) D_p(s) = D_p(t+s)
both = channels.dephasing_channel(q2, pdeph(t1))(channels.amplitude_damping_channel(q2, gam(t1))(r2))
coh_t = (both.matrix[0, 0, 0, 1] / r2.matrix[0, 0, 0, 1]).item()
assert abs(coh_t - math.exp(-t1 * (1 / (2 * T1) + 1 / Tphi))) < 1e-12       # 1/T_2 = 1/(2T_1) + 1/T_phi
```

<a id="s-6-5"></a>

## 6.5 Canaux aléatoires

*Notes de cours :* [section 6.5](../../lecture_notes/fr/06_quantum_channels.md#s-6-5) : [lemme 6.13](../../lecture_notes/fr/06_quantum_channels.md#r-6-13) et [théorème 6.15](../../lecture_notes/fr/06_quantum_channels.md#r-6-15) (pureté de sortie d'un canal aléatoire).

Une **isométrie aléatoire de Haar** $V:\mathbb C^k\to\mathbb C^N$ est, comme dans la bibliothèque
(`random_isometry`), le facteur $Q$ de $G=QR$, où $G$ est une matrice de Ginibre $N\times k$
(coefficients indépendants de densité $\pi^{-1}e^{-|z|^2}$), $Q^*Q=1$, et $R$ est triangulaire
supérieure à diagonale positive (sur le CPU, la bibliothèque prend une factorisation LAPACK
$G=Q'R'$ et multiplie la $i$-ième colonne de $Q'$ par $R'_{ii}/|R'_{ii}|$, ce qui donne cette
factorisation, comme dans le théorème 5.3, dès que tous les $R'_{ii}\ne0$, c'est-à-dire dès que $G$
est de rang $k$, ce qui a lieu presque sûrement d'après la démonstration du lemme 6.13 ; pour de
grands lots sur un GPU, elle utilise Gram–Schmidt, qui la donne directement).
`random_channel(alg, r)` tire de telles
$V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}$ indépendamment dans chaque bloc et utilise
les opérateurs de Kraus du théorème 6.5.

D'après le [théorème 6.15 des notes de cours](../../lecture_notes/fr/06_quantum_channels.md#r-6-15), pour un état pur $\psi$ dans un secteur de taille $k$ et $\Phi$
tiré par `random_channel` avec un rang de Kraus $r$, $\mathbb E\mathrm{Tr}\,\Phi(\psi)^2=(k+r)/(kr+1)$ et
$\mathbb E\,\Phi(\psi)=1_k/k$. Nous estimons ces deux quantités sur $20\,000$ échantillons, avec l'erreur type de la moyenne :

```python
B = 20_000
for r in (1, 2, 4):
    Phi_r = channels.random_channel(alg, kraus_rank=r, batch_size=B)
    for c, k in enumerate(alg.k_factors):
        psi = torch.randn(k, 1, dtype=dt)
        psi = psi / psi.norm()
        out_r = Phi_r(alg.from_blocks([psi @ psi.conj().T if d == c else torch.zeros(kk, kk, dtype=dt)
                                       for d, kk in enumerate(alg.k_factors)]))     # a batch of B outputs
        purity = (out_r @ out_r).trace.real
        mean, err = purity.mean().item(), purity.std().item() / math.sqrt(B)
        exact = (k + r) / (k * r + 1)
        print(f"k = {k}, r = {r}: E Tr Phi(psi)^2 = {mean:.4f} +- {err:.4f}   (exact {exact:.4f})")
        assert abs(mean - exact) < 5 * err + 1e-12
        avg = out_r.matrix.mean(0)[c, :k, :k]
        assert (avg - torch.eye(k, dtype=dt) / k).abs().max() < 0.02          # E Phi(psi) = 1/k
```

<a id="s-6-6"></a>

## 6.6 Points de comparaison

QuTiP (`to_choi`, `to_kraus`, `to_super`, `to_stinespring`), Qiskit (`qiskit.quantum_info` : les
classes `Kraus`, `Choi`, `SuperOp`, `Stinespring` et `random_quantum_channel` ; des modèles de bruit
comme `amplitude_damping_error` dans Qiskit Aer), toqito (`toqito.channel_ops`,
`toqito.channel_props`) et QETLAB (`ChoiMatrix`, `KrausOperators`, `ApplyMap`,
`RandomSuperoperator`) convertissent entre les représentations des sections 6.1–6.2, testent la
complète positivité et (pour la plupart) tirent des canaux aléatoires, pour des applications
$M_n\to M_m$ entre algèbres de matrices pleines.
Leurs conventions diffèrent (empilement par colonnes ou par lignes, ordre des facteurs de la matrice de
Choi) : comparez donc les conventions avant de comparer les nombres. Un canal qui préserve les
secteurs sur $\mathcal M$ peut y être traité comme un canal sur $M_D$, $D=\sum_ck_c$, avec des
superopérateurs de taille $D^2\times D^2$ au lieu de $k_c^2\times k_c^2$ par secteur. La bibliothèque
ajoute la structure par blocs elle-même : opérateurs de Kraus stockés par secteur et par lots,
applications entre algèbres différentes (`InterSectorChannel`), duales pour chacune des traces du
chapitre 3, et matrice de transition entre secteurs.

<a id="x-6"></a>

## Exercices

[Exercices](../../lecture_notes/fr/06_quantum_channels.md#x-6) des notes de cours, numériquement :

1. (Werner–Holevo) Vérifier avec `Channel.from_superoperator` que $\Phi(\rho)=\big(\mathrm{Tr}(\rho)1-\rho^T\big)/(k-1)$
   sur $M_3$ est un canal unifère de rang de Kraus $3$.
2. Vérifier que le canal de déphasage sur $M_2\oplus M_3$ a pour rang de Kraus $k_c$ dans le bloc $c$, alors que la
   bibliothèque stocke $1+k_{\max}$ opérateurs de Kraus.
3. Vérifier que le canal unitaire mixte d'opérateurs de Kraus $X^aZ^b/k$ est `center_expectation` sur $M_3$.
4. Vérifier $(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp$ pour les trois traces de la bibliothèque.

<details>
<summary>Solutions</summary>

*Exercice 1.*

```python
k = 3
M3 = TypeIAlgebra([k], [k], precision='double')
J_WH = choi(lambda X: (X.trace() * torch.eye(k, dtype=dt) - X.T) / (k - 1), k)
w = torch.linalg.eigvalsh(J_WH)
assert w.min() > -1e-12 and int((w > 1e-10).sum()) == k * (k - 1) // 2
WH = Channel.from_superoperator(M3, [to_superoperator(J_WH, k)])
assert WH.is_trace_preserving() and WH.is_unital()
r3 = states.random_density_matrix(M3).matrix[0, 0]
assert close(WH(M3.operator(r3[None, None])).matrix[0, 0], (torch.eye(k, dtype=dt) - r3.T) / (k - 1))
```

*Exercice 2.*

```python
ranks = [int((torch.linalg.eigvalsh(J_c[0]) > 1e-10).sum()) for J_c in zoo['dephasing'].choi()]
print(ranks, zoo['dephasing'].kraus_rank)                 # [2, 3] versus 4 stored operators
assert ranks == alg.k_factors and zoo['dephasing'].kraus_rank == 1 + max(alg.k_factors)
```

*Exercice 3.*

```python
Xs = torch.roll(torch.eye(k, dtype=dt), 1, dims=0)
Zs = torch.diag(torch.exp(2j * math.pi * torch.arange(k, dtype=torch.float64) / k)).to(dt)
W = [torch.linalg.matrix_power(Xs, a) @ torch.linalg.matrix_power(Zs, b) / k for a in range(k) for b in range(k)]
twirl = channels.kraus_channel([M3.operator(Wab[None, None]) for Wab in W])
r3op = states.random_density_matrix(M3)
assert close(twirl(r3op), channels.center_expectation(M3)(r3op))
```

*Exercice 4.*

```python
Psi2 = channels.random_inter_sector_channel(alg_out, alg, kraus_rank=2)
for t in trace:
    assert close((Psi2 @ Psi).adjoint(trace=t)(x), Psi.adjoint(trace=t)(Psi2.adjoint(trace=t)(x)))
```

</details>

**Suite :** [Tutoriel 7](07_tensor_products.md) : produits tensoriels, traces partielles et intrication dans la bibliothèque.
