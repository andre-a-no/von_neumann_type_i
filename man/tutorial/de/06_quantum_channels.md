[English](../06_quantum_channels.md) · [Español](../es/06_quantum_channels.md) · [Français](../fr/06_quantum_channels.md) · **Deutsch** · [中文](../zh/06_quantum_channels.md) · [日本語](../ja/06_quantum_channels.md) · [Русский](../ru/06_quantum_channels.md)

# Tutorial 6. Quantenkanäle

Dieses Tutorial begleitet [Kapitel 6](../../lecture_notes/de/06_quantum_channels.md) der Vorlesungsnotizen, das die Beschreibungen vollständig
positiver Abbildungen zwischen direkten Summen von Matrixalgebren beweist – Choi-Matrizen,
Kraus-Operatoren, die Stinespring-Dilatation – und duale Abbildungen, Fixpunkte, Übergangsmatrizen und
zufällige Kanäle untersucht. Hier sehen wir, wie die Bibliothek Kanäle speichert, und überprüfen diese
Ergebnisse numerisch.

Im Folgenden seien $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ und
$\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)$ mit minimalen zentralen Projektionen $P_c$ und
$Q_d$; $A=\bigoplus_cA_c$; $M_{k'\times k}$ ist der Raum der komplexen $k'\times k$-Matrizen;
$E_{ij}=|i\rangle\langle j|$ sind die Matrixeinheiten; $\mathrm{Tr}$ ist die ungewichtete Spur
$\sum_c\mathrm{Tr}\,A_c$ und $\mathrm{tr}_w=\sum_cw_c\mathrm{Tr}\,A_c$ eine treue Spur
(Satz 3.1). Blockmatrizen wirken auf $\mathbb C^n\otimes\mathbb C^k$, der Blockindex zuerst.

<a id="s-6-1"></a>

## 6.1 Vollständig positive Abbildungen und der Satz von Choi

*Vorlesungsnotizen:* [Abschnitt 6.1](../../lecture_notes/de/06_quantum_channels.md#s-6-1): [Proposition 6.1](../../lecture_notes/de/06_quantum_channels.md#r-6-1) (die Transposition), [Satz 6.2](../../lecture_notes/de/06_quantum_channels.md#r-6-2) (Choi), [Korollar 6.3](../../lecture_notes/de/06_quantum_channels.md#r-6-3) (Kraus-Rang) und [Proposition 6.4](../../lecture_notes/de/06_quantum_channels.md#r-6-4) (Abbildungen zwischen direkten Summen).

Die Bibliothek hat zwei Klassen. `Channel` speichert **sektorerhaltende** Abbildungen ($C'=C$, $k'=k$,
$\Phi_{dc}=0$ für $d\ne c$), deren Kraus-Operatoren $K_l=\bigoplus_cK^{cc}_l$ in $\mathcal M$ liegen,
als Tensor der Form `(batch, r, C, k_max, k_max)`; `InterSectorChannel` speichert die allgemeine
Abbildung aus [Proposition 6.4](../../lecture_notes/de/06_quantum_channels.md#r-6-4) mit der Form `(batch, r, C', C, k'_max, k_max)`. `Channel.choi()` liefert
$J(\Phi_{cc})$ für jeden Sektor und `Channel.superoperator()` die Matrix $S_c$ von $\Phi_{cc}$ auf
zeilenweise angeordneten Vektoren, $\mathrm{vec}(X)_{(a,b)}=X_{ab}$. Beide enthalten die Zahlen
$\Phi(E_{ij})_{ab}=J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$ in unterschiedlicher Anordnung.
`Channel.from_superoperator` löst `ValueError` aus, wenn eine Choi-Matrix nicht hermitesch ist,
$\|J-J^*\|>2\cdot10^{-10}s$ (in doppelter Genauigkeit), wobei $s=\max(1,\|J\|)$ mit der Spektralnorm;
andernfalls ersetzt es $J$ durch seinen hermiteschen Teil $(J+J^*)/2$, diagonalisiert diesen wie im
Beweis von [Satz 6.2](../../lecture_notes/de/06_quantum_channels.md#r-6-2), löst `ValueError` aus, wenn ein Eigenwert kleiner als $-10^{-10}s$ ist, und setzt
die übrigen negativen Eigenwerte auf $0$. Da $J^*=J$ genau dann gilt, wenn $\Phi(X^*)=\Phi(X)^*$ für alle $X$
([Vorlesungsnotizen, Abschnitt 6.1](../../lecture_notes/de/06_quantum_channels.md#s-6-1)), weist es bis auf diese Rundungstoleranzen genau die Abbildungen
zurück, die Adjungierte nicht erhalten, und unter den übrigen genau die, die nicht CP sind; zum Beispiel
werden $X\mapsto iX$ und die Transposition zurückgewiesen.

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

## 6.2 Stinespring-Dilatation und die Freiheit in den Kraus-Operatoren

*Vorlesungsnotizen:* [Abschnitt 6.2](../../lecture_notes/de/06_quantum_channels.md#s-6-2): [Satz 6.5](../../lecture_notes/de/06_quantum_channels.md#r-6-5) (Stinespring) und [Satz 6.7](../../lecture_notes/de/06_quantum_channels.md#r-6-7) (unitäre Freiheit). Für den Block $c=1$ des Kanals aus Abschnitt 6.1 prüft der Code, dass $V$ eine Isometrie ist, dass $\Phi(\rho)=\mathrm{Tr}_E(V\rho V^*)$ und $\Phi^*(B)=V^*(1\otimes B)V$ gelten und dass die Kraus-Operatoren $\sum_lu_{jl}K_l$ denselben Kanal ergeben.

Im Code ist der Umgebungsindex von $V$ der erste Index des Kraus-Tensors:

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

## 6.3 Duale Abbildungen bezüglich einer Spur

*Vorlesungsnotizen:* [Abschnitt 6.3](../../lecture_notes/de/06_quantum_channels.md#s-6-3): die duale Abbildung $\Phi^\sharp$ bezüglich $(w,v)$, [Proposition 6.8](../../lecture_notes/de/06_quantum_channels.md#r-6-8) und [Korollar 6.9](../../lecture_notes/de/06_quantum_channels.md#r-6-9).

In der Bibliothek ist `Channel.adjoint()` gleich $\Phi^*$ und `InterSectorChannel.adjoint(trace=...)`
gleich $\Phi^\sharp$ für dieselbe benannte Spur auf beiden Algebren (`'blunt'`, `'norm'`, `'tau_vN'`;
die Gewichte $1/(Ck_c)$ und $1/(C'k'_d)$ von `tau_vN` berücksichtigen auch die Anzahl der Sektoren).
Wir prüfen [Proposition 6.8](../../lecture_notes/de/06_quantum_channels.md#r-6-8) für eine zufällige Abbildung $M_2\oplus M_3\to M_1\oplus M_2\oplus M_2$ und
die drei Spuren sowie [Korollar 6.9(a)](../../lecture_notes/de/06_quantum_channels.md#r-6-9) für den Teilchenverlust $\Phi(x)=(x_0+\mathrm{Tr}\,x_1)\oplus0$ auf
$M_1\oplus M_2$, der TP für die ungewichteten Spuren ist, aber nicht für `Tr_norm`
([Vorlesungsnotizen, Abschnitt 6.3](../../lecture_notes/de/06_quantum_channels.md#s-6-3)):

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

## 6.4 Dekohärenz

*Vorlesungsnotizen:* [Abschnitt 6.4](../../lecture_notes/de/06_quantum_channels.md#s-6-4): [Satz 6.10](../../lecture_notes/de/06_quantum_channels.md#r-6-10) (Fixpunkte unitaler Kanäle), [Satz 6.11](../../lecture_notes/de/06_quantum_channels.md#r-6-11) (unitale Kanäle erhöhen die Entropie) und [Proposition 6.12](../../lecture_notes/de/06_quantum_channels.md#r-6-12) (die Übergangsmatrix).

Die Standardkanäle der Bibliothek wirken innerhalb jedes Blocks ($p,\gamma\in[0,1]$):

| Kanal | $\Phi(\rho)_c$ | Kraus-Operatoren im Block $c$ |
|---|---|---|
| `dephasing_channel` | $(1-p)\rho_c+p\mathrm{diag}\,\rho_c$ | $\sqrt{1-p}\,1$, $\sqrt p\,E_{ii}$ |
| `depolarizing_channel` | $(1-p)\rho_c+p\mathrm{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1-p}\,1$, $\sqrt{p/k_c}\,E_{ij}$ |
| `center_expectation` | $\mathrm{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1/k_c}\,E_{ij}$ |
| `amplitude_damping_channel` | Niveaus $i>0$ zerfallen nach $0$ | $E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}$, $\sqrt\gamma\,E_{0i}$ |
| `random_mixed_unitary_channel` | $\sum_jp_jU_j\rho U_j^*$ | $\sqrt{p_j}\,U_j$ |

Für jeden Kanal der Tabelle auf $M_2\oplus M_3$ berechnet der Code die Dimension seiner Fixpunkte (des
Kerns von $S_c-1$) und die der Kommutante seiner Kraus-Operatoren: Für die unitalen Kanäle stimmen sie
überein ([Satz 6.10](../../lecture_notes/de/06_quantum_channels.md#r-6-10)), für die Amplitudendämpfung haben sie gleiche Dimension, sind aber verschiedene
Räume. Dann prüft er [Satz 6.11](../../lecture_notes/de/06_quantum_channels.md#r-6-11) an $100$ zufälligen Zuständen und [Proposition 6.12](../../lecture_notes/de/06_quantum_channels.md#r-6-12) für die
Abbildung `Psi` aus Abschnitt 6.3 und für den Teilchenverlust.

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

Unter vollständiger Amplitudendämpfung sinkt die Entropie des maximal gemischten Zustands von
$M_2\oplus M_3$ von $\log5\approx1.609$ auf $H(\tfrac25,\tfrac35)\approx0.673$: Ein nicht unitaler
Kanal kann reinigen.

Auf $M_2$ kommutieren die Dephasierung $D_p$ und die Amplitudendämpfung $A_\gamma$, und ihre
Hintereinanderausführung multipliziert die Kohärenzen mit $(1-p)\sqrt{1-\gamma}$; mit $1-p=e^{-t/T_\varphi}$ und
$1-\gamma=e^{-t/T_1}$ ergibt das $1/T_2=1/(2T_1)+1/T_\varphi$ ([Vorlesungsnotizen, Abschnitt 6.4](../../lecture_notes/de/06_quantum_channels.md#s-6-4),
*Kohärenzen eines Qubits*). Der Code prüft diese Identitäten und die Kompositionsgesetze
$D_{p(t)}D_{p(s)}=D_{p(t+s)}$, $A_{\gamma(t)}A_{\gamma(s)}=A_{\gamma(t+s)}$:

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

## 6.5 Zufällige Kanäle

*Vorlesungsnotizen:* [Abschnitt 6.5](../../lecture_notes/de/06_quantum_channels.md#s-6-5): [Lemma 6.13](../../lecture_notes/de/06_quantum_channels.md#r-6-13) und [Satz 6.15](../../lecture_notes/de/06_quantum_channels.md#r-6-15) (Ausgangsreinheit eines zufälligen Kanals).

Eine **Haar-zufällige Isometrie** $V:\mathbb C^k\to\mathbb C^N$ ist, wie in der Bibliothek
(`random_isometry`), der Faktor $Q$ von $G=QR$, wobei $G$ eine $N\times k$-Ginibre-Matrix ist
(unabhängige Einträge mit der Dichte $\pi^{-1}e^{-|z|^2}$), $Q^*Q=1$ und $R$ eine obere
Dreiecksmatrix mit positiver Diagonale (auf der CPU nimmt die Bibliothek eine LAPACK-Zerlegung
$G=Q'R'$ und multipliziert die $i$-te Spalte von $Q'$ mit $R'_{ii}/|R'_{ii}|$; das ergibt diese
Zerlegung, wie in Satz 5.3, sobald alle $R'_{ii}\ne0$ sind, d. h. sobald $G$ Rang $k$ hat, was nach
dem Beweis von Lemma 6.13 fast sicher gilt; für große Stapel auf einer GPU verwendet sie
Gram–Schmidt, das sie direkt liefert). `random_channel(alg, r)` zieht solche
$V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}$ unabhängig in jedem Block und verwendet die
Kraus-Operatoren aus Satz 6.5.

Nach [Satz 6.15 der Vorlesungsnotizen](../../lecture_notes/de/06_quantum_channels.md#r-6-15) gilt für einen reinen Zustand $\psi$ in einem Sektor der Größe $k$
und einen von `random_channel` mit Kraus-Rang $r$ gezogenen Kanal $\Phi$: $\mathbb E\mathrm{Tr}\,\Phi(\psi)^2=(k+r)/(kr+1)$ und
$\mathbb E\,\Phi(\psi)=1_k/k$. Wir schätzen beides aus $20\,000$ Stichproben, mit dem Standardfehler des Mittelwerts:

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

## 6.6 Womit man vergleichen sollte

QuTiP (`to_choi`, `to_kraus`, `to_super`, `to_stinespring`), Qiskit (`qiskit.quantum_info`: die
Klassen `Kraus`, `Choi`, `SuperOp`, `Stinespring` und `random_quantum_channel`; Rauschmodelle wie
`amplitude_damping_error` in Qiskit Aer), toqito (`toqito.channel_ops`, `toqito.channel_props`) und
QETLAB (`ChoiMatrix`, `KrausOperators`, `ApplyMap`, `RandomSuperoperator`) rechnen zwischen den
Darstellungen der Abschnitte 6.1–6.2 um, testen vollständige Positivität und ziehen (die meisten von
ihnen) zufällige Kanäle, und zwar für Abbildungen $M_n\to M_m$ zwischen vollen Matrixalgebren.
Ihre Konventionen unterscheiden sich (spalten- oder zeilenweises Stapeln, die Reihenfolge der
Faktoren der Choi-Matrix); man vergleiche also die Konventionen, bevor man Zahlen vergleicht. Ein
sektorerhaltender Kanal auf $\mathcal M$ lässt sich dort als Kanal auf $M_D$, $D=\sum_ck_c$,
behandeln, mit Superoperatoren der Größe $D^2\times D^2$ statt $k_c^2\times k_c^2$ pro Sektor. Die
Bibliothek fügt die Blockstruktur selbst hinzu: nach Sektoren und gestapelt gespeicherte
Kraus-Operatoren, Abbildungen zwischen verschiedenen Algebren (`InterSectorChannel`), die dualen
Abbildungen für jede der Spuren aus Kapitel 3 und die Übergangsmatrix zwischen den Sektoren.

<a id="x-6"></a>

## Aufgaben

[Aufgaben](../../lecture_notes/de/06_quantum_channels.md#x-6) der Vorlesungsnotizen, numerisch:

1. (Werner–Holevo) Man prüfe mit `Channel.from_superoperator`, dass $\Phi(\rho)=\big(\mathrm{Tr}(\rho)1-\rho^T\big)/(k-1)$
   auf $M_3$ ein unitaler Kanal mit Kraus-Rang $3$ ist.
2. Man prüfe, dass der Dephasierungskanal auf $M_2\oplus M_3$ im Block $c$ den Kraus-Rang $k_c$ hat, während die
   Bibliothek $1+k_{\max}$ Kraus-Operatoren speichert.
3. Man prüfe, dass der gemischt-unitäre Kanal mit den Kraus-Operatoren $X^aZ^b/k$ auf $M_3$ gleich `center_expectation` ist.
4. Man prüfe $(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp$ für die drei Spuren der Bibliothek.

<details>
<summary>Lösungen</summary>

*Aufgabe 1.*

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

*Aufgabe 2.*

```python
ranks = [int((torch.linalg.eigvalsh(J_c[0]) > 1e-10).sum()) for J_c in zoo['dephasing'].choi()]
print(ranks, zoo['dephasing'].kraus_rank)                 # [2, 3] versus 4 stored operators
assert ranks == alg.k_factors and zoo['dephasing'].kraus_rank == 1 + max(alg.k_factors)
```

*Aufgabe 3.*

```python
Xs = torch.roll(torch.eye(k, dtype=dt), 1, dims=0)
Zs = torch.diag(torch.exp(2j * math.pi * torch.arange(k, dtype=torch.float64) / k)).to(dt)
W = [torch.linalg.matrix_power(Xs, a) @ torch.linalg.matrix_power(Zs, b) / k for a in range(k) for b in range(k)]
twirl = channels.kraus_channel([M3.operator(Wab[None, None]) for Wab in W])
r3op = states.random_density_matrix(M3)
assert close(twirl(r3op), channels.center_expectation(M3)(r3op))
```

*Aufgabe 4.*

```python
Psi2 = channels.random_inter_sector_channel(alg_out, alg, kraus_rank=2)
for t in trace:
    assert close((Psi2 @ Psi).adjoint(trace=t)(x), Psi.adjoint(trace=t)(Psi2.adjoint(trace=t)(x)))
```

</details>

**Weiter:** [Tutorial 7](07_tensor_products.md): Tensorprodukte, partielle Spuren und Verschränkung in der Bibliothek.
