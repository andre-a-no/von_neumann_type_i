[English](../08_dynamics.md) · [Español](../es/08_dynamics.md) · [Français](../fr/08_dynamics.md) · **Deutsch** · [中文](../zh/08_dynamics.md) · [日本語](../ja/08_dynamics.md) · [Русский](../ru/08_dynamics.md)

# Tutorial 8. Dynamik

Dieses Tutorial begleitet [Kapitel 8](../../lecture_notes/de/08_dynamics.md) der Vorlesungsnotizen, das beweist, dass ein Hamiltonoperator
in $`\mathcal M`$ die Sektorgewichte erhält ([Satz 8.3](../../lecture_notes/de/08_dynamics.md#r-8-3)) und dass Generatoren in GKSL-Form
quantendynamische Halbgruppen erzeugen ([Satz 8.6](../../lecture_notes/de/08_dynamics.md#r-8-6)), und das Sprünge zwischen Sektoren, Postselektion und
stationäre Zustände behandelt. Wir verwenden die exakten Löser und die Runge-Kutta-Löser von `dynamics` für die
Schrödinger-, die von-Neumann- und die Lindblad-Gleichung und prüfen die Ergebnisse gegen SciPy und gegen
GKSL-Generatoren auf $`M_D`$.

<a id="s-8-1"></a>

## 8.1 Abgeschlossene Systeme

*Vorlesungsnotizen:* [Abschnitt 8.1](../../lecture_notes/de/08_dynamics.md#s-8-1); [Proposition 8.2](../../lecture_notes/de/08_dynamics.md#r-8-2) (der Propagator $`U(t)=e^{-itH}`$), [Satz 8.3](../../lecture_notes/de/08_dynamics.md#r-8-3) (Ladungserhaltung).

Nach Proposition 8.2(1) zerfällt die Zeitentwicklung in $`C`$ unabhängige Probleme der Größen $`k_c`$,
eines pro Block. Die Bibliothek berechnet $`U(t)`$ mit einer Eigenzerlegung pro Blockgröße (`dynamics.propagator`), und
`dynamics.von_neumann` entwickelt eine Dichte in der Eigenbasis zu vielen Zeitpunkten auf einmal,
wobei $`\rho(t)_{ij}=e^{-it(\lambda_i-\lambda_j)}\rho_{ij}`$. Wir prüfen gegen `scipy.linalg.expm`:

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

Ein Einheitsvektor $`\psi=\bigoplus_c\psi_c`$ definiert auf $`\mathcal M`$ den Zustand $`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$, der
sich mit $`U(t)`$ entwickelt ([Vorlesungsnotizen, Abschnitt 8.1](../../lecture_notes/de/08_dynamics.md#s-8-1)); daher müssen `dynamics.schrodinger` (Vektoren im Blockformat `(batch, C, k_max)`)
und `von_neumann` darauf übereinstimmen. Für ein zeitabhängiges $`H(t)`$ verwendet die Bibliothek das
klassische Runge-Kutta-Verfahren mit fester Schrittweite (`dynamics.rk4`, `schrodinger_rk4`): Ein
Schritt multipliziert $`\psi`$ mit einem Polynom in Werten von $`H(t)\in\mathcal M`$, verschiebt also
niemals Amplitude zwischen Sektoren, ist aber nicht unitär: Für konstantes $`H`$ multipliziert ein
Schritt der Weite $`h`$ einen Eigenvektor von $`H`$ zum Eigenwert $`\lambda`$ mit $`T(-ih\lambda)`$,
$`T(z)=\sum_{k=0}^4z^k/k!`$, und $`|T(ix)|^2=1-\frac{x^6}{72}+\frac{x^8}{576}<1`$ für $`0<|x|<2\sqrt2`$
(ausmultiplizieren). RK4 erhält die Sektorgewichte also nicht exakt: Für konstantes $`H`$ verliert eine
Sektorkomponente mit einem von null verschiedenen Koeffizienten bei einem Eigenvektor mit
$`0<|h\lambda|<2\sqrt2`$ in jedem Schritt Gewicht. Wie stark sie sich in einem Lauf ändern, ist eine
numerische Beobachtung (unten). (Für einen einzelnen großen Sektor benötigt
`krylov.evolve` nur Matrix-Vektor-Produkte; siehe [Tutorial 9, Abschnitt 9.5](../../tutorial/de/09_spin_chains.md#s-9-5).)

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

Numerische Beobachtung: In diesem Lauf ändern sich die Sektorgewichte unter RK4 um $`1.1\cdot10^{-9}`$
(der ausgegebene Wert) und unter den exakten Lösern, für die Satz 8.3 exakte Erhaltung liefert, um
weniger als die Toleranz $`10^{-12}`$ der `assert`-Anweisungen.

<a id="s-8-2"></a>

## 8.2 Quantendynamische Halbgruppen

*Vorlesungsnotizen:* [Abschnitt 8.2](../../lecture_notes/de/08_dynamics.md#s-8-2); [Satz 8.6](../../lecture_notes/de/08_dynamics.md#r-8-6) (GKSL-Generatoren und die Dyson-Entwicklung), [Satz 8.7](../../lecture_notes/de/08_dynamics.md#r-8-7).

Die Bibliothek hat vier Lindblad-Routinen: `lindblad_rhs` gibt die rechte Seite als Funktion für
`rk4` zurück; `lindblad_evolve` integriert mit RK4; `lindblad_superoperator` gibt für jeden Sektor
die $`k_c^2\times k_c^2`$-Matrix von $`\mathcal L`$ auf dem zeilenweise gebildeten
$`\mathrm{vec}\,\rho_c`$ zurück; und `lindblad_channel` berechnet $`e^{t\mathcal L}`$ exakt, Sektor
für Sektor, als `Channel`. Wir prüfen Satz 8.6 für beliebige Sprünge auf $`M_3`$: Die Choi-Matrix von
$`e^{t\mathcal L}`$ ist positiv (Satz 6.2), die Abbildung ist TP, das Halbgruppengesetz gilt, RK4
stimmt mit der Exponentialfunktion überein, und der Sprunganteil
$`e^{t\mathcal L}-S_0(t)=\sum_{n\ge1}S_n(t)`$ ist selbst CP:

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

## 8.3 Lindblad-Gleichungen in einer Algebra mit Sektoren

*Vorlesungsnotizen:* [Abschnitt 8.3](../../lecture_notes/de/08_dynamics.md#s-8-3): der Generator $`\mathcal L`$ mit Sprüngen innerhalb der Sektoren und zwischen ihnen, den `dynamics.lindblad_rhs` implementiert; [Proposition 8.8](../../lecture_notes/de/08_dynamics.md#r-8-8), [Proposition 8.9](../../lecture_notes/de/08_dynamics.md#r-8-9) (Sektorgewichte).

**Beispiel: Teilchenverlust** ([Vorlesungsnotizen, Abschnitt 8.3](../../lecture_notes/de/08_dynamics.md#s-8-3)). Zwei bosonische Moden mit höchstens zwei
Teilchen, $`\mathcal M=M_1\oplus M_2\oplus M_3`$ (Sektoren $`N=0,1,2`$), der Bose-Hubbard-Dimer
$`H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)`$, und jede Mode verliert Teilchen mit der Rate $`\gamma`$: ein
`InterSectorChannel` mit den Kraus-Operatoren $`a_1,a_2`$, eingeschränkt auf $`N\to N-1`$. Nach
[Proposition 8.9(2)](../../lecture_notes/de/08_dynamics.md#r-8-9) genügen die Sektorgewichte $`p'=Qp`$ mit $`Q=R-\mathrm{diag}(\sum_dR_{dc})_c`$ für die Rate $`1`$,
$`R`$ = `transition_matrix()`, und ausgehend von $`N=2`$ folgen sie der Binomialverteilung $`p_2=q^2`$, $`p_1=2q(1-q)`$, $`p_0=(1-q)^2`$,
$`q=e^{-\gamma t}`$, was auch immer das Tunneln (Hopping) und die Wechselwirkung bewirken. Der Code prüft
$`\Phi^*(1)=0\oplus1\oplus2`$, d. h. $`a_1^*a_1+a_2^*a_2=N`$ auf den verwendeten Blöcken, und $`p(t)`$ gegen beide Formeln.

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

RK4 erhält die Positivität nicht garantiert: Für $`\rho'=\mathcal L\rho`$ ist ein Schritt der Weite
$`h`$ die Abbildung $`\rho\mapsto T(h\mathcal L)\rho`$ mit dem Polynom $`T`$ aus Abschnitt 8.1 (man setze
die vier Stufen ein), die spurerhaltend ist (da $`\mathrm{Tr}\,\mathcal L(X)=0`$ für alle $`X`$), aber
nicht positiv sein muss. Zum Beispiel bildet $`\mathcal L`$ auf $`M_2`$ mit der Basis
$`|g\rangle,|e\rangle`$, $`H=0`$ und einem Sprung $`|g\rangle\langle e|`$ mit der Rate $`\gamma`$ die Matrix
$`\mathrm{diag}(a,b)`$ auf $`\gamma b\mathrm{diag}(1,-1)`$ ab, also bildet ein Schritt mit
$`h\gamma=3`$ die Matrix $`\mathrm{diag}(0,1)`$ auf
$`\mathrm{diag}(1-T(-3),T(-3))=\mathrm{diag}(-\frac38,\frac{11}8)`$ ab. Numerische
Beobachtung: In diesem Lauf ist der kleinste Eigenwert von $`\rho(t)`$ gleich $`\approx-2\cdot10^{-8}`$
(der ausgegebene Wert, mit Schrittweiten bis $`0.02`$). `lindblad_channel` ist exakt, behandelt aber nur
Sprünge innerhalb der Sektoren. Schließlich prüfen wir Proposition 8.8(2) für eine zufällige
Abbildung zwischen den Sektoren und Proposition 8.8(3) für den Teilchenverlust, im Vergleich mit
GKSL-Generatoren auf $`M_6`$, die mit numpy aufgebaut werden:

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

## 8.4 Postselektion und die sprungfreie Entwicklung

*Vorlesungsnotizen:* [Abschnitt 8.4](../../lecture_notes/de/08_dynamics.md#s-8-4); [Proposition 8.10](../../lecture_notes/de/08_dynamics.md#r-8-10) (die sprungfreie Wahrscheinlichkeit $`w_0`$), [Proposition 8.11](../../lecture_notes/de/08_dynamics.md#r-8-11) (Postselektion auf „kein Teilchen verloren“).

Die Funktion `states.lueders_update(rho, P)` der Bibliothek gibt
$`P\rho P/\mathrm{Tr}(P\rho P)`$ und $`\mathrm{Tr}(\rho P)`$ zurück. Nun mögen die beiden
Moden mit verschiedenen Raten $`\gamma_1>\gamma_2`$ Teilchen verlieren; dann ist
$`A=\gamma_1n_1+\gamma_2n_2`$ nicht zentral (im Sektor $`N=2`$ ist es
$`\mathrm{diag}(2\gamma_1,\gamma_1+\gamma_2,2\gamma_2)`$ in der Basis
$`|20\rangle,|11\rangle,|02\rangle`$), Proposition 8.9(2) ist nicht anwendbar, und der postselektierte
Zustand auf $`N=2`$ ist die normierte Entwicklung mit $`H_{\mathrm{eff}}`$, dessen antihermitescher Anteil
$`-\frac i2A`$ in diesem Sektor kein Vielfaches von $`1`$ ist. Wir prüfen Proposition 8.11
für $`n=0`$ und – indem wir das Integral für $`S_1`$ mit der Simpson-Regel berechnen – für $`n=1`$:

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

Numerische Beobachtung: Bei $`t=2`$ beträgt das Gewicht $`p_2(2)=w_0(2)`$ des Sektors $`N=2`$ („kein
Verlust“) $`0.230`$ (oben ausgegeben). Ohne Tunneln ($`J=0`$) wäre es exakt
$`e^{-2\gamma_1t}=e^{-2.4}=0.091`$, da dann $`|20\rangle`$ ein Eigenvektor von $`H_{\mathrm{eff}}`$ zum
Eigenwert $`U-i\gamma_1`$ ist und $`w_0(t)=|e^{-it(U-i\gamma_1)}|^2`$ gilt (Proposition 8.10); hier hängen
die Sektorgewichte also sehr wohl vom Hamiltonoperator ab. Der postselektierte Zustand hat vom
Zustand des abgeschlossenen Systems den Spurabstand $`0.135`$ (oben ausgegeben): Das Bedingen auf $`N=2`$
ändert den Zustand, über die unitäre Entwicklung hinaus.

<a id="s-8-5"></a>

## 8.5 Stationäre Zustände und Relaxation

*Vorlesungsnotizen:* [Abschnitt 8.5](../../lecture_notes/de/08_dynamics.md#s-8-5); [Korollar 8.13](../../lecture_notes/de/08_dynamics.md#r-8-13), [Proposition 8.14](../../lecture_notes/de/08_dynamics.md#r-8-14) (spontane Emission).

Nach [Proposition 8.14](../../lecture_notes/de/08_dynamics.md#r-8-14) hat die spontane Emission auf $`M_2`$ ($`H=\omega|e\rangle\langle e|`$, ein Sprung
$`|g\rangle\langle e|`$ mit der Rate $`\gamma`$) $`e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}}`$, das Spektrum
$`\{0,-\gamma,-\gamma/2\pm i\omega\}`$ und die Spektrallücke $`\gamma/2`$. Wir prüfen das Spektrum mit `lindblad_superoperator`, den Kanal
und den Spurabstand mit `lindblad_channel` sowie
[Korollar 8.13](../../lecture_notes/de/08_dynamics.md#r-8-13) für zufällige Sprünge innerhalb der Sektoren von $`M_2\oplus M_3`$:

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

Numerische Beobachtung: Für diese zufälligen Sprünge findet der Code genau einen stationären Zustand
pro Sektor. Nach Satz D.6, sektorweise angewandt ([Vorlesungsnotizen, Abschnitt 8.5](../../lecture_notes/de/08_dynamics.md#s-8-5)), gilt dies, wenn jeder Sektor
einen treuen stationären Zustand hat und die einzigen Matrizen, die mit $`H_c`$ und allen $`L_{j,c}`$, $`L_{j,c}^*`$ kommutieren,
die Vielfachen von $`1`$ sind; der Code prüft diese beiden Voraussetzungen (in Gleitkommaarithmetik, für die oben
verwendeten Raten $`g_j=1`$) und dann die Folgerung.

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

## 8.6 Womit man vergleichen sollte

QuTiP (`sesolve`, `mesolve` für die Lindblad-Gleichung mit adaptiven ODE-Lösern, `mcsolve` für
Quantentrajektorien, `steadystate`, `liouvillian`), dynamiqs (JAX; gebatchte und differenzierbare
Schrödinger- und Lindblad-Löser auf GPUs), QuantumOptics.jl und SciPy (`scipy.linalg.expm` sowie
`scipy.sparse.linalg.expm_multiply` für $`e^{t\mathcal L}\rho`$ mit einem dünnbesetzten vektorisierten
Generator) lösen dieselben Gleichungen auf einem einzelnen Hilbertraum, und für die obigen Probleme
müssen sie dieselben Zahlen liefern. Auf $`\mathcal M`$ übergibt man ihnen blockdiagonale
$`D\times D`$-Matrizen; der vektorisierte Generator lebt dann auf einem $`D^2`$-dimensionalen Raum
(dünnbesetzt, aber mit allen Sektoren in einer gemeinsamen Rechnung gekoppelt), sofern der Benutzer
das Problem nicht von Hand in Sektoren zerlegt. Die Bibliothek hat die Blockstruktur eingebaut:
exakte Propagatoren und Lindblad-Kanäle Sektor für Sektor ($`k_c^2\times k_c^2`$-Superoperatoren),
Sprünge zwischen Sektoren als `InterSectorChannel`s mit der Kompression aus Proposition 8.8, Batches
von Hamiltonoperatoren und Sprungoperatoren (für Mittelungen über Unordnung) und automatische
Differentiation (die letzten beiden auch in dynamiqs). Sie hat keine adaptive Schrittweitensteuerung
(RK4 arbeitet mit fester Schrittweite und ist nur näherungsweise positiv), keinen Trajektorienlöser,
und `lindblad_channel` akzeptiert keine Sprünge zwischen Sektoren; dafür sind QuTiP oder dynamiqs
auf $`M_D`$ die Referenz. Neben RK4 hat `lindblad_evolve` zwei Zeitschrittverfahren, deren Schritte für jede Schrittweite Kanäle sind
(Abschnitt 8.7): das Splitting zweiter Ordnung `method='kraus'` und `method='exact'`, das auch Sprünge zwischen Sektoren
behandelt, über die Matrix von $`\mathcal L`$ auf $`\mathcal M`$ (für kleine Sektoren).

<a id="s-8-7"></a>

## 8.7 Zeitschritte, die Kanäle sind

*Vorlesungsnotizen:* [Abschnitt 8.7](../../lecture_notes/de/08_dynamics.md#s-8-7); [Proposition 8.15](../../lecture_notes/de/08_dynamics.md#r-8-15) (RK4 ist nicht positiv), [Satz 8.16](../../lecture_notes/de/08_dynamics.md#r-8-16) (ein Splitting-Verfahren, dessen Schritte Kanäle sind).

`lindblad_evolve` hat drei Zeitschrittverfahren, die mit `method` gewählt werden: `'rk4'` (die unveränderte Voreinstellung)
ist der Schritt $`R_h=T_4(h\mathcal L)`$; `'kraus'` ist das Splitting $`\tilde\Psi_h`$ aus Satz 8.16, gebaut aus
$`E_h=e^{-ihH_{\mathrm{eff}}/2}`$, dem Sprunganteil $`\Phi_{\mathrm{tot}}`$ (Sprünge innerhalb der Sektoren und zwischen ihnen)
und der Korrektur $`M_h^{-1/2}`$; `'exact'` ist $`e^{h\mathcal L}`$, berechnet aus der Matrix von $`\mathcal L`$ auf
$`\bigoplus_cM_{k_c}`$ (aufgefüllt auf $`C k_{\max}^2`$ Koordinaten; für einen konstanten Generator eine Exponentialfunktion pro
Intervalllänge, und `substeps` wird nicht verwendet). Zuerst prüfen wir Proposition 8.15 am Zweiniveausystem mit
$`\omega=\gamma=1`$ und $`h=1`$: Ein RK4-Schritt ist $`T_4(h\mathcal L)`$ mit der Matrix aus `lindblad_superoperator`,
$`\tau-|a|^2<0`$, und der reine Zustand $`\cos\frac\theta2|g\rangle+\sin\frac\theta2|e\rangle`$, $`\theta=0.2746`$, wird auf
eine Matrix mit negativem Eigenwert abgebildet; die anderen beiden Verfahren liefern Zustände.

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

Numerische Beobachtung: Hier ist $`\tau-|a|^2=-0.0091`$, und der kleinste Eigenwert nach einem RK4-Schritt ist
$`-8.5\cdot10^{-5}`$ (die ausgegebenen Werte); $`\theta=0.2746`$ ist auf vier Stellen der Winkel, der ihn minimiert.
Proposition 8.15 handelt von einem Qubit; in einer Algebra mit Sektoren zeigt sich derselbe Mangel in den Choi-Matrizen der
Sektorkomponenten $`R_{dc}`$, die genau dann alle positiv sind, wenn der Schritt CP ist (Proposition 6.4). Die Funktion
`step_map` unten rekonstruiert einen Schritt aus `lindblad_evolve` auf hermiteschen Matrizen (die Ausgabe wird
selbstadjungiert gemacht, daher wird der Schritt auf $`E_{ij}+E_{ji}`$ und $`i(E_{ij}-E_{ji})`$ ausgewertet und linear
fortgesetzt). Für den Teilchenverlust aus Abschnitt 8.3 prüfen wir Satz 8.16(1) und Proposition 8.8(1): Die Schritte von
`'kraus'` und `'exact'` sind CP und TP für Schrittweiten von $`0.1`$ bis $`1000`$, während der RK4-Schritt schon für
$`h=0.1`$ nicht CP ist.

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

Numerische Beobachtung: Für $`h=0.1`$ ist der kleinste Eigenwert einer Choi-Matrix des RK4-Schritts $`-9.7\cdot10^{-5}`$;
für $`h=3`$, wo einige Eigenwerte $`z`$ von $`h\mathcal L`$ die Ungleichung $`|T_4(z)|>1`$ erfüllen (RK4 ist instabil), ist er
$`-9.6\cdot10^{2}`$ (oben ausgegeben); die anderen beiden Verfahren bleiben bis auf Rundung Kanäle. Für $`h=0.1`$ erfüllen
alle Eigenwerte $`|T_4(z)|\le1`$: Der Verlust der Positivität ist keine Instabilität. Die Renormierung mit $`M_h^{-1/2}`$ wird
aus den Eigenwerten von $`M_h`$ berechnet, deren Konditionszahl exponentiell mit $`h`$ wachsen kann, wenn sich die
Zerfallsraten in einem Sektor unterscheiden; übersteigt sie $`\varepsilon^{-1/4}`$ ($`\varepsilon`$ die Maschinengenauigkeit),
so führt die Bibliothek den Schritt rekursiv als zwei Schritte der Weite $`h/2`$ aus, so dass der berechnete Schritt bis auf
Rundung spurerhaltend bleibt.

Schließlich die Genauigkeit. `'exact'` stimmt mit $`e^{t\hat{\mathcal L}}`$ des Generators $`\hat{\mathcal L}`$ auf $`M_6`$
mit den Sprüngen $`\sqrt\gamma a_1,\sqrt\gamma a_2`$ überein (Proposition 8.8(1) und (3), mit SciPy aus den Matrizen
`L_loss` aus Abschnitt 8.3 berechnet), und der Spurabstand der anderen beiden Verfahren dazu fällt wie $`h^2`$ (Satz 8.16(3))
beziehungsweise wie $`h^4`$:

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

Numerische Beobachtung: Halbieren von $`h`$ teilt den Fehler von `'kraus'` durch etwa $`4`$ und den von RK4 durch etwa
$`16`$ (die ausgegebenen Verhältnisse): Das Splitting erkauft Positivität für jede Schrittweite mit einer niedrigeren
Ordnung. Für einen zeitabhängigen Hamiltonoperator $`H`$ oder zeitabhängige Sprünge frieren `'kraus'` und `'exact'` den
Generator in der Mitte jedes Schritts ein, so dass jeder Schritt nach Satz 8.16(1) ein Kanal bleibt; die Ordnung dieses
Verfahrens mit eingefrorenem Mittelpunkt wird in den Vorlesungsnotizen nicht bewiesen.

<a id="x-8"></a>

## Aufgaben

[Aufgaben 2 und 4](../../lecture_notes/de/08_dynamics.md#x-8) der Vorlesungsnotizen: die Eigenwerte
$`-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}`$, $`\delta=(\gamma_1-\gamma_2)/4`$, von $`H_{\mathrm{eff}}`$ im Sektor $`N=1`$ aus
Abschnitt 8.4, mit einem einzigen Jordan-Block der Größe $`2`$ bei $`J=|\delta|`$; und das Spektrum $`\{0,0,-2\gamma\pm i\omega\}`$ der
reinen Dephasierung auf $`M_2`$ ($`L=\sigma_z`$). Man überprüfe beides numerisch.

<details>
<summary>Lösungen</summary>

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

**Weiter:** [Tutorial 9. Spinketten](09_spin_chains.md).
