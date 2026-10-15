[English](../11_optimisation.md) · [Español](../es/11_optimisation.md) · [Français](../fr/11_optimisation.md) · **Deutsch** · [中文](../zh/11_optimisation.md) · [日本語](../ja/11_optimisation.md) · [Русский](../ru/11_optimisation.md)

# Tutorial 11. Optimierung und Spurungleichungen

Dieses Tutorial begleitet [Kapitel 11](../../lecture_notes/de/11_optimisation.md) der Vorlesungsnotizen, das die Extrema dreier Spurfunktionale beweist
und erklärt, warum Stichproben unitärer Matrizen sie nicht finden. Hier werden die Parametrisierungen von
Nebenbedingungsmengen (`optimize.PositiveParam`, `UnitaryParam`, `SelfAdjointParam`) und der Multistart-Optimierer
(`optimize.extremize`) der Bibliothek an diesen Sätzen geprüft. Mit dem Michelson-Kontrast $`\Delta`$ und
$`\ell(\delta)=(1-\delta)/(1+\delta)`$ lauten die drei Funktionale:

| | Funktional | Variablen |
|---|---|---|
| Exp. 1 | $`z=\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)`$ | $`X,Y\ge0`$, $`U\in\mathcal M`$ unitär |
| Exp. 2 | $`z=\mathrm{Tr}\lvert XY\rvert-\mathrm{Tr}(\lvert X\rvert\,\lvert Y\rvert)`$ | $`X=UX_0`$, $`Y=VY_0`$, $`X_0,Y_0\ge0`$, $`U,V`$ unitär |
| Exp. 3 | $`z=\mathrm{Tr}(Y\lvert X\rvert Y)-\mathrm{Tr}\lvert YXY\rvert`$ | $`X=X^*`$, $`Y\ge0`$ |

wobei $`|X|`$, $`|Y|`$ auf $`\lambda_{\max}=1`$ normiert sind und ihre Kontraste als Parameter dienen. Die
**Einhüllenden** sind $`\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}`$ und
$`\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}`$, wobei alle übrigen Variablen frei sind.

<a id="s-11-1"></a>

## 11.1 Nebenbedingungsmengen und ihre Parametrisierungen

*Vorlesungsnotizen:* [Abschnitt 11.1](../../lecture_notes/de/11_optimisation.md#s-11-1); [Lemma 11.1](../../lecture_notes/de/11_optimisation.md#r-11-1), [Proposition 11.2](../../lecture_notes/de/11_optimisation.md#r-11-2).

`optimize.PositiveParam` mit vorgegebenem Kontrast $`\delta`$ implementiert die Parametrisierung aus
[Abschnitt 11.1](../../lecture_notes/de/11_optimisation.md#s-11-1) der Vorlesungsnotizen: Pro Start zieht es zwei verschiedene Plätze $`s_\pm`$ (`spectrum.pin_max`,
`spectrum.pin_min`; jedes geordnete Paar hat die Wahrscheinlichkeit $`1/(D(D-1))`$) und eine Haar-verteilte unitäre Matrix
$`U_0`$ (`unitary.base`, gezogen von `random_unitary_operator`, [Lemma 11.10](../../lecture_notes/de/11_optimisation.md#r-11-10)) und bildet die Parameter $`\theta`$
(`spectrum.theta`) und $`G`$ (`unitary.G`) ab auf
```math
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell(\delta),\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
```
`UnitaryParam` ist die Abbildung $`G\mapsto U_0e^{G-G^*}`$ allein, und `SelfAdjointParam` fügt Vorzeichen hinzu,
$`X=U\mathrm{diag}(\sigma\lambda)U^*`$. Nach [Proposition 11.2](../../lecture_notes/de/11_optimisation.md#r-11-2) ist das Bild die Menge der Operatoren der
Kontrastschale $`\mathcal P_\delta`$ mit $`1`$ und $`\ell(\delta)`$ auf den gezogenen Plätzen. Ein freier Kontrast in $`[a,b]`$ ist ein
weiterer Parameter; die Voreinstellung der Bibliothek ist $`[a,b]=[0,0.999]`$, und der Code unten behält sie bei.

Der Code prüft die Nebenbedingung und, konstruktiv, Proposition 11.2: Die Parameter aus dem Beweis reproduzieren
eine Zielmatrix.

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

## 11.2 Die von Neumannsche Spurungleichung

*Vorlesungsnotizen:* [Abschnitt 11.2](../../lecture_notes/de/11_optimisation.md#s-11-2); [Satz 11.4](../../lecture_notes/de/11_optimisation.md#r-11-4), [Korollar 11.5](../../lecture_notes/de/11_optimisation.md#r-11-5), [Proposition 11.6](../../lecture_notes/de/11_optimisation.md#r-11-6), [Korollar 11.7](../../lecture_notes/de/11_optimisation.md#r-11-7).

Nach [Vorlesungsnotizen, Satz 11.4](../../lecture_notes/de/11_optimisation.md#r-11-4) (von Neumannsche Spurungleichung) ist
$`|\mathrm{Tr}(AB)|\le\sum_is_i(A)s_i(B)`$, und nach [Korollar 11.5](../../lecture_notes/de/11_optimisation.md#r-11-5) gilt $`\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i`$ für selbstadjungierte $`P,Q`$ mit
Eigenwerten $`p,q`$ sowie $`\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i`$ für $`P,Q\ge0`$. Für Exp. 2 bei festen Spektren liefert
[Proposition 11.6](../../lecture_notes/de/11_optimisation.md#r-11-6) $`\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}(y^\downarrow_{c,i}-y^\uparrow_{c,i})`$, angenommen für $`X_0=\mathrm{diag}(x^\downarrow)`$, $`V=J`$ die umkehrende Permutation und $`Y_0=\mathrm{diag}(y^\uparrow)`$
(bzw. $`\mathrm{diag}(y^\downarrow)`$), und [Korollar 11.7](../../lecture_notes/de/11_optimisation.md#r-11-7) liefert die Einhüllende $`\sup z=-\inf z=\sum_c\lfloor k_c/2\rfloor(1-\ell_X)(1-\ell_Y)`$.

Hier ist $`1-\ell_X=2\delta/(1+\delta)`$ für $`\Delta(|X|)=\delta`$; für ungerades $`k_c`$ spielt der mittlere Eigenwert keine
Rolle. Der Code prüft Satz 11.4, Korollar 11.5 und die extremalen Konfigurationen aus Proposition 11.6 für $`k=3`$,
$`C=2`$.

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

## 11.3 Eine Spurungleichung mit Gleichheitsbedingung

*Vorlesungsnotizen:* [Abschnitt 11.3](../../lecture_notes/de/11_optimisation.md#s-11-3); [Lemma 11.8](../../lecture_notes/de/11_optimisation.md#r-11-8), [Proposition 11.9](../../lecture_notes/de/11_optimisation.md#r-11-9).

Nach [Vorlesungsnotizen, Proposition 11.9](../../lecture_notes/de/11_optimisation.md#r-11-9) ist $`z=\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\ge0`$ in Exp. 3, mit Gleichheit genau dann,
wenn $`YX_+Y^2X_-Y=0`$, insbesondere wenn $`XY=YX`$. Für $`X=\mathrm{diag}(1,-1,0)`$ und $`Y`$ mit den Zeilen $`(1,0,a)`$, $`(0,1,0)`$, $`(a,0,1)`$, $`-1<a<1`$, gilt Gleichheit,
obwohl $`XY\ne YX`$ für $`a\ne0`$. Der Code prüft beides (mit $`a=0.5`$); `X.abs()` ist $`|X|`$, und `trace_norm()` ist
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

## 11.4 Stichproben oder Optimierung

*Vorlesungsnotizen:* [Abschnitt 11.4](../../lecture_notes/de/11_optimisation.md#s-11-4); [Lemma 11.10](../../lecture_notes/de/11_optimisation.md#r-11-10), [Proposition 11.11](../../lecture_notes/de/11_optimisation.md#r-11-11), [Korollar 11.12](../../lecture_notes/de/11_optimisation.md#r-11-12), [Lemma 11.13](../../lecture_notes/de/11_optimisation.md#r-11-13).

Wie groß ist $`|\mathrm{Tr}(AU)|`$ typischerweise, verglichen mit seinem Maximum $`\|A\|_1`$? Die Bibliothek zieht $`U`$
als den Faktor $`Q`$ der QR-Zerlegung $`Z=QR`$, $`R_{ii}>0`$, einer **Ginibre-Matrix** $`Z`$ (unabhängige Einträge mit der
Dichte $`\pi^{-1}e^{-|z|^2}`$; Mezzadris Rezept), und zwar in jedem Block unabhängig.

Nach [Vorlesungsnotizen, Proposition 11.11](../../lecture_notes/de/11_optimisation.md#r-11-11) ist $`\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k`$ für Haar-verteiltes $`U\in U(k)`$, und nach
[Korollar 11.12](../../lecture_notes/de/11_optimisation.md#r-11-12) ist für $`X,Y\ge0`$ mit Spektren in $`[a,b]`$, $`a>0`$, das quadratische Mittel von $`|\mathrm{Tr}(XUY)|`$
höchstens das $`\min(1/\sqrt k,b^2/(a^2k))`$-Fache seines Maximums über $`U`$, und $`\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le b^4/(a^4k^2)`$.

$`N`$ Haar-Stichproben finden eine Verletzung von $`|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)`$ also höchstens mit
Wahrscheinlichkeit $`Nb^4/(a^4k^2)`$ (eine Vereinigung von $`N`$ Ereignissen), obwohl es für jedes nichtkommutierende
Paar eine gibt (Proposition 3.9). Die Schranke ist erst für $`k>\sqrt N\,b^2/a^2`$ aussagekräftig (für $`b/a=6`$ und
$`N=2000`$: $`k>1609`$); numerisch übersehen die Stichproben die Verletzung schon für viel kleinere $`k`$ (siehe unten).
Der Code prüft Proposition 11.11 per Monte-Carlo (der Stichprobenmittelwert liegt bei diesem Lauf mit festem Seed
innerhalb von fünf Standardfehlern vom exakten Wert) und gibt für in $`[0.2,1.2]`$ gleichverteilte Spektren ($`b/a=6`$) den
Quotienten aus Korollar 11.12, $`\sup z`$ und das beste $`z`$ aus $`2000`$ Stichproben aus, beide relativ zu
$`\mathrm{Tr}(XY)`$.

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

*Numerische Beobachtung* (die ausgegebene Tabelle; nur für diesen Lauf mit festem Seed geprüft): Der Quotient bleibt
unter der bewiesenen Schranke $`\min(1/\sqrt k,b^2/(a^2k))`$ und halbiert sich ungefähr, wann immer sich $`k`$ verdoppelt,
d. h. er fällt wie $`1/k`$, und ab $`k=4`$ erreicht die beste von $`2000`$ Stichproben nicht einmal $`z=0`$. Eine
Gradientensuche kommt dagegen im Lauf unten beiden Extrema nahe (wiederum eine numerische Beobachtung). Für das Paar
`Xf, Yf` sind beide Extrema über $`U`$ bekannt:

Nach [Vorlesungsnotizen, Lemma 11.13](../../lecture_notes/de/11_optimisation.md#r-11-13) ist $`\min_U|\mathrm{Tr}(AU)|=0`$ für jedes $`A\in M_k(\mathbb C)`$, $`k\ge2`$.

Also ist $`\inf_Uz=-\mathrm{Tr}(XY)`$ und $`\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)`$ (Proposition 3.9). Ein einziger
Aufruf von `extremize` nähert beide an, mit $`12`$ Starts, die $`z`$ maximieren, und $`12`$, die $`-z`$ maximieren (in diesem
Lauf bis auf eine relative Genauigkeit von $`10^{-3}`$, unten geprüft). Nach Proposition 11.11 ist das quadratische
Mittel von $`|\mathrm{Tr}(XUY)|`$ gleich $`\|YX\|_2/\sqrt k`$; mit den Schranken aus dem Beweis von Korollar 11.12
(Spektren in $`[a,b]`$) ist $`\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4`$ und $`\mathrm{Tr}(XY)\ge ka^2`$, sodass die
Markov-Ungleichung $`\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)`$ für
$`\varepsilon>0`$ liefert: Für großes $`k`$ liegt $`z`$ mit hoher Wahrscheinlichkeit nahe bei seinem Infimum.
*Numerische Beobachtung* (für diesen Lauf mit festem Seed geprüft): Stichproben kommen dem Infimum bis auf $`1\,\%`$
nahe, erreichen aber nicht einmal $`z=0`$, während das Supremum positiv ist:

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

## 11.5 Multistart-Suche nach Einhüllenden

*Vorlesungsnotizen:* [Abschnitt 11.5](../../lecture_notes/de/11_optimisation.md#s-11-5); [Proposition 11.14](../../lecture_notes/de/11_optimisation.md#r-11-14), [Korollar 11.7](../../lecture_notes/de/11_optimisation.md#r-11-7), [Proposition 11.9](../../lecture_notes/de/11_optimisation.md#r-11-9), [Proposition 11.15](../../lecture_notes/de/11_optimisation.md#r-11-15), [Korollar 11.16](../../lecture_notes/de/11_optimisation.md#r-11-16).

`optimize.extremize(objective, params, maximize)` führt Adam in `rounds` Runden zu je `steps` Schritten aus
(Lernrate `lr * decay**round`) auf dem Verlust $`\mp\sum_b\text{objective}_b`$ aus, zentriert die unitären Karten nach
jeder Runde neu und gibt **pro Start** den größten (für `maximize=False` den kleinsten) Wert der Zielfunktion unter
denen zurück, die an den Gradienteniterierten und einmal am Endpunkt berechnet wurden. Jedes Batch-Element ist ein
unabhängiger Start mit eigenen Plätzen, Vorzeichen, $`U_0`$ und eigenem Kontrast: Adam arbeitet koordinatenweise; hängt
also der $`b`$-te Eintrag der Zielfunktion nur von den $`b`$-ten Parametern ab (wie bei allen Zielfunktionen unten), so
beeinflussen sich die Starts nicht. Die Vorzeichen von `SelfAdjointParam` sind diskret: Zwischen den Runden probiert
`discrete_step` pro Start acht einzelne Vorzeichenwechsel (jeweils an einem gleichverteilt zufälligen Platz)
nacheinander aus und behält jeden Wechsel, der die Zielfunktion verbessert.

In exakter Arithmetik ist jede Iterierte zulässig (Proposition 11.2 und die Bemerkungen danach), also ist jeder
zurückgegebene Wert der Wert von $`z`$ an einem zulässigen Punkt; in Gleitkommaarithmetik ist er dieser Wert bis auf
Rundungsfehler (die Prüfungen unten erlauben $`10^{-12}`$). Daher liefert eine Maximierung eine **untere Schranke** für
das Supremum, eine Minimierung eine obere Schranke für das Infimum, und nicht mehr – eine lokale Suche kann in einem
lokalen Extremum stehen bleiben. Die andere Richtung erfordert einen Beweis, wie Korollar 11.7 und
Proposition 11.9, oder eine konvexe Relaxierung (Abschnitt 11.6). Eine berechnete Einhüllende wird geprüft (i) an
geschlossenen Formeln; (ii) an von Hand gebauten Konfigurationen, etwa direkten Summen kleinerer Lösungen, unter die
sie nicht fallen darf; (iii) an strukturierten Stichproben (Zwei-Niveau-Spektren, Permutationsbasen), deren Werte
weitere angenommene Werte sind.

Nach [Vorlesungsnotizen, Proposition 11.14](../../lecture_notes/de/11_optimisation.md#r-11-14) ist, wenn alle $`k_c\ge2`$ sind, die untere Einhüllende von Exp. 1 gleich $`\underline z(\delta)=-\big(D-2\delta/(1+\delta)\big)`$.

Der Code prüft Proposition 11.14 auf $`M_2\oplus M_2`$, Korollar 11.7 auf $`M_3\oplus M_3`$ (beide Kontraste
vorgegeben, $`k=3`$ ungerade) und Proposition 11.9; Suprema und Infima teilen sich wie oben einen Aufruf. Die
einseitigen Prüfungen (kein zurückgegebener Wert jenseits eines bewiesenen Extremums) sind Folgerungen aus diesen
Sätzen; die Prüfungen, dass die Suche dem Extremum bis auf die angegebene Toleranz nahekommt, sind numerische
Beobachtungen dieses Laufs mit festem Seed. Das Schlüsselwort `is_self_adjoint` wählt den Algorithmus für die
Spurnorm (Singulärwerte oder, für $`YXY`$, Eigenwerte), ohne in jedem Schritt zu testen. Die letzte Prüfung ist vom
Typ (ii): $`z`$ aus Exp. 3 ist additiv über direkte Summen, und der Kontrast einer direkten Summe ist der größte
Kontrast ihrer Summanden (Aufgabe 3), also gilt für die Einhüllenden von Exp. 3 mit $`\Delta(Y)`$ frei in $`[0,0.999]`$
wie im Code
```math
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
```
insbesondere ist die Summe zweier auf $`M_2`$ bei den Kontrasten $`\delta`$ und $`\delta'\le\delta`$ angenommener Werte ein
auf $`M_2\oplus M_2`$ beim Kontrast $`\delta`$ angenommener Wert.

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

*Numerische Beobachtungen* (die ausgegebenen Werte dieses Laufs). Die Suche auf $`M_2\oplus M_2`$ bleibt um etwa
$`0.2\,\%`$ hinter der direkten Summe zurück: Sie ist nicht vollständig konvergiert, und die beste untere Schranke für
$`\overline z_{M_2\oplus M_2}(0.5)`$ ist die größere Zahl – dies zu erkennen ist der Zweck von Prüfung (ii). Die Werte
für $`M_2`$ stimmen bis auf $`10^{-4}`$ mit $`(1-\delta)\big(2d/(1+d)\big)^2=(1-\delta)(1-\ell(d))^2`$ überein, wobei
$`d=0.999`$ der größte zugelassene Kontrast von $`Y`$ ist. Diese geschlossene Form, zuerst numerisch in der Begleitarbeit
beobachtet, ist die obere Einhüllende aus [Korollar 11.16](../../lecture_notes/de/11_optimisation.md#r-11-16) (unten); der
Code prüft nur die bewiesene Richtung (kein Wert darüber). Nach demselben Korollar ist die direkte Summe optimal:
$`\overline z_{M_2\oplus M_2}(0.5)=1.5\,(2d/(1+d))^2`$. Für größere Blöcke kennen wir keine geschlossene Form der oberen
Einhüllenden von Exp. 1 und Exp. 3; berechnete Kurven sind dann untere Schranken (angenommene Werte), nicht mehr.

Nach [Vorlesungsnotizen, Proposition 11.15](../../lecture_notes/de/11_optimisation.md#r-11-15) ist auf $`M_2`$ für Spektren
$`x_1\ge x_2`$ von $`X`$ (von $`|X|`$ in Exp. 3) und $`y_1\ge y_2`$ von $`Y`$ das Supremum $`\sup z`$ gleich
$`(x_1-x_2)^2(y_1-y_2)^2/(4(x_1+x_2)(y_1+y_2))`$ in Exp. 1 und $`\frac{2x_1x_2}{x_1+x_2}(y_1-y_2)^2`$ in Exp. 3. Ist
$`e_1`$ ein Eigenvektor von $`X`$ zum Eigenwert vom Betrag $`x_1`$ und $`f_1`$ einer von $`Y`$ zu $`y_1`$, so wird das
Supremum angenommen, wenn $`\gamma=|\langle e_1,f_1\rangle|^2`$ gleich $`\gamma_1=\frac14+\frac S{2q}`$ (Exp. 1, mit
$`|\mathrm{Tr}(XUY)|=\|YX\|_1`$) bzw. $`\gamma_3=\frac Sq`$ (Exp. 3, $`X=\mathrm{diag}(x_1,-x_2)`$) ist, wobei $`S=x_1y_2+x_2y_1`$
und $`q=(x_1+x_2)(y_1+y_2)`$. Nach
[Korollar 11.16](../../lecture_notes/de/11_optimisation.md#r-11-16) hat für $`C`$ Blöcke der Größe $`2`$ mit $`\Delta(|X|)=\delta`$ und
$`\Delta(Y)=\delta'`$ Exp. 1 die Werte $`\sup z=C\delta^2\delta'^2/((1+\delta)(1+\delta'))`$ und
$`\inf z=-(2C-1+\ell(\delta)\ell(\delta'))`$ und Exp. 3 den Wert $`\sup z=(C-\delta)\big(2\delta'/(1+\delta')\big)^2`$. Beide
Suprema wachsen mit $`\delta'`$: Mit $`\Delta(Y)`$ frei in $`[0,d]`$ sind die oberen Einhüllenden die Werte bei $`\delta'=d`$,
und für $`d\to1`$ streben sie gegen $`C\delta^2/(2(1+\delta))`$ und $`C-\delta`$, die Werte für ein singuläres $`Y`$
($`\Delta(Y)=1`$).

Der Code baut die reellen Maximierer aus dem Beweis für $`C=1,2`$ (`from_blocks`), prüft die Kontraste und die Werte,
prüft die Schranken aus Proposition 11.15 an zufälligen komplexen Paaren, erreicht die Grenzwerte $`d\to1`$ mit
Rang-1-Projektionen in jedem Block und führt eine Multistart-Suche für Exp. 1 auf $`M_2`$ entlang $`\Delta(X)`$ mit
$`\Delta(Y)`$ frei in $`[0,0.999]`$ aus, die unter der geschlossenen Form (einem Satz) bleiben muss und ihr in diesem Lauf
bis auf $`10^{-3}`$ nahekommt.

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

## 11.6 Womit man vergleichen sollte

**Manopt**, **pymanopt** und **geoopt** stellen Matrixmannigfaltigkeiten (Stiefel-Mannigfaltigkeiten, orthogonale
und unitäre Gruppen) mit Retraktionen wie QR, polar oder Cayley bereit sowie Riemannsche Löser (konjugierte
Gradienten, Trust-Region-Verfahren, Riemannsches Adam), eine Alternative zu Adam in einer neu zentrierten Karte (einen
Leistungsvergleich stellen wir hier nicht an). **cvxpy** übergibt konvexe Probleme in Operatorvariablen an SDP-Löser,
die ein Optimum zusammen mit einem dualen Zertifikat liefern, beides bis auf die Lösertoleranz; eine Relaxierung – ein
konvexes Problem über einer größeren zulässigen Menge – beschränkt das Supremum nach oben, in der Richtung, die eine
lokale Suche nicht liefern kann. **NCTSSOS**, **ncpol2sdpa** und **NCSOStools** beschränken Eigenwerte und
Spuren nichtkommutativer Polynome in allen Dimensionen zugleich (Abschnitt 0.3). Die Bibliothek ergänzt exakte
Parametrisierungen der Spektral- und Kontrastnebenbedingungen von $`\bigoplus_cM_{k_c}`$ (die Extrema können in jedem
Sektor liegen), unitäre Matrizen, die die Blöcke respektieren, und eine einzige gebündelte Rechnung für viele Starts
und Werte der Nebenbedingung, auf CPU oder GPU. Sie liefert Kandidaten für Extremalstellen und untere Schranken für
Suprema; den Rest müssen Beweise oder Relaxierungen beisteuern.

<a id="x-11"></a>

## Aufgaben

[Aufgaben 1 und 2](../../lecture_notes/de/11_optimisation.md#x-11) der Vorlesungsnotizen, numerisch: Die Cayley-Transformation $`C(S)=(1-S)(1+S)^{-1}`$ eines
schiefadjungierten $`S`$ ist unitär, hat keinen Eigenwert $`-1`$, und es gilt $`S=(1-C)(1+C)^{-1}`$ (Aufgabe 1); für Haar-verteiltes $`O`$
auf $`O(k)`$, gezogen mit `random_unitary_operator` einer reellen Algebra (`complex_valued=False`), gilt
$`\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k`$ (Aufgabe 2).

<details>
<summary>Lösungen</summary>

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

**Weiter:** [Tutorial 12](12_equivariant_learning.md): äquivariante Schichten als Elemente einer Kommutante, in der Bibliothek.
