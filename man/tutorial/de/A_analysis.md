[English](../A_analysis.md) · [Español](../es/A_analysis.md) · [Français](../fr/A_analysis.md) · **Deutsch** · [中文](../zh/A_analysis.md) · [日本語](../ja/A_analysis.md) · [Русский](../ru/A_analysis.md)

# Tutorial, Anhang A. Tatsachen aus Analysis und Wahrscheinlichkeitstheorie

Dieses Tutorial begleitet [Anhang A](../../lecture_notes/de/A_analysis.md) der Vorlesungsnotizen, der die Tatsachen aus Analysis und
Wahrscheinlichkeitstheorie beweist, die die Kapitel mit ihrer Nummer zitieren: Nullstellenmengen von Polynomen,
Eigenwertungleichungen, Gauß-Vektoren und ihre Richtungen, Gaußsche Momente und gleichverteilte Abstände. Die Prüfungen
verwenden nur NumPy, nicht die Bibliothek: Jeder Block zieht zufällige Vektoren oder Matrizen und vergleicht sie mit der
exakten Aussage. Die im ersten Block definierte Hilfsfunktion `mc` stellt sicher, dass ein Monte-Carlo-Mittelwert
höchstens $5$ Standardfehler vom exakten Wert entfernt liegt.

<a id="s-A-2"></a>

## A.2 Nullstellenmengen analytischer Funktionen

*Vorlesungsnotizen:* [Abschnitt A.2](../../lecture_notes/de/A_analysis.md#s-A-2); [Korollar A.7](../../lecture_notes/de/A_analysis.md#r-A-7).

Nach [Korollar A.7](../../lecture_notes/de/A_analysis.md#r-A-7) verschwindet ein nicht identisch verschwindendes Polynom an einem zufälligen Vektor mit
Dichte mit Wahrscheinlichkeit $0$. Für $\det Z$, wobei $Z$ eine reelle $3\times3$-Matrix mit unabhängigen
$N(0,1)$-Einträgen ist, ergibt das $\mathbb P(|\det Z|<\varepsilon)\to0$ für $\varepsilon\to0$.

```python
import math
import numpy as np

rng = np.random.default_rng(0)


def mc(samples, exact, label, kind="exact"):
    """Monte Carlo mean with standard error; the reference value must lie within 5 standard errors."""
    samples = np.asarray(samples, dtype=float)
    mean, err = samples.mean(), samples.std() / math.sqrt(samples.size)
    print(f"{label}: {mean:.4f} +- {err:.4f}   ({kind} {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12, label


# Corollary A.7: det Z = 0 has probability 0 for a real Gaussian 3 x 3 matrix, so by continuity of the measure
# P(|det Z| < eps) -> 0 as eps -> 0. Numerical observation (not a theorem of this appendix): the fractions
# below shrink by about a factor 10 per factor 10 in eps.
d = np.abs(np.linalg.det(rng.standard_normal((200_000, 3, 3))))
fractions = [(d < eps).mean() for eps in (1e-1, 1e-2, 1e-3)]
print("P(|det Z| < eps) for eps = 0.1, 0.01, 0.001:", fractions)
assert fractions[1] < 0.2 * fractions[0] and fractions[2] < 0.2 * fractions[1]
```

<a id="s-A-3"></a>

## A.3 Eigenwerte hermitescher Matrizen

*Vorlesungsnotizen:* [Abschnitt A.3](../../lecture_notes/de/A_analysis.md#s-A-3); [Satz A.8](../../lecture_notes/de/A_analysis.md#r-A-8), [Korollar A.9](../../lecture_notes/de/A_analysis.md#r-A-9), [Satz A.10](../../lecture_notes/de/A_analysis.md#r-A-10), [Lemma A.12](../../lecture_notes/de/A_analysis.md#r-A-12).

Der Block prüft an zufälligen hermiteschen Matrizen die Weylschen Ungleichungen
$\lambda_j(A)+\lambda_1(E)\le\lambda_j(A+E)\le\lambda_j(A)+\lambda_n(E)$ ([Satz A.10](../../lecture_notes/de/A_analysis.md#r-A-10)), die
Courant–Fischer-Schranke $\max_{x\in U,|x|=1}\langle x,Ax\rangle\ge\lambda_j(A)$ für $j$-dimensionales $U$
([Satz A.8](../../lecture_notes/de/A_analysis.md#r-A-8)), die Verschränkung $\lambda_j(A)\le\lambda_j(Q^*AQ)\le\lambda_{j+n-m}(A)$
([Korollar A.9](../../lecture_notes/de/A_analysis.md#r-A-9)) und die Lagrange-Formel für die Spektralprojektionen ([Lemma A.12](../../lecture_notes/de/A_analysis.md#r-A-12)).

```python
def herm(n, batch):
    G = rng.standard_normal((batch, n, n)) + 1j * rng.standard_normal((batch, n, n))
    return (G + G.conj().transpose(0, 2, 1)) / 2


n = 6
A, E = herm(n, 2000), 0.1 * herm(n, 2000)
lam_A, lam_AE, lam_E = np.linalg.eigvalsh(A), np.linalg.eigvalsh(A + E), np.linalg.eigvalsh(E)   # ascending
normE = np.linalg.norm(E, ord=2, axis=(1, 2))                                    # operator norm
assert (np.abs(lam_AE - lam_A).max(axis=1) <= normE + 1e-12).all()               # Theorem A.10
assert (lam_AE >= lam_A + lam_E[:, :1] - 1e-12).all() and (lam_AE <= lam_A + lam_E[:, -1:] + 1e-12).all()

a0, j = A[0], 3                          # Theorem A.8: every 3-dim subspace has max <x, Ax> >= lambda_3
for _ in range(200):
    Q, _ = np.linalg.qr(rng.standard_normal((n, j)) + 1j * rng.standard_normal((n, j)))
    assert np.linalg.eigvalsh(Q.conj().T @ a0 @ Q)[-1] >= lam_A[0, j - 1] - 1e-12
m = 4                                    # Corollary A.9: interlacing for a compression to 4 dimensions
Q, _ = np.linalg.qr(rng.standard_normal((n, m)) + 1j * rng.standard_normal((n, m)))
theta = np.linalg.eigvalsh(Q.conj().T @ a0 @ Q)
assert (theta >= lam_A[0, :m] - 1e-12).all() and (theta <= lam_A[0, n - m:] + 1e-12).all()

U, _ = np.linalg.qr(rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
mu, mult = np.array([-1.0, 0.5, 2.0]), [2, 3, 1]          # Lemma A.12 with multiplicities 2, 3, 1
spec = np.repeat(mu, mult)
A1 = U @ np.diag(spec) @ U.conj().T
for i, m_i in enumerate(mult):
    P = np.eye(n, dtype=complex)
    for l in range(len(mu)):
        if l != i:
            P = P @ (A1 - mu[l] * np.eye(n)) / (mu[i] - mu[l])
    Vi = U[:, spec == mu[i]]
    assert np.abs(P - Vi @ Vi.conj().T).max() < 1e-10 and abs(np.trace(P).real - m_i) < 1e-10
```

<a id="s-A-4"></a>

## A.4 Gauß-Vektoren

*Vorlesungsnotizen:* [Abschnitt A.4](../../lecture_notes/de/A_analysis.md#s-A-4); [Satz A.15](../../lecture_notes/de/A_analysis.md#r-A-15), [Proposition A.17](../../lecture_notes/de/A_analysis.md#r-A-17).

[Satz A.15](../../lecture_notes/de/A_analysis.md#r-A-15) liefert die Dichte $f_X$ eines nicht ausgearteten zentrierten Gauß-Vektors $X$; der Block
vergleicht $\mathbb P(X\in[0,1]^2)$ mit dem nach der Mittelpunktsregel berechneten Integral von $f_X$. Nach
[Proposition A.17](../../lecture_notes/de/A_analysis.md#r-A-17) ist $\mathbb E\mathrm{Tr}\,H^2$ gleich $n(n+1)/2$ für das GOE und $n^2$ für das GUE,
und $OHO^T$ hat die Verteilung der GOE-Matrix $H$; der Block prüft zwei Momente ihrer Einträge.

```python
Amat = rng.standard_normal((2, 4))                       # Theorem A.15: X = A g in R^2, A of rank 2
C = Amat @ Amat.T
X = rng.standard_normal((400_000, 4)) @ Amat.T
h = 1 / 400                                              # midpoint rule for the density on [0, 1]^2
grid = (np.arange(400) + 0.5) * h
xx = np.stack(np.meshgrid(grid, grid, indexing='ij'), -1).reshape(-1, 2)
dens = np.exp(-0.5 * np.einsum('ni,ij,nj->n', xx, np.linalg.inv(C), xx)) / (2 * math.pi * math.sqrt(np.linalg.det(C)))
mc(((X >= 0) & (X <= 1)).all(axis=1), dens.sum() * h * h, "P(X in [0,1]^2) against the integral of f_X",
   kind="midpoint rule")

nn, B = 4, 100_000                                       # Proposition A.17
G = rng.standard_normal((B, nn, nn))
goe = (G + G.transpose(0, 2, 1)) / 2
Gc = (rng.standard_normal((B, nn, nn)) + 1j * rng.standard_normal((B, nn, nn))) / math.sqrt(2)
gue = (Gc + Gc.conj().transpose(0, 2, 1)) / math.sqrt(2)
mc(np.einsum('bij,bji->b', goe, goe), nn * (nn + 1) / 2, "E Tr H^2, GOE")      # n + 2 * n(n-1)/2 * 1/2
mc(np.einsum('bij,bji->b', gue, gue).real, nn ** 2, "E Tr H^2, GUE")           # n + 2 * n(n-1)/2 * 1
O, _ = np.linalg.qr(rng.standard_normal((nn, nn)))
rot = O @ goe @ O.T                                                             # has the law of the GOE
mc(rot[:, 0, 0] ** 4, 3.0, "E (O H O^T)_11^4")
mc(rot[:, 0, 1] ** 2, 0.5, "E (O H O^T)_12^2")
```

<a id="s-A-5"></a>

## A.5 Die Richtung eines Gauß-Vektors

*Vorlesungsnotizen:* [Abschnitt A.5](../../lecture_notes/de/A_analysis.md#s-A-5); [Korollar A.19](../../lecture_notes/de/A_analysis.md#r-A-19), [Satz A.20](../../lecture_notes/de/A_analysis.md#r-A-20), [Korollar A.21](../../lecture_notes/de/A_analysis.md#r-A-21).

Für einen komplexen Standard-Gauß-Vektor $z\in\mathbb C^3$ besagt [Satz A.20](../../lecture_notes/de/A_analysis.md#r-A-20)(b), dass $|z|^2$ die
Dichte $t^{N-1}e^{-t}/(N-1)!$ hat und $z/|z|$ gleichverteilt auf der Sphäre und unabhängig von $|z|$ ist;
[Korollar A.21](../../lecture_notes/de/A_analysis.md#r-A-21) liefert die Momente eines gleichverteilten Einheitsvektors. Der Block prüft dies, das Volumen
$V_3=4\pi/3$ der Einheitskugel ([Korollar A.19](../../lecture_notes/de/A_analysis.md#r-A-19)) und dass die Phase von $z_1$ gleichverteilt ist
([Abschnitt A.5](../../lecture_notes/de/A_analysis.md#s-A-5), nach Satz A.20).

```python
N, B = 3, 400_000                                   # Theorem A.20 (b), Corollary A.21, in C^3
z = (rng.standard_normal((B, N)) + 1j * rng.standard_normal((B, N))) / math.sqrt(2)
R2 = (np.abs(z) ** 2).sum(axis=1)
u = z / np.sqrt(R2)[:, None]
t = 2.0
p_gamma = 1 - math.exp(-t) * sum(t ** k / math.factorial(k) for k in range(N))   # P(|z|^2 <= t) for Gamma(N)
mc(R2 <= t, p_gamma, "P(|z|^2 <= 2)")
mc(R2 ** 2, N * (N + 1), "E|z|^4")
mc(np.abs(u[:, 0]) ** 4, 2 / (N * (N + 1)), "E|u_1|^4")
mc(np.abs(u[:, 0] * u[:, 1]) ** 2, 1 / (N * (N + 1)), "E|u_1|^2 |u_2|^2")
mc((R2 <= t) * np.abs(u[:, 0]) ** 4, p_gamma * 2 / (N * (N + 1)), "E[1{|z|^2 <= 2} |u_1|^4], independence")
mc((u[:, 0] * u[:, 1].conj()).real, 0.0, "E Re u_1 conj(u_2)")

x3 = rng.uniform(-1, 1, (B, 3))                     # Corollary A.19: V_3 = s_3 / 3 = 4 pi / 3
mc((x3 ** 2).sum(axis=1) <= 1, math.pi / 6, "P(uniform point of [-1,1]^3 in the unit ball)")
g3 = rng.standard_normal((B, 3))
u3 = g3 / np.linalg.norm(g3, axis=1, keepdims=True)
mc(u3[:, 0] ** 4, 3 / 15, "E u_1^4 on S^2")
phase = np.angle(z[:, 0])                           # sigma_1 is normalised arc length
mc((phase >= 0) & (phase < 1), 1 / (2 * math.pi), "P(arg z_1 in [0, 1))")
```

<a id="s-A-6"></a>

## A.6 Gaußsche Momente: die Isserlis-Wick-Formel

*Vorlesungsnotizen:* [Abschnitt A.6](../../lecture_notes/de/A_analysis.md#s-A-6); [Satz A.24](../../lecture_notes/de/A_analysis.md#r-A-24).

Die Isserlis-Wick-Formel ([Satz A.24](../../lecture_notes/de/A_analysis.md#r-A-24)) schreibt ein Moment eines zentrierten Gauß-Vektors als Summe
über die $(2p-1)!!$ Paarungen von Produkten von Kovarianzen. Der Block zählt die Paarungen auf, vergleicht die Formel
mit Monte-Carlo-Momenten eines Gauß-Vektors in $\mathbb R^3$ und gewinnt $\mathbb E|z|^4=2$ für eine komplexe
Standard-Gauß-Variable $z$ aus der komplexen Form des Satzes zurück.

```python
def pairings(idx):
    """All pairings of the list idx, as lists of pairs."""
    if not idx:
        yield []
        return
    for i in range(1, len(idx)):
        for rest in pairings(idx[1:i] + idx[i + 1:]):
            yield [(idx[0], idx[i])] + rest


def wick(cov, idx):                                       # sum over pairings of the positions of idx
    return sum(math.prod(cov[idx[a], idx[b]] for a, b in p) for p in pairings(list(range(len(idx)))))


assert [len(list(pairings(list(range(2 * p))))) for p in range(1, 5)] == [1, 3, 15, 105]    # (2p-1)!!
Amat = rng.standard_normal((3, 3)) / math.sqrt(3)
cov = Amat @ Amat.T
Y = rng.standard_normal((1_000_000, 3)) @ Amat.T
for idx in ((0, 1, 2, 0), (0, 0, 0, 0), (0, 1, 1, 2, 2, 0)):
    mc(np.prod(Y[:, list(idx)], axis=1), wick(cov, idx), f"E Y{idx}")
mc(np.prod(Y[:, [0, 1, 2]], axis=1), 0.0, "E Y(0, 1, 2), odd")
M = np.array([[1, 1j], [1, -1j]])                  # (z, conj z) = M (x, y) with x, y independent N(0, 1/2)
cz = M @ M.T / 2                                     # E[Z_a Z_b] = [[0, 1], [1, 0]]
assert abs(wick(cz, (0, 1, 0, 1)) - 2) < 1e-12                           # E|z|^4 = 2
```

<a id="s-A-7"></a>

## A.7 Exponentialverteilte Variablen und gleichverteilte Abstände

*Vorlesungsnotizen:* [Abschnitt A.7](../../lecture_notes/de/A_analysis.md#s-A-7); [Satz A.25](../../lecture_notes/de/A_analysis.md#r-A-25), [Korollar A.26](../../lecture_notes/de/A_analysis.md#r-A-26).

Nach [Korollar A.26](../../lecture_notes/de/A_analysis.md#r-A-26) haben die $n+1$ Abstände von $n$ sortierten gleichverteilten Punkten die Verteilung
der normierten exponentialverteilten Variablen $(e_0,\dots,e_n)/T$, mit den Momenten aus Teil (c),
$\mathbb P(s_i/s_j\le t)=t/(1+t)$ und dem mittleren Abstandsverhältnis $2\ln2-1$. Der Block prüft dies für beide
Familien mit $n=5$ sowie die Unabhängigkeit von $(e_0,\dots,e_n)/T$ und $T$ ([Satz A.25](../../lecture_notes/de/A_analysis.md#r-A-25)) über
$\mathbb E[w_0T^2]=\mathbb Ew_0\,\mathbb ET^2=n+2$.

```python
n, B = 5, 200_000                                       # Corollary A.26
levels = np.sort(rng.uniform(0, 1, (B, n)), axis=1)
s = np.diff(np.concatenate([np.zeros((B, 1)), levels, np.ones((B, 1))], axis=1), axis=1)    # s_0, ..., s_n
e = rng.exponential(1.0, (B, n + 1))
w = e / e.sum(axis=1, keepdims=True)                    # normalised exponentials
for name, sp in (("uniform spacings", s), ("normalised exponentials", w)):
    mc(sp[:, 2], 1 / (n + 1), f"{name}: E s_2")
    mc(sp[:, 2] ** 2, 2 / ((n + 1) * (n + 2)), f"{name}: E s_2^2")
    mc(sp[:, 0] * sp[:, n], 1 / ((n + 1) * (n + 2)), f"{name}: E s_0 s_n")
    for tt in (0.5, 2.0):
        mc(sp[:, 1] / sp[:, 2] <= tt, tt / (1 + tt), f"{name}: P(s_1 / s_2 <= {tt})")
r = np.minimum(s[:, 1:-2], s[:, 2:-1]) / np.maximum(s[:, 1:-2], s[:, 2:-1])    # ratios of interior spacings
mc(r.mean(axis=1), 2 * math.log(2) - 1, "<r>, uniform levels")
T = e.sum(axis=1)
mc(w[:, 0] * T ** 2, n + 2, "E[w_0 T^2] = E w_0 E T^2 (independence)")
```

**Weiter:** [Tutorial, Anhang B. Schur-Weyl-Dualität und der Weingarten-Kalkül](B_schur_weyl.md).
