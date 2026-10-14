[English](../09_spin_chains.md) · [Español](../es/09_spin_chains.md) · [Français](../fr/09_spin_chains.md) · **Deutsch** · [中文](../zh/09_spin_chains.md) · [日本語](../ja/09_spin_chains.md) · [Русский](../ru/09_spin_chains.md)

# 9. Spinketten

Eine Kette aus $L$ Spins $\frac12$ lebt auf $(\mathbb C^2)^{\otimes L}$, und dichte lineare Algebra in Dimension
$2^L$ wird jenseits von $L\approx14$ teuer (eine dichte Matrix hat $4^L$ Einträge). Symmetrien verschieben diese Grenze: Die Magnetisierung zerlegt den Raum in $L+1$
Sektoren (Kapitel 4), Translationen zerlegen jeden Sektor erneut, und $SU(2)$ verkleinert die Blöcke noch einmal.
Dieses Kapitel konstruiert diese Sektoren für die XXZ-Kette, beweist, was sie über Grundzustände aussagen, prüft
die Bibliothek an freien Fermionen und am Bethe-Ansatz und geht über dichte Matrizen (Krylov-Verfahren) und
saubere Ketten (Unordnung) hinaus.

**Sie lernen**, warum die Algebra einer Kette $\bigoplus_NM_{\binom LN}$ ist und warum der Turm dieser Algebren
das Pascalsche Dreieck als Bratteli-Diagramm und den Index $4\cos^2\frac\pi{2L+4}$ hat, wie man Impulssektoren
abzählt, die Marshallsche Vorzeichenregel und ihre Folgen für Impuls und Spin von Grundzuständen (Lieb–Mattis),
die Jordan-Wigner-Lösung der XX-Kette, die Bethe-Gleichungen für zwei Magnonen, warum Lanczos monoton konvergiert
und wie man Niveaustatistik und die Verschränkung freier Fermionen misst.

Im ganzen Kapitel sind die Plätze $0,\dots,L-1$, die Basis von $(\mathbb C^2)^{\otimes L}$ besteht aus den
Konfigurationen $|b\rangle=|b_0\dots b_{L-1}\rangle$ mit $b_i=1$ für Spin nach oben (Konvention der Bibliothek: Bit
$i$ eines Labels ist Platz $i$), $s^a_i$ ($a=x,y,z$) sind die Spinmatrizen $\frac12\sigma_a$ auf Platz $i$,
$s^\pm_i=s^x_i\pm is^y_i$, $n_i=s^z_i+\frac12$, $S^a=\sum_is^a_i$, und $N=S^z+\frac L2$ ist die Anzahl der Spins nach
oben. Eine Kette ist **offen** (Bindungen $(i,i+1)$, $i<L-1$) oder **periodisch** (zusätzlich die Bindung
$(L-1,0)$; dann ist $L\ge3$, und Indizes werden modulo $L$ genommen). $J>0$ (antiferromagnetisch), $\hbar=1$, und
$\log$ ist der natürliche Logarithmus.

## 9.1 Die XXZ-Kette und ihre $U(1)$-Sektoren

Die **XXZ-Kette** mit Anisotropie $\Delta\in\mathbb R$ und Feldern $h_i\in\mathbb R$ ist
$$
H=J\sum_{(i,j)}\Big[\tfrac12\big(s^+_is^-_j+s^-_is^+_j\big)+\Delta\,s^z_is^z_j\Big]+\sum_ih_is^z_i ,
$$
summiert über die Bindungen; $\frac12(s^+_is^-_j+s^-_is^+_j)=s^x_is^x_j+s^y_is^y_j$ bewegt einen Spin nach oben entlang
der Bindung. Für $\Delta=1$ ist dies die Heisenberg-Kette $J\sum\mathbf s_i\cdot\mathbf s_j$, für $\Delta=0$ die
XX-Kette.

**Proposition 9.1.** $[H,S^z]=0$. Folglich ist $H\in U(1)'\cong\bigoplus_{N=0}^LM_{\binom LN}(\mathbb C)$, wobei der
Block $N$ auf dem Spann der Konfigurationen mit $N$ Spins nach oben wirkt.

*Beweis.* Die Diagonalterme kommutieren mit dem diagonalen $S^z$. Für $i\neq j$ kommutiert $s^+_is^-_j$ mit $s^z_l$,
$l\ne i,j$, und $[s^z_i+s^z_j,s^+_is^-_j]=[s^z_i,s^+_i]s^-_j+s^+_i[s^z_j,s^-_j]=s^+_is^-_j-s^+_is^-_j=0$ wegen
$[s^z,s^\pm]=\pm s^\pm$ (Abschnitt 4.5); für $s^-_is^+_j$ geht man zu den Adjungierten über. Folglich
kommutiert $H$ mit $N=S^z+\frac L2$, das die ganzzahligen Eigenwerte $0,\dots,L$ hat, mit Eigenräumen, die von den
Konfigurationen mit $N$ Spins nach oben aufgespannt werden und die Dimension $\binom LN$ haben; der Rest ist
Proposition 4.10(2) mit $Q=N$. $\square$

`SpinChain(L).xxz` baut die Blöcke direkt aus den Konfigurationen jedes Sektors (`chain.labels[N]`). Wir
vergleichen sie mit der vollen $2^L\times2^L$-Matrix, die aus Kronecker-Produkten gebildet wird:

```python
import functools, itertools, math, time
import numpy as np
import torch
from scipy.optimize import brentq
from torch_vn_algebra import SpinChain, SU2Chain, DensityMatrix, inclusions, krylov

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
close = lambda a, b, tol=1e-10: bool((torch.as_tensor(a) - torch.as_tensor(b)).abs().max() < tol)
def op(L, ops):                             # product of one-site operators {site: 2 x 2}; site 0 = lowest bit
    return functools.reduce(torch.kron, [ops.get(j, torch.eye(2, dtype=f64)) for j in reversed(range(L))])

sp = torch.tensor([[0., 0.], [1., 0.]], dtype=f64)                    # s^+ : |0> = down -> |1> = up
sm, sz = sp.T.contiguous(), torch.diag(torch.tensor([-0.5, 0.5], dtype=f64))
L = 6
chain = SpinChain(L, boundary='periodic')
h = 0.3 * torch.randn(L, dtype=f64)                                  # random fields
H_full = sum(0.5 * (op(L, {i: sp, j: sm}) + op(L, {i: sm, j: sp})) + 0.7 * op(L, {i: sz, j: sz})
             for i, j in chain.bonds) + sum(h[i] * op(L, {i: sz}) for i in range(L))   # J = 1, Delta = 0.7
Sz = sum(op(L, {i: sz}) for i in range(L))
assert close(H_full @ Sz, Sz @ H_full)                               # Proposition 9.1
H = chain.xxz(J=1.0, Delta=0.7, h=h)                                 # an element of (+)_N M_binom(6, N)
for N, k in enumerate(chain.algebra.k_factors):
    lab = chain.labels[N]                                            # configurations of sector N
    assert close(H.matrix[0, N, :k, :k], H_full[lab][:, lab])
```

**Der Turm der Ketten.** Sei $\mathcal A_L=\bigoplus_NM_{\binom LN}$, wirkend auf $L$ Plätzen. Das Hinzufügen eines
Platzes, $x\mapsto x\otimes1_2$, bildet $\mathcal A_L$ in $\mathcal A_{L+1}$ ab, denn $x\otimes1$ kommutiert mit
$S^z\otimes1+1\otimes s^z$. Das ergibt einen Turm $\mathcal A_0\subset\mathcal A_1\subset\mathcal A_2\subset\cdots$.

**Proposition 9.2.** Die Inklusion $\mathcal A_L\subset\mathcal A_{L+1}$ hat die $(L+2)\times(L+1)$-Inklusionsmatrix
$\Lambda_L$ mit $\lambda_{N,N}=\lambda_{N,N-1}=1$ und allen übrigen Einträgen $0$; die Blockgrößen gehorchen der
Pascalschen Regel $\binom{L+1}N=\binom LN+\binom L{N-1}$, das Bratteli-Diagramm des Turms ist das Pascalsche Dreieck,
und der bipartite Graph von $\Lambda_L$ ist ein Pfad mit $2L+3$ Ecken.

*Beweis.* Die Konfigurationen mit $N$ Spins nach oben auf $L+1$ Plätzen sind diejenigen mit $N$ Spins nach oben auf
den ersten $L$ Plätzen und dem neuen Spin nach unten sowie diejenigen mit $N-1$ und dem neuen Spin nach oben. Auf
diesen beiden Unterräumen wirkt $x\otimes1$ als $x_N$ bzw. $x_{N-1}$, also enthält Block $N$ von $\mathcal A_{L+1}$
die Blöcke $N$ und $N-1$ von $\mathcal A_L$ je einmal; die Dimensionen addieren sich gemäß der Pascalschen Regel. Im
bipartiten Graphen ist der Sektor $N$ von $\mathcal A_{L+1}$ mit $N-1$ und $N$ von $\mathcal A_L$ verbunden, was den
Pfad $0'-0-1'-1-\cdots-L'-L-(L+1)'$ ergibt (Striche für $\mathcal A_{L+1}$). $\square$

Ein Pfad von der Spitze des Bratteli-Diagramms zum Sektor $N$ auf Stufe $L$ ist eine Konfiguration $b$; davon gibt
es $\binom LN$, die Größe des Blocks. Für eine zusammenhängende Inklusion mit Inklusionsmatrix $\Lambda$
**definieren** wir den **Index** als $\|\Lambda\|^2$ (Operatornorm); für $M_k\otimes1\subset M_k\otimes M_m$ ist die
Inklusionsmatrix die $1\times1$-Matrix $(m)$, also ist der Index das $m^2$ aus Abschnitt 7.2. `Inclusion.index`
berechnet ihn. (Der Name ist der Indextheorie von V. Jones entlehnt, *Invent. Math.* 72 (1983) 1; in diesem Kurs ist
der Index einer Inklusion endlichdimensionaler Algebren diese Definition und nichts weiter.)

**Lemma 9.3.** Die Adjazenzmatrix $A_n$ des Pfades mit den Ecken $1,\dots,n$ hat die Eigenwerte
$2\cos\frac{\pi k}{n+1}$, $k=1,\dots,n$, mit Eigenvektoren $v_k(j)=\sin\frac{\pi jk}{n+1}$; folglich ist
$\|A_n\|=2\cos\frac\pi{n+1}$.

*Beweis.* Man setze $v_k(0)=v_k(n+1)=0$. Dann ist $(A_nv_k)(j)=v_k(j-1)+v_k(j+1)=2\cos\frac{\pi k}{n+1}\,v_k(j)$ wegen
$\sin(a-b)+\sin(a+b)=2\sin a\cos b$, und $v_k\ne0$, da $v_k(1)=\sin\frac{\pi k}{n+1}\neq0$. Die $n$ Eigenwerte sind
verschieden ($\cos$ ist auf $(0,\pi)$ streng monoton fallend), also sind es alle, und die Norm der symmetrischen
Matrix $A_n$ ist der größte Betrag, angenommen bei $k=1$ und $k=n$. $\square$

**Korollar 9.4.** $[\mathcal A_{L+1}:\mathcal A_L]=\|\Lambda_L\|^2=4\cos^2\frac\pi{2L+4}$, wachsend gegen $4$.

*Beweis.* Ordnet man die Ecken entlang des Pfades aus Proposition 9.2, so wird die Adjazenzmatrix
$B=\begin{pmatrix}0&\Lambda_L\\\Lambda_L^T&0\end{pmatrix}$ zu $A_{2L+3}$ (eine Ähnlichkeit mittels Permutation). Nach
Lemma 1.1 ist $\|B\|^2=\|B^2\|$ und $\|\Lambda_L\Lambda_L^T\|=\|\Lambda_L^T\Lambda_L\|=\|\Lambda_L\|^2$, und
$B^2=\Lambda_L\Lambda_L^T\oplus\Lambda_L^T\Lambda_L$, also $\|B\|=\|\Lambda_L\|$; Lemma 9.3 mit $n=2L+3$ liefert den Wert. $\square$

```python
tower = inclusions.chain_tower(8)                                    # A_0 c A_1 c ... c A_8
print(inclusions.bratteli_diagram(tower).levels[:6])                 # Pascal's triangle
for L_, inc in enumerate(tower):
    sizes = [math.comb(L_ + 1, N) for N in range(L_ + 2)]
    assert inc.M.k_factors == sizes and inc.is_connected
    assert list(inc.matrix @ np.array(inc.N.k_factors)) == sizes      # unitality = Pascal's rule
    n = 2 * L_ + 3
    A_path = np.eye(n, k=1) + np.eye(n, k=-1)
    assert abs(np.linalg.eigvalsh(A_path).max() - 2 * math.cos(math.pi / (n + 1))) < 1e-12   # Lemma 9.3
    assert abs(inc.index() - 4 * math.cos(math.pi / (2 * L_ + 4)) ** 2) < 1e-12           # Corollary 9.4
```

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). $S^z$ ist erhalten, sobald die Wechselwirkung invariant unter Drehungen um die
> $z$-Achse ist; $N$ zählt die Spins nach oben, oder die Teilchen nach der Jordan-Wigner-Abbildung aus
> Abschnitt 9.4. Im Turm wird jeder Ladungssektor der längeren Kette von zwei Sektoren der kürzeren gespeist.

## 9.2 Translationen, Impuls und die Marshallsche Vorzeichenregel

Auf einer periodischen Kette sei $T$ die Translation $T|b_0b_1\dots b_{L-1}\rangle=|b_{L-1}b_0\dots b_{L-2}\rangle$, die
Platz $i$ auf Platz $i+1$ verschiebt. $T$ ist eine Permutationsmatrix, also unitär, $T^L=1$ und
$Ts^a_iT^*=s^a_{i+1}$. Folglich ist $TS^zT^*=S^z$ und $THT^*=H$ für ein homogenes Feld $h_i=h$: $T$ bildet Bindungen
auf Bindungen ab. Seine Eigenwerte haben daher die Form $e^{ik}$ mit **Impulsen** $k=2\pi m/L$, $m=0,\dots,L-1$.

**Proposition 9.5 (Impulssektoren).** Sei $N$ fest. Für eine Konfiguration $a$ mit $N$ Spins nach oben sei $R_a$
die Größe ihrer Bahn $\{T^ja\}$. (a) $R_a$ teilt $L$; für $k=2\pi m/L$ mit $mR_a\equiv0\pmod L$ ist der Vektor
$|a,k\rangle=\sum_{j=0}^{R_a-1}e^{-ikj}T^j|a\rangle$ von null verschieden mit $T|a,k\rangle=e^{ik}|a,k\rangle$, und diese $R_a$
Vektoren bilden eine Orthogonalbasis des Spanns der Bahn. (b) Der Eigenwert $e^{2\pi im/L}$ von $T$ auf dem Sektor
$N$ hat die Vielfachheit $d_{N,m}=\#\{\text{Bahnen }:mR_a\equiv0\pmod L\}$, und jeder Operator, der mit $S^z$ und $T$
kommutiert – etwa die translationsinvariante XXZ-Kette –, ist blockdiagonal, liegt also in
$\bigoplus_N\bigoplus_mM_{d_{N,m}}$.

*Beweis.* (a) Der Stabilisator von $a$ in der zyklischen Gruppe $\{T^j\}\cong\mathbb Z_L$ ist eine Untergruppe, und
$R_a$ ist ihr Index, ein Teiler von $L$; die Vektoren $T^j|a\rangle$, $0\le j<R_a$, sind orthonormal, und
$T^{R_a}|a\rangle=|a\rangle$. Ist $e^{ikR_a}=1$, so gilt $T|a,k\rangle=e^{ik}\sum_{j=1}^{R_a}e^{-ikj}T^j|a\rangle=e^{ik}|a,k\rangle$,
da der Term $j=R_a$ gleich dem Term $j=0$ ist; der Koeffizient von $|a\rangle$ ist $1$. Die zulässigen
$m\in\{0,\dots,L-1\}$ sind die Vielfachen von $L/R_a$, insgesamt $R_a$ Stück, mit verschiedenen Eigenwerten des
unitären $T$; also sind die Vektoren orthogonal und spannen den $R_a$-dimensionalen Spann der Bahn auf. (b) Der
Sektor $N$ ist die orthogonale Summe der Spanne seiner Bahnen, jeder davon invariant unter $T$, und nach (a) tritt
der Eigenwert $e^{2\pi im/L}$ im Spann der Bahn von $a$ einmal auf, wenn $mR_a\equiv0\pmod L$, und sonst gar nicht.
Ein $X$, das mit $S^z$ und $T$ kommutiert, kommutiert
mit deren Spektralprojektionen $P_N$, $Q_m$ (Satz 1.5(d)), also ist $X=\sum_{N,m}P_NQ_mXP_NQ_m$. $\square$

`xxz_momentum` baut die Blöcke in der normierten Basis $|a,k\rangle$, mit einem Repräsentanten $a$ pro Bahn. Wir
zählen die Bahnen unabhängig ab und vergleichen die Spektren mit dem dichten Sektor:

```python
L = 12
ring = SpinChain(L, boundary='periodic')
rot = lambda x: ((x << 1) | (x >> (L - 1))) & ((1 << L) - 1)        # T on labels: site i -> i + 1
for N in range(L + 1):
    orbits = {frozenset(itertools.accumulate(range(L - 1), lambda x, _: rot(x), initial=a))
              for a in ring.labels[N].tolist()}                       # {a, Ta, ..., T^{L-1} a}
    d = [sum((m * len(o)) % L == 0 for o in orbits) for m in range(L)]
    assert d == ring.momentum_dims(N)                                   # Proposition 9.5(b)
print(ring.momentum_dims(6))
Hk = ring.xxz_momentum(J=1.0, Delta=1.0, h=0.2, sector=6)               # 12 blocks instead of one of size 924
E_k = torch.sort(torch.cat([e[0] for e in Hk.eigenvalues()])).values
assert close(E_k, torch.linalg.eigvalsh(ring.xxz(1.0, 1.0, h=0.2, sector=6).matrix[0, 0]))
```

Der Sektor $N=6$ mit $924$ Zuständen zerfällt in Blöcke der Größen $80,75,78,76,78,75,80,75,78,76,78,75$.

**Lemma 9.6 (Perron-Frobenius, symmetrischer Fall).** Sei $M$ eine reelle symmetrische Matrix mit $M_{ij}\ge0$ für
$i\ne j$, deren Graph ($i\sim j$ genau dann, wenn $M_{ij}>0$) zusammenhängend ist. Dann ist ihr größter Eigenwert
$\lambda$ einfach, mit einem Eigenvektor, dessen Einträge alle $>0$ sind.

*Beweis.* $\lambda=\max_{\|x\|=1}\langle x,Mx\rangle$, und das Maximum wird genau auf den Einheitsvektoren des
$\lambda$-Eigenraums angenommen (Entwicklung in einer Orthonormalbasis aus Eigenvektoren). Für ein solches $x$ erfüllt
der Vektor $|x|$ der Beträge $\langle|x|,M|x|\rangle\ge\operatorname{Re}\langle x,Mx\rangle=\lambda$, denn die Diagonalterme
stimmen überein und $M_{ij}|x_i||x_j|\ge\operatorname{Re}M_{ij}\bar x_ix_j$ für $i\ne j$; also ist $|x|$ ein Eigenvektor.
Ist $|x|_i=0$, so erzwingt $0=(\lambda-M_{ii})|x|_i=\sum_{j\ne i}M_{ij}|x|_j$, dass $|x|_j=0$ für alle Nachbarn $j$ von
$i$ gilt, und wegen des Zusammenhangs ist $x=0$. Also hat kein von null verschiedener Eigenvektor einen
Nulleintrag; ein Eigenraum der Dimension $\ge2$ enthält aber einen von null verschiedenen Vektor mit $x_1=0$. Daher
ist $\lambda$ einfach, und $|x|>0$ ist ein Eigenvektor. $\square$

Wir beweisen die Vorzeichenregel für jede bipartite Anordnung der Kopplungen; die Kette ist ein Spezialfall.

**Satz 9.7 (Marshallsche Vorzeichenregel).** Die Plätze $\{0,\dots,L-1\}$ seien in zwei Mengen $A\sqcup B$ zerlegt,
und es sei
$$
H=\sum_{i<j}\Big[\tfrac{J_{ij}}2\big(s^+_is^-_j+s^-_is^+_j\big)+K_{ij}\,s^z_is^z_j\Big]+\sum_ih_is^z_i
$$
mit reellen $K_{ij},h_i$ und reellen $J_{ij}$, sodass $J_{ij}\ge0$ gilt, wenn $i,j$ in verschiedenen Mengen liegen,
$J_{ij}\le0$, wenn sie in derselben Menge liegen, und der Graph $G$ auf den Plätzen mit einer Kante $i\sim j$ genau
dann, wenn $J_{ij}\ne0$, zusammenhängend ist. Sei $N_B(b)$ die Anzahl der Spins nach oben von $b$ auf $B$. In jedem
Sektor $N$ ist der kleinste Eigenwert von $H_N$ einfach, und sein Eigenvektor kann so gewählt werden, dass
$(-1)^{N_B(b)}\psi(b)>0$ für jede Konfiguration $b$ mit $N$ Spins nach oben gilt.

Insbesondere gilt dies für die XXZ-Kette mit $J>0$ und beliebigen $\Delta,h_i$, wenn die Kette offen oder periodisch
mit geradem $L$ ist, mit $A$ den geraden und $B$ den ungeraden Plätzen: Dann ist $J_{ij}=J$, $K_{ij}=J\Delta$ auf den
Bindungen und $0$ sonst, jede Bindung verbindet einen geraden mit einem ungeraden Platz (beim Ring, weil $L$ gerade
ist), und $G$ enthält den Pfad $0-1-\cdots-(L-1)$.

*Beweis.* Sei $V|b\rangle=(-1)^{N_B(b)}|b\rangle$; $V$ ist diagonal, $V=V^*=V^{-1}$. Der Operator $s^+_is^-_j$ bewegt einen
Spin nach oben von $j$ nach $i$; er ändert $N_B$ um $\pm1$, wenn $i,j$ in verschiedenen Mengen liegen, und lässt es
sonst unverändert, also ist entsprechend $Vs^+_is^-_jV=\mp s^+_is^-_j$, und ebenso für $s^-_is^+_j$. Zwei verschiedene
Konfigurationen $b,b'$ sind genau dann durch einen solchen Term verbunden, wenn $b'$ aus $b$ durch Vertauschen der
entgegengesetzten Spins auf zwei Plätzen $i,j$ entsteht, und dann sind $i,j$ durch $b,b'$ bestimmt (die Plätze, an
denen sie sich unterscheiden). Daher hat $M=-VH_NV$ die Nichtdiagonaleinträge $J_{ij}/2$ (verschiedene Mengen) bzw.
$-J_{ij}/2$ (dieselbe Menge), also in beiden Fällen $|J_{ij}|/2\ge0$, positiv genau dann, wenn $i\sim j$ in $G$, und
$0$ zwischen Konfigurationen, die sich anders unterscheiden. Der Graph von $M$ ist zusammenhängend. Denn man
schreibe $(i\,j)$ für die Vertauschung der Spins auf den Plätzen $i,j$; für verschiedene $a,b,c$ ist
$(a\,c)=(a\,b)(b\,c)(a\,b)$, also ist durch Induktion entlang eines Pfades in $G$ von $i$ nach $j$ jedes $(i\,j)$ eine
Komposition von Vertauschungen entlang Kanten von $G$, und jede solche Vertauschung fixiert entweder eine
Konfiguration (gleiche Spins) oder läuft entlang einer Kante des Graphen von $M$. Haben $b\ne b'$ je $N$ Spins nach
oben, so gibt es Plätze $i$, oben in $b$ und unten in $b'$, und $j$, unten in $b$ und oben in $b'$; $(i\,j)b$
unterscheidet sich von $b'$ an zwei Plätzen weniger als $b$, also führen endlich viele solche Schritte von $b$ zu
$b'$. Nach Lemma 9.6 ist der größte Eigenwert von $M$ – das Negative des kleinsten von $H_N$ – einfach, mit einem
Eigenvektor $\varphi>0$; man nehme $\psi=V\varphi$. $\square$

**Korollar 9.8 (Impuls des Grundzustands).** Auf einer periodischen Kette mit geradem $L$, $J>0$ und homogenem Feld
hat der Grundzustand des Sektors $N$ den Impuls $k=0$, wenn $N$ gerade ist, und $k=\pi$, wenn $N$ ungerade ist.

*Beweis.* $T$ kommutiert mit $H_N$, und der Grundzustand $\psi=V\varphi$ ist einfach, also gilt $T\psi=e^{ik}\psi$.
Nun ist $(T\psi)(b)=\psi(T^{-1}b)$, und $T^{-1}b$ hat auf Platz $i$ den Spin $b_{i+1}$; da $L$ gerade ist, ist $i+1$
für ungerades $i$ gerade, also $N_B(T^{-1}b)=N-N_B(b)$ und
$(VT\psi)(b)=(-1)^{N_B(b)}(-1)^{N-N_B(b)}\varphi(T^{-1}b)=(-1)^N\varphi(T^{-1}b)$. Außerdem ist $VT\psi=e^{ik}\varphi$.
Sowohl $\varphi$ als auch $\varphi\circ T^{-1}$ sind positiv, also $e^{ik}=(-1)^N$. $\square$

```python
for bc in ('open', 'periodic'):
    c = SpinChain(10, bc, complex_valued=False)
    for N in range(11):
        w, V = torch.linalg.eigh(c.xxz(1.0, 0.7, h=0.3 * torch.randn(10, dtype=f64), sector=N).matrix[0, 0])
        x = (-1.0) ** c.bits(N)[:, 1::2].sum(-1) * V[:, 0]                  # (-1)^{N_B(b)} psi(b)
        assert (x * x.sum().sign()).min() > 0 and (len(w) == 1 or w[1] - w[0] > 1e-8)   # Theorem 9.7
for N in range(1, L):
    lowest = [e[0, 0].item() for e in ring.xxz_momentum(1.0, 0.5, sector=N).eigenvalues()]
    m0 = int(np.argmin(lowest))
    assert m0 == (0 if N % 2 == 0 else L // 2)                              # Corollary 9.8: k = 0 or pi
    assert sorted(lowest)[1] - lowest[m0] > 1e-8
```

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). In einem Antiferromagneten legt
> das Gitter die Vorzeichen des Grundzustands fest: Einen Spin über eine Bindung zu bewegen, ändert das Vorzeichen.
> Physiker nutzen das, um Quanten-Monte-Carlo für bipartite Antiferromagneten ohne Vorzeichenproblem zu rechnen;
> bei frustrierten Ketten, in denen eine antiferromagnetische Kopplung zwei Plätze derselben Menge verbindet, ist
> die Voraussetzung von Satz 9.7 verletzt.

## 9.3 $SU(2)$ bei $\Delta=1$

Für $\Delta=1$ und $h=0$ kommutiert $H=J\sum\mathbf s_i\cdot\mathbf s_j$ mit $S^x,S^y,S^z$ (Proposition 4.15), also
ist $H\cong\bigoplus_S1_{2S+1}\otimes H_S$ mit Blöcken $H_S\in M_{m_S}$, $m_S=\binom L{L/2-S}-\binom L{L/2-S-1}$
(Proposition 4.13), und das Modul `su2` baut $H_S$ auf den Vektoren höchsten Gewichts (Proposition 4.15(b)).
Dasselbe gilt für jedes $H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j$ mit reellen $J_{ij}$: Nach Proposition 4.15(a) liegt
jedes $\mathbf s_i\cdot\mathbf s_j$ in $\{S_a\}'$. Sei $E_0(S)$ der kleinste Eigenwert von $H_S$, für die auftretenden Spins $S$.

**Satz 9.9 (Lieb-Mattis).** Sei $H=\sum_{i<j}J_{ij}\,\mathbf s_i\cdot\mathbf s_j$ mit $A,B,J_{ij}$ wie in Satz 9.7 (also
$K_{ij}=J_{ij}$, $h_i=0$), $n_A=|A|\ge n_B=|B|$ und $S_*=\frac{n_A-n_B}2$.
(a) Für $S_*\le S\le\frac L2$ (mit ganzzahligem $S-\frac L2$) hat der Grundzustand des Sektors $S^z=S$ den Gesamtspin
$S$; für $0\le M\le S_*$ (mit ganzzahligem $M-\frac L2$) hat der Grundzustand des Sektors $S^z=M$ den Gesamtspin $S_*$.
(b) Folglich ist $E_0(S)$ für $S\ge S_*$ ein einfacher Eigenwert von $H_S$, $E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)$,
und $E_0(S)>E_0(S_*)$ für $S<S_*$: Der Grundraum von $H$ ist ein einziges Multiplett, vom Spin $S_*$.
Für die Heisenberg-Kette mit $J>0$, offen oder periodisch mit geradem $L$, sind $A$ und $B$ die geraden und
ungeraden Plätze, $S_*=0$ für gerades $L$ und $S_*=\frac12$ für ungerades $L$ (offene Kette), und
$E_0(S_*)<E_0(S_*+1)<\dots<E_0(\frac L2)$.

*Beweis.* (a) Sei $M\ge0$ mit $M-\frac L2\in\mathbb Z$, $s=\max(M,S_*)$, und sei $\psi$ der Grundzustand des Sektors
$S^z=M$, einfach nach Satz 9.7: $\psi=V\varphi$ mit $\varphi>0$. $\mathbf S^2$ kommutiert mit $H$ und mit $S^z$, erhält also
den eindimensionalen Grundraum des Sektors: $\psi$ hat einen Gesamtspin $S'$. Wir konstruieren einen Vektor $\chi$
des Sektors mit $\mathbf S^2\chi=s(s+1)\chi$ und $\langle\psi,\chi\rangle\ne0$; da Eigenvektoren von $\mathbf S^2$ zu
verschiedenen Eigenwerten orthogonal sind, folgt $S'=s$. Man setze $p=\frac L2-s$; dann ist $0\le p\le n_B$, denn
$s\ge S_*$ bedeutet $p\le\frac{n_A+n_B}2-\frac{n_A-n_B}2$. Man wähle verschiedene $a_1,\dots,a_p\in A$ und verschiedene
$b_1,\dots,b_p\in B$, sei $R$ die Menge der übrigen $2s$ Plätze, und sei
$$
\chi_0=\bigotimes_{l=1}^p\tfrac1{\sqrt2}\big(|{\uparrow}_{a_l}{\downarrow}_{b_l}\rangle-|{\downarrow}_{a_l}{\uparrow}_{b_l}\rangle\big)\otimes|{\uparrow\cdots\uparrow}\rangle_R ,
\qquad \chi=(S^-)^{s-M}\chi_0 .
$$
Jedes Singulett wird von $s^\pm_{a_l}+s^\pm_{b_l}$ annihiliert und die Spins nach oben auf $R$ von $s^+_i$, also ist
$S^+\chi_0=0$; mit $S^z\chi_0=s\chi_0$ ergibt das $\mathbf S^2\chi_0=(S^-S^++(S^z)^2+S^z)\chi_0=s(s+1)\chi_0$
(Abschnitt 4.5), und $\mathbf S^2\chi=s(s+1)\chi$, da $\mathbf S^2$ mit $S^-$ kommutiert. Schreibt man
$S^-=\sum_l(s^-_{a_l}+s^-_{b_l})+S^-_R$, so kommutieren die Summanden, und die ersten annihilieren die Singuletts;
also ist $\chi$ das Produkt der Singuletts mit $(S^-_R)^{s-M}|{\uparrow\cdots\uparrow}\rangle_R=(s-M)!\sum_D|D\rangle$, der
Summe über die Konfigurationen von $R$ mit genau $s-M$ Spins nach unten (diese existieren, da $s-M\le2s=|R|$);
somit ist $\chi\ne0$, $S^z\chi=M\chi$, und die Koeffizienten von $\chi$ sind reell. Sein Vorzeichen: Ein Singulett
trägt $+$ bei, wenn sein Spin nach oben auf $A$ sitzt, und $-$, wenn er auf $B$ sitzt, und der $R$-Anteil ist
positiv. Ist $M=s$, so zeigen alle Plätze von $R$ nach oben; ist $M<s$, so ist $s=S_*$, $p=n_B$ und $R\subseteq A$. In
beiden Fällen ist $N_B(b)$ die Anzahl der Singuletts mit dem Spin nach oben auf $B$ plus eine Konstante, $|R\cap B|$
(falls $M=s$) bzw. $0$, also hat $(-1)^{N_B(b)}\chi(b)$ für alle $b$ mit $\chi(b)\ne0$ dasselbe Vorzeichen. Daher ist
$\langle\psi,\chi\rangle=\langle\varphi,V\chi\rangle=\pm\sum_b\varphi(b)|\chi(b)|\ne0$.

(b) Nach Proposition 4.15(a), Korollar 4.12 und Satz 4.11 trifft ein Multiplett vom Spin $S'$ den Sektor $S^z=M$ in
einem Vektor, wenn $S'\ge|M|$, und sonst gar nicht; also ist das Spektrum von $H$ auf dem Sektor $S^z=M\ge0$ die
Vereinigung der Spektren der $H_{S'}$, $S'\ge M$, mit Vielfachheiten. Nach (a) ist sein kleinster Eigenwert einfach
und gehört zum Spin $s=\max(M,S_*)$: Also ist $E_0(s)$ einfach in $H_s$, und $E_0(s)<E_0(S')$ für alle $S'\ge M$,
$S'\ne s$. Für $M=S\ge S_*$ ergibt das die Einfachheit und $E_0(S)<E_0(S')$ für $S'>S$; für $M=S<S_*$ ergibt es
$E_0(S_*)<E_0(S)$. Also ist $E_0(S_*)$ der kleinste Eigenwert von $H$, angenommen nur in $H_{S_*}$, wo er einfach ist.
Bei der Kette übertrifft die Anzahl der geraden Plätze die der ungeraden um $0$ ($L$ gerade) oder $1$ ($L$
ungerade). $\square$

Satz 9.7 geht auf W. Marshall (*Proc. Roy. Soc. A* 232 (1955) 48) zurück und Satz 9.9 auf E. Lieb und D. Mattis
(*J. Math. Phys.* 3 (1962) 749). Die **Singulett-Triplett-Lücke** $E_0(1)-E_0(0)$ des geraden Heisenberg-Rings ist
nach Satz 9.9 positiv.

Wir vergleichen die Blockgrößen für zwölf Spins, prüfen Satz 9.9 an den $U(1)$-Sektoren, berechnen die Lücke und
prüfen die Sätze 9.7 und 9.9 auf einem bipartiten Graphen, der keine Kette ist ($n_A=5$, $n_B=3$, $S_*=1$, zufällige
Kopplungen):

```python
sym = SU2Chain(12, boundary='periodic', complex_valued=False)
H_su2 = sym.heisenberg(J=1.0)                                       # (+)_S M_{m_S}, S = 0, ..., 6
u1_sizes = [math.comb(12, N) for N in range(13)]
print(sym.multiplicities, max(u1_sizes), sum(k ** 3 for k in u1_sizes), sum(m ** 3 for m in sym.multiplicities))
E0_S = sym.lowest_energies(H_su2)[0]                                # E_0(S), S = 0, ..., 6
E0_N = [torch.linalg.eigvalsh(ring.xxz(1.0, 1.0, sector=6 + S).matrix[0, 0])[0] for S in range(7)]
assert close(E0_S, torch.stack(E0_N))                               # Theorem 9.9: sector S^z = S has spin S
assert (E0_S[1:] > E0_S[:-1] + 1e-8).all()                          # E_0(S) strictly increasing
gap = {}
for n in (6, 8, 10, 12):
    s = SU2Chain(n, 'periodic', complex_valued=False)
    e = [torch.linalg.eigvalsh(s.heisenberg(spin=S).matrix[0, 0])[0].item() for S in (0, 1)]
    gap[n] = n * (e[1] - e[0])
print({n: round(g, 4) for n, g in gap.items()}, round(math.pi ** 2 / 2, 4))
assert 0 < gap[6] < gap[8] < gap[10] < gap[12] < math.pi ** 2 / 2
```

Die $SU(2)$-Blöcke haben die Größen $m_S=132,297,275,154,54,11,1$ gegenüber $924$ für den größten $U(1)$-Block, und
der Aufwand $\sum k^3$ der Diagonalisierung sinkt von $2.0\cdot10^9$ auf $5.3\cdot10^7$. $L$ mal die Lücke wächst
langsam, $4.108,4.181,4.232,4.270$ für $L=6,\dots,12$, gegen $\pi^2/2\approx4.935$.

**Korrelationen.** Die Spins der Plätze bilden einen Tensoroperator vom Rang eins, und nach dem
**Wigner-Eckart-Theorem** sind ihre Matrixelemente zwischen Multipletts ein Clebsch-Gordan-Koeffizient mal ein
**reduziertes Matrixelement** $\langle S',a\|\mathbf s_i\|S,b\rangle$, das nicht von den magnetischen Quantenzahlen
abhängt (A. R. Edmonds, *Angular Momentum in Quantum Mechanics*, Princeton 1957, Kap. 5; ohne Beweis).
`spin_correlation` nutzt dies, um $\mathbf s_i\cdot\mathbf s_j$ auf jedem $H_S$ allein aus Vektoren höchsten
Gewichts zu bauen; die Methode `'direct'` projiziert stattdessen den Bindungsoperator. Wir vergleichen beide mit dem
im $U(1)$-Sektor ausgeschriebenen Grundzustand:

```python
gs = sym.ground_state(H_su2)
S0, coef = gs.spin.item(), gs.coefficients[0, :sym.multiplicity(gs.spin.item())]
pairs = [(0, r) for r in range(12)]
corr_we = sym.expectation(sym.spin_correlation(pairs), S0, coef)
corr_dir = sym.expectation(sym.spin_correlation(pairs, method='direct'), S0, coef)
psi = sym.embed(S0, coef)                                            # the ground state in the sector N = 6
lab, z = sym.u1.labels[6], sym.u1.bits(6).to(f64) - 0.5
def s_dot_s(i, j):                                                   # <psi| s_i . s_j |psi>, from the definition
    if i == j:
        return 0.75
    flip = ((lab >> i) & 1) != ((lab >> j) & 1)
    partner = sym.u1.index_of[lab[flip] ^ ((1 << i) | (1 << j))]
    return ((psi ** 2 * z[:, i] * z[:, j]).sum() + 0.5 * (psi[flip] * psi[partner]).sum()).item()
assert S0 == 0 and close(corr_we, corr_dir)
assert close(corr_we, torch.tensor([s_dot_s(0, r) for r in range(12)], dtype=f64))
assert ((-1.0) ** torch.arange(1, 12) * corr_we[1:] > 0).all()       # signs (-1)^r
print([round(x, 4) for x in corr_we.tolist()])
```

Die Korrelationen $\langle\mathbf s_0\cdot\mathbf s_r\rangle=0.75,-0.4489,0.1878,-0.1658,0.1209,-0.1224,0.1070,\dots$
wechseln das Vorzeichen, wie es die Sätze 9.7 und 9.9 erzwingen: Das Singulett ist invariant unter Drehungen, die
die Achsen permutieren, und $\langle s^+_0s^+_r\rangle=0$ im Sektor $S^z=0$, also
$\langle\mathbf s_0\cdot\mathbf s_r\rangle=3\langle s^x_0s^x_r\rangle=
\frac34\langle s^+_0s^-_r+s^-_0s^+_r\rangle=\frac34\sum_{b_0\ne b_r}\psi(b)\psi(b')$, wobei $b'$ aus $b$ durch Vertauschen
der Spins auf den Plätzen $0$ und $r$ entsteht; da $N_B(b')-N_B(b)$ für ungerades $r$ gleich $\pm1$ und für gerades
$r$ gleich $0$ ist, hat jeder Term das Vorzeichen $(-1)^r$. Der Abfall ist langsam, wie $(-1)^r\sqrt{\log r}/r$ für
$L\to\infty$ (Affleck et al., *loc. cit.*; ohne Beweis).

## 9.4 Exakte Prüfungen: freie Fermionen und der Bethe-Ansatz

**Jordan-Wigner.** Auf einer offenen Kette setze man $c_j=\big(\prod_{l<j}(-1)^{n_l}\big)s^-_j$.

**Satz 9.10 (die XX-Kette ist frei).** (a) $\{c_i,c_j^*\}=\delta_{ij}$, $\{c_i,c_j\}=0$ und $c_j^*c_j=n_j$.
(b) Die offene XX-Kette ist $H=\frac J2\sum_{j<L-1}(s^+_js^-_{j+1}+s^-_js^+_{j+1})=\sum_{i,j}h_{ij}c_i^*c_j$ mit
$h=\frac J2A_L$. (c) Seien $u_k$ orthonormale reelle Eigenvektoren von $h$ zu den Eigenwerten
$\varepsilon_k=J\cos\frac{\pi k}{L+1}$ (Lemma 9.3) und $b_k=\sum_ju_k(j)c_j$. Für $K=\{k_1<\dots<k_N\}$ bilden die
Vektoren $b_K^*|0\rangle=b_{k_1}^*\cdots b_{k_N}^*|0\rangle$ ($|0\rangle$: alle Spins nach unten) eine orthonormale
Eigenbasis des Sektors $N$, mit Energien $\sum_{k\in K}\varepsilon_k$.

*Beweis.* (a) Operatoren auf verschiedenen Plätzen kommutieren, und auf einem Platz antikommutiert $(-1)^n$ mit
$s^\pm$. Für $i<j$ treffen $c_i$ und $c_j^*$ (oder $c_j$) auf Platz $i$ als $s^-_i$ und $(-1)^{n_i}$ aufeinander –
eine Antivertauschung – und kommutieren ansonsten, also antikommutieren sie. Für $i=j$ quadrieren sich die Strings
zu $1$, und $c_jc_j^*+c_j^*c_j=s^-_js^+_j+s^+_js^-_j=1$, $c_j^*c_j=s^+_js^-_j=n_j$. (b)
$c_j^*c_{j+1}=s^+_j(-1)^{n_j}s^-_{j+1}=s^+_js^-_{j+1}$, da $s^+_j\neq0$ nur auf $n_j=0$; dann geht man zu den
Adjungierten über. (c) $\{b_k,b_l^*\}=\sum_ju_k(j)u_l(j)=\delta_{kl}$ und $\{b_k,b_l\}=0$, und
$H=\sum_k\varepsilon_kb_k^*b_k$ nach der Spektralzerlegung von $h$. Aus den CAR folgt $[b_l^*b_l,b_k^*]=\delta_{lk}b_k^*$,
also $[H,b_k^*]=\varepsilon_kb_k^*$; mit $H|0\rangle=0$ ergibt das $Hb_K^*|0\rangle=(\sum_K\varepsilon_k)b_K^*|0\rangle$. Da
$b_k|0\rangle=0$ (jedes $c_j$ senkt einen Spin ab), liefert das Durchschieben der Vernichter nach rechts mit den CAR
$\langle b_{K'}^*0,b_K^*0\rangle=\delta_{KK'}$ durch Induktion über $N$. Jedes $c_j^*$ erhöht $N$ um eins, also liegt
$b_K^*|0\rangle$ im Sektor $|K|$: $\binom LN$ orthonormale Vektoren in einem Raum dieser Dimension. $\square$

```python
L = 10
xx = SpinChain(L)
eps = [math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1)]
for N in range(L + 1):
    E = torch.linalg.eigvalsh(xx.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
    E_ff = sorted(sum(K) for K in itertools.combinations(eps, N))     # sums over N-element subsets
    assert close(E, torch.tensor(E_ff, dtype=f64))                    # all 1024 levels, Theorem 9.10
    assert abs(E[0].item() - sum(sorted(eps)[:N])) < 1e-12           # ground state: the N lowest modes
```

**Magnonen.** Auf dem Ring mit $h=0$ hat der vollständig polarisierte Zustand $|0\rangle$ die Energie
$E_F=J\Delta L/4$. Ein Spin nach oben bei $x$ bricht zwei Bindungen auf:
$H|x\rangle=(E_F-J\Delta)|x\rangle+\frac J2(|x-1\rangle+|x+1\rangle)$ ($L\ge3$).

**Proposition 9.11 (ein Magnon).** Im Sektor $N=1$ des Rings erfüllt $|k\rangle=\sum_xe^{-ikx}|x\rangle$, $k=2\pi m/L$,
die Gleichungen $T|k\rangle=e^{ik}|k\rangle$ und $H|k\rangle=\big(E_F+J(\cos k-\Delta)\big)|k\rangle$.

*Beweis.* $\sum_xe^{-ikx}(|x-1\rangle+|x+1\rangle)=\sum_y(e^{-ik(y+1)}+e^{-ik(y-1)})|y\rangle=2\cos k\,|k\rangle$ unter
Verwendung von $e^{ikL}=1$; dieselbe Umindizierung liefert $T|k\rangle=e^{ik}|k\rangle$. $\square$

Zwei Magnonen streuen aneinander. Wir schreiben $|x_1,x_2\rangle$, $x_1<x_2$, für die Spins nach oben bei $x_1,x_2$; ein
Paar von Nachbarn bricht zwei Bindungen auf, ein getrenntes Paar vier.

**Proposition 9.12 (zwei Magnonen; Bethe-Gleichungen).** Seien $k_1,k_2\in\mathbb C$, $K=k_1+k_2$ und
$(A_{12},A_{21})\neq0$ mit (i) $A_{12}(1+e^{iK}-2\Delta e^{ik_2})+A_{21}(1+e^{iK}-2\Delta e^{ik_1})=0$ und
(ii) $A_{21}=e^{ik_2L}A_{12}$, $A_{12}=e^{ik_1L}A_{21}$. Seien
$\psi(x_1,x_2)=A_{12}e^{i(k_1x_1+k_2x_2)}+A_{21}e^{i(k_2x_1+k_1x_2)}$ auf $\mathbb Z^2$ und
$\Psi=\sum_{0\le x_1<x_2<L}\psi(x_1,x_2)|x_1,x_2\rangle$. Ist $\Psi\neq0$, so gilt
$H\Psi=\big(E_F+J(\cos k_1+\cos k_2-2\Delta)\big)\Psi$. Verschwinden die Nenner nicht, so besagen (i) und (ii)
$$
e^{ik_1L}=-\frac{1+e^{iK}-2\Delta e^{ik_1}}{1+e^{iK}-2\Delta e^{ik_2}},\qquad
e^{ik_2L}=-\frac{1+e^{iK}-2\Delta e^{ik_2}}{1+e^{iK}-2\Delta e^{ik_1}} .
$$

*Beweis.* Sei $E=E_F+J(\cos k_1+\cos k_2-2\Delta)$. (1) Für alle $(x_1,x_2)\in\mathbb Z^2$ liefert jede ebene Welle Werte
$\psi(x_1\pm1,x_2)$ und $\psi(x_1,x_2\pm1)$ mit der Summe $2(\cos k_1+\cos k_2)\psi(x_1,x_2)$. (2) Durch Einsetzen ist
$\psi(x,x)+\psi(x+1,x+1)-2\Delta\psi(x,x+1)$ gleich $e^{iKx}$ mal der linken Seite von (i), verschwindet also für alle
$x$. (3) Nach (ii) ist
$\psi(x_2,x_1+L)=A_{12}e^{ik_2L}e^{i(k_2x_1+k_1x_2)}+A_{21}e^{ik_1L}e^{i(k_1x_1+k_2x_2)}=\psi(x_1,x_2)$, und zweimalige
Anwendung ergibt $\psi(x_1+L,x_2+L)=\psi(x_1,x_2)$. Die Repräsentanten einer Konfiguration $\{y_1<y_2\}$ in
$D=\{x_1<x_2<x_1+L\}$ sind $(y_1+nL,y_2+nL)$ und $(y_2+nL,y_1+(n+1)L)$, also ist $\psi$ eine Funktion der
Konfigurationen. Sei $(x_1,x_2)\in D$. Ist $2\le x_2-x_1\le L-2$, so bleiben alle vier Sprünge in $D$, und der
Diagonalterm ist $E_F-2J\Delta$, also $(H\Psi)(x_1,x_2)=(E_F-2J\Delta)\psi+\frac J2\cdot2(\cos k_1+\cos k_2)\psi=E\psi$
nach (1). Ist $x_2=x_1+1$, so sind die Sprünge $(x_1+1,x_2)$, $(x_1,x_2-1)$ blockiert, und der Diagonalterm ist
$E_F-J\Delta$; nach (1) und (2) ist
$(H\Psi)(x_1,x_2)=(E_F-J\Delta)\psi+\frac J2\big(2(\cos k_1+\cos k_2)\psi-\psi(x_1,x_1)-\psi(x_1+1,x_1+1)\big)=E\psi$.
Ist $x_2=x_1+L-1$, so verwende man den Repräsentanten $(x_2,x_1+L)$, der vom vorigen Typ ist ($L\ge3$). $\square$

Proposition 9.12 sagt nur, dass Lösungen Eigenvektoren liefern; sie sagt nicht, dass alle Eigenvektoren so
entstehen. Wir schreiben $k_{1,2}=K/2\pm q$ und $c=\cos\frac K2$. Dann ist
$1+e^{iK}-2\Delta e^{ik_{1,2}}=2e^{iK/2}(c-\Delta e^{\pm iq})$, also lösen $A_{12}=c-\Delta e^{iq}$,
$A_{21}=-(c-\Delta e^{-iq})$ die Gleichung (i). Multiplikation der beiden Gleichungen (ii) ergibt $e^{iKL}=1$,
$K=2\pi n/L$, und dann ist (ii) äquivalent zu $F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}=0$, einer Gleichung
in einer Variablen. Wir suchen ihre Nullstellen $q\in(0,\pi)$, wo $e^{-iqL/2}F(q)$ reell oder rein imaginär ist,
weil $c-\Delta e^{\pm iq}$ zueinander komplex konjugiert sind (Streuzustände), und auf den Geraden $q=i\kappa$,
$\pi+i\kappa$, $\kappa>0$, wo $F$ reell ist (gebundene Zustände); für $L=9$, $\Delta=0.6$ bauen wir $\Psi$ zu jeder
Nullstelle und vergleichen mit dem Sektor $N=2$ mit $36$ Zuständen:

```python
L, Delta = 9, 0.6
ring9 = SpinChain(L, boundary='periodic')
H2 = ring9.xxz(J=1.0, Delta=Delta, sector=2).matrix[0, 0]
x1, x2 = np.triu_indices(L, 1)                                       # all pairs x1 < x2
rows = ring9.index_of[torch.tensor((1 << x1) | (1 << x2))]
def roots(g, grid):                                                  # sign changes of g on the grid, refined
    v = np.array([g(t) for t in grid])
    return [brentq(g, grid[i], grid[i + 1], xtol=1e-15) for i in np.nonzero(v[:-1] * v[1:] < 0)[0]]
E_bethe, vecs, n_bound = [], [], 0
for n in range(L):
    K, c = 2 * np.pi * n / L, math.cos(np.pi * n / L)
    F = lambda q: (-1) ** n * np.exp(1j * q * L) * (c - Delta * np.exp(-1j * q)) + c - Delta * np.exp(1j * q)
    G = lambda q: (lambda z: z.real + z.imag)(np.exp(-0.5j * q * L) * F(q))    # real, or imaginary
    bound = [s + 1j * t for s in (0, np.pi)                          # q = i kappa and pi + i kappa
             for t in roots(lambda t: F(s + 1j * t).real, np.linspace(1e-3, 20, 4000))]
    n_bound += len(bound)
    for q in roots(G, (np.arange(1, 2000) + 0.5 ** 0.5) * np.pi / 2000) + bound:
        psi = np.exp(0.5j * K * (x1 + x2)) * ((c - Delta * np.exp(1j * q)) * np.exp(1j * q * (x1 - x2))
                                              - (c - Delta * np.exp(-1j * q)) * np.exp(-1j * q * (x1 - x2)))
        v = torch.zeros(36, dtype=c128)
        v[rows] = torch.tensor(psi / np.linalg.norm(psi))
        E = 9 * Delta / 4 + 2 * (c * np.cos(q) - Delta)                 # E_F + J (cos k1 + cos k2 - 2 Delta)
        assert abs(E.imag) < 1e-12 and close(H2 @ v, E.real * v, 1e-9)   # Proposition 9.12
        E_bethe.append(E.real)
        vecs.append(v)
assert torch.linalg.matrix_rank(torch.stack(vecs)) == 36             # a complete eigenbasis
assert close(torch.sort(torch.tensor(E_bethe)).values, torch.linalg.eigvalsh(H2))
print(len(E_bethe), n_bound)
ring8, Phi = SpinChain(8, 'periodic'), torch.zeros(28, dtype=c128)     # even L: Phi = sum_x (-1)^x |x, x+1>
for x in range(8):
    Phi[ring8.index_of[(1 << x) | (1 << (x + 1) % 8)]] = (-1) ** x
assert close(ring8.xxz(1.0, Delta, sector=2).matrix[0, 0] @ Phi, (8 * Delta / 4 - Delta) * Phi)   # E_F - J Delta
```

Für $L=9$ stammen alle $36$ Eigenzustände aus den Bethe-Gleichungen: $34$ Streuzustände und $2$ gebundene
Zustände. Die Vollständigkeit wird hier überprüft, nicht bewiesen, und für gerades $L$ ist sie verletzt. Ist
$\Delta\ne0$, so erfüllt $\Phi=\sum_x(-1)^x|x,x+1\rangle$ die Gleichungen $T\Phi=-\Phi$ und $H\Phi=(E_F-J\Delta)\Phi$ (die
Sprünge nach $|x-1,x+1\rangle$ und $|x,x+2\rangle$ heben sich paarweise auf). Ein Bethe-Vektor erfüllt
$T\Psi=e^{-iK}\Psi$, da $\psi(x_1-1,x_2-1)=e^{-iK}\psi(x_1,x_2)$; ist $e^{iK}=-1$, so gilt
$\cos k_1+\cos k_2=2\cos\frac K2\cos\frac{k_1-k_2}2=0$, und seine Energie ist $E_F-2J\Delta$. Also ist $\Phi$ zu allen
Bethe-Vektoren orthogonal (die letzten Zeilen des Codes prüfen $H\Phi$ für $L=8$). Mit $\Delta=0.6$ und $4\le L\le13$
findet die obige Suche für ungerades $L$ alle Zustände und für gerades $L$ alle außer $\Phi$. Für $N>2$ Magnonen hat
der Ansatz $N!$ Terme und $N$ gekoppelte Bethe-Gleichungen (H. Bethe, *Z. Phys.* 71 (1931) 205; M. Karbach,
G. Müller, *Computers in Physics* 11 (1997) 36).

> **Physikalische Bedeutung.** Ein Magnon ist eine Spinwelle: ein umgeklappter Spin, delokalisiert über den Ring,
> mit Impuls $k$. Zwei Magnonen streuen mit einer Phasenverschiebung (reelles $q$) oder binden aneinander
> (komplexes $q$). Bethes Einsicht war, dass für $N$ Magnonen die Streuung in Zweikörperphasen faktorisiert; seine
> Gleichungen liefern zum Beispiel die Grundzustandsenergie pro Platz $\frac14-\log2$ der unendlichen
> Heisenberg-Kette (L. Hulthén, *Ark. Mat. Astron. Fys.* 26A (1938) Nr. 11; ohne Beweis).

## 9.5 Lanczos- und Krylov-Verfahren

Der Sektor $N=L/2$ hat die Dimension $12870$ für $L=16$ und $184756$ für $L=20$, aber $H$ ist dünn besetzt ($O(L)$
von null verschiedene Einträge pro Zeile, `xxz_sparse`), und Krylov-Verfahren brauchen nur Produkte $Hv$. Sei $H$
hermitesch von der Größe $n$ und $v\ne0$. Der **Krylov-Raum** ist
$\mathcal K_m=\operatorname{span}\{v,Hv,\dots,H^{m-1}v\}$.

**Proposition 9.13 (Lanczos).** Sei $\dim\mathcal K_m=m$. Sei $q_1=v/\|v\|$, und für $j<m$ sei $q_{j+1}$ die normierte
Komponente von $Hq_j$ orthogonal zu $\mathcal K_j$. Dann ist $q_1,\dots,q_m$ eine Orthonormalbasis von $\mathcal K_m$;
mit $Q_m=(q_1\cdots q_m)$ ist $T_m=Q_m^*HQ_m$ reell symmetrisch und tridiagonal mit Diagonale
$\alpha_j=\langle q_j,Hq_j\rangle$ und Nebendiagonale $\beta_j=\langle q_{j+1},Hq_j\rangle>0$; und $HQ_m=Q_mT_m+r_me_m^T$, wobei
$r_m\perp\mathcal K_m$ genau dann verschwindet, wenn $\mathcal K_{m+1}=\mathcal K_m$ (andernfalls $r_m=\beta_mq_{m+1}$,
$q_{m+1}=r_m/\|r_m\|$).

*Beweis.* Es gelte $\operatorname{span}\{q_1,\dots,q_j\}=\mathcal K_j$ (wahr für $j=1$). Dann ist
$H\mathcal K_j\subseteq\mathcal K_j+\mathbb CHq_j$, da $Hq_i\in\mathcal K_{i+1}$; also
$\mathcal K_{j+1}=\mathcal K_j+\mathbb CH^jv\subseteq\mathcal K_j+\mathbb CHq_j\subseteq\mathcal K_{j+1}$ und
$\mathcal K_{j+1}=\mathcal K_j+\mathbb CHq_j$. Für $j<m$ wächst die Dimension, also ist die zu $\mathcal K_j$ orthogonale
Komponente von $Hq_j$ von null verschieden, und die Induktion läuft weiter; für $j=m$ verschwindet diese
Komponente, $r_m$, genau dann, wenn $\mathcal K_{m+1}=\mathcal K_m$. Für $i>j+1$ ist $q_i\perp\mathcal K_{j+1}\ni Hq_j$, also
$(T_m)_{ij}=0$, und $(T_m)_{ji}=0$, da $T_m$ hermitesch ist. $\alpha_j$ ist reell, und $\beta_j=\langle q_{j+1},Hq_j\rangle$
ist die Norm der Komponente, $>0$. Also ist $Hq_j=\beta_{j-1}q_{j-1}+\alpha_jq_j+\beta_jq_{j+1}$ für $j<m$ und
$Hq_m=\beta_{m-1}q_{m-1}+\alpha_mq_m+r_m$. $\square$

Die Eigenwerte $\theta^{(m)}_1\le\dots\le\theta^{(m)}_m$ von $T_m$ sind die **Ritz-Werte**.

**Satz 9.14 (monotone Konvergenz).** Seien $\lambda_1\le\dots\le\lambda_n$ die Eigenwerte von $H$ und
$\dim\mathcal K_m=m$. Dann gilt $\theta^{(m)}_j\ge\lambda_j$, und $\theta^{(m)}_j\ge\theta^{(m+1)}_j$, falls
$\dim\mathcal K_{m+1}=m+1$ ($j\le m$). Hat $v$ eine von null verschiedene Komponente im Eigenraum von $\lambda_1$ und
ist $\mathcal K_{m+1}=\mathcal K_m$, so ist $\theta^{(m)}_1=\lambda_1$.

*Beweis.* Nach Courant-Fischer (Standard der linearen Algebra) ist
$\theta^{(m)}_j=\min_{U}\max_{x\in U,\|x\|=1}\langle x,Hx\rangle$ über die $j$-dimensionalen $U\subseteq\mathcal K_m$, da
$Q_m$ eine Isometrie auf $\mathcal K_m$ ist und $\langle y,T_my\rangle=\langle Q_my,HQ_my\rangle$. Dieselbe Formel über
$U\subseteq\mathcal K_{m+1}$ liefert $\theta^{(m+1)}_j$ und über $U\subseteq\mathbb C^n$ liefert $\lambda_j$; Minimieren über
mehr Unterräume ergibt einen kleineren Wert. Ist $\mathcal K_{m+1}=\mathcal K_m$, so ist $\mathcal K_m$ $H$-invariant und
enthält $p(H)v$ für jedes Polynom $p$, insbesondere $P_1v\ne0$ ($P_1$ die Spektralprojektion zu $\lambda_1$, nach
Lagrange-Interpolation ein Polynom in $H$, wie in Satz 1.5); also ist $\theta^{(m)}_1\le\langle x,Hx\rangle=\lambda_1$ bei
$x=P_1v/\|P_1v\|\in\mathcal K_m$. $\square$

Wie schnell $\theta^{(m)}_1\to\lambda_1$ konvergiert, hängt von der Lücke $\lambda_2-\lambda_1$ im Verhältnis zur Breite des
Spektrums ab (Kaniel-Paige-Saad-Schranken über Tschebyschow-Polynome; Y. Saad, *Numerical Methods for Large
Eigenvalue Problems*, SIAM 2011, Kap. 6; ohne Beweis). In Gleitkommaarithmetik verlieren die $q_j$ ihre
Orthogonalität; `krylov.lanczos` reorthogonalisiert vollständig. Für die XX-Kette ist die exakte
Grundzustandsenergie bekannt (Satz 9.10):

```python
L, N = 12, 6
xx = SpinChain(L, complex_valued=False)
Hs = xx.xxz_sparse(J=1.0, Delta=0.0, sector=N)                      # sparse: diagonal + hopping
alpha, beta, Q = krylov.lanczos(Hs.matvec, torch.randn(1, xx.sector_dim(N), dtype=f64), 40)
T = torch.diag(alpha[0]) + torch.diag(beta[0], 1) + torch.diag(beta[0], -1)
assert close(Q[0] @ Q[0].T, torch.eye(40, dtype=f64)) and close(Q[0] @ Hs.to_dense()[0] @ Q[0].T, T)
ritz = torch.tensor([torch.linalg.eigvalsh(T[:m, :m])[0].item() for m in range(1, 41)])   # T_m = leading block
E_xx = sum(sorted(math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1))[:N])
assert (ritz[1:] <= ritz[:-1] + 1e-12).all() and (ritz >= E_xx - 1e-12).all()             # Theorem 9.14
print([f"{(ritz[m - 1] - E_xx).item():.1e}" for m in (5, 10, 20, 40)])
```

Der Fehler des kleinsten Ritz-Werts fällt monoton, $0.42$, $3.1\cdot10^{-3}$ und $3.1\cdot10^{-6}$ nach $5$, $10$ und
$20$ Schritten, und erreicht nach $40$ Schritten das Rundungsniveau, in Dimension $924$. `krylov.ground_state`
startet Lanczos neu, bis das Residuum $\|H\psi-E\psi\|$ unter $10^{-10}$ liegt. Wir prüfen es an den
Heisenberg-Ring-Energien, die in `scripts/chains.py` auf neun Stellen tabelliert sind (die Referenzwerte des
Begleitartikels), dicht bis $L=12$ und mit Lanczos für $L=16$:

```python
RING = {4: -2.0, 6: -2.802775638, 8: -3.651093409, 10: -4.515446354, 12: -5.387390917,
        14: -6.263549534, 16: -7.142296361}
for n in (4, 6, 8, 10, 12):
    e0 = torch.linalg.eigvalsh(SpinChain(n, 'periodic').heisenberg(sector=n // 2).matrix[0, 0])[0].item()
    assert abs(e0 - RING[n]) < 5e-10
    print(n, round(e0 / n, 6))
t0 = time.time()
E16, psi16, res = krylov.ground_state(SpinChain(16, 'periodic', complex_valued=False).xxz_sparse(sector=8),
                                      generator=torch.Generator().manual_seed(0))
print(f"L = 16: E0 = {E16.item():.10f}, residual {res.item():.1e}, {time.time() - t0:.2f} s; "
      f"infinite chain: {0.25 - math.log(2):.6f}")
assert abs(E16.item() - RING[16]) < 5e-10 and res.item() < 1e-9
```

Für $L=16$ (Dimension $12870$) braucht Lanczos auf einem CPU-Thread deutlich weniger als eine Sekunde. Die Energie
pro Platz, $-0.5,-0.467129,-0.456387,-0.451545,-0.448949$ für $L=4,\dots,12$, nähert sich von unten dem Wert
$\frac14-\log2=-0.443147$ der unendlichen Kette.

**Zeitentwicklung.** Für $e^{-i\tau H}v$ projiziert man auf $\mathcal K_m$: $w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1$.

**Proposition 9.15.** Sei $\dim\mathcal K_m=m$. Für jedes $\tau$ gilt $\|w_\tau\|=\|v\|$ und
$\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle$. Ist $\mathcal K_m$ $H$-invariant, so ist $w_\tau=e^{-i\tau H}v$.

*Beweis.* $Q_m$ ist eine Isometrie und $e^{-i\tau T_m}$ unitär ($T_m$ ist reell symmetrisch), das liefert die Norm. Da
$Q_m^*HQ_m=T_m$ mit $e^{-i\tau T_m}$ kommutiert, ist
$\langle w_\tau,Hw_\tau\rangle=\|v\|^2\langle e^{-i\tau T_m}e_1,T_me^{-i\tau T_m}e_1\rangle=\|v\|^2\langle e_1,T_me_1\rangle=\langle v,Hv\rangle$. Ist
$\mathcal K_m$ invariant, so gilt $HQ_m=Q_mT_m$ (Proposition 9.13, $r_m=0$), also $H^jQ_m=Q_mT_m^j$ und
$e^{-i\tau H}Q_m=Q_me^{-i\tau T_m}$; man wende dies auf $e_1$ an, mit $Q_me_1=v/\|v\|$. $\square$

Jeder Krylov-Schritt erhält also Norm und Energie, ob er genau ist oder nicht; Genauigkeit erfordert ein im
Vergleich zu $\tau\|H\|$ großes $m$, und `krylov.evolve` unterteilt das Zeitintervall und schätzt den Fehler a
posteriori (M. Hochbruck, C. Lubich, *SIAM J. Numer. Anal.* 34 (1997) 1911; ohne Beweis). Wir führen einen
Quench aus dem Néel-Zustand $|{\uparrow\downarrow\uparrow\downarrow}\cdots\rangle$ auf dem Heisenberg-Ring durch und
verfolgen die alternierende Magnetisierung $m_s=\frac2L\sum_i(-1)^i\langle s^z_i\rangle$, im Vergleich mit dem dichten
Propagator für $L=12$:

```python
times, m_s = [0.0, 0.5, 1.0, 2.0, 4.0], {}
stag = lambda c, N, out: [((p[0].abs() ** 2) @ (c.bits(N).to(f64) - 0.5) @ torch.tensor(
    [(-1.0) ** i for i in range(c.L)], dtype=f64) * 2 / c.L).item() for p in out]
for n in (12, 16):
    c = SpinChain(n, 'periodic')
    Hn = c.xxz_sparse(J=1.0, Delta=1.0, sector=n // 2)
    psi0 = c.vector_in_sector('10' * (n // 2))[None]
    out = krylov.evolve(Hn, psi0, times)                               # (times, batch, dim)
    energy = torch.stack([(p.conj() * Hn.matvec(p)).sum(-1).real for p in out])
    assert close(out.norm(dim=-1), 1.0) and close(energy, energy[0])   # Proposition 9.15
    m_s[n] = stag(c, n // 2, out)
    print(n, [round(x, 4) for x in m_s[n]])
    if n == 12:
        w, U = torch.linalg.eigh(c.xxz(1.0, 1.0, sector=6).matrix[0, 0])
        for t, p in zip(times, out):
            assert close(p[0], U @ (torch.exp(-1j * t * w) * (U.mH @ psi0[0])))
print(f"{max(abs(a - b) for a, b in zip(m_s[12][:4], m_s[16][:4])):.1e}")   # L = 12 vs 16 up to t = 2
```

$m_s$ fällt von $1$ über $0.770$ und $0.279$ auf $-0.183$ bei $t=2$, für beide Größen bis auf $1.7\cdot10^{-5}$
übereinstimmend; bei $t=4$ unterscheiden sie sich ($0.026$ gegenüber $0.035$): Die endliche Größe des kleineren
Rings wird sichtbar.

## 9.6 Unordnung: Niveaustatistik und Verschränkung

Die **Heisenberg-Kette im Zufallsfeld** hat unabhängige, auf $[-W,W]$ gleichverteilte Felder $h_i$. Sie erhält
$S^z$ (Proposition 9.1) und generisch nichts weiter: Die Felder brechen die Translationen, die Spiegelungen, die
Spinumkehr und $SU(2)$. Daher wird die Niveaustatistik pro Sektor $N$ berechnet (Satz 5.15 und Abschnitt 5.4: Das
Mischen von Sektoren täuscht unabhängige Niveaus vor). Bei schwacher Unordnung stoßen sich die Niveaus wie im GOE
ab, $\langle r\rangle\approx0.5307$ (Abschnitt 5.3); bei starker Unordnung lokalisieren die Eigenzustände
(**Vielteilchenlokalisierung**), und die Niveaus werden unabhängig, $\langle r\rangle=2\log2-1\approx0.3863$
(Proposition 5.14). Wir berechnen $\langle r\rangle$ in der mittleren Hälfte des Sektors $N=L/2$, mit einem
Standardfehler über die unabhängigen Realisierungen (jede liefert ein mittleres Verhältnis):

```python
gen = torch.Generator().manual_seed(1)
r_mean = {}
for n, R in ((8, 200), (10, 100)):
    c = SpinChain(n, 'periodic', complex_valued=False)
    k = c.sector_dim(n // 2)
    for W in (0.5, 1, 2, 4, 8, 16):
        E = torch.linalg.eigvalsh(c.random_field_heisenberg(W, R, sector=n // 2, generator=gen).matrix[:, 0])
        r = SpinChain.level_spacing_ratio(E[:, k // 4: 3 * k // 4])     # one value per realisation
        r_mean[n, W] = (r.mean().item(), (r.std() / math.sqrt(R)).item())
    print(n, '  '.join(f"W={W}: {m:.3f}+-{s:.3f}" for (n_, W), (m, s) in r_mean.items() if n_ == n))
assert abs(r_mean[10, 0.5][0] - 0.5307) < 0.015 and abs(r_mean[10, 16][0] - 0.3863) < 0.015
assert r_mean[10, 1][0] > r_mean[10, 2][0] > r_mean[10, 4][0] > r_mean[10, 8][0]
```

Für $L=10$ geht das mittlere Verhältnis von $0.530$ ($W=0.5$) über $0.484$ ($W=2$) und $0.410$ ($W=4$) bis $0.384$
($W=16$), mit Standardfehlern von $0.003$–$0.004$. Für $L=8$ liegt der Wert bei $W=0.5$ ($0.507$) unter dem bei
$W=1$ ($0.525$): Bei $W=0$ ist die Kette integrabel und besitzt alle obigen Symmetrien, und ein GOE-Wert erfordert
eine Unordnung, die groß gegenüber dem Niveauabstand ist. Der Übergangsbereich ist breit und driftet mit $L$; eine
Finite-Size-Skalierung verortet ihn bei $W\approx3.7$ (D. Luitz, N. Laflorencie, F. Alet, *Phys. Rev. B* 91 (2015)
081103), und ob eine lokalisierte Phase für $L\to\infty$ überlebt, ist umstritten (J. Šuntajs, J. Bonča, T. Prosen,
L. Vidmar, *Phys. Rev. E* 102 (2020) 062144).

**Verschränkung freier Fermionen.** Für die Eigenzustände der offenen XX-Kette ist die Verschränkungsentropie der
ersten $n$ Plätze durch eine $n\times n$-Matrix bestimmt.

**Satz 9.16 (Peschels Formel).** Sei $\Phi=b_{k_1}^*\cdots b_{k_N}^*|0\rangle$ wie in Satz 9.10,
$C_{ij}=\langle\Phi,c_i^*c_j\Phi\rangle$ und $C_A=(C_{ij})_{i,j<n}$. Dann gilt $0\le C_A\le1$, und die
Verschränkungsentropie der Plätze $A=\{0,\dots,n-1\}$ ist, wobei $\nu$ die Eigenwerte von $C_A$ mit Vielfachheit
durchläuft,
$$
S_A=\sum_\nu \eta(\nu),\qquad \eta(\nu)=-\nu\log\nu-(1-\nu)\log(1-\nu)\quad(0\log0=0).
$$

*Beweis.* Für einen reellen Vektor $x$ setze man $b(x)=\sum_jx_jc_j$; dann ist $b(x)^*b(y)^*=-b(y)^*b(x)^*$ und
$\{b(x),b(y)^*\}=\langle x,y\rangle$. Sei $U\subseteq\mathbb R^L$ der Spann der $u_{k_a}$, und seien $P_A$, $P_B=1-P_A$ die
Koordinatenprojektionen auf $\{0,\dots,n-1\}$ und dessen Komplement. Diagonalisieren der Form $\langle P_Ax,P_Ay\rangle$
auf $U$ liefert eine Orthonormalbasis $\phi_1,\dots,\phi_N$ von $U$ mit $\langle P_A\phi_a,P_A\phi_b\rangle=\nu_a\delta_{ab}$,
also $\langle P_B\phi_a,P_B\phi_b\rangle=(1-\nu_a)\delta_{ab}$. Da $(x_1,\dots,x_N)\mapsto b(x_1)^*\cdots b(x_N)^*|0\rangle$
multilinear und alternierend ist, gilt $\Phi=\pm\prod_a\big(\sqrt{\nu_a}f_a^*+\sqrt{1-\nu_a}\,g_a^*\big)|0\rangle$ mit
$f_a=b(P_A\phi_a/\sqrt{\nu_a})$, $g_a=b(P_B\phi_a/\sqrt{1-\nu_a})$ (Terme mit Koeffizient null weggelassen).
Ausmultiplizieren ergibt $\Phi=\sum_I\pm\prod_{a\in I}\sqrt{\nu_a}\prod_{a\notin I}\sqrt{1-\nu_a}\,F_I^*G_{I^c}^*|0\rangle$,
summiert über die Teilmengen $I$, mit den geordneten Produkten $F_I^*$, $G_{I^c}^*$. Für $j\ge n$ enthält der String
von $c_j^*$ ganz $A$: $c_j^*=(-1)^{N_A}\otimes\tilde c_j^*$ mit Jordan-Wigner-Operatoren $\tilde c_j$ der Plätze $\ge n$,
während für $i<n$ gilt $c_i^*=\hat c_i^*\otimes1$. Da $(-1)^{N_A}$ auf $|0_A\rangle$ gleich $1$ ist, folgt
$F_I^*G_{I^c}^*|0\rangle=\xi_I\otimes\chi_{I^c}$ mit $\xi_I=\hat F_I^*|0_A\rangle$, $\chi_{I^c}=\tilde G_{I^c}^*|0_B\rangle$, und
diese bilden nach der CAR-Rechnung aus Satz 9.10 Orthonormalsysteme. Das ist eine Schmidt-Zerlegung (Lemma 7.7) mit
Gewichten $p_I=\prod_{a\in I}\nu_a\prod_{a\notin I}(1-\nu_a)$, der Verteilung unabhängiger Bernoulli-Variablen, also
$S_A=-\sum_Ip_I\log p_I=\sum_a\eta(\nu_a)$. Schließlich ist $C_{ij}=\sum_a\phi_a(i)\phi_a(j)$ (man entwickle $c_i^*c_j$
nach den $b(\phi_a)$), also hat $C_A=\sum_a(P_A\phi_a)(P_A\phi_a)^T$ die Eigenwerte $\nu_a$ auf den orthogonalen Vektoren
$P_A\phi_a$ und sonst $0$, und $\eta(0)=0$. $\square$

`entanglement_entropy` berechnet den reduzierten Zustand Sektor für Sektor (er kommutiert mit der Magnetisierung
des Blocks, Proposition 7.8) und seine Entropie. Der Grundzustand des Sektors $N$ ist $\Phi$ mit den $N$
niedrigsten Moden, einfach, weil die $\varepsilon_k$ verschieden sind:

```python
L, N = 12, 6
xx = SpinChain(L)
w, V = torch.linalg.eigh(xx.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
rho = DensityMatrix(xx.sector_algebra(N), matrix=torch.outer(V[:, 0], V[:, 0].conj())[None, None])
eps_k, u = np.linalg.eigh(np.eye(L, k=1) / 2 + np.eye(L, k=-1) / 2)  # h = A_L / 2, ascending
C = u[:, :N] @ u[:, :N].T                                            # correlation matrix of the N lowest modes
bern = lambda nu: -(nu * np.log(nu) + (1 - nu) * np.log(1 - nu))
S_peschel = [bern(np.clip(np.linalg.eigvalsh(C[:n, :n]), 1e-15, 1 - 1e-15)).sum() for n in range(1, L)]
S_lib = [xx.entanglement_entropy(rho, n, sector=N).item() for n in range(1, L)]
assert np.allclose(S_lib, S_peschel, atol=1e-9)                       # Theorem 9.16
print([round(s, 4) for s in S_lib[:6]])
```

Die Entropien $0.6931,0.5199,0.7456,0.6283,0.7719,0.6560$ für $n=1,\dots,6$ oszillieren mit der Parität von $n$ und
wachsen bei festem $n/L$ wie $\frac16\log L$, das Gesetz der konformen Feldtheorie für eine offene Kette
(I. Peschel, *J. Phys. A* 36 (2003) L205; P. Calabrese, J. Cardy, *J. Stat. Mech.* (2004) P06002; ohne Beweis).

## 9.7 Womit man vergleichen sollte

QuSpin (Python) baut Spinketten-Hamiltonoperatoren in Sektoren von Magnetisierung, Impuls, Parität und Spinumkehr
und verwendet die dünn besetzten Eigenwertlöser und Matrixexponentiale von SciPy; spezialisierte Codes für exakte
Diagonalisierung (XDiag, SPINPACK, HΦ) sind für große dünn besetzte Probleme optimiert, und mit Gittersymmetrien und
großen Rechnern haben solche Codes etwa 40–50 Spins erreicht; ITensor, TeNPy und TensorKit.jl (mit
$SU(2)$-Fusionsbäumen) behandeln lange Ketten mit Matrixproduktzuständen. Was die Bibliothek hinzufügt, ist
bescheiden und spezifisch: die Kette als `TypeIAlgebra` $\bigoplus_NM_{\binom LN}$, sodass die Zustände, Kanäle,
Entropien und Dynamiken der Kapitel 3–8 Sektor für Sektor gelten; den Inklusionsturm mit seinem Index;
$SU(2)$-Blöcke aus Vektoren höchsten Gewichts mit Wigner-Eckart-Korrelationen; Stapel von
Unordnungsrealisierungen in einem einzigen Aufruf; und reduzierte Zustände, die die Sektoren nie verlassen. Sie hat
keine Spiegelungs- oder Spinumkehrsektoren, Impulssektoren nur für homogene XXZ-Ketten, dichte Verfahren bis
$L\approx16$, Krylov-Verfahren bis $L\approx24$ (etwa $26$ auf einer großen GPU) und keine Tensornetzwerke.

## Aufgaben

1. (Impulssektoren von Hand) Man bestimme für $L=4$ und $N=2$ die Bahnen der sechs Konfigurationen unter $T$ und
   die Dimensionen $d_{2,m}$, $m=0,\dots,3$.
2. (Summenregeln) Sei $\psi$ der Grundzustand des Heisenberg-Rings mit geradem $L$. Man zeige
   $E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$ und $\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0$.
3. (Zu viel Symmetrie) Man zeige, dass im Sektor $N=L/2$ der Vektor mit lauter Koeffizienten $1$ ein Vielfaches
   von $(S^-)^{L/2}|{\uparrow\cdots\uparrow}\rangle$ und ein Eigenvektor des Heisenberg-Rings zur Energie $JL/4$ ist und
   dass jeder Krylov-Raum $\mathcal K_m(H,v)$ in jedem $H$-invarianten Unterraum liegt, der $v$ enthält. Was liefert
   `krylov.ground_state`, wenn man von diesem Vektor aus startet?

<details>
<summary>Lösungen</summary>

1. $0011\to1001\to1100\to0110\to0011$ ist eine Bahn mit $R=4$, und $0101\to1010\to0101$ eine mit $R=2$. Nach
   Proposition 9.5 trägt die erste zu jedem $m$ bei, die zweite zu den $m$ mit $2m\equiv0\pmod4$, d. h. $m=0,2$:
   $d_{2,m}=2,1,2,1$.

```python
assert SpinChain(4, 'periodic').momentum_dims(2) == [2, 1, 2, 1]
```

2. Der Grundzustand ist ein einfaches Singulett (Satz 9.9), und $T$ kommutiert mit $H$, also $T\psi=e^{ik}\psi$. Da
   $\mathbf s_i\cdot\mathbf s_{i+1}=T^i(\mathbf s_0\cdot\mathbf s_1)T^{*i}$, ist $\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle$ für alle
   $L$ Bindungen gleich, und $E_0=J\sum_i\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$. Da
   $\mathbf S^2\psi=0$, ist
   $0=\langle\mathbf S^2\rangle=\sum_{i,j}\langle\mathbf s_i\cdot\mathbf s_j\rangle=L\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle$, wiederum
   wegen der Translationsinvarianz.

```python
assert abs(12 * corr_we[1].item() - gs.energy.item()) < 1e-10 and abs(corr_we.sum().item()) < 1e-10
```

3. $S^-$ senkt einen der Spins nach oben ab; $L/2$-malige Anwendung auf den vollständig polarisierten Zustand
   erreicht jede Konfiguration mit $L/2$ Spins nach oben in $(L/2)!$ Reihenfolgen, alle mit Koeffizient $1$.
   $|{\uparrow\cdots\uparrow}\rangle$ hat alle $L$ Bindungen parallel, Energie $JL/4$, und $S^-$ kommutiert mit $H$,
   also ist der Vektor ein Eigenvektor zur selben Energie. Ist $W$ invariant und $v\in W$, so ist $H^jv\in W$ für
   alle $j$. Hier ist $\mathcal K_m=\mathbb Cv$ für alle $m$, und `ground_state` liefert $E=JL/4$, das obere Ende des
   Spektrums (jedes $\mathbf s_i\cdot\mathbf s_j\le\frac14$); ein zufälliger Startvektor vermeidet dies mit
   Wahrscheinlichkeit eins.

```python
c10 = SpinChain(10, 'periodic', complex_valued=False)
E_sym, _, _ = krylov.ground_state(c10.xxz_sparse(1.0, 1.0, sector=5), v0=torch.ones(1, 252, dtype=f64))
assert abs(E_sym.item() - 2.5) < 1e-10 and E_sym.item() > RING[10] + 7     # far from E_0
```

</details>

**Weiter:** [Kapitel 10](10_quantum_computing.md) behandelt Quantencomputing mit Sektoren: dekohärenzfreie
Unterräume und rauschfreie Teilsysteme, Fehlerkorrektur in Operatoralgebren und Mixer, die ein erhaltenes
Hamming-Gewicht respektieren – die Algebra $\bigoplus_wM_{\binom nw}$ dieses Kapitels in einer neuen Rolle.
