[English](../05_random_matrices.md) · [Español](../es/05_random_matrices.md) · [Français](../fr/05_random_matrices.md) · **Deutsch** · [中文](../zh/05_random_matrices.md) · [日本語](../ja/05_random_matrices.md) · [Русский](../ru/05_random_matrices.md)

# 5. Zufallsmatrizen in einer Algebra

Die Theorie der Zufallsmatrizen fragt, wie ein „typischer“ Operator aussieht. In einer Algebra
$\bigoplus_cM_{k_c}$ muss die Antwort die Sektoren respektieren: Die natürlichen zufälligen
unitären Elemente sind die der Algebra, nicht die der großen Matrixalgebra $M_D$, und die
natürlichen zufälligen Hamiltonoperatoren haben unabhängige Blöcke. Dieses Kapitel konstruiert
das Haar-Maß auf der unitären Gruppe der Algebra, berechnet seine ersten und zweiten Momente
exakt und zeigt dann, was in der Spektralstatistik schiefgeht, wenn man die Sektoren ignoriert.

**Sie lernen**, warum die QR-Zerlegung einer Gaußschen Matrix eine Haar-verteilte unitäre Matrix
liefert (und warum die Phasen festgelegt werden müssen), dass die Mittelung über die unitäre Gruppe
der Algebra die bedingte Erwartung auf das Zentrum ist, die Weingarten-Formel zweiter Ordnung, die
Wignersche Vermutung, das Abstandsverhältnis und seinen exakten Wert für Gaußsche
$3\times3$-Matrizen, und warum eine Niveaustatistik, die an einem Spektrum mit vermischten Sektoren
berechnet wird, wie die unabhängiger Niveaus aussieht.

Im Folgenden sei $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ mit $D=\sum_ck_c$, $P_c$ seien die
minimalen zentralen Projektionen und $E$ die bedingte Erwartung auf das Zentrum aus
Proposition 2.8. $U(k)$ ist die Gruppe der unitären $k\times k$-Matrizen und
$U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}$ die unitäre Gruppe von $\mathcal M$. Zufallsvariablen
leben auf einem Wahrscheinlichkeitsraum mit Erwartungswert $\mathbb E$; der Erwartungswert einer
Zufallsmatrix wird eintragsweise genommen. Eine *standardnormalverteilte komplexe* Zufallsvariable
hat unabhängige Real- und Imaginärteile $\sim N(0,\tfrac12)$, sodass $\mathbb E|z|^2=1$ (das zieht
`torch.randn` für einen komplexen dtype laut seiner Dokumentation; der erste Codeblock überprüft die
zweiten Momente); eine *Ginibre-Matrix* hat unabhängige standardnormalverteilte komplexe Einträge.

## 5.1 Haar-Maß

Sei $G$ eine kompakte Gruppe von Matrizen (eine abgeschlossene Untergruppe eines $U(n)$). Ein
Borelsches Wahrscheinlichkeitsmaß $\mu$ auf $G$ heißt **linksinvariant**, wenn $\mu(gB)=\mu(B)$ für
alle $g\in G$ und alle Borelmengen $B$ gilt, äquivalent $\int f(gh)\,d\mu(h)=\int f\,d\mu$ für jede
beschränkte Borel-Funktion $f$; **rechtsinvariant** wird mit $Bg$ und $f(hg)$ definiert.

**Lemma 5.1 (Eindeutigkeit).** Ist $\mu$ ein linksinvariantes und $\nu$ ein rechtsinvariantes
Borelsches Wahrscheinlichkeitsmaß auf $G$, so ist $\mu=\nu$. Folglich trägt $G$ höchstens ein
linksinvariantes Borelsches Wahrscheinlichkeitsmaß, und falls es existiert, ist es auch
rechtsinvariant und invariant unter $g\mapsto g^{-1}$.

*Beweis.* Sei $f:G\to\mathbb R$ stetig (also beschränkt). Nach dem Satz von Fubini ist
$\int\!\big(\int f(gh)\,d\mu(h)\big)d\nu(g)=\int\!\big(\int f(gh)\,d\nu(g)\big)d\mu(h)$. Das innere
Integral links ist wegen der Linksinvarianz gleich $\int f\,d\mu$, das rechts wegen der
Rechtsinvarianz gleich $\int f\,d\nu$; also ist $\int f\,d\mu=\int f\,d\nu$. Zwei Borelsche
Wahrscheinlichkeitsmaße auf einem metrischen Raum mit denselben Integralen stetiger Funktionen
sind gleich (Korollar A.2: Die stetigen Funktionen $\max(0,1-n\operatorname{dist}(\cdot,K))$ fallen
gegen die Indikatorfunktion einer abgeschlossenen Menge $K$, und die abgeschlossenen Mengen
bestimmen das Maß nach dem Lemma von Dynkin, Lemma A.1). Sind $\mu_1,\mu_2$ linksinvariant, so ist
$\nu(B)=\mu_2(B^{-1})$ rechtsinvariant, da $(Bg)^{-1}=g^{-1}B^{-1}$; also ist $\mu_1=\nu=\mu_2$, und
$\mu_1$ ist rechtsinvariant und inversionsinvariant. $\square$

**Definition.** Das **Haar-Maß** von $G$ ist sein linksinvariantes Borelsches
Wahrscheinlichkeitsmaß. Ein zufälliges Element mit dieser Verteilung heißt **Haar-verteilt**.

Lemma 5.1 sagt nichts über die Existenz. Der Kurs braucht nur die Haar-Maße von $U(k)$, $U(1)$,
$U(\mathcal M)$, $SU(k)$, $O(k)$ und der Diagonalgruppen $T$ aus Lemma C.2, und jedes davon wird
explizit konstruiert: $U(k)$ in Satz 5.3, $U(\mathcal M)$ in Korollar 5.4, $SU(k)$ im Anschluss
daran, $O(k)$ in Satz B.17 und $T$ in Lemma C.2 (siehe Abschnitt B.7).

**Lemma 5.2 (QR-Zerlegung).** Jede invertierbare Matrix $Z\in M_k(\mathbb C)$ lässt sich als
$Z=QR$ schreiben, mit $Q$ unitär und $R$ eine obere Dreiecksmatrix mit positiver Diagonale. $Q$ und
$R$ sind eindeutig und hängen stetig von $Z$ ab.

*Beweis.* Man wende das Gram-Schmidt-Verfahren auf die Spalten $z_1,\dots,z_k$ von $Z$ an, die
linear unabhängig sind: $q_j=w_j/\|w_j\|$ mit $w_j=z_j-\sum_{i<j}\langle q_i,z_j\rangle q_i$; dabei
ist $w_j\neq0$, weil
$z_j\notin\operatorname{span}\{z_1,\dots,z_{j-1}\}=\operatorname{span}\{q_1,\dots,q_{j-1}\}$. Dann
ist $z_j=\sum_{i<j}\langle q_i,z_j\rangle q_i+\|w_j\|q_j$, d. h. $Z=QR$ mit $R_{jj}=\|w_j\|>0$, und
jeder Schritt ist auf den invertierbaren Matrizen stetig. Ist $Q_1R_1=Q_2R_2$, so ist
$T=Q_2^*Q_1=R_2R_1^{-1}$ unitär und eine obere Dreiecksmatrix mit positiver Diagonale. Ihre Inverse
$T^*$ ist eine untere und zugleich eine obere Dreiecksmatrix (die Inverse einer oberen
Dreiecksmatrix ist eine obere Dreiecksmatrix), also ist $T$ diagonal; eine unitäre
Diagonalmatrix mit positiver Diagonale ist $1$. $\square$

**Satz 5.3 (Haar-verteilte unitäre Matrizen aus Gaußschen Matrizen).** Sei $Z$ eine
$k\times k$-Ginibre-Matrix. Dann ist $Z$ mit Wahrscheinlichkeit eins invertierbar, und der Faktor
$Q$ aus Lemma 5.2 ist Haar-verteilt auf $U(k)$. Ist überdies $Z=Q'R'$ irgendeine QR-Zerlegung mit
$Q'$ unitär und $R'$ obere Dreiecksmatrix (etwa die von `torch.linalg.qr` gelieferte), so ist
$Q=Q'\Lambda$ mit $\Lambda=\operatorname{diag}(R'_{ii}/|R'_{ii}|)$ Haar-verteilt.

*Beweis.* *Invertierbarkeit.* $|\det Z|^2$ ist ein reelles Polynom in den $2k^2$ reellen
Koordinaten von $Z$, das nicht identisch null ist (es ist gleich $1$ in $Z=1$), und $Z$ hat die
Dichte $\pi^{-k^2}e^{-\|Z\|_2^2}$ (Proposition A.16(b)); also ist $P(\det Z=0)=0$ nach
Korollar A.7(b).

*Invarianz der Gaußschen Verteilung.* Für $V\in U(k)$ ist die Abbildung $Z\mapsto VZ$ reell-linear
und erhält $\|Z\|_2$ (Abschnitt 1.4); sie ist also eine orthogonale Abbildung von
$M_k(\mathbb C)\cong\mathbb R^{2k^2}$, die das Lebesgue-Maß (Lemma A.3(b)) und die Dichte erhält:
$VZ$ hat dieselbe Verteilung wie $Z$ (Proposition A.16(b)).

*Äquivarianz.* Ist $Z=QR$, so ist $VZ=(VQ)R$ die Zerlegung von $VZ$ aus Lemma 5.2, also
$Q(VZ)=VQ(Z)$. Da $Q$ auf den invertierbaren Matrizen stetig ist (Lemma 5.2), ist $Q(Z)$ ein
zufälliges Element von $U(k)$, und für jede Borelmenge $B\subseteq U(k)$ gilt
$P(VQ(Z)\in B)=P(Q(VZ)\in B)=P(Q(Z)\in B)$: Die Verteilung von $Q(Z)$ ist linksinvariant und nach
Lemma 5.1 das Haar-Maß (das beweist auch, dass das Haar-Maß von $U(k)$ existiert). Schließlich ist
$Q'R'=(Q'\Lambda)(\Lambda^*R')$, und $\Lambda^*R'$ ist eine obere Dreiecksmatrix mit Diagonale
$|R'_{ii}|>0$, also $Q'\Lambda=Q(Z)$ wegen der Eindeutigkeit. $\square$

Derselbe Beweis mit reellen Gaußschen Matrizen liefert das Haar-Maß der orthogonalen Gruppe
$O(k)$, mit der Korrektur $\operatorname{diag}(\operatorname{sgn}R'_{ii})$; das wird in Satz B.17
ausgeführt. Ohne die Korrektur $\Lambda$ hängt die Verteilung von $Q'$ von der
Vorzeichenkonvention der QR-Routine ab. Für die von `torch.linalg.qr` aufgerufene Routine (LAPACK)
findet der folgende Code $\mathbb E\operatorname{Re}\operatorname{Tr}Q'$ viele Standardfehler von
$0$ entfernt, während eine Haar-verteilte unitäre Matrix $\mathbb EU=0$ hat (Korollar 5.5): Dieses
$Q'$ ist also nicht Haar-verteilt. Das ist eine Aussage über diese Routine, vom Code überprüft,
kein Satz (das Rezept mit $\Lambda$ stammt von F. Mezzadri, *How to generate random matrices from
the classical compact groups*, Notices AMS 54 (2007), 592–604).

**Korollar 5.4 (Haar-Maß der Algebra).** $U(\mathcal M)=\{\bigoplus_cU_c:U_c\in U(k_c)\}\cong\prod_cU(k_c)$,
und $U=\bigoplus_cU_c$ ist genau dann Haar-verteilt auf $U(\mathcal M)$, wenn die Blöcke $U_c$
unabhängig sind und jedes $U_c$ Haar-verteilt auf $U(k_c)$ ist.

*Beweis.* $U^*U=\bigoplus_cU_c^*U_c$, also ist $U$ genau dann unitär, wenn jeder Block es ist. Seien
die $U_c$ unabhängig und Haar-verteilt. Für $V=\bigoplus_cV_c\in U(\mathcal M)$ sind die Blöcke
$V_cU_c$ von $VU$ unabhängig (Funktionen unabhängiger Zufallsvariablen), und $V_cU_c$ hat die
Verteilung von $U_c$; also hat $VU$ die Verteilung von $U$. Die Verteilung von $U$ ist daher
linksinvariant und nach Lemma 5.1 das Haar-Maß der kompakten Gruppe $U(\mathcal M)$. Umgekehrt hat
wegen der Eindeutigkeit jedes Haar-verteilte $U$ diese Verteilung. $\square$

Die Gruppe $U(\mathcal M)$ ist viel kleiner als $U(D)$: Für $C\ge2$ liegt eine Haar-verteilte
unitäre Matrix $U$ aus $U(D)$ mit Wahrscheinlichkeit null in $\mathcal M$. Denn nach Satz 5.3 und
Lemma 5.1 hat $U$ die Verteilung des Faktors $Q$ einer $D\times D$-Ginibre-Matrix $Z$, und nach der
Gram-Schmidt-Konstruktion in Lemma 5.2 ist die erste Spalte von $Q$ gleich $z/\|z\|$ für die erste
Spalte $z$ von $Z$, einen standardnormalverteilten komplexen Gauß-Vektor in $\mathbb C^D$ mit einer
Dichte (Proposition A.16(a)). Für jedes $i$ ist $\{z_i=0\}$ ein echter Unterraum von $\mathbb C^D$,
also $P(z_i=0)=0$ (Korollar A.7(c)), und alle Einträge der ersten Spalte sind fast sicher von null
verschieden. Ein Element von $\mathcal M$ hat in der ersten Spalte außerhalb der Zeilen des ersten
Blocks Nulleinträge, und diese Zeilen existieren für $C\ge2$.

Das Haar-Maß von $SU(k)$ erhält man auf dieselbe Weise: Ist $U$ Haar-verteilt auf $U(k)$ und
$\theta=\arg\det U\in(-\pi,\pi]$ (eine Borel-Funktion von $U$), so ist $V=e^{-i\theta/k}U\in SU(k)$,
und für $W\in SU(k)$ liefert $\det(WU)=\det U$ die Beziehung $V(WU)=WV(U)$; da $WU$ die Verteilung
von $U$ hat, ist die Verteilung von $V$ linksinvariant auf $SU(k)$, nach Lemma 5.1 also das
Haar-Maß von $SU(k)$.

**Korollar 5.5 (erste Momente der Einträge).** Ist $U$ Haar-verteilt auf $U(k)$, so gilt
$\mathbb EU_{ij}=0$ und $\mathbb E|U_{ij}|^2=1/k$.

*Beweis.* $-1\in U(k)$, also hat $-U$ die Verteilung von $U$, und $\mathbb EU=-\mathbb EU$. Für eine
Permutationsmatrix $P$ hat $PU$ die Verteilung von $U$, also hängt $\mathbb E|U_{ij}|^2$ nicht von
$i$ ab; wegen $\sum_i|U_{ij}|^2=1$ ist es gleich $1/k$. $\square$

Die Bibliothek zieht Haar-verteilte unitäre Elemente von $U(\mathcal M)$ mit
`random_unitary_operator`, mit einer einzigen gestapelten QR-Zerlegung pro Blockgröße für alle
Sektoren und Stichproben zugleich (auf einer GPU mit Gram-Schmidt, was wegen der Eindeutigkeit in
Lemma 5.2 dasselbe $Q$ liefert). Das Maß `'haar'` ist Satz 5.3 in jedem Block; für eine komplexe
Algebra multipliziert `'haar_su'` wie oben mit $e^{-i\theta/k}$; `'coe'` ist das Ensemble
$S=W^TW$ aus Aufgabe 3 (die übrigen Maße der Bibliothek werden in diesem Kurs nicht verwendet). Wir
überprüfen die Korollare 5.4 und 5.5 sowie, dass der Sampler der Bibliothek genau die korrigierte
QR-Zerlegung aus Satz 5.3 ist:

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinChain, channels

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def mc(samples, exact, name):
    """Monte Carlo mean of real samples with its standard error; must agree within 5 standard errors."""
    mean, err = samples.mean().item(), samples.std().item() / math.sqrt(samples.numel())
    print(f"{name}: {mean:.4f} +- {err:.4f}   (theory {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12
    return mean

alg = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')        # M_1 (+) M_2 (+) M_3, D = 6
B = 100_000
U = alg.random_unitary_operator(B)                                   # Haar on U(1) x U(2) x U(3)
assert close(U.adjoint() @ U, alg.identity())                        # unitary elements of the algebra
u = [U.matrix[:, c, :k, :k] for c, k in enumerate(alg.k_factors)]   # the blocks U_c
for c, k in enumerate(alg.k_factors):                                # Corollary 5.5
    mc(u[c][:, 0, 0].real, 0.0, f"E Re (U_{c})_11")
    mc(u[c][:, 0, 0].abs() ** 2, 1 / k, f"E |(U_{c})_11|^2")
mc(u[1][:, 0, 0].abs() ** 2 * u[2][:, 0, 0].abs() ** 2, 1 / 6, "E |(U_1)_11|^2 |(U_2)_11|^2")   # independent blocks

z = torch.randn(B, dtype=dt)                                         # torch.randn draws standard complex Gaussians:
mc(z.real ** 2, 0.5, "E (Re z)^2")                                   # Re z, Im z ~ N(0, 1/2), uncorrelated
mc(z.imag ** 2, 0.5, "E (Im z)^2")
mc(z.real * z.imag, 0.0, "E Re z Im z")

SU = alg.random_unitary_operator(4, measure='haar_su').matrix
assert torch.allclose(torch.linalg.det(SU[:, 2, :3, :3]), torch.ones(4, dtype=dt))

torch.manual_seed(1)
Z = torch.randn(B, 2, 2, dtype=dt)                                   # Ginibre matrices
torch.manual_seed(1)
V = alg.random_unitary(2, batch_size=B)                              # the library, same random numbers
Q, R = torch.linalg.qr(Z)                                            # LAPACK: diag(R) real, either sign
d = R.diagonal(dim1=-2, dim2=-1)
assert close(V, Q * (d / d.abs()).unsqueeze(-2))                     # Theorem 5.3: V = Q Lambda
raw = Q.diagonal(dim1=-2, dim2=-1).sum(-1).real                      # Re Tr Q without the correction
print(f"E Re Tr Q without the phase correction: {raw.mean():.4f} +- {raw.std() / math.sqrt(B):.4f}")
assert abs(raw.mean()) > 20 * raw.std() / math.sqrt(B)               # far from 0: not Haar
```

*Numerische Beobachtung.* Ohne die Korrektur findet der Code
$\mathbb E\operatorname{Re}\operatorname{Tr}Q'\approx-0.84$ für $k=2$, während eine Haar-verteilte
unitäre Matrix $\mathbb EU=0$ hat.

## 5.2 Momente und der Weingarten-Kalkül

**Lemma 5.6.** Kommutiert $T\in M_k(\mathbb C)$ mit allen unitären Diagonalmatrizen und allen
Permutationsmatrizen, so ist $T\in\mathbb C1$. Insbesondere ist die Kommutante von $U(k)$ gleich
$\mathbb C1$.

*Beweis.* Für $\Delta=\operatorname{diag}(e^{i\theta_1},\dots,e^{i\theta_k})$ mit paarweise
verschiedenen $e^{i\theta_j}$ liefert $\Delta T=T\Delta$ die Gleichung
$(e^{i\theta_i}-e^{i\theta_j})T_{ij}=0$, also ist $T$ diagonal. Für die Permutationsmatrix
$P_\sigma e_j=e_{\sigma(j)}$ liefert $P_\sigma TP_\sigma^{-1}=T$ die Gleichung
$T_{\sigma(j)\sigma(j)}=T_{jj}$, also sind alle Diagonaleinträge gleich. $\square$

**Proposition 5.7 (der Haar-Twirl ist die bedingte Erwartung auf das Zentrum).** Sei $U$
Haar-verteilt auf $U(\mathcal M)$. Für jedes $X\in\mathcal M$ gilt
$$
\mathbb E\,[UXU^*]=E(X)=\bigoplus_c\frac{\operatorname{Tr}X_c}{k_c}\,1_{k_c}.
$$

*Beweis.* $(UXU^*)_c=U_cX_cU_c^*$, und $U_c$ ist Haar-verteilt auf $U(k_c)$ (Korollar 5.4). Sei
$T_c=\mathbb E[U_cX_cU_c^*]$ (die Einträge sind beschränkt). Für $V\in U(k_c)$ hat $VU_c$ die
Verteilung von $U_c$, also ist $VT_cV^*=\mathbb E[(VU_c)X_c(VU_c)^*]=T_c$. Nach Lemma 5.6 ist
$T_c=\lambda1$, und
$k_c\lambda=\operatorname{Tr}T_c=\mathbb E\operatorname{Tr}(U_cX_cU_c^*)=\operatorname{Tr}X_c$. $\square$

Die Normierung ist die von Proposition 2.8 und von `channels.center_expectation`. Kommutiert
$X\in\mathcal M$ mit $U(\mathcal M)$, so ist $UXU^*=X$, also ist $X=E(X)$ zentral; das Haar-Mittel
ist das kontinuierliche Gegenstück zur Gruppenmittelung aus Kapitel 2, Aufgabe 3. Ist $U$
stattdessen Haar-verteilt auf $U(D)$, so liefert derselbe Beweis (mit $C=1$)
$\mathbb E[UXU^*]=\frac{\operatorname{Tr}X}{D}1$: Die Sektorgewichte gehen verloren, und $UXU^*$
muss nicht einmal in $\mathcal M$ liegen (für $C\ge2$, $X=P_1$ und die Drehung $U$ um $\pi/4$ in
der Ebene von $e_1$ und $e_{k_1+1}$ hat $UP_1U^*$ an der Position $(1,k_1+1)$ den Eintrag
$\frac12$).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Ein Haar-zufälliges unitäres Element der Algebra ist eine zufällige
> Zeitentwicklung, die die Superauswahlregel respektiert. Die Mittelung darüber verwürfelt alles
> innerhalb jedes Sektors und bewahrt genau den klassischen Teil: das Gewicht jeder Ladung. Das
> ist der „Twirl“ der Quanteninformation, und es ist die bedingte Erwartung auf das Zentrum aus
> Kapitel 2.

Für die zweiten Momente brauchen wir die Kommutante von $\{V\otimes V\}$ auf
$\mathbb C^k\otimes\mathbb C^k$. Sei $F$ die **Vertauschung** (*swap*), $F(x\otimes y)=y\otimes x$;
sie ist selbstadjungiert und unitär, und es gilt
$\operatorname{Tr}1=k^2$, $\operatorname{Tr}F=\sum_{i,j}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k$, $F^2=1$.

**Lemma 5.8.** Sei $k\ge2$. Eine Matrix $X$ auf $\mathbb C^k\otimes\mathbb C^k$ kommutiert genau
dann mit $V\otimes V$ für alle $V\in U(k)$, wenn $X=a1+bF$ mit $a,b\in\mathbb C$; und $1,F$ sind
linear unabhängig.

*Beweis.* $F(V\otimes V)(x\otimes y)=Vy\otimes Vx=(V\otimes V)F(x\otimes y)$, also kommutieren $1$
und $F$ mit allen $V\otimes V$. Umgekehrt kommutiere $X$ mit ihnen; man schreibe
$X_{ij,pq}=\langle e_i\otimes e_j,X(e_p\otimes e_q)\rangle$.
*Unitäre Diagonalmatrizen:* $\Delta\otimes\Delta$ multipliziert $e_i\otimes e_j$ mit
$e^{i(\theta_i+\theta_j)}$, also gilt $(e^{i(\theta_i+\theta_j)}-e^{i(\theta_p+\theta_q)})X_{ij,pq}=0$.
Sind die Multimengen $\{i,j\}$ und $\{p,q\}$ verschieden, so kommt ein Index $m$ in der ersten
$n\in\{\pm1,\pm2\}$-mal häufiger vor; mit $\theta_m=\pi/2$ und allen anderen $\theta=0$ ist die
Klammer gleich $e^{i(\theta_p+\theta_q)}(e^{in\pi/2}-1)\neq0$, also $X_{ij,pq}=0$. Folglich ist
$X(e_i\otimes e_i)=a_ie_i\otimes e_i$ und für $i\neq j$
$X(e_i\otimes e_j)=b_{ij}\,e_i\otimes e_j+c_{ij}\,e_j\otimes e_i$.
*Permutationen:* Wendet man $X(P_\sigma\otimes P_\sigma)=(P_\sigma\otimes P_\sigma)X$ auf
$e_i\otimes e_j$ an, so erhält man $a_{\sigma(i)}=a_i$, $b_{\sigma(i)\sigma(j)}=b_{ij}$,
$c_{\sigma(i)\sigma(j)}=c_{ij}$. Die Permutationen wirken transitiv auf geordneten Paaren
verschiedener Indizes, also ist $a_i=a'$, $b_{ij}=b'$, $c_{ij}=c'$; somit ist $X=b'1+c'F$ auf dem
Spann der $e_i\otimes e_j$ mit $i\ne j$.
*Eine Drehung:* Sei $R$ die unitäre Matrix mit $Re_1=(e_1+e_2)/\sqrt2$, $Re_2=(-e_1+e_2)/\sqrt2$, $Re_m=e_m$ für $m\ge3$.
Dann ist $(R\otimes R)(e_1\otimes e_1)=\frac12(e_1\otimes e_1+e_1\otimes e_2+e_2\otimes e_1+e_2\otimes e_2)$, und
$X(R\otimes R)(e_1\otimes e_1)=\frac12\big(a'e_1\otimes e_1+(b'+c')(e_1\otimes e_2+e_2\otimes e_1)+a'e_2\otimes e_2\big)$
muss gleich $(R\otimes R)X(e_1\otimes e_1)=a'(R\otimes R)(e_1\otimes e_1)$ sein. Der Vergleich der
Koeffizienten von $e_1\otimes e_2$ ergibt $a'=b'+c'$, und das ist der Wert von $b'1+c'F$ auf
$e_i\otimes e_i$. Also ist $X=b'1+c'F$. Schließlich bilden $1$ und $F$ den Vektor $e_1\otimes e_2$
auf die linear unabhängigen Vektoren $e_1\otimes e_2$ und $e_2\otimes e_1$ ab. $\square$

**Korollar 5.9 (Schur-Weyl-Dualität zweiter Ordnung).** Sei $k\ge2$ und $P_\pm=\frac12(1\pm F)$.
Dies sind die Projektionen auf den symmetrischen und den antisymmetrischen Unterraum
$\mathrm{Sym}=\{x:Fx=x\}$ und $\mathrm{Alt}=\{x:Fx=-x\}$ der Dimensionen $k(k+1)/2$ und
$k(k-1)/2$. Die Darstellung $V\mapsto V\otimes V$ von $U(k)$ ist auf $\mathrm{Sym}$ und auf
$\mathrm{Alt}$ irreduzibel, und die beiden Darstellungen sind nicht äquivalent.

*Beweis.* Wegen $F=F^*$ und $F^2=1$ sind $P_\pm$ Projektionen mit $P_++P_-=1$ und Bildern
$\mathrm{Sym}$, $\mathrm{Alt}$; ihre Dimensionen sind $\operatorname{Tr}P_\pm=(k^2\pm k)/2$. Der
Spann $\mathcal A$ von $\{V\otimes V\}$ ist eine $*$-Algebra ($1\otimes1$, Produkte
$VW\otimes VW$, Adjungierte $V^*\otimes V^*$) mit Kommutante
$\mathcal A'=\operatorname{span}\{1,F\}=\operatorname{span}\{P_+,P_-\}$ nach Lemma 5.8. Sei
$K\subseteq\mathrm{Sym}$ invariant unter allen $V\otimes V$. Nach Lemma 2.1 liegt seine Projektion
in $\mathcal A'$, ist also $xP_++yP_-$ mit $x,y\in\{0,1\}$ (sie ist eine Projektion); sie
verschwindet auf $\mathrm{Alt}$, also ist $y=0$, und $K=0$ oder $K=\mathrm{Sym}$. Dasselbe Argument
gilt für $\mathrm{Alt}$. Verflechtet $T:\mathrm{Sym}\to\mathrm{Alt}$ die beiden Darstellungen, so
kommutiert $TP_+$ (aufgefasst als Operator auf $\mathbb C^k\otimes\mathbb C^k$) mit $\mathcal A$,
also ist $TP_+=xP_++yP_-$; es verschwindet auf $\mathrm{Alt}$, also ist $y=0$, und sein Bild liegt
in $\mathrm{Alt}\cap\mathrm{Sym}=0$, also ist $T=0$. $\square$

**Satz 5.10 (Weingarten-Formel zweiter Ordnung).** Sei $k\ge2$, $U$ Haar-verteilt auf $U(k)$ und
$X$ eine Matrix auf $\mathbb C^k\otimes\mathbb C^k$. Dann gilt
$$
\mathbb E\,[(U\otimes U)X(U\otimes U)^*]=\alpha\,1+\beta F,\qquad
\alpha=\frac{\operatorname{Tr}X-\operatorname{Tr}(FX)/k}{k^2-1},\quad
\beta=\frac{\operatorname{Tr}(FX)-\operatorname{Tr}X/k}{k^2-1}.
$$
Äquivalent gilt mit $\mathrm{Wg}(e)=\frac1{k^2-1}$, $\mathrm{Wg}(\tau)=-\frac1{k(k^2-1)}$ auf der
Gruppe $S_2=\{e,\tau\}$
$$
\mathbb E\big[U_{i_1j_1}U_{i_2j_2}\overline{U_{i'_1j'_1}U_{i'_2j'_2}}\big]
=\sum_{\sigma,\pi\in S_2}\delta_{i_1i'_{\sigma(1)}}\delta_{i_2i'_{\sigma(2)}}\,
\delta_{j_1j'_{\pi(1)}}\delta_{j_2j'_{\pi(2)}}\,\mathrm{Wg}(\sigma\pi^{-1}).
$$

*Beweis.* Sei $T=\mathbb E[(U\otimes U)X(U\otimes U)^*]$. Wie in Proposition 5.7 kommutiert $T$
mit allen $V\otimes V$, also ist $T=\alpha1+\beta F$ nach Lemma 5.8. Bildet man Spuren, so ist
$\operatorname{Tr}T=\operatorname{Tr}X$ und, da $F$ mit $U\otimes U$ kommutiert,
$\operatorname{Tr}(FT)=\mathbb E\operatorname{Tr}((U\otimes U)^*F(U\otimes U)X)=\operatorname{Tr}(FX)$.
Das ergibt $k^2\alpha+k\beta=\operatorname{Tr}X$ und $k\alpha+k^2\beta=\operatorname{Tr}(FX)$, ein
lineares Gleichungssystem mit Determinante $k^2(k^2-1)\neq0$, dessen Lösung die angegebene ist.
Für die Einträge nehme man $X=|e_{j_1}\otimes e_{j_2}\rangle\langle e_{j'_1}\otimes e_{j'_2}|$: Dann
ist $\langle e_{i_1}\otimes e_{i_2},T\,e_{i'_1}\otimes e_{i'_2}\rangle$ die linke Seite,
$\operatorname{Tr}X=\delta_{j_1j'_1}\delta_{j_2j'_2}$ (der Term $\pi=e$),
$\operatorname{Tr}(FX)=\delta_{j_1j'_2}\delta_{j_2j'_1}$ (der Term $\pi=\tau$), und $1$, $F$
tragen $\delta_{i_1i'_1}\delta_{i_2i'_2}$ (der Term $\sigma=e$) und
$\delta_{i_1i'_2}\delta_{i_2i'_1}$ ($\sigma=\tau$) bei. Fasst man die Terme zusammen, so ist
$\alpha=\sum_\pi\mathrm{Wg}(e\pi^{-1})[\pi]$ und $\beta=\sum_\pi\mathrm{Wg}(\tau\pi^{-1})[\pi]$,
wobei $[\pi]$ die $j$-Deltas von $\pi$ bezeichnet. $\square$

Die Matrix $\big(\mathrm{Wg}(\sigma^{-1}\pi)\big)_{\sigma,\pi}$ ist die Inverse der Gram-Matrix
$\begin{pmatrix}k^2&k\\k&k^2\end{pmatrix}$ von $1$ und $F$ bezüglich
$\langle A,B\rangle=\operatorname{Tr}(A^*B)$; in dieser Form lässt sich das Ergebnis
verallgemeinern.

**Korollar 5.11.** Für $k\ge2$ gilt:
$$
\mathbb E|U_{11}|^4=\frac2{k(k+1)},\qquad
\mathbb E|U_{11}|^2|U_{12}|^2=\frac1{k(k+1)},\qquad
\mathbb E\big[U_{11}U_{22}\overline{U_{12}U_{21}}\big]=-\frac1{(k-1)k(k+1)}.
$$

*Beweis.* In Satz 5.10 seien (i) alle Indizes gleich $1$: Jedes Produkt von Deltas ist $1$, also
ist das Moment $2(\mathrm{Wg}(e)+\mathrm{Wg}(\tau))=2\frac{k-1}{k(k^2-1)}$. (ii) $i=(1,1)$,
$j=(1,2)$, $i'=(1,1)$, $j'=(1,2)$: Nur $\pi=e$ bleibt übrig, was $\mathrm{Wg}(e)+\mathrm{Wg}(\tau)$
ergibt. (iii) $i=(1,2)$, $j=(1,2)$, $i'=(1,2)$, $j'=(2,1)$: Nur $\sigma=e$, $\pi=\tau$ bleiben
übrig, was $\mathrm{Wg}(\tau)$ ergibt. $\square$

Das letzte Moment ist negativ: Die Einträge einer Haar-verteilten unitären Matrix sind nicht
unabhängig, und die Korrelationen sind von relativer Ordnung $1/k$. Für die Algebra sind Einträge
verschiedener Blöcke unabhängig (Korollar 5.4), sodass gemischte Momente über Sektoren hinweg in
Momente einzelner Blöcke faktorisieren.

Für Momente beliebiger Ordnung greift dieselbe Überlegung, sobald die Kommutante von
$\{V^{\otimes n}\}$ bekannt ist; sie wird von den Permutationen der $n$ Tensorfaktoren aufgespannt
(Schur-Weyl-Dualität, Satz B.6(c)). Anhang B führt das aus.

**Satz 5.12 (Weingarten-Formel beliebiger Ordnung).** Sei $n\ge1$ und $U$ Haar-verteilt auf
$U(k)$. Dann gilt
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\pi\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\pi(m)}}\;\mathrm{Wg}(\sigma\pi^{-1},k),
$$
wobei $\mathrm{Wg}(\cdot,k)$ eine Funktion auf $S_n$ ist, die nur vom Zyklentyp abhängt. Für
$n\le k$ ist die Matrix $\big(\mathrm{Wg}(\sigma^{-1}\pi,k)\big)_{\sigma,\pi}$ die Inverse der
Gram-Matrix $\big(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)}\big)_{\sigma,\pi}$; für $n>k$ ist diese
Matrix singulär, und man nimmt ihre Pseudoinverse. Momente mit verschieden vielen Faktoren $U$ und
$\overline U$ verschwinden.

*Beweis.* Das ist Satz B.13 mit $d=k$: Teil (a) liefert $W_{\sigma\pi}=\mathrm{Wg}(\sigma^{-1}\pi)$
mit einer Funktion $\mathrm{Wg}$ des Zyklentyps, Teil (b) die Invertierbarkeit der Gram-Matrix genau
für $n\le k$, Teil (d) die Formel für die Momente und Teil (e) das Verschwinden unbalancierter
Momente. $\square$

(Die Formel stammt von D. Weingarten, J. Math. Phys. 19 (1978), 999–1001, und B. Collins, Int.
Math. Res. Not. 2003, 953–982; siehe auch B. Collins, P. Śniady, Comm. Math. Phys. 264 (2006),
773–795.)

Wir überprüfen Proposition 5.7 an `channels.center_expectation`, Lemma 5.8 mit linearer Algebra
(die Kommutante von drei gezogenen $V\otimes V$, als Kern berechnet, enthält $F$ und ist
$2$-dimensional, nach Lemma 5.8 also bereits $\operatorname{span}\{1,F\}$) sowie Satz 5.10 und
Korollar 5.11 per Monte-Carlo mit dem $U(3)$-Block der obigen Stichproben:

```python
X = alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])      # a fixed element
T = (U @ X @ U.adjoint()).matrix                                    # B samples of U X U*
mean, err = T.mean(0), T.std(0) / math.sqrt(B)
EX = channels.center_expectation(alg)(X).matrix[0]
random = err > 1e-12                                                 # entries that are not constant
z = ((mean - EX).abs()[random] / err[random]).max().item()
print(f"twirl over U(M) vs E(X): largest deviation {z:.2f} standard errors")
assert ((mean - EX).abs() <= 5 * err + 1e-12).all()                 # Proposition 5.7

def commutant_dim(ops, tol=1e-9):
    """Dimension of {X : XA = AX for all A in ops}, as the null space of X -> (AX - XA)_A."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    ops = [A.contiguous() for A in ops]                              # torch.kron needs contiguous input
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return int((torch.linalg.svdvals(L) < tol).sum())

k = 3
F = torch.zeros(k * k, k * k, dtype=dt)
for i in range(k):
    for j in range(k):
        F[j * k + i, i * k + j] = 1                                  # F (e_i (x) e_j) = e_j (x) e_i
VV = [torch.kron(W, W) for W in alg.random_unitary(k, batch_size=3).contiguous()]
assert all(close(F @ W, W @ F) for W in VV)
assert commutant_dim(VV) == 2                                        # Lemma 5.8: span{1, F}

u3 = u[2]                                                            # Haar on U(3)
mc(u3[:, 0, 0].abs() ** 4, 2 / (k * (k + 1)), "E|U_11|^4")
mc(u3[:, 0, 0].abs() ** 2 * u3[:, 0, 1].abs() ** 2, 1 / (k * (k + 1)), "E|U_11|^2 |U_12|^2")
mc((u3[:, 0, 0] * u3[:, 1, 1] * (u3[:, 0, 1] * u3[:, 1, 0]).conj()).real,
   -1 / ((k - 1) * k * (k + 1)), "E U_11 U_22 conj(U_12 U_21)")

n_s = 20_000                                                         # Theorem 5.10 for a random X
X2 = torch.randn(k * k, k * k, dtype=dt)
trX, trFX = X2.diagonal().sum(), (F @ X2).diagonal().sum()
alpha, beta = (trX - trFX / k) / (k ** 2 - 1), (trFX - trX / k) / (k ** 2 - 1)
UU = torch.einsum('bij,bkl->bikjl', u3[:n_s], u3[:n_s]).reshape(n_s, k * k, k * k)   # U (x) U
T2 = UU @ X2 @ UU.mH
mean, err = T2.mean(0), T2.std(0) / math.sqrt(n_s)
assert ((mean - alpha * torch.eye(k * k, dtype=dt) - beta * F).abs() <= 5 * err).all()
```

## 5.3 Niveaustatistik

Sei $H$ eine zufällige selbstadjungierte $n\times n$-Matrix mit Eigenwerten
$\lambda_1\le\dots\le\lambda_n$. Ihre **Abstände** sind $s_i=\lambda_{i+1}-\lambda_i$. Zwei
Ensembles sind Standard: das **GOE**, $H=(G+G^T)/2$ mit reellem $G$ mit unabhängigen
$N(0,1)$-Einträgen, und das **GUE**, $H=(G+G^*)/\sqrt2$ mit einer Ginibre-Matrix $G$
(Abschnitt 1.7). In beiden Fällen sind die Diagonaleinträge $N(0,1)$-verteilt, die Einträge oberhalb
der Diagonale im GOE $N(0,\tfrac12)$-verteilt und haben im GUE unabhängige Real- und Imaginärteile
mit Verteilung $N(0,\tfrac12)$, die Einträge $h_{ij}$, $i\le j$, sind unabhängig, und die Dichte
von $H$ ist proportional zu $e^{-\operatorname{Tr}H^2/2}$ (Proposition A.17; die Proportionalität
folgt aus $\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$).

**Proposition 5.13 (Wignersche Vermutung).** Sei $n=2$, $s=\lambda_2-\lambda_1$ und
$S=s/\mathbb Es$. Für das GOE ($\beta=1$) und das GUE ($\beta=2$) hat $S$ die Dichten
$$
p_1(S)=\frac\pi2\,S\,e^{-\pi S^2/4},\qquad p_2(S)=\frac{32}{\pi^2}\,S^2\,e^{-4S^2/\pi},
$$
mit den Verteilungsfunktionen $P(S\le t)=1-e^{-\pi t^2/4}$ und
$P(S\le t)=\operatorname{erf}(2t/\sqrt\pi)-\frac{4t}\pi e^{-4t^2/\pi}$. Insbesondere ist
$P(S\le\varepsilon)$ für $\varepsilon\to0$ von der Ordnung $\varepsilon^{\beta+1}$.

*Beweis.* Für $H=\begin{pmatrix}a&b\\\bar b&c\end{pmatrix}$ sind die Eigenwerte
$\frac{a+c}2\pm\frac12\sqrt{(a-c)^2+4|b|^2}$, also $s=\sqrt{(a-c)^2+4|b|^2}$. Nach Proposition A.17
sind $a$, $c$ und die reellen Koordinaten von $b$ unabhängig, $a,c\sim N(0,1)$, und jede reelle
Koordinate von $b$ ist $N(0,\frac12)$-verteilt; also ist $a-c\sim N(0,2)$ (Proposition A.14(c)), und
$2\operatorname{Re}b$, $2\operatorname{Im}b$ sind $N(0,2)$-verteilt, alle unabhängig. *GOE:* $x=a-c$
und $y=2b$ sind unabhängige $N(0,2)$-Variablen, und in Polarkoordinaten (Korollar A.19 mit $n=2$,
$s_2=2\pi$) ist
$P(s\le t)=\int_{x^2+y^2\le t^2}\frac{e^{-(x^2+y^2)/4}}{4\pi}\,dx\,dy=\int_0^t\frac r2e^{-r^2/4}dr=1-e^{-t^2/4}$.
Dann ist $\mathbb Es=\int_0^\infty P(s>t)\,dt=\int_0^\infty e^{-t^2/4}dt=\sqrt\pi$ und
$P(S\le t)=P(s\le\sqrt\pi t)=1-e^{-\pi t^2/4}$; man differenziere. *GUE:* $x=a-c$,
$y=2\operatorname{Re}b$, $y'=2\operatorname{Im}b$ sind unabhängig $N(0,2)$-verteilt, und nach
Satz A.20 ($n=3$, $s_3=4\pi$) ist die Dichte von $s=\sqrt{x^2+y^2+y'^2}$ gleich
$4\pi r^2(4\pi)^{-3/2}e^{-r^2/4}=\frac{r^2e^{-r^2/4}}{2\sqrt\pi}$. Ihr Mittelwert ist
$\frac1{2\sqrt\pi}\int_0^\infty r^3e^{-r^2/4}dr=\frac{8}{2\sqrt\pi}=\frac4{\sqrt\pi}$. Die
Substitution $r=4S/\sqrt\pi$ liefert $p_2$, und partielle Integration von $p_2$
($\int x^2e^{-ax^2}dx=-\frac{x}{2a}e^{-ax^2}+\frac1{2a}\int e^{-ax^2}dx$ mit $a=4/\pi$) liefert die
Verteilungsfunktion. Das Verhalten für kleines $\varepsilon$ folgt aus
$p_\beta(S)\sim cS^\beta$. $\square$

Der Exponent zählt Dimensionen: Eine reelle symmetrische $2\times2$-Matrix ist genau dann
entartet, wenn zwei reelle Bedingungen gelten ($a=c$, $b=0$), eine komplexe hermitesche genau dann,
wenn drei gelten. Für jedes $n$ haben die geordneten Eigenwerte $\lambda_1<\dots<\lambda_n$ des GOE
($\beta=1$) und des GUE ($\beta=2$) die gemeinsame Dichte
$$
Z_{n,\beta}^{-1}\prod_{i<j}(\lambda_j-\lambda_i)^\beta\,e^{-\sum_i\lambda_i^2/2}
$$
(Satz C.7, in Anhang C mit der Weylschen Integrationsformel bewiesen), und der kleinste Abstand
$\delta$ erfüllt $c\,\varepsilon^{\beta+1}\le P(\delta<\varepsilon)\le C\varepsilon^{\beta+1}$ für
$0<\varepsilon\le1$, mit Konstanten $0<c\le C$, die von $n$ und $\beta$ abhängen (Proposition C.9).
Der Limes $n\to\infty$ wird in diesem Kurs nicht behandelt (siehe M. L. Mehta, *Random Matrices*,
3. Aufl., Elsevier, 2004).

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Wigner schlug vor,
> dass sich die Energieniveaus komplexer Quantensysteme (schwerer Kerne) statistisch wie die
> Eigenwerte eines zufälligen Hamiltonoperators mit der richtigen Symmetrie verhalten: reell
> symmetrisch, wenn es eine Zeitumkehrsymmetrie gibt ($\beta=1$), andernfalls komplex hermitesch
> ($\beta=2$). Die Niveaus *stoßen* sich dann *ab*: Zwei Niveaus derselben Symmetrie kommen einander
> fast nie nahe (von Neumann und Wigner, 1929). Abstoßung ist die typische Signatur von
> Quantenchaos, unabhängige Niveaus (ohne Abstoßung) sind die integrabler oder lokalisierter
> Systeme; das sind die Vermutungen von Bohigas-Giannoni-Schmit und von Berry-Tabor aus der Physik,
> die in diesem Kurs weder bewiesen noch verwendet werden.

Die Abstände hängen von der lokalen Niveaudichte ab, die entlang des Spektrums variiert; um sie mit
$p_\beta$ zu vergleichen, braucht man ein **Unfolding** (Entfaltung), die Variablentransformation
$\lambda\mapsto\bar N(\lambda)$ mit einer glatten Näherung $\bar N$ der Zählfunktion, die den
mittleren Abstand zu $1$ macht. Das Ergebnis hängt von der Wahl von $\bar N$ ab. Das Verhältnis
aufeinanderfolgender Abstände vermeidet diese Wahl.

**Definition.** Das **Abstandsverhältnis** ist
$r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})\in[0,1]$ (V. Oganesyan, D. A. Huse, Phys. Rev. B 75
(2007), 155111). Es ändert sich nicht unter $\lambda\mapsto a\lambda+b$, $a>0$, hängt also nicht
vom lokalen mittleren Abstand ab; deshalb verwendet man es ohne Unfolding.

Ein Verhältnis braucht drei Niveaus, also sind die kleinsten Gaußschen Matrizen mit einem
Abstandsverhältnis die $3\times3$-Matrizen; für sie lässt sich das mittlere Verhältnis exakt
berechnen, so wie die Wignersche Vermutung für $2\times2$.

**Proposition 5.14 (mittlere Abstandsverhältnisse).** (a) *Unabhängige Niveaus.* Seien $s,s'$
unabhängig mit derselben Exponentialverteilung (wie zwei aufeinanderfolgende Abstände der Niveaus
$e_0+\dots+e_{k-1}$, gebildet aus unabhängigen exponentialverteilten $e_j$, Abschnitt A.7). Dann
gilt für $r=\min(s,s')/\max(s,s')$: $P(r\le t)=\frac{2t}{1+t}$ auf $[0,1]$, die Dichte ist
$\frac2{(1+t)^2}$, und
$$
\mathbb Er=2\ln2-1\approx0.3863 .
$$
(b) *$3\times3$-GOE und -GUE.* Sei $H$ eine $3\times3$-Matrix aus dem GOE ($\beta=1$) oder dem GUE
($\beta=2$) mit Eigenwerten $\lambda_1<\lambda_2<\lambda_3$, Abständen $s_1=\lambda_2-\lambda_1$,
$s_2=\lambda_3-\lambda_2$ und $r=\min(s_1,s_2)/\max(s_1,s_2)$. Dann hat $x=s_2/s_1$ eine Dichte
proportional zu
$$
f_\beta(x)=\frac{(x+x^2)^\beta}{(1+x+x^2)^{1+3\beta/2}}\qquad\text{on }(0,\infty),
$$
$r$ hat die Dichte $f_\beta/I_0$ auf $(0,1)$, wobei $I_0=\int_0^1f_\beta$, und
$$
\mathbb Er=4-2\sqrt3\approx0.5359\quad(\text{GOE}),\qquad
\mathbb Er=\frac{2\sqrt3}\pi-\frac12\approx0.6027\quad(\text{GUE}).
$$

*Beweis.* (a) Wegen der Skaleninvarianz sei die Rate $1$. Für $x=s/s'$ ist wegen der
Unabhängigkeit und nach dem Satz von Tonelli
$P(x\le t)=\int_0^\infty e^{-y}(1-e^{-ty})\,dy=1-\frac1{1+t}=\frac t{1+t}$.
Für $t<1$ sind die Ereignisse $\{x\le t\}$ und $\{1/x\le t\}$ disjunkt, $x$ und $1/x$ haben
dieselbe Verteilung ($s,s'$ sind austauschbar), und $\{r\le t\}$ ist ihre Vereinigung, also ist
$P(r\le t)=\frac{2t}{1+t}$. Dann ist
$\mathbb Er=\int_0^1\frac{2t}{(1+t)^2}dt=2\Big[\ln(1+t)+\frac1{1+t}\Big]_0^1=2\ln2-1$.

(b) *Die Verteilung von $x$.* Nach Satz C.7 hat der Vektor $\lambda=(\lambda_1,\lambda_2,\lambda_3)$
die Dichte $Z^{-1}\Delta(\lambda)^\beta e^{-|\lambda|^2/2}$ auf $C=\{\lambda_1<\lambda_2<\lambda_3\}$,
mit $\Delta(\lambda)=\prod_{i<j}(\lambda_j-\lambda_i)$. Die lineare Abbildung
$$
(m,s_1,s_2)\mapsto\Big(m-\tfrac{2s_1+s_2}3,\ m+\tfrac{s_1-s_2}3,\ m+\tfrac{s_1+2s_2}3\Big)
$$
hat Determinante $1$ und bildet $\mathbb R\times(0,\infty)^2$ auf $C$ ab; ihre Inverse ist
$m=\frac13\sum_i\lambda_i$, $s_1=\lambda_2-\lambda_1$, $s_2=\lambda_3-\lambda_2$. In diesen Variablen
ist $\Delta=s_1s_2(s_1+s_2)$ und $|\lambda|^2=3m^2+\frac23(s_1^2+s_1s_2+s_2^2)$ (ausmultiplizieren).
Nach Lemma A.3(c) und dem Satz von Fubini (Integration über $m$) hat $(s_1,s_2)$ die Dichte
$c\,\big(s_1s_2(s_1+s_2)\big)^\beta e^{-(s_1^2+s_1s_2+s_2^2)/3}$ auf $(0,\infty)^2$ mit einer
Konstanten $c>0$. Sei $q(x)=1+x+x^2$. Für jede Borelmenge $A\subseteq(0,\infty)$ ergibt die
Substitution $s_2=xs_1$ im inneren Integral zusammen mit dem Satz von Tonelli
$$
P(x\in A)=c\int_0^\infty\!\!\int_A s_1^{3\beta+1}(x+x^2)^\beta e^{-s_1^2q(x)/3}\,dx\,ds_1
=c\,\gamma_\beta\int_A(x+x^2)^\beta\Big(\frac3{q(x)}\Big)^{1+3\beta/2}dx ,
$$
wobei die Substitution $s_1=(3/q(x))^{1/2}u$ den Faktor liefert, mit
$\gamma_\beta=\int_0^\infty u^{3\beta+1}e^{-u^2}du\in(0,\infty)$. Also hat $x$ die Dichte
$f_\beta/\int_0^\infty f_\beta$.

*Reduktion auf $[0,1]$.* $f_\beta(1/x)=x^2f_\beta(x)$ (Zähler und Nenner mit Potenzen von $x$
multiplizieren). Mit $x=1/y$ ist
$\int_{1/A}f_\beta(x)\,dx=\int_Af_\beta(1/y)y^{-2}dy=\int_Af_\beta$ für jede Borelmenge
$A\subseteq(0,1)$, wobei $1/A=\{1/y:y\in A\}$; insbesondere ist $\int_0^\infty f_\beta=2I_0$ mit
$I_j=\int_0^1x^jf_\beta(x)\,dx$. Da $r=\min(x,1/x)$ und $P(x=1)=0$, ist
$\{r\in A\}=\{x\in A\}\cup\{x\in1/A\}$ (disjunkt), also
$$
P(r\in A)=\frac{2\int_Af_\beta}{2I_0}=\int_A\frac{f_\beta}{I_0},\qquad \mathbb Er=\frac{I_1}{I_0} .
$$
Wegen $x+x^2=q-1$ ist $f_\beta=g_\beta(q)$ mit $g_1(u)=u^{-3/2}-u^{-5/2}$ und $g_2(u)=(u-1)^2u^{-4}$.
Man schreibe $x=\frac12(q'(x)-1)$ mit $q'=2x+1$ und substituiere $u=q(x)$ (von $q(0)=1$ bis
$q(1)=3$):
$$
I_1=\frac12\int_1^3g_\beta(u)\,du-\frac12I_0,\qquad \mathbb Er=\frac{\int_1^3g_\beta}{2I_0}-\frac12 .
$$
Dabei ist $\int_1^3g_1=\big[-2u^{-1/2}+\frac23u^{-3/2}\big]_1^3=\frac43-\frac{16}{9\sqrt3}$ und
$\int_1^3g_2=\big[-u^{-1}+u^{-2}-\frac13u^{-3}\big]_1^3=\frac8{81}$.

*Die Integrale $I_0$.* Sei $L_a=\int_0^1q^{-a}dx$. Aus $(2x+1)^2=4q-3$ folgt
$\frac d{dx}\big[(2x+1)q^{-a}\big]=2q^{-a}-a(2x+1)^2q^{-a-1}=(2-4a)q^{-a}+3aq^{-a-1}$, und Integration
über $[0,1]$ ergibt
$$
3a\,L_{a+1}=3^{1-a}-1+(4a-2)L_a .
$$
GOE: $f_1=q^{-3/2}-q^{-5/2}$. Mit $a=\frac12$ ist $L_{3/2}=\frac23(\sqrt3-1)$; mit $a=\frac32$ ist
$L_{5/2}=\frac29\big(\frac1{\sqrt3}-1+4L_{3/2}\big)=\frac{18\sqrt3-22}{27}$. Also ist
$I_0=L_{3/2}-L_{5/2}=\frac4{27}$ und
$\mathbb Er=\frac{27}8\big(\frac43-\frac{16}{9\sqrt3}\big)-\frac12=\frac92-2\sqrt3-\frac12=4-2\sqrt3$.
GUE: $f_2=q^{-2}-2q^{-3}+q^{-4}$. Wegen $q=(x+\frac12)^2+\frac34$ ist
$L_1=\frac2{\sqrt3}\big[\arctan\frac{2x+1}{\sqrt3}\big]_0^1=\frac2{\sqrt3}\big(\frac\pi3-\frac\pi6\big)=\frac\pi{3\sqrt3}$, und die
Rekursion mit $a=1,2,3$ liefert $L_2=\frac23L_1=\frac{2\pi}{9\sqrt3}$, $L_3=L_2-\frac19$, $L_4=\frac19\big(10L_3-\frac89\big)$. Also ist
$I_0=L_2-2L_3+L_4=\frac{2\pi}{81\sqrt3}$ und $\mathbb Er=\frac{81\sqrt3}{4\pi}\cdot\frac8{81}-\frac12=\frac{2\sqrt3}\pi-\frac12$. $\square$

Die Werte in (b) wurden von Y. Y. Atas, E. Bogomolny, O. Giraud und G. Roux gefunden (Phys. Rev.
Lett. 110 (2013), 084101). Für größere $n$ wird das mittlere Verhältnis in diesem Kurs nicht
berechnet; der folgende Code *beobachtet* es numerisch für das GOE mit $n=64$.

Teil (a) gilt exakt für $n$ unabhängige gleichverteilte Niveaus auf einem Intervall: Ihre $n+1$
Abstände (einschließlich der beiden Randlücken) haben die Verteilung von
$(e_0,\dots,e_n)/\sum_je_j$ für unabhängige exponentialverteilte $e_j$, sodass jedes Verhältnis
$s_i/s_{i+1}$ aufeinanderfolgender Abstände wie $e_i/e_{i+1}$ verteilt ist und jedes
Abstandsverhältnis $r_i$ die Verteilung von $r$ aus (a) hat (Satz A.25 und Korollar A.26(b),(d)).

Die Funktion `SpinChain.level_spacing_ratio` (in `chains`) gibt das mittlere Verhältnis jeder
Zeile von Niveaus zurück. Wir überprüfen die Propositionen 5.13 und 5.14 und halten den GOE-Wert
für $n=64$ fest:

```python
ratio = SpinChain.level_spacing_ratio
n2 = 200_000
alg_r = TypeIAlgebra([2], [2], complex_valued=False, precision='double')     # real symmetric 2 x 2
alg_c = TypeIAlgebra([2], [2], precision='double')
G = torch.randn(n2, 1, 2, 2, dtype=torch.float64)
Gc = torch.randn(n2, 1, 2, 2, dtype=dt)
goe2 = alg_r.operator((G + G.mT) / 2, is_self_adjoint=True)
gue2 = alg_c.operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
cdf = {1: lambda t: 1 - math.exp(-math.pi * t ** 2 / 4),
       2: lambda t: math.erf(2 * t / math.sqrt(math.pi)) - 4 * t / math.pi * math.exp(-4 * t ** 2 / math.pi)}
for beta, H, mean_s in ((1, goe2, math.sqrt(math.pi)), (2, gue2, 4 / math.sqrt(math.pi))):
    lam = H.eigenvalues()[0]                                         # (n2, 2), ascending
    S = (lam[:, 1] - lam[:, 0]) / mean_s
    mc(S, 1.0, f"beta={beta}: E S")
    for t in (0.25, 1.0, 2.0):
        mc((S <= t).double(), cdf[beta](t), f"beta={beta}: P(S <= {t})")

levels = torch.rand(4000, 50, dtype=torch.float64)                   # independent levels
mc(ratio(levels), 2 * math.log(2) - 1, "<r>, independent uniform levels")      # Proposition 5.14(a)

n3 = 100_000                                                         # Proposition 5.14(b): 3 x 3
G = torch.randn(n3, 1, 3, 3, dtype=torch.float64)
Gc = torch.randn(n3, 1, 3, 3, dtype=dt)
goe3 = TypeIAlgebra([3], [3], complex_valued=False, precision='double').operator((G + G.mT) / 2, is_self_adjoint=True)
gue3 = TypeIAlgebra([3], [3], precision='double').operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
h = 1e-4
xs = (torch.arange(int(1 / h), dtype=torch.float64) + 0.5) * h      # midpoint rule on [0, 1]
for beta, H, exact in ((1, goe3, 4 - 2 * math.sqrt(3)), (2, gue3, 2 * math.sqrt(3) / math.pi - 0.5)):
    f = (xs + xs ** 2) ** beta / (1 + xs + xs ** 2) ** (1 + 1.5 * beta)
    assert abs((xs * f).sum() / f.sum() - exact) < 1e-7              # E r = I_1 / I_0
    r3 = ratio(H.eigenvalues()[0])                                   # one ratio per matrix
    mc(r3, exact, f"beta={beta}: <r>, 3 x 3")
    mc((r3 <= 0.5).double(), (f * (xs <= 0.5)).sum().item() / f.sum().item(), f"beta={beta}: P(r <= 1/2), 3 x 3")

n = 64                                                               # numerical observation, no theorem
alg64 = TypeIAlgebra([n], [n], complex_valued=False, precision='double')
G = torch.randn(800, 1, n, n, dtype=torch.float64)
lam = alg64.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()[0]
r_goe = ratio(lam[:, n // 4: 3 * n // 4])                            # 800 samples of the mean ratio
print(f"<r>, GOE(64), middle half of the spectrum: {r_goe.mean():.4f} +- {r_goe.std() / math.sqrt(800):.4f}")

def same_mean(a, b, name):
    """Two independent Monte Carlo samples of the same law: means agree within 5 combined standard errors."""
    ea, eb = a.std().item() / math.sqrt(a.numel()), b.std().item() / math.sqrt(b.numel())
    print(f"{name}: {a.mean():.4f} vs {b.mean():.4f} +- {math.hypot(ea, eb):.4f}")
    assert abs(a.mean().item() - b.mean().item()) <= 5 * math.hypot(ea, eb)
```

*Numerische Beobachtung.* Der Code findet $\langle r\rangle\approx0.53$ für die mittlere Hälfte des
Spektrums des GOE mit $n=64$, unterhalb des exakten $3\times3$-Werts $4-2\sqrt3\approx0.536$ und weit
oberhalb des Werts $2\ln2-1\approx0.386$ unabhängiger Niveaus. Diese Zahl ist ein Ergebnis der
Rechnung, kein Satz; unten dient sie als Referenz für GOE-Blöcke derselben Größe $64$, deren
Abstandsverhältnisse genau dieselbe Verteilung haben.

## 5.4 Warum Sektoren nicht vermischt werden dürfen

Sei $H=\bigoplus_cH_c\in\mathcal M$ ein zufälliges selbstadjungiertes Element mit **unabhängigen
Blöcken**, zum Beispiel GOE-Blöcken. Als Operator auf $\mathbb C^D$ hat es als Spektrum die
Vereinigung der Spektren der Blöcke. Innerhalb eines Blocks stoßen sich die Niveaus ab; wir zeigen,
dass Niveaus verschiedener Blöcke dies nicht tun.

Wir nehmen an, dass jeder Block eine beschränkte stetige **Niveaudichte** $\rho_c$ hat:
$\mathbb E\sum_i\varphi(\lambda_i(H_c))=\int\varphi\,\rho_c$ für beschränkte Borel-Funktionen
$\varphi$ (also $\int\rho_c=k_c$). Für GOE- und GUE-Blöcke gilt dies nach Korollar C.8, das
$\rho_c$ aus der gemeinsamen Dichte von Satz C.7 herleitet.

**Satz 5.15 (keine Abstoßung zwischen Sektoren).** Für $c\neq d$ und jede beschränkte
Borel-Funktion $\psi:\mathbb R^2\to\mathbb R$ gilt
$$
\mathbb E\sum_{i,j}\psi\big(\lambda_i(H_c),\lambda_j(H_d)\big)=\iint\psi(x,y)\,\rho_c(x)\rho_d(y)\,dx\,dy .
$$
Ist folglich $N_\varepsilon$ die Anzahl der Paare von Niveaus aus verschiedenen Blöcken mit
Abstand $<\varepsilon$, so gilt $\mathbb EN_\varepsilon/(2\varepsilon)\to\sum_{c<d}\int\rho_c\rho_d$
für $\varepsilon\to0$, und der kleinste Abstand $\delta$ des vereinigten Spektrums erfüllt
$$
P(\delta<\varepsilon)\ \ge\ \frac{\mathbb EN_\varepsilon}{\sum_{c<d}k_ck_d}
=\frac{2\varepsilon\sum_{c<d}\int\rho_c\rho_d}{\sum_{c<d}k_ck_d}+o(\varepsilon).
$$

*Beweis.* Die geordneten Eigenwerte sind stetige Funktionen der Matrix (Korollar A.11), also sind
$\Lambda_c=(\lambda_i(H_c))_i$ und $\Lambda_d$ unabhängige Zufallsvektoren. Nach dem Satz von
Fubini (ihre gemeinsame Verteilung ist das Produkt ihrer Verteilungen) ist
$\mathbb E\sum_{i,j}\psi(\Lambda_{c,i},\Lambda_{d,j})=\mathbb E\sum_j g(\Lambda_{d,j})$ mit
$g(y)=\mathbb E\sum_i\psi(\Lambda_{c,i},y)=\int\psi(x,y)\rho_c(x)\,dx$, was beschränkt und
Borel-messbar ist; man wende die Definition von $\rho_d$ auf $g$ an. Mit
$\psi(x,y)=1\{|x-y|<\varepsilon\}$ ist
$$
\frac{\mathbb EN_\varepsilon}{2\varepsilon}=\sum_{c<d}\int\rho_c(x)\Big(\frac1{2\varepsilon}\int_{x-\varepsilon}^{x+\varepsilon}\rho_d(y)\,dy\Big)dx ;
$$
die Klammer strebt wegen der Stetigkeit gegen $\rho_d(x)$ und ist durch $\sup\rho_d$ beschränkt,
und $\rho_c$ ist integrierbar, also liefert majorisierte Konvergenz den Grenzwert. Schließlich ist
$N_\varepsilon\le\sum_{c<d}k_ck_d$, also
$\mathbb EN_\varepsilon\le\big(\sum_{c<d}k_ck_d\big)P(N_\varepsilon\ge1)$, und aus
$N_\varepsilon\ge1$ folgt $\delta<\varepsilon$. $\square$

Für zwei GOE- oder GUE-Blöcke ist $P(\delta<\varepsilon)$ also mindestens von der Ordnung
$\varepsilon$: Ihre Niveaudichten sind überall positiv (in Korollar C.8 ist $q(x,\cdot)>0$ außerhalb
der Nullmenge, auf der zwei Koordinaten von $(x,y)$ übereinstimmen), also ist $\int\rho_c\rho_d>0$.
Innerhalb eines $2\times2$-GOE-Blocks ist sie gleich $1-e^{-\varepsilon^2/4}\le\varepsilon^2/4$
(Proposition 5.13, nicht normierter Abstand), und innerhalb eines $n\times n$-GOE-Blocks höchstens
eine Konstante mal $\varepsilon^2$ (Proposition C.9). Der Test auf Niveauabstoßung schlägt für das
vereinigte Spektrum fehl, obwohl jeder Block ein perfektes GOE ist.

Mit mehr Sektoren ist der Effekt stärker. Wir beweisen keinen Grenzwertsatz für viele Sektoren;
das folgende Experiment ist eine numerische Beobachtung. In der Bibliothek sind die Blöcke von
vornherein getrennt: `eigenvalues()` gibt ein Spektrum pro Sektor zurück, und sie zu vermischen
erfordert ein explizites `torch.cat`. Wir legen $C$ unabhängige GOE-Blöcke der Größe $64$ in eine
Algebra und vergleichen das mittlere Verhältnis pro Sektor mit dem des vereinigten Spektrums. Jeder
Sektor ist eine GOE-Matrix der Größe $64$, seine Verhältnisse haben also genau die Verteilung von
`r_goe` oben, und der Code testet das; das Verhalten des vereinigten Spektrums wird nur beobachtet:

```python
k, merged_r = 64, {}
for C in (1, 2, 4, 8, 16):
    alg_C = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double')
    G = torch.randn(800 // C, C, k, k, dtype=torch.float64)
    lam = alg_C.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()    # one spectrum per sector
    per_sector = torch.stack([ratio(l[:, k // 4: 3 * k // 4]) for l in lam])
    merged = torch.cat(lam, dim=-1).sort(dim=-1).values                       # the spectrum on C^D
    r_m = ratio(merged[:, k * C // 4: 3 * k * C // 4])
    merged_r[C] = r_m.mean().item()
    print(f"C = {C:2d}: <r> per sector {per_sector.mean():.4f}, merged {merged_r[C]:.4f}"
          f" +- {r_m.std() / math.sqrt(len(r_m)):.4f}")
    same_mean(per_sector.flatten(), r_goe, f"C = {C:2d}: sectors vs GOE(64)")   # every sector is GOE(64)
assert merged_r[1] > merged_r[2] > merged_r[4] > merged_r[16]       # numerical observation
assert abs(merged_r[16] - (2 * math.log(2) - 1)) < 0.01              # numerical observation
```

*Numerische Beobachtung.* Jeder Sektor hat $\langle r\rangle\approx0.53$, den GOE(64)-Wert aus
Abschnitt 5.3 (diese Übereinstimmung ist der Test einer Tatsache: Die Sektoren sind
GOE(64)-Matrizen). Das vereinigte Spektrum fällt bei zwei Sektoren auf $0.42$ und bei $8$ und $16$
Sektoren auf etwa $0.39$, weniger als $0.01$ vom Wert $2\ln2-1$ unabhängiger Niveaus entfernt
(Proposition 5.14(a)). Dass sich das vereinigte Spektrum mit wachsendem $C$ unabhängigen Niveaus
nähert, beobachtet der Code; in diesem Kurs wird es nicht bewiesen.

Dasselbe geschieht in einem physikalischen Modell. Die Heisenberg-Kette mit Zufallsfeldern
(Kapitel 9) erhält die Magnetisierung, also liegt ihr Hamiltonoperator in
$\bigoplus_{N=0}^LM_{\binom LN}$ (Proposition 9.1). Ihre Blöcke sind *nicht* unabhängig – sie
hängen alle von denselben Zufallsfeldern ab –, sodass Satz 5.15 in der angegebenen Form nicht
anwendbar ist; das folgende Experiment zeigt denselben Effekt, als numerische Beobachtung.

```python
chain = SpinChain(10, complex_valued=False)                         # sectors N = 0, ..., 10
lam = chain.random_field_heisenberg(W=0.5, batch_size=20).eigenvalues()
k5 = chain.sector_dim(5)                                            # 252
r_sector = ratio(lam[5][:, k5 // 4: 3 * k5 // 4]).mean().item()
spectrum = torch.cat(lam, dim=-1).sort(dim=-1).values               # all 1024 levels mixed
r_all = ratio(spectrum[:, 256:768]).mean().item()
print(f"Heisenberg chain, L = 10: <r> in sector N = 5: {r_sector:.4f}, all sectors mixed: {r_all:.4f}")
assert r_sector > 0.49 and r_all < 0.42                             # numerical observation
```

*Numerische Beobachtung.* Innerhalb des Sektors $N=5$ findet der Code
$\langle r\rangle\approx0.53$, nahe am obigen GOE(64)-Wert, und das Spektrum mit allen Sektoren
vermischt ergibt $0.39$, nahe an $2\ln2-1$.

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Niveaus mit
> verschiedenen erhaltenen Quantenzahlen wechselwirken nicht, also kreuzen sie sich ungehindert;
> nur Niveaus derselben Symmetrie stoßen sich ab. Ein Spektrum, das Symmetriesektoren vermischt,
> zeigt daher Poisson-artige Statistik, selbst wenn die Dynamik in jedem Sektor vollständig
> chaotisch ist, und eine nicht aufgelöste Symmetrie kann mit Integrabilität oder Lokalisierung
> verwechselt werden. Die Faustregel des Fachs – vor der Berechnung einer Niveaustatistik jede
> Symmetrie auflösen – besagt, dass man in der Algebra $\bigoplus_cM_{k_c}$ rechnen soll, nicht in
> $M_D$.

## 5.5 Womit man vergleichen sollte

Haar-verteilte unitäre Matrizen: `scipy.stats.unitary_group` und `ortho_group`, die
Haar-zufälligen unitären Matrizen von QuTiP sowie jede QR-Zerlegung mit `numpy`/`torch`, gefolgt
von der Phasenkorrektur. Weingarten-Kalkül: Die symbolischen Pakete IntU (Mathematica;
Z. Puchała, J. A. Miszczak), RTNI (Mathematica und Python; M. Fukuda, R. König, I. Nechita,
J. Phys. A 52 (2019), 425303) und das Python-Paket Haarpy berechnen Momente jeder Ordnung
symbolisch; Satz 5.10 ist ihr Fall zweiter Ordnung und Satz 5.12 (in Anhang B bewiesen) die
allgemeine Formel, die sie implementieren. Niveaustatistiken werden meist mit wenigen Zeilen
`numpy` berechnet, und Codes für exakte Diagonalisierung wie QuSpin liefern Spektren Sektor für
Sektor. Was die Bibliothek hinzufügt, ist die unitäre Gruppe einer Algebra als ein Objekt:
unabhängige Haar-Blöcke aller Größen, in gestapelten Aufrufen gezogen, der Twirl, identifiziert mit
der bedingten Erwartung `channels.center_expectation`, und Spektren, die pro Sektor zurückgegeben
werden, sodass das Vermischen von Sektoren nie der Standardfall ist.

## Aufgaben

1. Sei $U$ Haar-verteilt auf $U(\mathcal M)$ und $\operatorname{Tr}U=\sum_c\operatorname{Tr}U_c$.
   Man zeige $\mathbb E|\operatorname{Tr}U|^2=C$ und
   $\mathbb E|\operatorname{Tr}U|^4=\sum_cm(k_c)+2C(C-1)$ mit $m(1)=1$, $m(k)=2$ für $k\ge2$, und
   vergleiche mit $1$ und $2$ für eine Haar-verteilte unitäre Matrix aus $U(D)$, $D\ge2$. Man
   überprüfe dies für $\mathcal M=M_1\oplus M_2\oplus M_3$.
2. Man zeige mit Satz 5.12, dass $\mathbb E|U_{11}|^{2n}=\frac{n!}{k(k+1)\cdots(k+n-1)}$ für
   $n\le k$ (Korollar B.14 erweitert dies auf $n>k$). Man berechne die Weingarten-Matrix der
   Ordnung $3$ für $k=3$ numerisch und überprüfe $\mathbb E|U_{11}|^6=\frac1{10}$.
3. Das **zirkuläre orthogonale Ensemble** (`measure='coe'`) ist $S=W^TW$ mit $W$ Haar-verteilt auf
   $U(k)$. Man zeige, dass $S$ unitär und symmetrisch ist und dass seine Verteilung unter
   $S\mapsto V^TSV$ für jedes $V\in U(k)$ invariant ist. Seien $\lambda_1,\dots,\lambda_k$ die
   Eigenwerte von $S$ (mit Vielfachheit) und $N_B(S)=\#\{j:\lambda_j\in B\}$ für Borelmengen $B$ des
   Einheitskreises $\mathbb T$. Man zeige, dass $N_B(S)$ eine Zufallsvariable ist, dass
   $\nu(B)=\frac1k\mathbb EN_B(S)$ die normierte Bogenlänge ist (die mittlere Dichte der Eigenphasen
   ist gleichförmig) und dass $\mathbb E\operatorname{Tr}S^m=0$ für jede ganze Zahl $m\neq0$ gilt. Man
   berechne $\langle r\rangle$ für die Eigenphasen des COE mit $k=32$ und vergleiche es mit dem
   GOE-Wert aus Abschnitt 5.3 (eine numerische Beobachtung).

<details>
<summary>Lösungen</summary>

1. Man schreibe $t_c=\operatorname{Tr}U_c$; die $t_c$ sind unabhängig (Korollar 5.4). Die
   Multiplikation von $U_c$ mit $e^{i\theta}$ ändert seine Verteilung nicht und multipliziert
   $t_c$ mit $e^{i\theta}$; daher verschwindet ein Moment $\mathbb E[t_at_b\bar t_c\bar t_d]$, außer
   wenn jeder Sektor in $(a,b)$ ebenso oft vorkommt wie in $(c,d)$, d. h. $\{a,b\}=\{c,d\}$ als
   Multimengen, und ebenso ist $\mathbb E[t_a\bar t_b]=0$ für $a\neq b$. Nach Proposition 5.7 ist
   $\mathbb E[(U_c)_{ii}\overline{(U_c)_{jj}}]=\langle e_i,\mathbb E[U_ce_ie_j^*U_c^*]e_j\rangle=\delta_{ij}/k_c$,
   also $\mathbb E|t_c|^2=1$ und $\mathbb E|\operatorname{Tr}U|^2=C$. In $\mathbb E|\sum_at_a|^4$
   liefern die Terme mit $a=b$ den Beitrag $\sum_a\mathbb E|t_a|^4$, die mit $a\neq b$ den Beitrag
   $2\sum_{a\neq b}\mathbb E|t_a|^2\,\mathbb E|t_b|^2=2C(C-1)$ (zwei Wahlen von $(c,d)$). Für
   $k=1$ ist $|t|=1$. Für $k\ge2$ ist
   $\mathbb E|t|^4=\sum_{i,j,i',j'}\mathbb E[U_{ii}U_{jj}\overline{U_{i'i'}U_{j'j'}}]$; in Satz 5.10
   erzwingen die beiden Terme mit $\sigma=\pi$ auf jeweils eine Weise $\{i',j'\}=\{i,j\}$ und
   liefern nach Summation je $k^2\mathrm{Wg}(e)$, die beiden mit $\sigma\neq\pi$ erzwingen
   $i=j=i'=j'$ und liefern je $k\,\mathrm{Wg}(\tau)$: insgesamt
   $2\big(\frac{k^2}{k^2-1}-\frac1{k^2-1}\big)=2$. Der Fall $C=1$, $k=D$ liefert $1$ und $2$ für
   $U(D)$. Für $M_1\oplus M_2\oplus M_3$: $\mathbb E|\operatorname{Tr}U|^4=1+2+2+12=17$, gegenüber
   $2$ für eine Haar-verteilte unitäre Matrix aus $U(6)$.

```python
tr = U.trace                                                        # Tr U = sum_c Tr U_c, shape (B,)
mc(tr.abs() ** 2, 3.0, "E|Tr U|^2")
mc(tr.abs() ** 4, 17.0, "E|Tr U|^4")
```

2. Sind alle Indizes gleich $1$, so ist jedes Delta $1$, also
   $\mathbb E|U_{11}|^{2n}=\sum_{\sigma,\pi}\mathrm{Wg}(\sigma\pi^{-1})=\mathbf 1^TW\mathbf 1$, wobei
   $W=\Gamma^{-1}=(\mathrm{Wg}(\sigma^{-1}\pi))$ und $\Gamma=(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)})$.
   Dabei wurde benutzt, dass $\mathrm{Wg}$ nur vom Zyklentyp abhängt und dass
   $\sigma\pi^{-1}=\sigma(\pi^{-1}\sigma)\sigma^{-1}$ zu $\pi^{-1}\sigma=(\sigma^{-1}\pi)^{-1}$
   konjugiert ist, das den Zyklentyp von $\sigma^{-1}\pi$ hat. Jede Zeilensumme von $\Gamma$ ist
   $c_n(k)=\sum_{\rho\in S_n}k^{\#\mathrm{cycles}(\rho)}$, also $\Gamma\mathbf 1=c_n(k)\mathbf 1$,
   $W\mathbf 1=\mathbf 1/c_n(k)$ und $\mathbf 1^TW\mathbf 1=n!/c_n(k)$. Schließlich ist
   $c_n(x)=x(x+1)\cdots(x+n-1)$ nach Induktion: Eine Permutation von $\{1,\dots,n\}$ entsteht aus
   einer von $\{1,\dots,n-1\}$ entweder, indem man $n$ als Fixpunkt hinzufügt (ein Zyklus mehr,
   Faktor $x$), oder, indem man $n$ in seinem Zyklus hinter einem der $n-1$ Elemente einfügt
   ($n-1$ Möglichkeiten, gleiche Zyklenzahl); dies ist eine Bijektion, also
   $c_n(x)=(x+n-1)c_{n-1}(x)$.

```python
import itertools

def n_cycles(p):                                                    # number of cycles of a permutation
    seen, count = set(), 0
    for i in range(len(p)):
        count += i not in seen
        while i not in seen:
            seen.add(i)
            i = p[i]
    return count

perms = list(itertools.permutations(range(3)))                      # S_3; s.index(j) = s^{-1}(j)
Gram = torch.tensor([[3.0 ** n_cycles([s.index(p[i]) for i in range(3)]) for p in perms] for s in perms],
                    dtype=torch.float64)                            # k^{#cycles(s^{-1} p)}, k = 3
Wg3 = torch.linalg.inv(Gram)                                        # the Weingarten matrix, n = k = 3
assert abs(Wg3.sum().item() - 6 / (3 * 4 * 5)) < 1e-12
mc(u3[:, 0, 0].abs() ** 6, 0.1, "E|U_11|^6")
```

3. $S^*S=W^*\bar WW^TW=1$, da $W^T$ unitär ist, und $S^T=W^TW=S$. Für $V\in U(k)$ hat $WV$ die
   Verteilung von $W$ (Rechtsinvarianz, Lemma 5.1), und $(WV)^T(WV)=V^TSV$. Mit
   $V=e^{i\alpha/2}1$: $e^{i\alpha}S$ hat für jedes reelle $\alpha$ die Verteilung von $S$.

   *Messbarkeit.* Nach Satz 1.3 ist $S=Y\operatorname{diag}(\lambda_1,\dots,\lambda_k)Y^*$ mit
   unitärem $Y$. Ist $v$ ein Einheitsvektor mit $Sv=\mu v$, so liefert die Darstellung
   $v=\sum_jc_jYe_j$ die Gleichungen $c_j(\lambda_j-\mu)=0$, also liegt $v$ in
   $E_K=\operatorname{span}\{Ye_j:\lambda_j\in K\}$, sobald $\mu\in K$; und $\dim E_K=N_K(S)$. Daher
   gilt für $K\subseteq\mathbb T$ und $m\ge1$: $N_K(S)\ge m$ genau dann, wenn $S$ $m$ orthonormale
   Eigenvektoren mit Eigenwerten in $K$ hat. Ist $K$ abgeschlossen, so ist die Menge $F_{K,m}$ der
   unitären Matrizen mit dieser Eigenschaft abgeschlossen in $U(k)$: Gilt $S_l\to S$ und hat $S_l$
   orthonormale Eigenvektoren $v_i^{(l)}$ mit Eigenwerten $\mu_i^{(l)}\in K$ ($i\le m$), so hat eine
   Teilfolge $v_i^{(l)}\to v_i$ und $\mu_i^{(l)}\to\mu_i\in K$ für alle $i\le m$
   (Bolzano-Weierstraß in der kompakten Menge $(\text{Einheitssphäre}\times K)^m$), und die
   Grenzwerte sind orthonormal mit $Sv_i=\mu_iv_i$.
   Seien $0\le\theta_1(S)\le\dots\le\theta_k(S)<2\pi$ die Eigenphasen, sodass die Eigenwerte die
   $e^{i\theta_j(S)}$ sind (Korollar 1.4: $|\lambda_j|=1$). Für $0\le t<2\pi$ ist der Bogen
   $K_t=\{e^{i\varphi}:0\le\varphi\le t\}$ abgeschlossen und
   $\{\theta_j\le t\}=\{N_{K_t}\ge j\}=F_{K_t,j}$, also ist jedes $\theta_j$ eine Borel-Funktion auf
   $U(k)$. Folglich ist $N_B(S)=\sum_j1_B(e^{i\theta_j(S)})$ eine Zufallsvariable, und
   $\nu(B)=\frac1k\sum_jP(e^{i\theta_j(S)}\in B)$ ist ein Borelsches Wahrscheinlichkeitsmaß auf
   $\mathbb T$, das Mittel der Verteilungen der $k$ Eigenwerte.

   *Gleichförmigkeit.* Die Eigenwerte von $e^{i\alpha}S$ sind die $e^{i\alpha}\lambda_j$, also ist
   $N_B(e^{i\alpha}S)=N_{e^{-i\alpha}B}(S)$, und da $e^{i\alpha}S$ die Verteilung von $S$ hat, ist
   $\nu(e^{-i\alpha}B)=\frac1k\mathbb EN_B(e^{i\alpha}S)=\nu(B)$. Also ist $\nu$ drehinvariant und
   nach Proposition A.22 (mit $N=1$) die normierte Bogenlänge: $\mathbb EN_B(S)=k\,|B|/2\pi$ für
   jeden Bogen der Länge $|B|$. Ebenso ist $\operatorname{Tr}S^m=\sum_j\lambda_j^m$ (auch für
   $m<0$, da $S^{-1}=S^*$) eine beschränkte stetige Funktion von $S$, und
   $\mathbb E\operatorname{Tr}S^m=\mathbb E\operatorname{Tr}(e^{i\alpha}S)^m=e^{im\alpha}\,\mathbb E\operatorname{Tr}S^m$
   für alle $\alpha$, sodass es für $m\neq0$ verschwindet: Das sind $k$-mal die
   Fourier-Koeffizienten von $\nu$. Insbesondere ist die mittlere Dichte der Eigenphasen konstant,
   und die Abstände entlang des Kreises brauchen kein Unfolding.

   *Numerische Beobachtung.* Der Code findet $\langle r\rangle\approx0.53$ für das COE mit $k=32$,
   weniger als $0.015$ vom GOE(64)-Wert aus Abschnitt 5.3 entfernt. Diese Übereinstimmung ist ein
   Ergebnis der Rechnung; kein Satz dieses Kurses setzt die Abstandsstatistik des COE mit der des
   GOE in Beziehung (zu den zirkulären Ensembles siehe F. J. Dyson, J. Math. Phys. 3 (1962),
   140–156).

```python
alg_s = TypeIAlgebra([32], [32], precision='double')
torch.manual_seed(2)
S = alg_s.random_unitary_operator(400, measure='coe').matrix[:, 0]
torch.manual_seed(2)
W = alg_s.random_unitary(32, batch_size=400)                         # same random numbers
assert close(S, W.mT @ W)                                            # the library's COE is W^T W
assert close(S, S.mT) and close(S @ S.mH, torch.eye(32, dtype=dt))
for m in (1, 2, 3):                                                  # E Tr S^m = 0
    trm = torch.linalg.matrix_power(S, m).diagonal(dim1=-2, dim2=-1).sum(-1)
    mc(trm.real, 0.0, f"E Re Tr S^{m}")
    mc(trm.imag, 0.0, f"E Im Tr S^{m}")
phases = torch.linalg.eigvals(S).angle().sort(dim=-1).values         # in (-pi, pi]
mc(((phases > 0) & (phases <= 1)).sum(-1).double(), 32 / (2 * math.pi), "E #eigenphases in (0, 1]")
s = torch.diff(torch.cat([phases, phases[:, :1] + 2 * math.pi], dim=-1), dim=-1)   # around the circle
r = torch.minimum(s, s.roll(-1, -1)) / torch.maximum(s, s.roll(-1, -1))
r_coe = r.mean(-1)
print(f"<r>, COE(32) eigenphases: {r_coe.mean():.4f}; GOE(64), Section 5.3: {r_goe.mean():.4f}")
assert abs(r_coe.mean() - r_goe.mean()) < 0.015                      # numerical observation
```

</details>

**Weiter:** [Kapitel 6](06_quantum_channels.md) behandelt Quantenkanäle: Kraus- und
Stinespring-Darstellungen, duale Abbildungen bezüglich einer Spur und Abbildungen zwischen
Sektoren. Die bedingte Erwartung $E$, die hier als Haar-Mittel auftritt, ist das erste Beispiel.
