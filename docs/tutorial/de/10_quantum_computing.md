[English](../10_quantum_computing.md) · [Español](../es/10_quantum_computing.md) · [Français](../fr/10_quantum_computing.md) · **Deutsch** · [中文](../zh/10_quantum_computing.md) · [日本語](../ja/10_quantum_computing.md) · [Русский](../ru/10_quantum_computing.md)

# 10. Quantencomputing

Rauschen auf einem Quantencomputer ist ein Kanal, und ein Kanal erzeugt eine Algebra: die
$*$-Algebra seiner Kraus-Operatoren. Ihre Struktur (Kapitel 2) zeigt, welche Information das
Rauschen unberührt übersteht; eine Bedingung an Produkte von Kraus-Operatoren, in die eine weitere
Algebra eingeht, zeigt, welche Information eine Wiederherstellung zurückgewinnen kann. Der letzte
Teil behandelt Algorithmen, die eine Ladung erhalten: Ihre Zustände leben in der Algebra der
Ladungssektoren, und die Amplitudendämpfung verschiebt sie zwischen den Sektoren.

**Sie lernen**, warum ein unitaler Kanal und sein dualer Kanal dieselben Fixpunkte haben (den
Kommutanten des Rauschens), den Satz von Knill–Laflamme mit vollständigem Beweis, wie man nur die
Bezeichnung eines Sektors schützt und wie der XY-Mixer des Quantum Alternating Operator Ansatz, die
Amplitudendämpfung und die Postselektion in der Algebra $\bigoplus_wM_{\binom nw}$ der
Hamming-Gewichts-Sektoren aussehen.

Im Folgenden ist ein **Kanal** auf $M_D(\mathbb C)$ eine Abbildung $\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*$
mit **Kraus-Operatoren**, die $\sum_iE_i^*E_i=1$ erfüllen (spurerhaltend); er heißt **unital**, wenn
außerdem $\sum_iE_iE_i^*=1$ gilt. Seine **duale Abbildung** $\Phi^*(T)=\sum_iE_i^*TE_i$ erfüllt
$\operatorname{Tr}(\Phi(\rho)T)=\operatorname{Tr}(\rho\,\Phi^*(T))$ (Abschnitt 6.2). (Nach Satz 6.2 hat jede
vollständig positive Abbildung diese Form; hier sind Kanäle durch Kraus-Operatoren *gegeben*.) Die
**Wechselwirkungsalgebra** $\mathcal A$ von $\Phi$ ist die von $1$ und den $E_i$ erzeugte
$*$-Algebra, der Spann aller Produkte der $E_i$ und $E_i^*$; eine Matrix kommutiert genau dann mit
$\mathcal A$, wenn sie mit allen $E_i$ und $E_i^*$ kommutiert. Ein Zustand $\rho$ hat seinen
**Träger** in einem Unterraum $K$ mit Projektion $P_K$, wenn $\rho=P_K\rho P_K$. Auf $n$ Qubits ist
$D=2^n$; $X,Y,Z$ sind die Pauli-Matrizen, $Z=\operatorname{diag}(1,-1)$, $X_j$ ist $X$ auf dem Qubit $j$ (nummeriert mit $1,\dots,n$;
Abschnitt 10.4 folgt der Bibliothek, Plätze $0,\dots,n-1$), $|b\rangle$, $b\in\{0,1\}^n$, ist die
Rechenbasis, $|b|$ das **Hamming-Gewicht** und $P_w$ die Projektion auf $H_w=\operatorname{span}\{|b\rangle:|b|=w\}$.

## 10.1 Rauschalgebren, rauschfreie Teilsysteme und Fixpunkte

Nach Satz 2.6 gibt es ein unitäres $V$ mit $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$ und
$V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$: Jeder Kraus-Operator wirkt nur auf den ersten
Tensorfaktor jedes Summanden.

**Proposition 10.1 (rauschfreie Teilsysteme; dekohärenzfreie Unterräume).** (a) Man schreibe
$V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}$, und sei $\rho=V(\sigma\otimes\tau)V^*$, im Summanden $c$
gelegen, für Zustände $\sigma$ auf $\mathbb C^{k_c}$ und $\tau$ auf $\mathbb C^{m_c}$. Dann ist
$\Phi(\rho)=V(\Phi_c(\sigma)\otimes\tau)V^*$ mit dem Kanal $\Phi_c(\sigma)=\sum_iE_{i,c}\sigma E_{i,c}^*$ auf
$M_{k_c}$: Der Zustand $\tau$ bleibt unberührt. Außerdem ist $\Phi^*(T)=T$ für jedes $T\in\mathcal A'$.
(b) Ist $K$ ein Unterraum mit $E_iP_K=e_iP_K$, $e_i\in\mathbb C$, für alle $i$, so ist $\Phi(\rho)=\rho$
für jeden Zustand $\rho$ mit Träger in $K$.

*Beweis.* (a) $\Phi(\rho)=V\sum_i(E_{i,c}\sigma E_{i,c}^*\otimes\tau)V^*$, und $\Phi_c$ ist spurerhaltend,
weil der Summand $c$ von $V^*(\sum_iE_i^*E_i)V=1$ gleich $(\sum_iE_{i,c}^*E_{i,c})\otimes1_{m_c}$ ist. Für
$T\in\mathcal A'$ ist $\Phi^*(T)=T\sum_iE_i^*E_i=T$. (b) $P_KE_i^*=(E_iP_K)^*=\bar e_iP_K$, also
$E_i\rho E_i^*=E_iP_K\rho P_KE_i^*=|e_i|^2\rho$, und $\sum_i|e_i|^2P_K=P_K\sum_iE_i^*E_iP_K=P_K$. $\square$

Ein Summand mit $k_c=1$ ist der Fall (b): ein **dekohärenzfreier Unterraum** (Zanardi–Rasetti 1997,
Lidar–Chuang–Whaley 1998). Für $k_c>1$ ist der Faktor $\mathbb C^{m_c}$ ein **rauschfreies
Teilsystem** (Knill–Laflamme–Viola 2000): Das Rauschen wirkt, aber nur auf den anderen Faktor.
Proposition 10.1 gilt für jeden Kanal; für unitale Kanäle sind die erhaltenen Observablen *genau*
die von $\mathcal A'$:

**Satz 10.2 (Fixpunkte unitaler Kanäle).** Sei $\Phi$ unital. Dann gilt
$$
\operatorname{Fix}(\Phi)=\{T\in M_D(\mathbb C):\Phi(T)=T\}=\mathcal A'=\operatorname{Fix}(\Phi^*).
$$
Insbesondere ist $\operatorname{Fix}(\Phi)$ eine von-Neumann-Algebra, $*$-isomorph zu $\bigoplus_cM_{m_c}$.

*Beweis.* Die erste Gleichheit ist Satz 6.10 für die Algebra $M_D(\mathbb C)$ (ein Block): $T$ ist
genau dann ein Fixpunkt, wenn es mit allen $E_i$ und $E_i^*$ kommutiert, d. h. genau dann, wenn
$T\in\mathcal A'$. Die duale Abbildung $\Phi^*$ ist wieder ein unitaler Kanal: Ihre Kraus-Operatoren
$E_i^*$ erfüllen $\sum_iE_iE_i^*=1$ (spurerhaltend) und $\sum_iE_i^*E_i=1$ (unital), und sie erzeugen
dieselbe $*$-Algebra $\mathcal A$. Also liefert Satz 6.10, angewandt auf $\Phi^*$,
$\operatorname{Fix}(\Phi^*)=\mathcal A'$. Die letzte Behauptung ist Satz 2.6. $\square$

Für jeden Kanal gilt $\mathcal A'\subseteq\operatorname{Fix}(\Phi^*)$ (Proposition 10.1(a)), aber ohne
Unitalität gilt der Satz nicht: Für die Amplitudendämpfung ist $\operatorname{Fix}(\Phi)=\mathbb C|0\rangle\langle0|$,
während $\mathcal A'=\operatorname{Fix}(\Phi^*)=\mathbb C1$ (Abschnitt 6.4 und Aufgabe 1). (Historisch: Kribs 2003,
Arias–Gheondea–Gudder 2002.)

**Korollar 10.3 (kollektive Dephasierung).** Sei $U=e^{-i\theta Z_{\mathrm{tot}}}$,
$Z_{\mathrm{tot}}=\sum_jZ_j$, mit $e^{2i\theta d}\neq1$ für $d=1,\dots,n$, und
$\Phi(\rho)=\tfrac12(U\rho U^*+U^*\rho U)$. Dann ist $\mathcal A=\operatorname{span}\{P_0,\dots,P_n\}$, jedes
$H_w$ ist ein dekohärenzfreier Unterraum, und $\operatorname{Fix}(\Phi)=\bigoplus_wB(H_w)\cong\bigoplus_wM_{\binom nw}$
hat die Dimension $\sum_w\binom nw^2$.

*Beweis.* $Z_{\mathrm{tot}}|b\rangle=(n-2|b|)|b\rangle$, also $U=\sum_we^{-i\theta(n-2w)}P_w$ mit
verschiedenen Eigenwerten (ihre Quotienten sind $e^{2i\theta(w-w')}$, $0<|w-w'|\le n$). Die
Lagrange-Polynome aus dem Beweis von Satz 1.5(a) liefern $P_w=p_w(U)\in\mathcal A$; umgekehrt liegen
$U,U^*$ in der Algebra $\operatorname{span}\{P_w\}$. Also ist $\mathcal A=\operatorname{span}\{P_w\}$
($k_c=1$, $m_c=\binom nw$), und Proposition 10.1(b) ist auf jedes $H_w$ anwendbar. Kommutiert $T$ mit
allen $P_w$, so ist $P_vTP_{v'}=TP_vP_{v'}=0$ für $v\neq v'$, d. h. $T$ ist blockdiagonal, und
umgekehrt; Satz 10.2 liefert den Rest. $\square$

Kommutanten sind Nullräume von $T\mapsto(AT-TA)_A$ (Kapitel 2), Fixpunkte der von
$\sum_iE_i\otimes\bar E_i-1$ (der Matrix von $\Phi-\mathrm{id}$ auf zeilenweise vektorisierten
Matrizen). Für $n=4$ haben die drei Räume des Korollars die Dimensionen $5$, $70=1+16+36+16+1$ und
$70$:

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

$\bigoplus_wM_{\binom nw}$ ist die Algebra einer `SpinChain` (Kapitel 9), deren Sektoren die Anzahlen
der Spins nach oben sind, d. h. die Hamming-Gewichte. Dort sind die Kraus-Operatoren zentral,
$U=\sum_we^{-i\theta(n-2w)}P_w$, und der Kanal ist die Identität: Die Bibliothek speichert den
rauschfreien Teil.

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

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Ein fluktuierendes Feld, das auf allen Qubits gleich ist, verschiebt
> die Phase von $|b\rangle$ proportional zu $n-2|b|$. Zustände gleichen Hamming-Gewichts erhalten eine
> gemeinsame Phase, sodass $|01\rangle,|10\rangle$ ein gegen dieses Rauschen immunes Qubit kodieren
> (gefangene Ionen: Kielpinski et al., *Science* 291, 2001), während Kohärenzen zwischen
> verschiedenen Gewichten zerfallen: Das Rauschen erzwingt eine Superauswahlregel (Abschnitt 3.3).

**Kollektive Dekohärenz.** Wir nennen einen Kanal **kollektiv**, wenn seine Kraus-Operatoren Vielfache
von Drehungen $e^{-i\theta\,\mathbf n\cdot\mathbf S}$ sind ($\theta\in\mathbb R$, $\mathbf n\in\mathbb R^3$ ein Einheitsvektor),
$S_\alpha=\tfrac12\sum_j\sigma_{\alpha,j}$, $\sigma\in\{X,Y,Z\}$; das modelliert ein Feld, das auch in seiner Richtung
fluktuiert (das Modell, kein Satz). Da $\mathbf n\cdot\mathbf S$ selbstadjungiert ist, ist $e^{-i\theta\,\mathbf n\cdot\mathbf S}$
ein Polynom darin (Beweis von Satz 1.5(a)), also liegen diese Kraus-Operatoren in der von $1,S_x,S_y,S_z$ erzeugten
$*$-Algebra $\mathcal B$.

**Proposition 10.4 (Singulett-Unterraum).** Sei $K=\ker S_x\cap\ker S_y\cap\ker S_z$. Jeder Kanal,
dessen Kraus-Operatoren in $\mathcal B$ liegen, lässt jeden Zustand mit Träger in $K$ fest.

*Beweis.* $\mathcal B$ wird von den Wörtern in $S_x,S_y,S_z$ (und $1$) aufgespannt; Wörter der Länge
$\ge1$ annullieren $K$, also wirkt jedes Element von $\mathcal B$ auf $K$ als Skalar. Man wende
Proposition 10.1(b) an. $\square$

Für $n=4$ liefert Proposition 4.13 ([Kapitel 4](04_symmetries.md), Fall $L=4$)
$\mathcal B\cong M_5\oplus M_3\oplus M_1$ mit den Vielfachheiten $1,3,2$ (Gesamtspin $2,1,0$), also ist $\dim\mathcal B=35$,
und $\mathcal B'\cong M_1\oplus M_3\oplus M_2$ hat die Dimension $14$. Außerdem ist $\dim K=m_0=2$: Die $S_a$ wirken
summandenweise, also ist $K$ die Summe seiner Durchschnitte mit den Summanden; auf dem Summanden zum Spin $S$ wirken sie
als $S^{(S)}_a\otimes1$ (Korollar 4.12), und
$\sum_bu_b\otimes f_b$ ($f_b$ eine Basis von $\mathbb C^{m_S}$) wird genau dann von allen annulliert, wenn jedes $u_b$ von
allen $S^{(S)}_a$ annulliert wird. Diese $u$ bilden einen invarianten Unterraum des irreduziblen Spins $S$ (Satz 4.11), der
für $S>0$ nicht der ganze Raum ist (da $S^{(S)}_z\neq0$), also $0$, während er für $S=0$ gleich $\mathbb C^1$ ist. Der Code
findet dieselben Zahlen. Proposition 10.1(a) und ihr Beweis gelten wörtlich für jede $*$-Algebra der Form
$V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*$, die die Kraus-Operatoren enthält, hier $\mathcal B$: Vier Qubits tragen ein
dekohärenzfreies Qubit ($k=1$, $m=2$) und ein rauschfreies Qutrit ($k=3$, $m=3$) für jeden Kanal, dessen Kraus-Operatoren
in $\mathcal B$ liegen.

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

## 10.2 Fehlerkorrektur: der Satz von Knill–Laflamme

Ein **Code** ist ein von null verschiedener Unterraum $\mathcal C\subseteq\mathbb C^D$ mit Projektion $P$. Er **korrigiert**
$\Phi$, wenn es einen Kanal $\mathcal R$ (die **Wiederherstellung**) mit $\mathcal R(\Phi(\rho))=\rho$
für jeden Zustand $\rho$ mit Träger in $\mathcal C$ gibt.

**Satz 10.5 (Knill–Laflamme).** $\mathcal C$ korrigiert $\Phi$ genau dann, wenn es Zahlen
$\lambda_{ij}$ gibt mit
$$
PE_i^*E_jP=\lambda_{ij}P\qquad\text{for all }i,j.
$$

*Beweis.* *Notwendigkeit.* $\mathcal R$ habe die Kraus-Operatoren $R_m$; dann hat
$\mathcal R\circ\Phi$ die Kraus-Operatoren $R_mE_i$. Für einen Einheitsvektor $\psi\in\mathcal C$ ist
$\sum_{m,i}R_mE_i|\psi\rangle\langle\psi|E_i^*R_m^*=|\psi\rangle\langle\psi|$, und für $\phi\perp\psi$ liefert
der Erwartungswert in $\phi$ die Gleichung $\sum_{m,i}|\langle\phi,R_mE_i\psi\rangle|^2=0$, also
$R_mE_i\psi\in\mathbb C\psi$. Somit ist $R_mE_iP=c_{mi}P$: Eine lineare Abbildung $T$ mit
$Tx\in\mathbb Cx$ für alle $x$ ist ein Vielfaches von $1$ (ist $Tx=ax$, $Ty=by$, $T(x+y)=c(x+y)$ mit
linear unabhängigen $x,y$, so ist $a=c=b$). Wegen $\sum_mR_m^*R_m=1$ ist
$PE_i^*E_jP=\sum_m(R_mE_iP)^*(R_mE_jP)=\sum_m\bar c_{mi}c_{mj}P$.

*Hinlänglichkeit.* $\Lambda=(\lambda_{ij})$ ist hermitesch, da
$\lambda_{ji}P=(PE_i^*E_jP)^*=\bar\lambda_{ij}P$, und positiv, da
$\sum_{ij}\bar c_ic_j\lambda_{ij}P=PF^*FP\ge0$ für $F=\sum_jc_jE_j$. Man schreibe $\Lambda=UDU^*$ mit
unitärem $U$, $D=\operatorname{diag}(d_k)$, $d_k\ge0$ (Satz 1.3), und $F_k=\sum_iu_{ik}E_i$. Dann ist
$\sum_kF_k\rho F_k^*=\sum_{ij}(UU^*)_{ij}E_i\rho E_j^*=\Phi(\rho)$ und $PF_k^*F_lP=(U^*\Lambda U)_{kl}P=d_k\delta_{kl}P$.
Ist $d_k=0$, so ist $\|F_kPx\|^2=0$ für alle $x$, also $F_kP=0$. Ist $d_k>0$, so sei
$V_k=F_kP/\sqrt{d_k}$; dann ist $V_k^*V_l=\delta_{kl}P$, also $V_kP=V_k$, die $Q_k=V_kV_k^*$ sind
paarweise orthogonale Projektionen, und $Q=\sum_kQ_k$ ist eine Projektion. Die Wiederherstellung
$$
\mathcal R(\sigma)=\sum_{k:d_k>0}V_k^*\sigma V_k+(1-Q)\sigma(1-Q)
$$
ist ein Kanal, da $\sum_kV_kV_k^*+(1-Q)^2=1$. Für $\rho=P\rho P$ ist $\Phi(\rho)=\sum_{d_k>0}d_kV_k\rho V_k^*$;
wegen $V_l^*V_k=\delta_{kl}P$ und $(1-Q)V_k=(1-Q)Q_kV_k=0$ ist $\mathcal R(\Phi(\rho))=\sum_kd_k\rho$, und
$\sum_kd_k=\operatorname{Tr}\Lambda=1$, weil $\sum_i\lambda_{ii}P=P\sum_iE_i^*E_iP=P$. $\square$

(Knill–Laflamme 1997; Nielsen–Chuang, Theorem 10.1.) $\mathcal R$ misst das **Syndrom** – welcher der
orthogonalen Unterräume $Q_k\mathbb C^D$ den Zustand enthält – und dreht den Zustand in den Code
zurück.

**Beispiel: der Bitflip-Code.** $\mathcal C=\operatorname{span}\{|000\rangle,|111\rangle\}$ und „höchstens ein
Bitflip“, $E_0=\sqrt{p_0}\,1$, $E_j=\sqrt{p_j}\,X_j$. Für $i\neq j$ kippt $X_iX_j$ (mit $X_0=1$) ein oder
zwei Bits und bildet $\mathcal C$ in $\mathcal C^\perp$ ab; also ist $\Lambda=\operatorname{diag}(p_0,\dots,p_3)$,
$V_k=X_kP$, $Q=1$ und $\mathcal R(\sigma)=\sum_{k=0}^3PX_k\sigma X_kP$. Kippt jedes Qubit unabhängig mit
Wahrscheinlichkeit $p$, so ist die Bedingung verletzt (unter den Kraus-Operatoren sind Vielfache von $X_1X_2$ und $X_3$, und
$PX_1X_2X_3P=|000\rangle\langle111|+|111\rangle\langle000|$ ist kein Vielfaches von $P$), und $\mathcal R$ bildet $|0_L\rangle=|000\rangle$ auf
$|1_L\rangle$ ab, sobald $\ge2$ Bitflips auftreten, und lässt es sonst an seinem Platz: Die Fidelität
von $|0_L\rangle$ nach der Wiederherstellung ist die Wahrscheinlichkeit für höchstens einen Bitflip,
$(1-p)^3+3p(1-p)^2=1-3p^2+2p^3=0.972$ für $p=0.1$. In der Bibliothek ist $M_8$ eine Algebra mit einem
Block, und `R @ Phi` ist die Komposition $\mathcal R\circ\Phi$:

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

## 10.3 Fehlerkorrektur mit Operatoralgebren

Das Heisenberg-Bild erlaubt es, weniger als den ganzen Code zu schützen. Sei $\mathfrak A$ eine
von-Neumann-Algebra auf $\mathcal C=P\mathbb C^D$ (eine $*$-Algebra von Operatoren auf $\mathcal C$,
die $P$ enthält), die zu schützenden Observablen. $\mathfrak A$ heißt **korrigierbar** für $\Phi$,
wenn es einen Kanal $\mathcal R$ gibt mit
$$
P\,(\mathcal R\circ\Phi)^*(A)\,P=A\qquad\text{for all }A\in\mathfrak A,
$$
d. h. $\operatorname{Tr}(\mathcal R(\Phi(\rho))A)=\operatorname{Tr}(\rho A)$ für alle Zustände $\rho$ mit Träger in
$\mathcal C$ und alle $A\in\mathfrak A$ (die linke Seite ist $\operatorname{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)$, und ein
Operator $Y\in B(\mathcal C)$ mit $\operatorname{Tr}(\rho Y)=0$ für alle diese $\rho$ verschwindet:
$\langle x,Yx\rangle=\langle Px,YPx\rangle=0$ für alle $x$, und Polarisation wie im Beweis von Korollar 1.4). Für
$\mathfrak A=B(\mathcal C)$ ist das die Korrektur aus Abschnitt 10.2: Erfüllt $\sigma=\mathcal R(\Phi(\rho))$ die Gleichung
$P\sigma P=\rho$ (man nehme alle $A\in B(\mathcal C)$), so ist, da $\mathcal R\circ\Phi$ spurerhaltend ist,
$\operatorname{Tr}(1-P)\sigma(1-P)=\operatorname{Tr}\sigma-\operatorname{Tr}P\sigma P=1-1=0$, also $\sigma^{1/2}(1-P)=0$ und
$\sigma=\rho$; umgekehrt liefert $\mathcal R(\Phi(\rho))=\rho$ die Gleichheit der Spuren. Das andere Extrem ist eine
kommutative Algebra: Nur eine klassische Bezeichnung wird geschützt.

**Proposition 10.6 (Korrektur einer Sektorbezeichnung).** Sei $P=\sum_{c=1}^CP_c$ mit paarweise
orthogonalen Projektionen $P_c\neq0$, und sei $\mathfrak A=\operatorname{span}\{P_1,\dots,P_C\}$. Äquivalent
sind:
1. Es gibt $N_1,\dots,N_C\ge0$, $\sum_cN_c=1$, mit $\operatorname{Tr}(N_c\Phi(\rho))=1$ für jedes $c$ und
   jeden Zustand $\rho$ mit Träger in $P_c\mathbb C^D$;
2. $P_cE_i^*E_jP_d=0$ für alle $i,j$ und $c\neq d$, d. h. $PE_i^*E_jP\in\mathfrak A'$, dem Kommutanten
   von $\mathfrak A$ in $B(\mathcal C)$;
3. $\mathfrak A$ ist korrigierbar für $\Phi$.

*Beweis.* $T\in B(\mathcal C)$ kommutiert genau dann mit allen $P_c$, wenn $P_cTP_d=0$ für $c\neq d$
(wie in Korollar 10.3).

(2)$\Rightarrow$(1). Die Unterräume $S_c=\operatorname{span}\{E_iP_cx\}$ sind paarweise orthogonal, da
$\langle E_iP_cx,E_jP_dy\rangle=\langle x,P_cE_i^*E_jP_dy\rangle=0$. Sei $N_c$ die Projektion auf $S_c$ für
$c<C$ und $N_C=1-\sum_{c<C}N_c$, ebenfalls eine Projektion. Für $\rho$ mit Träger in $P_c\mathbb C^D$
hat $\Phi(\rho)=\sum_iE_iP_c\rho P_cE_i^*$ sein Bild in $S_c$, und $N_c$ ist die Identität auf $S_c$
(für $c=C$, weil $S_C\perp S_{c'}$, $c'<C$), also $\operatorname{Tr}(N_c\Phi(\rho))=\operatorname{Tr}\Phi(\rho)=1$.
(1)$\Rightarrow$(2). Für einen Einheitsvektor $x\in P_c\mathbb C^D$ und $d\neq c$ haben die Zahlen
$\operatorname{Tr}(N_{d'}\Phi(|x\rangle\langle x|))\ge0$ die Summe $1$, und die für $d'=c$ ist $1$; also
$0=\operatorname{Tr}(N_d\Phi(|x\rangle\langle x|))=\sum_i\|N_d^{1/2}E_ix\|^2$, folglich $N_dE_iP_c=0$ und,
summiert über $d$, $N_cE_iP_c=E_iP_c$. Daher ist $P_cE_i^*E_jP_d=P_cE_i^*N_cE_jP_d=0$ für $c\neq d$.
(1),(2)$\Rightarrow$(3). Sei $\omega_c=P_c/\operatorname{Tr}P_c$ und $\mathcal R(\sigma)=\sum_c\operatorname{Tr}(N_c\sigma)\,\omega_c$
(„messen, dann präparieren“), ein Kanal mit den Kraus-Operatoren $\sqrt\mu\,|u\rangle\langle f|N_c^{1/2}$
für die Eigenpaare $(\mu,u)$ von $\omega_c$ und die Vektoren $f$ einer Orthonormalbasis. Da
$\operatorname{Tr}(\mathcal R(\sigma)P_d)=\operatorname{Tr}(N_d\sigma)$, ist $\mathcal R^*(P_d)=N_d$. Mit
$N_dE_iP_c=\delta_{cd}E_iP_c$ (aus dem vorigen Schritt) und dann $\sum_iE_i^*E_i=1$ folgt
$$
P(\mathcal R\circ\Phi)^*(P_d)P=P\,\Phi^*(N_d)\,P=\sum_i\sum_{c,c'}P_{c'}E_i^*N_dE_iP_c=\sum_i\sum_{c'}P_{c'}E_i^*E_iP_d=P_d\Big(\sum_iE_i^*E_i\Big)P_d=P_d .
$$
(3)$\Rightarrow$(1). Sei $N_c=\mathcal R^*(P_c)$ für $c<C$ und $N_C=\mathcal R^*(1-\sum_{c<C}P_c)$; sie sind
positiv, und ihre Summe ist $\mathcal R^*(1)=1$. Für $\rho$ mit Träger in $P_c\mathbb C^D$ und $d<C$ ist
$\operatorname{Tr}(N_d\Phi(\rho))=\operatorname{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(P_d)P)=\operatorname{Tr}(\rho P_d)=\delta_{cd}$,
also auch $\operatorname{Tr}(N_C\Phi(\rho))=\delta_{cC}$. $\square$

**Satz 10.7 (Bény–Kempf–Kribs).** Eine von-Neumann-Algebra
$\mathfrak A$ auf $\mathcal C$ ist genau dann korrigierbar für $\Phi$, wenn $PE_i^*E_jP\in\mathfrak A'$
für alle $i,j$, wobei $\mathfrak A'$ der Kommutant von $\mathfrak A$ in $B(\mathcal C)$ ist.

*Beweis.* Dies ist Satz D.9, bewiesen in [Anhang D](D_gksl.md), dessen Beweis auch einen Wiederherstellungskanal
konstruiert: Er misst den Sektor von $\mathfrak A\cong\bigoplus_cM_{k_c}$ und ein Syndrom und macht den Fehler auf jedem
Faktor $\mathbb C^{k_c}$ wie im Beweis von Satz 10.5 rückgängig. $\square$

Für $\mathfrak A=B(\mathcal C)$ ist $\mathfrak A'=\mathbb CP$, und dies ist Satz 10.5; Proposition 10.6
ist der Fall des kommutativen $\mathfrak A=\operatorname{span}\{P_c\}$ (Korollar D.10(a),(b)). Für
$\mathfrak A=M_k\otimes1_m$ auf $\mathcal C=\mathbb C^k\otimes\mathbb C^m$ erhält man die Operator-Codes
(Subsystem-Codes); die rauschfreien Teilsysteme aus Abschnitt 10.1 sind solche Codes, korrigiert durch
$\mathcal R=\mathrm{id}$ (Korollar D.10(c), mit vertauschten Tensorfaktoren). Nach Satz 2.6 ist jedes $\mathfrak A$
$*$-isomorph zu einer Algebra $\bigoplus_cM_{k_c}$; ist $C\ge2$ und $k_c\ge2$ für ein $c$, so schützt ein korrigierbares
$\mathfrak A$ die Sektorbezeichnung $c$ zusammen mit der Quanteninformation jedes $M_{k_c}$ – ein **hybrider** Code
(Bemerkung nach Korollar D.10, mit einem Beispiel, in dem die Knill–Laflamme-Bedingung verletzt ist, Satz 10.7 aber
greift). Historisch: C. Bény, A. Kempf, D. W. Kribs, *Phys. Rev. Lett.* 98, 100502 (2007) und *Phys. Rev. A* 76,
042303 (2007); Operator-Codes: Kribs–Laflamme–Poulin, *Phys. Rev. Lett.* 94, 180501 (2005).

**Beispiel.** Der Wiederholungscode unter unabhängigen Phasenflips (Wahrscheinlichkeit $q$ pro
Qubit), gefolgt von höchstens einem Bitflip. Die Kraus-Operatoren sind Vielfache von $X_kZ^s$,
$s\in\{0,1\}^3$, und jedes $E_i^*E_j$ kippt höchstens zwei Bits, also gilt Bedingung 2 für
$P_0=|000\rangle\langle000|$, $P_1=|111\rangle\langle111|$: Das klassische Bit überlebt und wird durch
Mehrheitsentscheid ausgelesen, $N_0=\sum_{|b|\le1}|b\rangle\langle b|$, $N_1=1-N_0$. Das Qubit überlebt
nicht: $PZ_1P$ ist kein Vielfaches von $P$, und nach der Wiederherstellung $\mathcal R$ aus Abschnitt
10.2 erleidet der Zustand $|{+_L}\rangle=(|000\rangle+|111\rangle)/\sqrt2$ den logischen Phasenflip
$Z_1Z_2Z_3$ genau dann, wenn eine ungerade Anzahl von Phasenflips aufgetreten ist; seine Fidelität
ist $(1+(1-2q)^3)/2=0.8645$ für $q=0.05$. (Der Kanal ist $\Phi_{\mathrm{one}}\circ\mathcal Z$, mit dem Kanal
$\Phi_{\mathrm{one}}$ „höchstens ein Bitflip“ aus dem Bitflip-Beispiel und den unabhängigen Phasenflips $\mathcal Z$;
auf $\mathcal C$ wirkt jedes $Z^s$ als $Z_1^{|s|}$, also hat $\mathcal Z(\rho)=\sum_sq^{|s|}(1-q)^{3-|s|}Z_1^{|s|}\rho Z_1^{|s|}$
wieder Träger in $\mathcal C$, und $\mathcal R\circ\Phi_{\mathrm{one}}$ ist nach Satz 10.5 die Identität auf solchen
Zuständen. $Z_1$ bildet $|{+_L}\rangle$ auf $|{-_L}\rangle\perp|{+_L}\rangle$ ab, und $|s|$ ist mit Wahrscheinlichkeit
$\frac12\big((1-q+q)^3+(1-q-q)^3\big)$ gerade.)

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

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Eine kommutative
> Algebra zu schützen heißt, klassische Information zu schützen (ein Bit, die Bezeichnung eines
> Superauswahlsektors), was viel billiger ist, als ein Qubit zu schützen: Rauschen, das mit den
> Sektorprojektionen kommutiert, ist harmlos (dann ist $P_cE_i^*E_jP_d=E_i^*E_jP_cP_d=0$ für $c\neq d$,
> Bedingung 2 von Proposition 10.6).

## 10.4 Mixer mit Nebenbedingungen

Der Quantum Alternating Operator Ansatz (Hadfield et al., *Algorithms* 12, 34, 2019) präpariert
$|\beta,\gamma\rangle=e^{-i\beta_pH_M}e^{-i\gamma_pH_C}\cdots e^{-i\beta_1H_M}e^{-i\gamma_1H_C}|b_0\rangle$ mit einem
diagonalen Kosten-Hamiltonoperator $H_C$ und einem **Mixer** $H_M$. Für eine Nebenbedingung $|b|=w$
(Graphbisektion, Portfolios aus $w$ Anlagen) nimmt man $|b_0|=w$ und den **XY-Mixer** auf einem
Ring, $H_M=\tfrac14\sum_{(i,j)}(X_iX_j+Y_iY_j)$.

**Proposition 10.8.** $H_M$, $H_C$, $e^{-i\beta H_M}$ und $e^{-i\gamma H_C}$ liegen in
$\mathcal M=\{P_0,\dots,P_n\}'=\bigoplus_wB(H_w)\cong\bigoplus_wM_{\binom nw}$, also ist $|\beta,\gamma\rangle\in H_w$:
Der Ansatz verlässt den zulässigen Unterraum nie.

*Beweis.* Auf den Qubits $i,j$ ist $X\otimes X+Y\otimes Y=2(|01\rangle\langle10|+|10\rangle\langle01|)$ (man
verwende $Y|0\rangle=i|1\rangle$, $Y|1\rangle=-i|0\rangle$): $X_iX_j+Y_iY_j$ vertauscht die Bits $i,j$, wenn sie
verschieden sind, und annulliert $|b\rangle$ sonst, bildet also jedes $H_w$ in sich ab. Ein diagonales
$H_C$ kommutiert mit den diagonalen $P_w$, und ebenso die Exponentiale (Satz 1.5(d)). $\square$

In der Bibliothek ist $\mathcal M$ die Algebra `SpinChain(n).algebra` (Bit $1$ = Spin nach oben,
also $S^z_i=-\tfrac12Z_i$), und $H_M$ ist ihr XXZ-Hamiltonoperator
$J\sum_{(i,j)}\big(\tfrac12(S^+_iS^-_j+S^-_iS^+_j)+\Delta S^z_iS^z_j\big)$ mit $J=1$, $\Delta=0$, da
$S^+_iS^-_j+S^-_iS^+_j=\tfrac12(X_iX_j+Y_iY_j)$; wir prüfen das an der dichten $64\times64$-Matrix. Für
$H_C$ nehmen wir die Graphbisektion der Dreiecke $\{0,1,2\}$, $\{3,4,5\}$, verbunden durch die Kanten
$(2,3)$, $(5,0)$: $H_C=\sum_{(i,j)\in E}\tfrac12(1-Z_iZ_j)$, $Z_iZ_j=4S^z_iS^z_j$, zählt die geschnittenen
Kanten. Die beste Bisektion schneidet $2$ Kanten (die beiden Dreiecke; jede andere Menge von drei Ecken
teilt beide Dreiecke, und ein geteiltes Dreieck hat zwei geschnittene Kanten); wir starten bei $|101010\rangle$, das $6$ schneidet.
Schließlich wird der ganze Ansatz mit dichten $64\times64$-Matrizen unabhängig von der Bibliothek
nachgerechnet.

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

Die **Amplitudendämpfung** des Qubits $i$ ($|1\rangle\to|0\rangle$ mit Wahrscheinlichkeit $g$) hat die
Kraus-Operatoren $K_0^{(i)}=|0\rangle\langle0|_i+\sqrt{1-g}\,|1\rangle\langle1|_i$ und
$K_1^{(i)}=\sqrt g\,|0\rangle\langle1|_i$, der das Gewicht senkt: eine Abbildung zwischen Sektoren,
`chain.site_amplitude_damping(i, g)` (ein `InterSectorChannel`).

**Proposition 10.9 (Leckage und Postselektion).** Sei $\rho_0$ mit Träger in $H_w$, seien $U_l$
Unitäre in $\mathcal M$, $\mathcal D$ die Dämpfung aller Qubits und
$\rho_l=\mathcal D(U_l\rho_{l-1}U_l^*)$. Dann ist jedes $\rho_l$ blockdiagonal mit Gewichten $\le w$, und
$$
P_w\rho_pP_w=(1-g)^{wp}\;U\rho_0U^*,\qquad U=U_p\cdots U_1 :
$$
Das Gewicht $w$ wird mit Wahrscheinlichkeit $(1-g)^{wp}$ gefunden, und der postselektierte Zustand ist
der rauschfreie.

*Beweis.* Die Kraus-Operatoren von $\mathcal D$ sind $K_s=\prod_iK^{(i)}_{s_i}$, $s\in\{0,1\}^n$; $K_s$
bildet $H_v$ in $H_{v-|s|}$ ab, d. h. $P_vK_s=K_sP_{v+|s|}$ ($P_v=0$ für $v>n$). Ist $\sigma$
blockdiagonal mit Gewichten $\le w$, so verschwindet
$P_vK_s\sigma K_s^*P_{v'}=K_sP_{v+|s|}\sigma P_{v'+|s|}K_s^*$ für $v\neq v'$ und für $v>w-|s|$: $\mathcal D(\sigma)$
ist blockdiagonal mit Gewichten $\le w$, und sein Block zum Gewicht $w$ stammt allein von $s=0$. Da
$K_0|b\rangle=(1-g)^{|b|/2}|b\rangle$, ist $P_w\mathcal D(\sigma)P_w=(1-g)^wP_w\sigma P_w$. Die $U_l$ kommutieren
mit den $P_v$; Induktion über $l$. $\square$

Ein `InterSectorChannel` behält nur die Diagonalblöcke seiner Ausgabe; hier gibt es keine anderen
(der Beweis von Proposition 10.9 gilt wörtlich für die Dämpfung eines einzelnen Qubits, mit $s\in\{0,1\}$),
also geht nichts verloren. Mit $g=0.05$ und $p=2$ ist die Akzeptanz $0.95^6=0.7351$:

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

Numerische Beobachtung (die ausgegebenen Werte): Der unzulässige Anteil senkt den erwarteten Schnitt
von $4.3789$ auf $4.3123$, eine scheinbare Verbesserung, die einen Optimierer der Winkel in die Irre
führen würde; nach der Postselektion ist der Wert exakt (Proposition 10.9).

**Kontinuierliche Zeit.** Dasselbe gilt für die Lindblad-Gleichung (Satz 8.6, [Kapitel 8](08_dynamics.md))
$$
\dot\rho=-i[H,\rho]+\gamma_1\sum_i\Big(S^-_i\rho S^+_i-\tfrac12\{S^+_iS^-_i,\rho\}\Big)+\gamma_2\sum_k\Big(L_k\rho L_k^*-\tfrac12\{L_k^*L_k,\rho\}\Big)
$$
mit $H\in\mathcal M$ und Sprüngen $L_k\in\mathcal M$ innerhalb der Sektoren, z. B. reiner Dephasierung
$L_k=S^z_k$. Hier ist $\sum_iS^+_iS^-_i=\sum_i|1\rangle\langle1|_i=N$ der **Teilchenzahloperator**,
$NP_v=vP_v$: Er ist zentral.

**Proposition 10.10.** Hat $\rho(0)$ seinen Träger in $H_w$, so ist $P_w\rho(t)P_w=e^{-\gamma_1wt}\sigma(t)$,
wobei $\sigma$ dieselbe Gleichung mit $\gamma_1=0$ und $\sigma(0)=\rho(0)$ löst. Folglich ist die Akzeptanz
$e^{-\gamma_1wt}$, und der postselektierte Zustand $\sigma(t)$ hängt nicht von $\gamma_1$ ab.

*Beweis.* Sei $B_{vv'}=P_v\rho P_{v'}$. Da $H,L_k,N$ mit den $P_v$ kommutieren und
$P_vS^-_i=S^-_iP_{v+1}$, ist
$\dot B_{vv'}=-i[H,B_{vv'}]+\gamma_1\sum_iS^-_iB_{v+1,v'+1}S^+_i-\tfrac{\gamma_1}2(v+v')B_{vv'}+\gamma_2\mathcal L_2(B_{vv'})$,
wobei $\mathcal L_2$ der Dephasierungsanteil ist. Die Blöcke mit $\max(v,v')>w$ erfüllen ein
abgeschlossenes lineares System mit Anfangswerten null, verschwinden also (Eindeutigkeit in Lemma
8.1). Dann ist $\dot B_{ww}=-i[H,B_{ww}]-\gamma_1wB_{ww}+\gamma_2\mathcal L_2(B_{ww})$, sodass
$e^{\gamma_1wt}B_{ww}$ die Gleichung von $\sigma$ mit demselben Anfangswert löst und daher gleich
$\sigma$ ist (Lemma 8.1). Diese Gleichung hat wieder GKSL-Form, also ist $\operatorname{Tr}\sigma(t)=1$
(Satz 8.6). $\square$

Man vergleiche Proposition 8.11, die ein nicht zentrales $\sum_jL_j^*L_j$, aber keine Sprünge
innerhalb der Sektoren zulässt: Hier ist der Verlustterm $\gamma_1N$ zentral, sodass der sprungfreie
Faktor der Skalar $e^{-\gamma_1wt}$ ist und Dephasierung innerhalb des Sektors erlaubt ist.

Für $\gamma_2=0$ ist $\sigma(t)=e^{-iHt}\rho(0)e^{iHt}$: Die Postselektion ist perfekt. Der
Begleitartikel (`examples/paper/ex9_qaoa_xy_mixer.py`) lässt $|111000\rangle$ unter $H_M$ bis $T=2$
evolvieren; wir reproduzieren seine Zeile $\gamma_1=\gamma_2=0.1$: $P(w=3)=e^{-0.6}=0.5488$ (Proposition
10.10) sowie, als numerische Beobachtung des folgenden Codes, die Fidelität $0.4805$ und die postselektierte
Fidelität $0.8755$; nach Proposition 10.10 ist letztere dieselbe wie für $\gamma_1=0$.

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

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Energierelaxation ($T_1$), einer der wichtigsten Fehler
> supraleitender Qubits, ist ein Sprung in einen niedrigeren Sektor. Die Messung des
> Hamming-Gewichts – eine zentrale Projektion, harmlos innerhalb eines Sektors (Abschnitt 3.3) –
> entdeckt ihn, und die Postselektion entfernt ihn um den Preis der Akzeptanzrate. Dephasierung
> kommutiert mit den Sektorprojektionen und wird nicht entdeckt.

## 10.5 Womit man vergleichen sollte

Qiskit (`qiskit.quantum_info`, Rauschmodelle von Aer, `QAOAAnsatz` mit einem eigenen Mixer), Cirq,
PennyLane (`qml.qaoa.xy_mixer`, `default.mixed`) und QuTiP/`qutip-qip` simulieren all dies auf $n$
vollen Qubits und sind die richtigen Werkzeuge für Schaltkreise; stim sampelt Clifford-Schaltkreise
mit Pauli-artigem Rauschen auf Tausenden von Qubits (Amplitudendämpfung ist kein Pauli-Kanal: Pauli-Kanäle sind unital, die Amplitudendämpfung nicht). Uns
ist in diesen Paketen keine Funktion bekannt, die die Wechselwirkungsalgebra eines Kanals, ihren
Kommutanten oder seine ganze Fixpunktalgebra berechnet (QuTiPs `steadystate` liefert einen
stationären Zustand eines Lindblad-Operators); zu numerischen Blockzerlegungen von
Matrix-$*$-Algebren siehe Murota et al., *Japan J. Indust. Appl. Math.* 27 (2010). Die Bibliothek
ergänzt die Speicherung in $\bigoplus_wM_{\binom nw}$, Abbildungen zwischen Sektoren und Postselektion
als zentrale Projektion, gestapelt. Ein Sektor braucht $\binom nw^2$ Einträge statt $4^n$
($1.7\cdot10^8$ gegenüber $4.3\cdot10^9$ für $n=16$, $w=8$), aber die ganze Algebra wird auf
$(n+1)\binom n{\lfloor n/2\rfloor}^2$ Einträge ($2.8\cdot10^9$) aufgefüllt, eine bescheidene Ersparnis. Sie
ist kein Schaltkreissimulator und hat keinen Stabilisatorformalismus.

## Aufgaben

1. (Satz 10.2 braucht Unitalität; vergleiche Abschnitt 6.4.) Man zeige für die Amplitudendämpfung
   eines Qubits, $E_0=|0\rangle\langle0|+\sqrt{1-g}\,|1\rangle\langle1|$, $E_1=\sqrt g\,|0\rangle\langle1|$, $0<g<1$,
   dass $\mathcal A'=\mathbb C1$, $\operatorname{Fix}(\Phi)=\mathbb C|0\rangle\langle0|$ und
   $\operatorname{Fix}(\Phi^*)=\mathbb C1$.
2. Man beweise $\sum_w\binom nw^2=\binom{2n}n$. Was sind für den Kanal aus Korollar 10.3 mit
   $\theta=\pi/2$ (der dessen Voraussetzung verletzt) $\mathcal A$ und
   $\operatorname{Fix}(\Phi)$? Man prüfe $\dim\operatorname{Fix}(\Phi)$ für $n=4$.
3. (Linearität der Fehlerkorrektur.) Es gelte $PE_a^*E_bP=\lambda_{ab}P$ für endlich viele Operatoren
   $E_a$ (nicht notwendig Kraus-Operatoren eines Kanals), und $\mathcal R$ werde aus $\Lambda$ wie im
   Beweis von Satz 10.5 gebildet. Man zeige, dass $\mathcal R$ jeden Kanal korrigiert, dessen
   Kraus-Operatoren in $\operatorname{span}\{E_a\}$ liegen. Man folgere, dass der Bitflip-Code den
   kohärenten Fehler $e^{-i\phi X_1}$ korrigiert.

<details>
<summary>Lösungen</summary>

1. $E_1$, $E_1^*$ und ihre Produkte ergeben alle Matrixeinheiten, also ist $\mathcal A=M_2$ und
   $\mathcal A'=\mathbb C1$. Für $T=(x_{kl})$ ist $\Phi(T)=\begin{pmatrix}x_{00}+gx_{11}&\sqrt{1-g}\,x_{01}\\\sqrt{1-g}\,x_{10}&(1-g)x_{11}\end{pmatrix}=T$
   genau dann, wenn $x_{11}=x_{01}=x_{10}=0$, und
   $\Phi^*(T)=\begin{pmatrix}x_{00}&\sqrt{1-g}\,x_{01}\\\sqrt{1-g}\,x_{10}&(1-g)x_{11}+gx_{00}\end{pmatrix}$
   ist genau dann gleich $T$, wenn $T=x_{00}1$. ($\Phi$ ist nicht unital, und $\Phi^*$ ist unital,
   aber nicht spurerhaltend.)

```python
E_ad = [torch.tensor([[1, 0], [0, math.sqrt(0.7)]], dtype=dt), torch.tensor([[0, math.sqrt(0.3)], [0, 0]], dtype=dt)]
fix_ad = fixed_points(E_ad)
assert len(fix_ad) == 1 and len(commutant_basis(E_ad + [E.mH for E in E_ad])) == 1
assert torch.allclose(fix_ad[0] / fix_ad[0][0, 0], torch.tensor([[1, 0], [0, 0]], dtype=dt))
assert len(fixed_points([E.mH for E in E_ad])) == 1                     # Fix(Phi^*) = C 1
```

2. Man zähle die $n$-elementigen Teilmengen von $\{1,\dots,2n\}$ nach der Anzahl $w$ ihrer Elemente
   in der ersten Hälfte: $\binom{2n}n=\sum_w\binom nw\binom n{n-w}=\sum_w\binom nw^2$. Für
   $\theta=\pi/2$ sieht $U=(-i)^n\sum_w(-1)^wP_w$ nur die Parität von $w$, also ist
   $\mathcal A=\operatorname{span}\{P_{\mathrm{even}},P_{\mathrm{odd}}\}$ mit $P_{\mathrm{even}}=\sum_{w\text{ even}}P_w$,
   und nach Satz 10.2 besteht $\operatorname{Fix}(\Phi)$ aus den bezüglich der Parität blockdiagonalen
   Operatoren: Dimension $2\cdot(2^{n-1})^2=128$ für $n=4$.

```python
U2 = torch.linalg.matrix_exp(-1j * (math.pi / 2) * sum(on(Z, j, 4) for j in range(4)))
assert len(fixed_points([U2 / math.sqrt(2), U2.mH / math.sqrt(2)])) == 128
assert sum(math.comb(4, v) ** 2 for v in range(5)) == math.comb(8, 4) == 70
```

3. Mit $F_k=\sum_au_{ak}E_a$ ist $\operatorname{span}\{E_a\}=\operatorname{span}\{F_k\}$, also ist jeder
   Kraus-Operator des Kanals $G_m=\sum_kh_{mk}F_k$, und $G_mP=\sum_{d_k>0}h_{mk}\sqrt{d_k}V_k$. Wegen
   $V_l^*V_k=\delta_{lk}P$ und $(1-Q)V_k=0$ ist $\mathcal R(G_m\rho G_m^*)=\sum_l|h_{ml}|^2d_l\,\rho$ für
   $\rho=P\rho P$, also $\mathcal R(\Phi(\rho))=c\rho$; Spurbildung ergibt $c=1$. Für
   $E_a\in\{1,X_1,X_2,X_3\}$ ist $\Lambda=1$, $V_k=X_kP$, und $\mathcal R$ ist die Wiederherstellung
   aus Abschnitt 10.2; $e^{-i\phi X_1}=\cos\phi\,1-i\sin\phi\,X_1$ liegt im Spann.

```python
for phi in (0.2, 1.0, 2.5):
    G1 = torch.linalg.matrix_exp(-1j * phi * on(X, 0, 3))
    out = (R @ channels.kraus_channel([op8(G1)]))(rho_L)
    assert (out.matrix - rho_L.matrix).abs().max() < 1e-12       # corrected exactly
```

</details>

**Weiter:** [Kapitel 11](11_optimisation.md): Optimierung über die unitären Gruppen einer Algebra und
Spurungleichungen.
