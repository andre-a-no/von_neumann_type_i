[English](../D_gksl.md) · [Español](../es/D_gksl.md) · [Français](../fr/D_gksl.md) · **Deutsch** · [中文](../zh/D_gksl.md) · [日本語](../ja/D_gksl.md) · [Русский](../ru/D_gksl.md)

# Anhang D. Generatoren quantendynamischer Halbgruppen

Kapitel 8 hat bewiesen, dass ein Generator in Gorini–Kossakowski–Sudarshan–Lindblad-Form (GKSL-Form) eine
quantendynamische Halbgruppe erzeugt (Satz 8.6), und hat die Umkehrung (Satz 8.7) sowie ihre Version für
$\mathcal M=\bigoplus_cM_{k_c}$ verwendet. Dieser Anhang beweist beides in endlicher Dimension, zusammen mit der
Freiheit in der GKSL-Form und einem Kriterium für die Eindeutigkeit des stationären Zustands. Das Werkzeug ist die
Choi-Matrix aus Kapitel 6: Ein Generator ist im Allgemeinen nicht vollständig positiv, aber seine Choi-Matrix ist positiv auf
dem orthogonalen Komplement des maximal verschränkten Vektors. Dieselbe Art von Argument – Positivität
einer Familie von Operatoren, dann eine Zerlegung – beweist im letzten Abschnitt den in Kapitel 10 formulierten Satz
über Fehlerkorrektur mit Operatoralgebren (Satz 10.7), mit einem expliziten Wiederherstellungskanal.

**Sie lernen**, warum der Generator einer vollständig positiven Halbgruppe bedingt vollständig
positiv ist, wie man $H$ und die Sprungoperatoren an der Choi-Matrix eines Generators abliest, in welchem Sinn
sie eindeutig sind, warum jeder Generator auf einer Algebra mit Sektoren die in Abschnitt 8.3 verwendete Form hat, wann
der stationäre Zustand eindeutig ist und welche Observablenalgebren ein Wiederherstellungskanal schützen kann.

Im Folgenden ist $M_n=M_n(\mathbb C)$, die $E_{ij}=|i\rangle\langle j|$ sind seine Matrixeinheiten, ein
**Superoperator** ist eine lineare Abbildung auf $M_n$ (oder auf $\mathcal M$), und
$$
\Omega=\sum_{i=1}^ne_i\otimes e_i\in\mathbb C^n\otimes\mathbb C^n,\qquad
J(\Psi)=\sum_{i,j=1}^nE_{ij}\otimes\Psi(E_{ij})=[\Psi(E_{ij})]_{i,j=1}^n
$$
sind der maximal verschränkte Vektor und die **Choi-Matrix** eines Superoperators $\Psi$ (Satz 6.2;
hier muss $\Psi$ nicht positiv sein). Für Vektoren $v,w$ ist $|v\rangle\langle w|$ der Operator
$x\mapsto\langle w,x\rangle v$. Ein Superoperator $\Psi$ heißt **hermitizitätserhaltend**, wenn
$\Psi(A^*)=\Psi(A)^*$ für alle $A$ gilt. Quantendynamische Halbgruppen und ihre Generatoren sind die aus
Abschnitt 8.2 (Lemma 8.5): $T_t=e^{t\mathcal L}$ mit Kanälen $T_t$.

## D.1 Choi-Matrizen der Bausteine

**Lemma D.1.** Sei $\Psi$ ein Superoperator auf $M_n$.

(a) $\Psi\mapsto J(\Psi)$ ist eine lineare Bijektion von den Superoperatoren auf $M_n$ auf $M_n(M_n)=M_{n^2}$.
(b) Jedes $v\in\mathbb C^n\otimes\mathbb C^n$ ist $v=(1\otimes K)\Omega$ für genau ein $K\in M_n$, und
$\langle(1\otimes K)\Omega,(1\otimes L)\Omega\rangle=\operatorname{Tr}(K^*L)$. Insbesondere ist
$\langle\Omega,(1\otimes K)\Omega\rangle=\operatorname{Tr}K$.
(c) Für $K,L\in M_n$ hat der Superoperator $A\mapsto KAL^*$ die Choi-Matrix
$|(1\otimes K)\Omega\rangle\langle(1\otimes L)\Omega|$. Insbesondere ist $J(\mathrm{id})=|\Omega\rangle\langle\Omega|$.
(d) $\Psi$ ist genau dann hermitizitätserhaltend, wenn $J(\Psi)$ selbstadjungiert ist. Jede CP-Abbildung ist hermitizitätserhaltend.

*Beweis.* (a) Der $(i,j)$-Block von $J(\Psi)$ ist $\Psi(E_{ij})$, und die $E_{ij}$ bilden eine Basis von $M_n$; also ist
$J$ linear und injektiv, und jede Blockmatrix $[X_{ij}]$ ist die Choi-Matrix des Superoperators
mit $E_{ij}\mapsto X_{ij}$. (b) $(1\otimes K)\Omega=\sum_ie_i\otimes Ke_i$ hat als $i$-ten Block $Ke_i$, die
$i$-te Spalte von $K$, die $K$ bestimmt, und jedes $v$ entsteht auf diese Weise. Ferner ist
$\langle\sum_ie_i\otimes Ke_i,\sum_je_j\otimes Le_j\rangle=\sum_i\langle Ke_i,Le_i\rangle=\operatorname{Tr}(K^*L)$;
für die letzte Aussage ersetze man $(K,L)$ durch $(1,K)$. (c) Der $(i,j)$-Block der rechten Seite ist
$Ke_i(Le_j)^*=KE_{ij}L^*$. (d) Der $(i,j)$-Block von $J(\Psi)^*$ ist $\Psi(E_{ji})^*$. Ist $\Psi$
hermitizitätserhaltend, so ist $\Psi(E_{ji})^*=\Psi(E_{ji}^*)=\Psi(E_{ij})$. Gilt umgekehrt
$\Psi(E_{ji})^*=\Psi(E_{ij})$ für alle $i,j$ und ist $A=\sum a_{ij}E_{ij}$, so ist
$\Psi(A)^*=\sum\bar a_{ij}\Psi(E_{ij})^*=\sum\bar a_{ij}\Psi(E_{ji})=\Psi(A^*)$. Eine CP-Abbildung hat eine Kraus-Form
(Satz 6.2), und $(\sum_lK_lAK_l^*)^*=\sum_lK_lA^*K_l^*$. $\square$

**Lemma D.2 (Zerlegung).** Seien $q_1,\dots,q_C$ orthonormale Vektoren in
$\mathbb C^n\otimes\mathbb C^n$, geschrieben als $q_c=(1\otimes Q_c)\Omega$, sei $W=\operatorname{span}\{q_c\}$, und sei
$\Psi$ ein hermitizitätserhaltender Superoperator mit $\langle v,J(\Psi)v\rangle\ge0$ für alle $v\perp W$.
Man setze $J=J(\Psi)$ und $\xi_c=Jq_c-\frac12\sum_d\langle q_d,Jq_c\rangle q_d$, und $K_c\in M_n$ sei definiert durch
$(1\otimes K_c)\Omega=\xi_c$. Dann gibt es $L_1,\dots,L_r\in M_n$, $r=\operatorname{rank}(P_W^\perp JP_W^\perp)$
mit der Projektion $P_W^\perp$ auf $W^\perp$, so dass
$$
\Psi(A)=\sum_{j=1}^rL_jAL_j^*+\sum_{c=1}^C\big(K_cAQ_c^*+Q_cAK_c^*\big)\qquad(A\in M_n),
$$
$\operatorname{Tr}(Q_c^*L_j)=0$ für alle $c,j$ und $\operatorname{Tr}(L_j^*L_k)=0$ für $j\neq k$.

*Beweis.* $J$ ist nach Lemma D.1(d) selbstadjungiert. Sei $P=\sum_c|q_c\rangle\langle q_c|$ die Projektion auf $W$
und $P^\perp=P_W^\perp=1-P$. Ausmultiplizieren von $J=(P+P^\perp)J(P+P^\perp)$ mit $PJP^\perp=PJ-PJP$ und
$P^\perp JP=JP-PJP$ ergibt
$$
J=P^\perp JP^\perp+\big(PJ+JP-PJP\big).
$$
Der erste Term ist selbstadjungiert, und $\langle v,P^\perp JP^\perp v\rangle=\langle P^\perp v,JP^\perp v\rangle\ge0$,
weil $P^\perp v\perp W$; also ist $P^\perp JP^\perp\ge0$. Nach Satz 1.3 ist
$P^\perp JP^\perp=\sum_{j=1}^r|l_j\rangle\langle l_j|$, wobei die $l_j$ die mit $\sqrt{\mu_j}$ skalierten Eigenvektoren zu den
von null verschiedenen Eigenwerten $\mu_j$ sind: Sie sind paarweise orthogonal und liegen im Bild von
$P^\perp JP^\perp$, das orthogonal zu $W$ ist. Man schreibe $l_j=(1\otimes L_j)\Omega$; nach Lemma D.1(b) ist
$\operatorname{Tr}(Q_c^*L_j)=\langle q_c,l_j\rangle=0$ und $\operatorname{Tr}(L_j^*L_k)=\langle l_j,l_k\rangle=0$ für $j\ne k$,
und nach Lemma D.1(c) ist $P^\perp JP^\perp$ die Choi-Matrix von $A\mapsto\sum_jL_jAL_j^*$. Für den zweiten Term gilt
$JP=\sum_c|Jq_c\rangle\langle q_c|$ und $PJP=\sum_{c,d}\langle q_d,Jq_c\rangle|q_d\rangle\langle q_c|$, also
$$
\sum_c|\xi_c\rangle\langle q_c|=JP-\tfrac12PJP,\qquad\sum_c|q_c\rangle\langle\xi_c|=PJ-\tfrac12PJP
$$
(die zweite Gleichung ist die adjungierte der ersten; $PJP$ ist selbstadjungiert), und ihre Summe ist $PJ+JP-PJP$. Nach
Lemma D.1(c) sind $|\xi_c\rangle\langle q_c|$ und $|q_c\rangle\langle\xi_c|$ die Choi-Matrizen von
$A\mapsto K_cAQ_c^*$ und $A\mapsto Q_cAK_c^*$. Also haben beide Seiten der behaupteten Identität die Choi-Matrix $J$,
und Lemma D.1(a) liefert die Identität. $\square$

## D.2 Bedingte vollständige Positivität

Ein Superoperator $\Psi$ auf $M_n$ heißt **bedingt vollständig positiv** (CCP), wenn
$\langle v,J(\Psi)v\rangle\ge0$ für jedes $v\perp\Omega$ gilt. Nach Lemma D.1(b) sind diese $v$ die Vektoren
$(1\otimes K)\Omega$ mit $\operatorname{Tr}K=0$.

**Satz D.3.** Für einen Superoperator $\mathcal L$ auf $M_n$ sind äquivalent:

1. $e^{t\mathcal L}$ ist CP für jedes $t\ge0$;
2. $\mathcal L$ ist hermitizitätserhaltend und bedingt vollständig positiv;
3. es gibt eine CP-Abbildung $\Phi$ auf $M_n$ und $G\in M_n$ mit $\mathcal L(A)=\Phi(A)+GA+AG^*$ für alle $A$.

In 3 kann man $\Phi(A)=\sum_{j=1}^rL_jAL_j^*$ mit $r\le n^2-1$, $\operatorname{Tr}L_j=0$ und
$\operatorname{Tr}(L_j^*L_k)=0$ für $j\neq k$ sowie $G$ mit $\operatorname{Tr}G\in\mathbb R$ wählen.

*Beweis.* (1)$\Rightarrow$(2). Nach Lemma 8.1 gilt $(e^{t\mathcal L}-\mathrm{id})/t\to\mathcal L$ für $t\to0$, und $J$ ist
linear auf einem endlichdimensionalen Raum, also stetig; daher ist
$J(e^{t\mathcal L})=J(\mathrm{id})+t\,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)$, wobei die letzte Choi-Matrix gegen
$J(\mathcal L)$ strebt. Für $v\perp\Omega$ ist $\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0$ (Lemma D.1(c)), also
gilt für $t>0$, da $J(e^{t\mathcal L})\ge0$ (Satz 6.2),
$$
0\le\frac1t\langle v,J(e^{t\mathcal L})v\rangle=\big\langle v,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)v\big\rangle
\;\xrightarrow[t\to0]{}\;\langle v,J(\mathcal L)v\rangle .
$$
Jedes $e^{t\mathcal L}$ ist CP, also hermitizitätserhaltend (Lemma D.1(d)); daher ist
$\mathcal L(A^*)=\lim_t(e^{t\mathcal L}(A)^*-A^*)/t=\mathcal L(A)^*$.

(2)$\Rightarrow$(3). Man wende Lemma D.2 mit $C=1$ und $q_1=\Omega/\sqrt n$ an, d. h. $Q_1=1/\sqrt n$:
$\mathcal L(A)=\sum_jL_jAL_j^*+GA+AG^*$ mit $G=K_1/\sqrt n$, und $\Phi=\sum_jL_j\cdot L_j^*$ ist CP
(Satz 6.2). Ferner ist $\operatorname{Tr}L_j=\sqrt n\operatorname{Tr}(Q_1^*L_j)=0$, $r\le\operatorname{rank}P_W^\perp=n^2-1$,
und nach Lemma D.1(b) ist
$\operatorname{Tr}G=\langle\Omega,\xi_1\rangle/\sqrt n=\langle q_1,\xi_1\rangle=\frac12\langle q_1,J(\mathcal L)q_1\rangle\in\mathbb R$,
da $J(\mathcal L)$ selbstadjungiert ist.

(3)$\Rightarrow$(1). Die Schritte 1–4 des Beweises von Satz 8.6 verwenden nur, dass der Generator die Form
$\mathcal L_0+\Phi$ hat, mit $\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)$ für eine Matrix
$H_{\mathrm{eff}}$ und CP-Abbildung $\Phi$. Mit $H_{\mathrm{eff}}=iG$ ist $\mathcal L_0(\rho)=G\rho+\rho G^*$, und diese
Schritte zeigen, dass $e^{t\mathcal L}=\sum_{m\ge0}S_m(t)$ CP ist. $\square$

Positivität allein würde nicht genügen. Sei $n\ge2$, $\theta$ die Transposition auf $M_n$ und $\mathcal L=\theta-\mathrm{id}$. Wegen $\theta^2=\mathrm{id}$ ist
$$
e^{t\mathcal L}=e^{-t}\big(\cosh t\,\mathrm{id}+\sinh t\,\theta\big)=\tfrac{1+e^{-2t}}2\,\mathrm{id}+\tfrac{1-e^{-2t}}2\,\theta ,
$$
eine Konvexkombination positiver spurerhaltender Abbildungen (Proposition 6.1). Aber $J(\theta)=[E_{ji}]=\sum_{i,j}E_{ij}\otimes E_{ji}$ ist der Swap-Operator $e_k\otimes e_l\mapsto e_l\otimes e_k$ und bildet daher
$v=e_1\otimes e_2-e_2\otimes e_1$, das orthogonal zu $\Omega$ ist, auf $-v$ ab (wie im Beweis von Proposition 6.1); mit
$\|v\|^2=2$ und $\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0$ folgt
$\langle v,J(\mathcal L)v\rangle=-2$: $\mathcal L$ ist nicht CCP, und tatsächlich ist
$\langle v,J(e^{t\mathcal L})v\rangle=-(1-e^{-2t})<0$ für jedes $t>0$. Eine Halbgruppe positiver spurerhaltender Abbildungen
besteht nicht notwendig aus Kanälen.

Wir prüfen Satz D.3 an einem GKSL-Generator mit beliebigen Sprüngen auf $M_3$, erzeugt mit `lindblad_superoperator`
(die Matrix $S$ von $\mathcal L$ auf dem zeilenweise gebildeten $\operatorname{vec}\rho$, so dass $J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$ wie in
Abschnitt 6.1), und an $\theta-\mathrm{id}$:

```python
import math
import numpy as np
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, channels, dynamics

torch.manual_seed(0)
rng = np.random.default_rng(0)
f64, c128 = torch.float64, torch.complex128

def choi_of(S, n):
    """J_{(i,a),(j,b)} = S_{(a,b),(i,j)}: the Choi matrix of the map with row-major matrix S (Section 6.1)."""
    return S.reshape(n, n, n, n).transpose(2, 0, 3, 1).reshape(n * n, n * n)

def superop(f, n):
    """Row-major matrix of a linear map f on n x n matrices."""
    E = np.eye(n * n).reshape(n * n, n, n)
    return np.stack([f(X).reshape(-1) for X in E], axis=1)

vec_of = lambda K: K.T.reshape(-1)                        # (1 (x) K) Omega = sum_i e_i (x) K e_i
op_of = lambda v, n: v.reshape(n, n).T                    # its inverse (Lemma D.1(b))
cplx = lambda *s: rng.standard_normal(s) + 1j * rng.standard_normal(s)

n = 3
M3 = TypeIAlgebra([n], [n], precision='double')
H = cplx(n, n)
H = (H + H.conj().T) / 2
Ls = [cplx(n, n) / 2 for _ in range(2)]                   # arbitrary jumps, not traceless
to_op = lambda X: M3.operator(torch.from_numpy(X).reshape(1, 1, n, n))
S = dynamics.lindblad_superoperator(M3, to_op(H), [to_op(L) for L in Ls])[0][0].numpy()
J = choi_of(S, n)
q = vec_of(np.eye(n)) / math.sqrt(n)                      # the unit vector Omega / sqrt(n)
Pperp = np.eye(n * n) - np.outer(q, q)
assert np.abs(J - J.conj().T).max() < 1e-12               # Hermiticity preserving (Lemma D.1(d))
print(f"min eigenvalue of J(L): {np.linalg.eigvalsh(J).min():.2f}; "
      f"on the complement of Omega: {np.linalg.eigvalsh(Pperp @ J @ Pperp).min():.1e}")
assert np.linalg.eigvalsh(J).min() < -0.1                 # L itself is not CP ...
assert np.linalg.eigvalsh(Pperp @ J @ Pperp).min() > -1e-12        # ... but CCP (Theorem D.3)
for t in (0.01, 0.1, 1.0):
    assert np.linalg.eigvalsh(choi_of(sla.expm(t * S), n)).min() > -1e-12

L_theta = superop(lambda X: X.T, n) - np.eye(n * n)     # L = transpose - id: Hermiticity preserving, not CCP
anti = vec_of(np.eye(n)[:, [1]] @ np.eye(n)[[0], :] - np.eye(n)[:, [0]] @ np.eye(n)[[1], :])   # e1(x)e2 - e2(x)e1
assert abs(anti @ q) < 1e-15 and abs(anti @ choi_of(L_theta, n) @ anti + 2) < 1e-12
X = cplx(n, n)
for t in (0.01, 0.1, 1.0):
    T_t = sla.expm(t * L_theta)
    b = (1 - math.exp(-2 * t)) / 2
    assert abs(anti @ choi_of(T_t, n) @ anti + (1 - math.exp(-2 * t))) < 1e-12        # <v, J(e^{tL}) v>
    assert np.abs(T_t - ((1 - b) * np.eye(n * n) + b * superop(lambda Y: Y.T, n))).max() < 1e-12
    assert abs(np.trace((T_t @ X.reshape(-1)).reshape(n, n)) - np.trace(X)) < 1e-12   # trace preserving
    assert np.linalg.eigvalsh(choi_of(T_t, n)).min() < -0.9 * b                     # but not CP
```

## D.3 Der GKSL-Satz

**Satz D.4 (Gorini–Kossakowski–Sudarshan; Lindblad; Satz 8.7).** Ein Superoperator $\mathcal L$ auf
$M_n$ ist genau dann der Generator einer quantendynamischen Halbgruppe auf $M_n$, wenn es $H=H^*$ und
$L_1,\dots,L_r$ in $M_n$ gibt mit
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_{j=1}^r\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
$$
Man kann $r\le n^2-1$, $\operatorname{Tr}H=0$, $\operatorname{Tr}L_j=0$ und $\operatorname{Tr}(L_j^*L_k)=0$ für $j\neq k$ wählen.

*Beweis.* „Wenn“ ist Satz 8.6. Sei umgekehrt $(T_t)$ eine quantendynamische Halbgruppe mit Generator
$\mathcal L$, also $T_t=e^{t\mathcal L}$ (Lemma 8.5). Da jedes $T_t$ CP ist, liefert Satz D.3
$\mathcal L(\rho)=\sum_jL_j\rho L_j^*+G\rho+\rho G^*$ mit $r\le n^2-1$ spurfreien, paarweise orthogonalen $L_j$ und
$\operatorname{Tr}G\in\mathbb R$. Man setze $H=\frac i2(G-G^*)$ und $R=\frac12(G+G^*)$, beide selbstadjungiert; dann ist
$G=-iH+R$, $G^*=iH+R$ und $G\rho+\rho G^*=-i[H,\rho]+\{R,\rho\}$. Außerdem ist
$\operatorname{Tr}H=\frac i2(\operatorname{Tr}G-\overline{\operatorname{Tr}G})=0$. Da jedes $T_t$ TP ist, gilt
$\operatorname{Tr}\mathcal L(\rho)=\frac d{dt}\operatorname{Tr}T_t(\rho)\big|_{t=0}=0$ (Lemma 8.1), und wegen der Zyklizität
$$
0=\operatorname{Tr}\mathcal L(\rho)=\operatorname{Tr}\Big(\rho\,\Big(\sum_jL_j^*L_j+2R\Big)\Big)\qquad\text{for all }\rho .
$$
Mit $\rho=X=\sum_jL_j^*L_j+2R$, das selbstadjungiert ist, folgt $\operatorname{Tr}(X^*X)=0$, also $X=0$ und
$R=-\frac12\sum_jL_j^*L_j$; das ist die GKSL-Form. $\square$

**Satz D.5 (Freiheit in der GKSL-Form).** Seien $(H;L_1,\dots,L_r)$ und $(H';L'_1,\dots,L'_s)$ GKSL-Daten
($H,H'$ selbstadjungiert) mit Generatoren $\mathcal L$, $\mathcal L'$ wie in Satz D.4.

(a) $\mathcal L'=\mathcal L$ in jedem der folgenden Fälle:
(i) $L'_j=L_j$ und $H'=H+c1$, $c\in\mathbb R$;
(ii) $H'=H$ und, nach Auffüllen der kürzeren Familie mit Nullen, $L'_j=\sum_ku_{jk}L_k$ mit einer unitären Matrix $u$;
(iii) $L'_j=L_j+a_j1$ mit $a_j\in\mathbb C$ und $H'=H-\frac i2\sum_j(\bar a_jL_j-a_jL_j^*)$, das selbstadjungiert ist.

(b) Ist $\operatorname{Tr}H=\operatorname{Tr}H'=0$, sind alle $L_j,L'_j$ spurfrei und ist $\mathcal L=\mathcal L'$, so ist $H=H'$ und,
nach Auffüllen, $L'_j=\sum_ku_{jk}L_k$ mit einer unitären Matrix $u$.

(c) Folglich haben zwei GKSL-Daten genau dann denselben Generator, wenn die einen aus den anderen durch endlich viele
Schritte (i)–(iii) hervorgehen. Insbesondere sind der Operator $H$ der normierten Form aus Satz D.4 und ihr
dissipativer Anteil $\Phi(\rho)=\sum_jL_j\rho L_j^*$ durch $\mathcal L$ bestimmt.

*Beweis.* (a) (i) $[c1,\rho]=0$. (ii) Nach Satz 6.7 ist $\Phi'=\sum_jL'_j\cdot L_j'^*$ gleich
$\Phi=\sum_jL_j\cdot L_j^*$; die duale Abbildung hängt nicht von den Kraus-Operatoren ab (Abschnitt 6.2), also ist
$\sum_jL_j'^*L'_j=\Phi^*(1)=\sum_jL_j^*L_j$. (iii) Für einen einzelnen Sprung $L$ und $a\in\mathbb C$ ergibt Ausmultiplizieren
$$
(L+a)\rho(L+a)^*-\tfrac12\{(L+a)^*(L+a),\rho\}-\Big(L\rho L^*-\tfrac12\{L^*L,\rho\}\Big)=\tfrac12[B,\rho],
\qquad B=\bar aL-aL^* ,
$$
denn die Terme $\bar aL\rho+a\rho L^*$ aus dem ersten Produkt verbinden sich mit
$-\frac12(aL^*\rho+a\rho L^*+\bar aL\rho+\bar a\rho L)$, und die Terme mit $|a|^2$ heben sich weg. Wegen $B^*=-B$ ist die
Matrix $\frac i2B$ selbstadjungiert, also ist $H'=H-\frac i2\sum_jB_j$ selbstadjungiert, und
$-i[H'-H,\rho]=-\frac12\sum_j[B_j,\rho]$ hebt die Änderung des Dissipators auf.

(b) Sei $\Phi=\sum_jL_j\cdot L_j^*$ und $G=-iH-\frac12\Phi^*(1)$, so dass $\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*$.
Nach Lemma D.1(c) ist
$J(\mathcal L)=\sum_j|l_j\rangle\langle l_j|+|g\rangle\langle\Omega|+|\Omega\rangle\langle g|$ mit $l_j=(1\otimes L_j)\Omega$ und
$g=(1\otimes G)\Omega$, und $l_j\perp\Omega$, weil $\operatorname{Tr}L_j=0$ (Lemma D.1(b)). Mit der Projektion $P^\perp$
auf $\Omega^\perp$ ist daher $P^\perp J(\mathcal L)P^\perp=\sum_j|l_j\rangle\langle l_j|=J(\Phi)$, und ebenso
$P^\perp J(\mathcal L')P^\perp=J(\Phi')$. Also ist $J(\Phi)=J(\Phi')$, $\Phi=\Phi'$ (Lemma D.1(a)), und Satz 6.7 liefert die
unitäre Matrix $u$. Dann ist $\Phi^*(1)=\Phi'^*(1)$, also $[H-H',\rho]=0$ für alle $\rho$. Eine Matrix $X$, die mit allen
$E_{ij}$ kommutiert, ist ein Vielfaches von $1$ (man vergleiche die Einträge von $XE_{ij}$ und $E_{ij}X$), also ist $H-H'=\lambda1$, und
$\lambda=0$ durch Spurbildung.

(c) Nach (a) ändern die Schritte (i)–(iii) den Generator nicht, und die Umkehrung jedes Schritts ist ein Schritt
derselben Art ((iii) mit $-a_j$, da $\bar a_jL'_j-a_jL_j'^*=\bar a_jL_j-a_jL_j^*$). Beliebige Daten werden durch (iii) mit
$a_j=-\operatorname{Tr}L_j/n$ und anschließend (i) mit $c=-\operatorname{Tr}H''/n$, wobei $H''$ der
Hamiltonoperator nach (iii) ist, in normierte Form gebracht; zwei normierte Daten mit demselben Generator hängen nach (b) über (ii) zusammen. $\square$

Der Beweis von Satz D.4 ist ein Algorithmus: $\Phi$ aus der Kompression von $J(\mathcal L)$ auf $\Omega^\perp$, dann
$G$ aus $\xi_1$, dann $H$. Wir führen ihn für den obigen Generator aus, bauen den Generator mit der Bibliothek wieder auf und
vergleichen mit der normierten Form der ursprünglichen Daten, berechnet nach Satz D.5(iii) und (i):

```python
def gksl_from_generator(S, n, tol=1e-10):
    """(H, [L_j]) with Tr H = Tr L_j = 0 from the matrix of a GKSL generator: the proof of Theorem D.4."""
    J = choi_of(S, n)
    q = vec_of(np.eye(n)) / math.sqrt(n)
    Pperp = np.eye(n * n) - np.outer(q, q)
    w, U = np.linalg.eigh(Pperp @ J @ Pperp)              # = J(Phi) >= 0
    jumps = [math.sqrt(x) * op_of(U[:, k], n) for k, x in enumerate(w) if x > tol]
    xi = J @ q - 0.5 * (q.conj() @ J @ q) * q
    G = op_of(xi, n) / math.sqrt(n)
    return 0.5j * (G - G.conj().T), jumps

H_rec, L_rec = gksl_from_generator(S, n)
S_rec = dynamics.lindblad_superoperator(M3, to_op(H_rec), [to_op(L) for L in L_rec])[0][0].numpy()
assert np.abs(S_rec - S).max() < 1e-12                    # the same generator
assert abs(np.trace(H_rec)) < 1e-12 and all(abs(np.trace(L)) < 1e-12 for L in L_rec)

a = [-np.trace(L) / n for L in Ls]                        # Theorem D.5(iii): make the jumps traceless
L_norm = [L + aj * np.eye(n) for L, aj in zip(Ls, a)]
H_norm = H - 0.5j * sum(np.conj(aj) * L - aj * L.conj().T for L, aj in zip(Ls, a))
H_norm = H_norm - np.trace(H_norm) / n * np.eye(n)        # Theorem D.5(i)
assert np.abs(H_norm - H_norm.conj().T).max() < 1e-12
assert np.abs(H_rec - H_norm).max() < 1e-12               # H is unique (Theorem D.5(b))
jump_choi = lambda Ks: sum(np.outer(vec_of(K), vec_of(K).conj()) for K in Ks)
assert np.abs(jump_choi(L_rec) - jump_choi(L_norm)).max() < 1e-12      # so is sum_j L_j . L_j^*
Wr, Wn = np.stack([vec_of(K) for K in L_rec], 1), np.stack([vec_of(K) for K in L_norm], 1)
u = np.linalg.lstsq(Wr, Wn, rcond=None)[0]                # L_norm_j = sum_k u_kj L_rec_k
assert np.abs(Wr @ u - Wn).max() < 1e-12 and np.abs(u.conj().T @ u - np.eye(2)).max() < 1e-12   # unitary
print("eigenvalues of the normalised H:", np.round(np.linalg.eigvalsh(H_rec), 4), " jumps:", len(L_rec))
```

**Satz D.6 (Eindeutigkeit des stationären Zustands).** Sei $\mathcal L$ der GKSL-Generator aus Satz D.4
mit den Daten $(H;L_1,\dots,L_r)$. Hat $\mathcal L$ einen treuen stationären Zustand $\sigma$ ($\sigma>0$) und sind die einzigen
Matrizen, die mit $H$ und allen $L_j$, $L_j^*$ kommutieren, die Vielfachen von $1$, so ist $\ker\mathcal L=\mathbb C\sigma$:
$\sigma$ ist der einzige stationäre Zustand. Außerdem ist $\ker\mathcal L^*=\mathbb C1$.

*Beweis.* Die duale Abbildung bezüglich der Spurpaarung, $\operatorname{Tr}(\mathcal L(\rho)A)=\operatorname{Tr}(\rho\,\mathcal L^*(A))$, ist
$\mathcal L^*(A)=i[H,A]+\sum_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})$: Wegen der Zyklizität ist
$\operatorname{Tr}(-i[H,\rho]A)=\operatorname{Tr}(\rho\,i[H,A])$, $\operatorname{Tr}(L\rho L^*A)=\operatorname{Tr}(\rho L^*AL)$ und
$\operatorname{Tr}(\{L^*L,\rho\}A)=\operatorname{Tr}(\rho\{L^*L,A\})$ (das ist Aufgabe 1 aus Kapitel 8).
Direktes Ausmultiplizieren, bei dem sich die Hamiltonschen Terme wegheben, weil $A\mapsto i[H,A]$ eine Derivation ist, ergibt
für jedes $X$
$$
\mathcal L^*(X^*X)-X^*\mathcal L^*(X)-\mathcal L^*(X^*)X=\sum_j[L_j,X]^*[L_j,X]:
$$
Für einen einzelnen Sprung heben sich von den sechs Antikommutatortermen vier paarweise weg, und die beiden anderen,
$\frac12X^*L^*LX+\frac12X^*L^*LX$, addieren sich, so dass die linke Seite $L^*X^*XL-X^*L^*XL-L^*X^*LX+X^*L^*LX$ ist, also
$(XL-LX)^*(XL-LX)$. Sei $X\in\ker\mathcal L^*$. Dann ist $\mathcal L^*(X^*)=\mathcal L^*(X)^*=0$, da $\mathcal L^*$
hermitizitätserhaltend ist, also $\mathcal L^*(X^*X)=\sum_j[L_j,X]^*[L_j,X]=:Y\ge0$, und
$\operatorname{Tr}(\sigma Y)=\operatorname{Tr}(\mathcal L(\sigma)X^*X)=0$. Wegen $\sigma>0$ verschwindet
$\operatorname{Tr}(\sigma Y)=\operatorname{Tr}(\sigma^{1/2}Y\sigma^{1/2})$ mit $\sigma^{1/2}Y\sigma^{1/2}\ge0$ nur, wenn
$\sigma^{1/2}Y\sigma^{1/2}=0$ ist, d. h. $Y=0$; Spurbildung ergibt $\sum_j\operatorname{Tr}([L_j,X]^*[L_j,X])=0$ mit nichtnegativen
Summanden, also $[L_j,X]=0$ für alle $j$. Dasselbe Argument für $X^*\in\ker\mathcal L^*$
liefert $[L_j,X^*]=0$, d. h. $[L_j^*,X]=0$. Dann kommutiert $X$ mit $L_j^*L_j$, der dissipative Anteil von $\mathcal L^*(X)$
verschwindet, und $0=\mathcal L^*(X)=i[H,X]$. Nach Voraussetzung ist $X\in\mathbb C1$; umgekehrt ist $\mathcal L^*(1)=0$. Also ist
$\ker\mathcal L^*=\mathbb C1$. In den Basen $(E_{ij})$ und $(E_{ji})$ von $M_n$, die bezüglich der Paarung
$(\rho,A)\mapsto\operatorname{Tr}(\rho A)$ dual zueinander sind, ist die Matrix von $\mathcal L^*$ die Transponierte derjenigen von $\mathcal L$; also haben beide
denselben Rang, und $\dim\ker\mathcal L=\dim\ker\mathcal L^*=1$. Ein Zustand in $\mathbb C\sigma$ mit Spur $1$ ist $\sigma$. $\square$

Dies ist der endlichdimensionale Fall eines Kriteriums von Frigerio [*Commun. Math. Phys.* 63, 269 (1978)],
das in Abschnitt 8.5 verwendet wird. Beide Voraussetzungen sind nötig. In Aufgabe 4 aus Kapitel 8 ist der Zustand $\tfrac12\,1$ treu
und stationär, aber $H$ und $\sigma_z$ kommutieren mit allen Diagonalmatrizen, und die stationären Zustände spannen den
zweidimensionalen Raum der Diagonalmatrizen auf. Umgekehrt seien auf $M_3$ $H=0$, $L_1=|e_1\rangle\langle e_3|$,
$L_2=|e_2\rangle\langle e_3|$. Die von $L_1,L_2$ erzeugte $*$-Algebra enthält $L_1L_2^*=E_{12}$, $L_1^*L_1=E_{33}$ und damit
alle Matrixeinheiten, also kommutieren nur die Vielfachen von $1$ mit $L_j,L_j^*$; aber jedes $\rho=P\rho P$ mit $P=E_{11}+E_{22}$
ist stationär, da $L_j\rho=0=\rho L_j^*$ und $\{E_{33},\rho\}=0$, also ist $\dim\ker\mathcal L\ge4$; nach Satz D.6 ist kein stationärer
Zustand treu. Für den obigen zufälligen Generator prüft der Code die Voraussetzungen und die Folgerung
(die Matrix von $\mathcal L^*$ auf den zeilenweise gebildeten Vektoren von $A^T$ ist $S^T$: $\operatorname{vec}(\mathcal L^*(A)^T)=S^T\operatorname{vec}(A^T)$) sowie das zweite Beispiel:

```python
def kernel_dim(S, tol=1e-9):
    return int((np.linalg.svd(S, compute_uv=False) < tol).sum())

def commutant_dim(ops, n):                    # dim {X : XA = AX for all A in ops}
    return kernel_dim(np.concatenate([np.kron(A, np.eye(n)) - np.kron(np.eye(n), A.T) for A in ops]))

w, V = np.linalg.eig(S)
sigma = V[:, np.argmin(abs(w))].reshape(n, n)             # row-major: vec(X)_{(a,b)} = X_ab
sigma = sigma / np.trace(sigma)
assert np.abs(S @ sigma.reshape(-1)).max() < 1e-12        # stationary
assert np.abs(sigma - sigma.conj().T).max() < 1e-10 and np.linalg.eigvalsh(sigma).min() > 1e-3   # faithful
assert commutant_dim([H] + Ls + [L.conj().T for L in Ls], n) == 1
assert kernel_dim(S) == 1 and kernel_dim(S.T) == 1        # Theorem D.6
assert np.abs(S.T @ np.eye(n).reshape(-1)).max() < 1e-12  # ker L^* = C 1

e3 = np.eye(3)
L_abs = [np.outer(e3[0], e3[2]), np.outer(e3[1], e3[2])]       # L_1 = |e_1><e_3|, L_2 = |e_2><e_3|
S_abs = superop(lambda R: sum(L @ R @ L.T - 0.5 * (L.T @ L @ R + R @ L.T @ L) for L in L_abs), 3)
assert commutant_dim(L_abs + [L.T for L in L_abs], 3) == 1
print("second example: dim ker L =", kernel_dim(S_abs))
assert kernel_dim(S_abs) == 4                             # no faithful stationary state, no uniqueness
```

## D.4 Algebren mit Sektoren

Sei $\mathcal M=\bigoplus_{c=1}^CM_{k_c}\subseteq M_D$ die Algebra der blockdiagonalen Matrizen, $D=\sum_ck_c$, mit
minimalen zentralen Projektionen $P_c$, und sei $E(X)=\sum_cP_cXP_c$ die bedingte Erwartung von $M_D$ auf
$\mathcal M$. Wie im Beweis von Proposition 8.8 gilt $E(B)=B$, $E(XB)=E(X)B$ und $E(BX)=BE(X)$ für $B\in\mathcal M$
sowie $E(X^*)=E(X)^*$. Für eine CP-Abbildung $\Phi:\mathcal M\to\mathcal M$ ist $\Phi^*$ ihre ungewichtete duale Abbildung (Proposition 6.8).

**Satz D.7.** Sei $\mathcal L$ der Generator einer quantendynamischen Halbgruppe $(T_t)$ auf $\mathcal M$. Dann gibt
es eine CP-Abbildung $\Phi:\mathcal M\to\mathcal M$ und $H=H^*\in\mathcal M$ mit
$$
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]\qquad(\rho\in\mathcal M).
$$
Nach Proposition 6.4 ist $\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l$ mit Kraus-Operatoren zwischen den Sektoren;
umgekehrt erzeugt jedes solche $\mathcal L$ eine quantendynamische Halbgruppe auf $\mathcal M$ (Proposition 8.8(1)).
Also sind die Generatoren aus Abschnitt 8.3 (mit $\Phi=\sum_jg_jL_j\cdot L_j^*+\sum_mg_m\Phi_m$) genau die
Generatoren quantendynamischer Halbgruppen auf $\mathcal M$.

*Beweis.* *Schritt 1: eine CP-Abbildung auf $M_D$.* Nach Proposition 6.4 ist
$T_t(\rho)_d=\sum_{c,l}K^{dc}_l(t)\rho_cK^{dc}_l(t)^*$. Sei $N^{dc}_l(t)\in M_D$ die Matrix mit dem Block $K^{dc}_l(t)$ an der
Position $(d,c)$ und Nullen sonst. Für $X\in M_D$ ist $N^{dc}_lXN^{dc*}_l$ der Block
$K^{dc}_lX_{cc}K^{dc*}_l$ an der Position $(d,d)$, wobei $X_{cc}$ der $c$-te Diagonalblock ist, d. h. der $c$-te Block
von $E(X)$. Daher ist $\hat T_t=T_t\circ E=\sum_{d,c,l}N^{dc}_l\cdot N^{dc*}_l$ CP auf $M_D$ (Satz 6.2), mit
$\hat T_0=E$. Nach Lemma 8.5 gilt $(\hat T_t-E)/t=\big((T_t-\mathrm{id})/t\big)\circ E\to\mathcal L\circ E$.

*Schritt 2: Positivität abseits der Sektorvektoren.* Sei $p_c=(1\otimes P_c)\Omega\in\mathbb C^D\otimes\mathbb C^D$ und
$W=\operatorname{span}\{p_c\}$. Nach Lemma D.1(c) ist $J(E)=\sum_c|p_c\rangle\langle p_c|$, also gilt für $v\perp W$ und $t>0$
$0\le\frac1t\langle v,J(\hat T_t)v\rangle=\langle v,J((\hat T_t-E)/t)v\rangle\to\langle v,J(\mathcal L\circ E)v\rangle$.
Jedes $T_t$ ist hermitizitätserhaltend (wegen seiner Kraus-Form), also auch $\mathcal L$ (wie in Satz D.3), und ebenso $E$
und $\mathcal L\circ E$.

*Schritt 3: Zerlegung.* Nach Lemma D.1(b) ist $\langle p_c,p_d\rangle=\operatorname{Tr}(P_cP_d)=\delta_{cd}k_c$, also sind die
$q_c=p_c/\sqrt{k_c}=(1\otimes Q_c)\Omega$ mit $Q_c=P_c/\sqrt{k_c}$ orthonormal. Lemma D.2 für
$\Psi=\mathcal L\circ E$ liefert eine CP-Abbildung $\hat\Phi=\sum_j\hat L_j\cdot\hat L_j^*$ auf $M_D$ und $K_c\in M_D$ mit
$$
\mathcal L(E(A))=\hat\Phi(A)+\sum_c\big(K_cAP_c+P_cAK_c^*\big)/\sqrt{k_c}\qquad(A\in M_D).
$$
Für $A\in\mathcal M$ ist $E(A)=A$ und $AP_c=P_cA$, also ist die Summe gleich $\tilde GA+A\tilde G^*$ mit
$\tilde G=\sum_cK_cP_c/\sqrt{k_c}$.

*Schritt 4: Kompression.* Für $A\in\mathcal M$ ist $\mathcal L(A)\in\mathcal M$, also
$\mathcal L(A)=E(\mathcal L(A))=E(\hat\Phi(A))+GA+AG^*$ mit $G=E(\tilde G)\in\mathcal M$. Die Abbildung
$\Phi=E\circ\hat\Phi$ auf $\mathcal M$ ist
$\Phi(A)_d=\sum_{c,j}(P_d\hat L_jP_c)A_c(P_d\hat L_jP_c)^*$: Sie ist CP nach Proposition 6.4, mit den Blöcken der $\hat L_j$ als
Kraus-Operatoren.

*Schritt 5: Spurerhaltung.* Wie in Satz D.4 ist $\operatorname{Tr}\mathcal L(A)=0$ für $A\in\mathcal M$, und
$\operatorname{Tr}\Phi(A)=\operatorname{Tr}(A\Phi^*(1))$, also $\operatorname{Tr}(AX)=0$ für alle $A\in\mathcal M$, wobei
$X=\Phi^*(1)+G+G^*\in\mathcal M$ selbstadjungiert ist. Mit $A=X$ folgt $X=0$ (die ungewichtete Spur ist treu, Satz 3.1).
Mit $H=\frac i2(G-G^*)\in\mathcal M$ ist $G=-iH-\frac12\Phi^*(1)$ und $G\rho+\rho G^*=-i[H,\rho]-\frac12\{\Phi^*(1),\rho\}$. $\square$

Die Vektoren $p_c$ treten an die Stelle von $\Omega=\sum_cp_c$, weil $\hat T_0=E$ ist und nicht $\mathrm{id}$: $J(\mathcal L\circ E)$ muss
auf $\Omega^\perp$ nicht positiv sein, und im folgenden Beispiel ist sie es nicht.

**Korollar D.8 (Heisenberg-Bild).** Sei $(S_t)_{t\ge0}$ eine stetige Familie unitaler CP-Abbildungen auf
$\mathcal M$ mit $S_0=\mathrm{id}$ und $S_{s+t}=S_sS_t$. Ihr Generator ist $A\mapsto\Psi(A)+K^*A+AK$ mit
einer CP-Abbildung $\Psi$ auf $\mathcal M$, $K\in\mathcal M$ und $\Psi(1)+K+K^*=0$.

*Beweis.* Sei $T_t=S_t^*$ (ungewichtete duale Abbildungen). Sie sind CP (Proposition 6.8) und TP (Korollar 6.9(b)),
$T_0=\mathrm{id}$, und $T_{s+t}=(S_tS_s)^*=S_s^*S_t^*=T_sT_t$, weil
$\operatorname{Tr}(S_tS_s(x)y)=\operatorname{Tr}(xS_s^*S_t^*(y))$ und $S_tS_s=S_{s+t}$. Da $\Psi\mapsto\Psi^*$ linear auf einem
endlichdimensionalen Raum ist, ist $t\mapsto T_t$ stetig und $(T_t-\mathrm{id})/t\to\mathcal L_S^*$, wobei $\mathcal L_S$ der
Generator von $(S_t)$ ist. Satz D.7 liefert $\mathcal L_S^*(\rho)=\Phi(\rho)+G\rho+\rho G^*$ mit $G+G^*=-\Phi^*(1)$, und wegen der
Zyklizität ist $\operatorname{Tr}((G\rho+\rho G^*)A)=\operatorname{Tr}(\rho(G^*A+AG))$, also nach Korollar 6.9(c)
$\mathcal L_S(A)=\Phi^*(A)+G^*A+AG$. Man nehme $\Psi=\Phi^*$ (CP nach Proposition 6.8) und $K=G$. $\square$

Dies ist der endlichdimensionale Fall des nach Satz 8.7 zitierten Satzes von Christensen und Evans. Wir
prüfen Satz D.7 auf $M_1\oplus M_2\subseteq M_3$ für einen Generator von `lindblad_rhs` mit einem Sprung innerhalb der
Sektoren und einem zufälligen Sprung zwischen ihnen und folgen dabei den fünf Schritten des Beweises:

```python
alg = TypeIAlgebra([1, 2], [1, 2], precision='double')    # M = M_1 (+) M_2 in M_3
D, ks = 3, alg.k_factors
P = [np.diag([1.0, 0, 0]), np.diag([0, 1.0, 1.0])]
E = lambda X: sum(Pc @ X @ Pc for Pc in P)               # the conditional expectation onto M
def blocks(X):                                            # dense block-diagonal 3 x 3 -> library layout
    out = torch.zeros(1, 2, 2, 2, dtype=c128)
    out[0, 0, :1, :1] = torch.from_numpy(X[:1, :1])
    out[0, 1] = torch.from_numpy(X[1:, 1:])
    return out
dense = lambda Y: sla.block_diag(Y[0, 0, :1, :1].numpy(), Y[0, 1].numpy())
Hs = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
inside = alg.operator(torch.randn(1, 2, 2, 2, dtype=c128))       # a jump inside the sectors
between = channels.random_inter_sector_channel(alg, alg, kraus_rank=2)
rhs = dynamics.lindblad_rhs(Hs, [inside, between], rates=[0.7, 1.3], algebra=alg)
L_M = lambda X: dense(rhs(0.0, blocks(X)))               # the generator on M (X block diagonal)
J = choi_of(superop(lambda X: L_M(E(X)), D), D)          # J(L o E) on M_3
p = [vec_of(Pc) for Pc in P]
Q = np.stack([pc / np.linalg.norm(pc) for pc in p], 1)    # orthonormal q_c = p_c / sqrt(k_c)
Pw = np.eye(D * D) - Q @ Q.T                              # projection onto the complement of W
assert np.linalg.eigvalsh(Pw @ J @ Pw).min() > -1e-12     # Step 2
w, U = np.linalg.eigh(Pw @ J @ Pw)                        # Step 3 (Lemma D.2)
jumps = [math.sqrt(x) * op_of(U[:, k], D) for k, x in enumerate(w) if x > 1e-10]
xis = [J @ Q[:, c] - 0.5 * sum((Q[:, d] @ J @ Q[:, c]) * Q[:, d] for d in range(2)) for c in range(2)]
G = E(sum(op_of(xi, D) @ Pc / math.sqrt(k) for xi, Pc, k in zip(xis, P, ks)))   # Step 4
Phi = lambda X: E(sum(K @ X @ K.conj().T for K in jumps))
Phi_dual_one = E(sum(K.conj().T @ K for K in jumps))      # Phi^*(1)
H_M = 0.5j * (G - G.conj().T)
assert np.abs(G + G.conj().T + Phi_dual_one).max() < 1e-12      # Step 5
for _ in range(3):
    X = E(cplx(D, D))
    gksl = Phi(X) - 0.5 * (Phi_dual_one @ X + X @ Phi_dual_one) - 1j * (H_M @ X - X @ H_M)
    assert np.abs(L_M(X) - gksl).max() < 1e-12            # Theorem D.7
    dual = lambda A: E(sum(K.conj().T @ A @ K for K in jumps)) + G.conj().T @ A + A @ G   # Corollary D.8
    assert abs(np.trace(L_M(X) @ (Y := E(cplx(D, D)))) - np.trace(X @ dual(Y))) < 1e-12
Pomega = np.eye(D * D) - np.outer(sum(p), sum(p)) / D     # projection onto the complement of Omega
print(f"numerical observation: min eigenvalue of J(L o E) on the complement of Omega: "
      f"{np.linalg.eigvalsh(Pomega @ J @ Pomega).min():.2f}")
assert np.linalg.eigvalsh(Pomega @ J @ Pomega).min() < -0.1
```

## D.5 Fehlerkorrektur mit Operatoralgebren

Wir verwenden den Rahmen aus Abschnitt 10.3. $\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*$ ist ein Kanal auf $M_D$, $\mathcal C=P\mathbb C^D$
ein Code, $B(\mathcal C)=\{X\in M_D:X=PXP\}$ die Operatoren auf $\mathcal C$, $\mathfrak A\subseteq B(\mathcal C)$ eine $*$-Algebra
mit Einselement $P$ (nach Satz 2.2 eine von-Neumann-Algebra auf $\mathcal C$) und $\mathfrak A'$ ihr Kommutant in
$B(\mathcal C)$. $\mathfrak A$ heißt **korrigierbar** für $\Phi$, wenn es einen Kanal $\mathcal R$ auf $M_D$ mit
$P(\mathcal R\circ\Phi)^*(A)P=A$ für alle $A\in\mathfrak A$ gibt. Hat $\mathcal R$ die Kraus-Operatoren $R_m$, so hat
$\mathcal R\circ\Phi$ die Kraus-Operatoren $R_mE_i$, und $\sum_mR_m^*R_m=1$.

**Satz D.9 (Bény–Kempf–Kribs; Satz 10.7).** $\mathfrak A$ ist genau dann korrigierbar für $\Phi$, wenn
$PE_i^*E_jP\in\mathfrak A'$ für alle $i,j$. In diesem Fall korrigiert der im Beweis konstruierte Kanal $\mathcal R$ die Algebra $\mathfrak A$.

*Beweis.* *Notwendigkeit.* $\mathcal R$ korrigiere $\mathfrak A$; man setze $W_{mi}=R_mE_iP$, so dass
$\sum_{m,i}W_{mi}^*XW_{mi}=P(\mathcal R\circ\Phi)^*(X)P$ für alle $X\in M_D$. Dies ist gleich $X$ für $X\in\mathfrak A$ und gleich $P$ für
$X=1$, weil $(\mathcal R\circ\Phi)^*$ unital ist (Korollar 6.9(b)). Sei $p\in\mathfrak A$ eine Projektion; dann ist $p\le P$
und $\sum W_{mi}^*(1-p)W_{mi}=P-p$. Für $x=px$ ist
$\sum_{m,i}\|(1-p)W_{mi}x\|^2=\langle x,(P-p)x\rangle=0$, und für $x=(P-p)x$ ist
$\sum_{m,i}\|pW_{mi}x\|^2=\langle x,px\rangle=0$. Daher ist $(1-p)W_{mi}p=0$ und $pW_{mi}(P-p)=0$, und wegen
$W_{mi}=W_{mi}P$ gilt
$$
W_{mi}p=pW_{mi}p=pW_{mi}P=pW_{mi}.
$$
Jedes selbstadjungierte $A\in\mathfrak A$ ist eine Linearkombination von Projektionen in $\mathfrak A$: Nach Satz 1.3 auf $\mathcal C$ ist
$A=\sum_k\lambda_kp_k$ mit verschiedenen $\lambda_k$ und orthogonalen Projektionen $p_k$, $\sum_kp_k=P$, und
$p_k=\prod_{l\neq k}(A-\lambda_lP)/(\lambda_k-\lambda_l)\in\mathfrak A$ (man wende beide Seiten auf die Eigenvektoren in $\mathcal C$ und auf $\mathcal C^\perp$ an; das leere Produkt ist $P$). Jedes Element
von $\mathfrak A$ ist $A_1+iA_2$ mit selbstadjungierten $A_1,A_2\in\mathfrak A$. Also ist $W_{mi}A=AW_{mi}$ für alle $A\in\mathfrak A$ und,
durch Übergang zu Adjungierten, $W_{mi}^*A=AW_{mi}^*$. Da $\sum_mW_{mi}^*W_{mj}=PE_i^*(\sum_mR_m^*R_m)E_jP=PE_i^*E_jP$, folgt
$$
A\,PE_i^*E_jP=\sum_mW_{mi}^*AW_{mj}=PE_i^*E_jP\,A\qquad(A\in\mathfrak A),
$$
d. h. $PE_i^*E_jP\in\mathfrak A'$.

*Hinlänglichkeit.* Nach Satz 2.6 für $\mathfrak A$ auf $\mathcal C$ gibt es eine Orthonormalbasis $|c;a,b\rangle$
($c=1,\dots,C$, $a=1,\dots,k_c$, $b=1,\dots,m_c$) von $\mathcal C$ – die Bilder der $e_a\otimes f_b$ im $c$-ten Summanden unter
dem unitären $V$ –, in der $\mathfrak A$ aus den Operatoren $\sum_{c,a,a',b}(A_c)_{aa'}|c;a,b\rangle\langle c;a',b|$,
$A_c\in M_{k_c}$, und $\mathfrak A'$ aus den Operatoren $\sum_{c,a,b,b'}(Y_c)_{bb'}|c;a,b\rangle\langle c;a,b'|$, $Y_c\in M_{m_c}$, besteht.
Die Voraussetzung liefert daher Matrizen $\Lambda^c_{ij}\in M_{m_c}$ mit
$$
\langle c;a,b|E_i^*E_j|d;a',b'\rangle=\delta_{cd}\,\delta_{aa'}\,(\Lambda^c_{ij})_{bb'} .
$$
*Schritt 1.* Für jedes $c$ und jedes Paar $\alpha=(i,b)$ sei $F^c_\alpha:\mathbb C^{k_c}\to\mathbb C^D$,
$F^c_\alpha e_a=E_i|c;a,b\rangle$, und $\Lambda^c_{\alpha\beta}=(\Lambda^c_{ij})_{bb'}$ für $\beta=(j,b')$. Die obige Gleichung besagt
$F^{c*}_\alpha F^d_\beta=\delta_{cd}\Lambda^c_{\alpha\beta}1_{k_c}$: eine Knill–Laflamme-Bedingung für die „Fehler“ $F^c_\alpha$ auf dem
Code $\mathbb C^{k_c}$, zusätzlich mit Orthogonalität zwischen den Sektoren.

*Schritt 2.* Wie im Beweis von Satz 10.5 ist $\Lambda^c$ selbstadjungiert, da
$\Lambda^c_{\beta\alpha}1=F^{c*}_\beta F^c_\alpha=(F^{c*}_\alpha F^c_\beta)^*=\bar\Lambda^c_{\alpha\beta}1$, und positiv, da
$\sum\bar z_\alpha z_\beta\Lambda^c_{\alpha\beta}1=(\sum z_\alpha F^c_\alpha)^*(\sum z_\beta F^c_\beta)\ge0$. Man schreibe
$\Lambda^c=U^c\operatorname{diag}(d^c_\kappa)U^{c*}$ mit unitärem $U^c$ und $d^c_\kappa\ge0$, und
$G^c_\kappa=\sum_\alpha u^c_{\alpha\kappa}F^c_\alpha$. Dann ist $G^{c*}_\kappa G^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}d^c_\kappa1$, also
$G^c_\kappa=0$ für $d^c_\kappa=0$, und andernfalls ist $V^c_\kappa=G^c_\kappa/\sqrt{d^c_\kappa}$ eine Isometrie mit
$V^{c*}_\kappa V^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}1$. Umkehrung ergibt $F^c_\alpha=\sum_\kappa\bar u^c_{\alpha\kappa}G^c_\kappa
=\sum_{\kappa:d^c_\kappa>0}\bar u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,V^c_\kappa$, also
$$
V^{c'*}_\kappa F^d_\beta=\delta_{c'd}\,\bar u^d_{\beta\kappa}\sqrt{d^d_\kappa}\,1,\qquad
F^{c*}_\alpha V^{c'}_\kappa=\delta_{cc'}\,u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,1 .
$$
*Schritt 3: die Wiederherstellung.* Die $V^c_\kappa V^{c*}_\kappa$ ($d^c_\kappa>0$) sind paarweise orthogonale Projektionen; sei $Q$ ihre
Summe. Sei $\iota_c:\mathbb C^{k_c}\to\mathbb C^D$, $\iota_ce_a=|c;a,1\rangle$, eine Isometrie, und
$$
\mathcal R(\sigma)=\sum_{c}\sum_{\kappa:d^c_\kappa>0}\iota_cV^{c*}_\kappa\,\sigma\,V^c_\kappa\iota_c^*+(1-Q)\sigma(1-Q).
$$
Dies ist ein Kanal: $\sum_{c,\kappa}V^c_\kappa\iota_c^*\iota_cV^{c*}_\kappa+(1-Q)^2=Q+1-Q=1$. Er misst den Sektor $c$ und
das Syndrom $\kappa$, macht den Fehler auf $\mathbb C^{k_c}$ rückgängig und schreibt das Ergebnis in die Kopie $b=1$.

*Schritt 4: Korrektheit.* $\mathcal R^*(X)=\sum_{c,\kappa}V^c_\kappa\iota_c^*X\iota_cV^{c*}_\kappa+(1-Q)X(1-Q)$, und für
$A\in\mathfrak A$ mit Blöcken $A_c$ ist $\iota_c^*A\iota_c=A_c$. Da jedes $F^d_\beta$ eine Kombination der $V^d_\kappa$ ist,
gilt $(1-Q)F^d_\beta=0$. Mit Schritt 2 folgt für $\alpha=(i,b)$, $\beta=(j,b')$
$$
\langle c;a,b|E_i^*\mathcal R^*(A)E_j|d;a',b'\rangle=\sum_{c',\kappa}\langle e_a,F^{c*}_\alpha V^{c'}_\kappa A_{c'}V^{c'*}_\kappa F^d_\beta e_{a'}\rangle
=\delta_{cd}\sum_\kappa u^c_{\alpha\kappa}d^c_\kappa\bar u^c_{\beta\kappa}(A_c)_{aa'}=\delta_{cd}(\Lambda^c_{ij})_{bb'}(A_c)_{aa'} .
$$
Man setze $j=i$ und summiere über $i$: Nach der Gleichung der Voraussetzung mit $\sum_iE_i^*E_i=1$ ist $\sum_i\Lambda^c_{ii}=1_{m_c}$, also
$$
\langle c;a,b|(\mathcal R\circ\Phi)^*(A)|d;a',b'\rangle=\delta_{cd}\,\delta_{bb'}(A_c)_{aa'}=\langle c;a,b|A|d;a',b'\rangle .
$$
Sowohl $P(\mathcal R\circ\Phi)^*(A)P$ als auch $A$ liegen in $B(\mathcal C)$ und haben in der Basis $|c;a,b\rangle$ dieselbe Matrix, also sind sie
gleich. $\square$

**Korollar D.10.** (a) Für $\mathfrak A=B(\mathcal C)$ ist $\mathfrak A'=\mathbb CP$, und Satz D.9 ist der Satz
10.5 von Knill–Laflamme ($C=1$, $m_1=1$; Schritt 3 ist die Wiederherstellung aus Satz 10.5). (b) Für $\mathfrak A=\operatorname{span}\{P_1,\dots,P_C\}$ mit
orthogonalen Projektionen $P_c\ne0$, $\sum_cP_c=P$, ist er die Äquivalenz von 2 und 3 in Proposition 10.6. (c) (*Rauschfreie
Teilsysteme.*) In der Notation von Proposition 10.1(a) sei $\mathcal C=V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})$ der $c$-te
Summand zur Wechselwirkungsalgebra $\mathcal A$ von $\Phi$ und $\mathfrak A=V(1_{k_c}\otimes M_{m_c})V^*$ auf $\mathcal C$, eine
Algebra der Form $M_{m}\otimes1_{k}$ nach Vertauschen der beiden Faktoren. Dann wird $\mathfrak A$ durch
$\mathcal R=\mathrm{id}$ korrigiert, und $PE_i^*E_jP=V(E_{i,c}^*E_{j,c}\otimes1_{m_c})V^*\in\mathfrak A'$ ($P=VV^*$).

*Beweis.* (a) $\mathfrak A'=\mathbb CP$ ist der Fall $C=1$, $m_1=1$ der obigen Beschreibung von $\mathfrak A'$, und
$PE_i^*E_jP\in\mathbb CP$ ist die Knill–Laflamme-Bedingung. Korrigierbarkeit von $B(\mathcal C)$ ist Korrektur im Sinne von
Abschnitt 10.2: Ist $\mathcal R(\Phi(\rho))=\rho$ für alle Zustände $\rho$ mit Träger in $\mathcal C$, so ist
$\operatorname{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)=\operatorname{Tr}(\rho A)$ für alle diese, und diese $\rho$ spannen $B(\mathcal C)$ auf; die
Umkehrung wird in Abschnitt 10.3 gezeigt. Mit $C=1$, $m_1=1$ sind die Kraus-Operatoren $\iota_1V^{1*}_\kappa$ aus Schritt 3 die
aus Satz 10.5, weil $\iota_1\iota_1^*=P$. (b) $T\in B(\mathcal C)$ kommutiert genau dann mit allen $P_c$, wenn
$P_cTP_d=0$ für $c\ne d$, und das ist Bedingung 2 aus Proposition 10.6. (c) Jedes $A\in\mathfrak A$, durch $0$ auf die
anderen Summanden fortgesetzt, liegt in $\mathcal A'=V(\bigoplus_c1_{k_c}\otimes M_{m_c})V^*$ (Satz 2.6), also ist $\Phi^*(A)=A$ nach Proposition
10.1(a), und $PAP=A$. Die Formel für $PE_i^*E_jP$ folgt aus $V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}$, und
$E_{i,c}^*E_{j,c}\otimes1_{m_c}$ kommutiert mit $1_{k_c}\otimes M_{m_c}$. $\square$

Nach Satz 2.6 ist jede $\mathfrak A$ $*$-isomorph zu einer Algebra $\bigoplus_cM_{k_c}$. Ist $C\ge2$ und $k_c\ge2$ für mindestens ein $c$, so
schützt eine korrigierbare $\mathfrak A$ die Sektorbezeichnung $c$ (klassische Information) zusammen mit der Quanteninformation
jedes $M_{k_c}$; solche Codes heißen **hybride** Codes. Als Beispiel nehmen wir $\mathfrak A=(M_2\otimes1_2)\oplus M_1$ auf einem zufälligen
fünfdimensionalen Code in $\mathbb C^8$ und ein Rauschen, das über $\sum_iB_{c,i}^*B_{c,i}=1_{m_c}$ auf die Eichfaktoren
$\mathbb C^{m_c}$ wirkt, die Sektoren auf orthogonale Unterräume abbildet, $\mathcal C^\perp$ durcheinanderwirft und dann durch eine
zufällige unitäre Matrix gemischt wird (Satz 6.7). Die Knill–Laflamme-Bedingung ist verletzt, die Bedingung von Satz D.9 erfüllt, und die
Wiederherstellung aus dem Beweis stellt jedes $A\in\mathfrak A$ wieder her, nicht aber einen generischen Operator auf $\mathcal C$:

```python
D, dims = 8, [(2, 2), (1, 1)]                             # A = (M_2 (x) 1_2) (+) (M_1 (x) 1_1), dim C = 5
U8 = np.linalg.qr(cplx(D, D))[0]
Vc = np.linalg.qr(cplx(D, 5))[0]                          # orthonormal basis |c; a, b> of the code
W = [Vc[:, :4], Vc[:, 4:]]                                # sector c: columns indexed by (a, b), row-major
Pcode = sum(Wc @ Wc.conj().T for Wc in W)
Y = [U8[:, :4], U8[:, 4:6]]                               # orthogonal images C^{k_c} (x) C^{s_c}, s_c = 2
r = 3
B = [np.linalg.qr(cplx(2 * r, m))[0].reshape(r, 2, m) for _, m in dims]   # sum_i B_i^* B_i = 1_{m_c}
E_code = [sum(Yc @ np.kron(np.eye(k), Bc[i]) @ Wc.conj().T for Yc, Bc, Wc, (k, m) in zip(Y, B, W, dims))
          for i in range(r)]
kraus = E_code + [np.linalg.qr(cplx(D, D))[0] @ (np.eye(D) - Pcode)]
mix = np.linalg.qr(cplx(r + 1, r + 1))[0]                 # Theorem 6.7: an equivalent family
kraus = [sum(mix[i, l] * K for l, K in enumerate(kraus)) for i in range(r + 1)]
assert np.abs(sum(K.conj().T @ K for K in kraus) - np.eye(D)).max() < 1e-12      # a channel

units = []                                                 # matrix units of A
for Wc, (k, m) in zip(W, dims):
    for a1 in range(k):
        for a2 in range(k):
            units.append(Wc @ np.kron(np.eye(k)[:, [a1]] @ np.eye(k)[[a2], :], np.eye(m)) @ Wc.conj().T)
comm = lambda X: max(np.abs(X @ A - A @ X).max() for A in units)
assert max(comm(Pcode @ Ki.conj().T @ Kj @ Pcode) for Ki in kraus for Kj in kraus) < 1e-12   # Theorem D.9
KL_dev = max(np.abs(Pcode @ Ki.conj().T @ Kj @ Pcode - np.trace(Pcode @ Ki.conj().T @ Kj) / 5 * Pcode).max()
             for Ki in kraus for Kj in kraus)
print(f"Knill-Laflamme condition violated by {KL_dev:.2f}")
assert KL_dev > 0.05

R_kraus, Qsum = [], np.zeros((D, D), complex)             # Steps 1-3 of the proof
for Wc, (k, m) in zip(W, dims):
    F = [Ki @ Wc[:, b::m] for Ki in kraus for b in range(m)]   # F_(i,b) e_a = E_i |c; a, b>
    Lam = np.array([[(Fa.conj().T @ Fb)[0, 0] for Fb in F] for Fa in F])
    assert max(np.abs(Fa.conj().T @ Fb - Lam[x, y] * np.eye(k)).max()
               for x, Fa in enumerate(F) for y, Fb in enumerate(F)) < 1e-12      # Step 1
    d, u = np.linalg.eigh(Lam)
    iota = Wc[:, 0::m]                                    # e_a -> |c; a, 1>
    for kk in range(len(d)):
        if d[kk] > 1e-12:
            Vk = sum(u[x, kk] * Fx for x, Fx in enumerate(F)) / math.sqrt(d[kk])
            R_kraus.append(iota @ Vk.conj().T)
            Qsum += Vk @ Vk.conj().T
R_kraus.append(np.eye(D) - Qsum)
assert np.abs(sum(R.conj().T @ R for R in R_kraus) - np.eye(D)).max() < 1e-10    # R is a channel
dual = lambda X: sum(Ki.conj().T @ R.conj().T @ X @ R @ Ki for R in R_kraus for Ki in kraus)
for _ in range(3):
    A = sum(complex(z) * Ue for z, Ue in zip(cplx(len(units)), units))
    assert np.abs(Pcode @ dual(A) @ Pcode - A).max() < 1e-10                       # A is corrected
B_gen = Pcode @ cplx(D, D) @ Pcode
print(f"a generic operator on the code: error {np.abs(Pcode @ dual(B_gen) @ Pcode - B_gen).max():.2f}")
assert np.abs(Pcode @ dual(B_gen) @ Pcode - B_gen).max() > 0.1
assert max(comm(R @ Ki @ Pcode) for R in R_kraus for Ki in kraus) < 1e-10   # proof of necessity: W_mi in A'
```

Die Ergebnisse dieses Anhangs gehen zurück auf V. Gorini, A. Kossakowski und E. C. G. Sudarshan
[*J. Math. Phys.* 17, 821 (1976)] und G. Lindblad [*Commun. Math. Phys.* 48, 119 (1976)] (Sätze D.3–D.5),
E. Christensen und D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)] (Korollar D.8 in allgemeinen von-Neumann-Algebren),
A. Frigerio (Satz D.6) sowie C. Bény, A. Kempf und D. W. Kribs [*Phys. Rev. Lett.* 98, 100502
(2007); *Phys. Rev. A* 76, 042303 (2007)] (Satz D.9); die obigen Beweise sind für die im Kurs verwendeten
endlichdimensionalen Aussagen vollständig.
