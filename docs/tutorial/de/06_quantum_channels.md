[English](../06_quantum_channels.md) · [Español](../es/06_quantum_channels.md) · [Français](../fr/06_quantum_channels.md) · **Deutsch** · [中文](../zh/06_quantum_channels.md) · [日本語](../ja/06_quantum_channels.md) · [Русский](../ru/06_quantum_channels.md)

# 6. Quantenkanäle

Ein Quantenkanal ist die allgemeinste Art, wie sich ein Zustand ändern kann: eine unitäre
Zeitentwicklung, eine Messung, deren Ergebnis vergessen wird, die Wechselwirkung mit einer Umgebung,
der Verlust eines Teilchens. Mathematisch ist er eine vollständig positive, spurerhaltende lineare
Abbildung. Dieses Kapitel beweist die Standardbeschreibungen solcher Abbildungen – Choi-Matrizen,
Kraus-Operatoren, die Stinespring-Dilatation – für Abbildungen zwischen endlichdimensionalen
von-Neumann-Algebren und untersucht dann, was an direkten Summen besonders ist: Abbildungen
zwischen Sektoren, duale Abbildungen bezüglich der Spuren aus Kapitel 3 und den klassischen Teil
eines Kanals. Die Kapitel 2 und 3 haben die Kanäle der Bibliothek bereits verwendet; hier folgt ihre
Theorie.

**Sie lernen**, warum Positivität nicht genügt, wie man zwischen Choi-Matrizen, Kraus-Operatoren,
Superoperatoren und Isometrien umrechnet, warum „spurerhaltend“ von der Spur abhängt, sobald eine
Abbildung den Sektor wechselt, welche Observablen die Dekohärenz überleben und was ein zufälliger
Kanal im Mittel mit einem reinen Zustand macht.

Im Folgenden seien $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ und
$\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)$ mit minimalen zentralen Projektionen $P_c$ und
$Q_d$; $A=\bigoplus_cA_c$; $M_{k'\times k}$ ist der Raum der komplexen $k'\times k$-Matrizen;
$E_{ij}=|i\rangle\langle j|$ sind die Matrixeinheiten; $\operatorname{Tr}$ ist die ungewichtete Spur
$\sum_c\operatorname{Tr}A_c$ und $\operatorname{tr}_w=\sum_cw_c\operatorname{Tr}A_c$ eine treue Spur
(Satz 3.1). Blockmatrizen wirken auf $\mathbb C^n\otimes\mathbb C^k$, der Blockindex zuerst.

## 6.1 Vollständig positive Abbildungen und der Satz von Choi

Eine lineare Abbildung $\Phi:\mathcal M\to\mathcal N$ heißt **positiv**, wenn aus $A\ge0$ stets
$\Phi(A)\ge0$ folgt. Die $n\times n$-Matrizen $[A_{ij}]$ mit Einträgen in $\mathcal M$ bilden wieder
eine direkte Summe von Matrixalgebren, $M_n(\mathcal M)\cong\bigoplus_cM_n(M_{k_c})$ vermöge
$[A_{ij}]\mapsto([(A_{ij})_c]_{ij})_c$; daher bedeutet $[A_{ij}]\ge0$, dass jede Blockmatrix
$[(A_{ij})_c]\in M_{nk_c}$ positiv ist. $\Phi$ heißt **$n$-positiv**, wenn
$\Phi_n:[A_{ij}]\mapsto[\Phi(A_{ij})]$ positiv ist, und **vollständig positiv** (CP), wenn es für jedes
$n$ $n$-positiv ist. Es heißt **spurerhaltend** (TP), wenn $\operatorname{Tr}\Phi(A)=\operatorname{Tr}A$
gilt, **unital**, wenn $\Phi(1)=1$, und ein **Kanal**, wenn es CP und TP ist. Wir verwenden
wiederholt, dass aus $X\ge0$ für jede rechteckige Matrix $T$ auch $TXT^*\ge0$ folgt, da
$\langle y,TXT^*y\rangle=\langle T^*y,XT^*y\rangle$.

**Proposition 6.1.** Für $k\ge2$ ist die Transposition $\theta(A)=A^T$ auf $M_k$ positiv, aber nicht
$2$-positiv.

*Beweis.* Für $A\ge0$ ist $\langle x,A^Tx\rangle=\sum_{ij}\bar x_iA_{ji}x_j=\langle\bar x,A\bar x\rangle\ge0$.
Mit $\Omega=e_1\otimes e_1+e_2\otimes e_2\in\mathbb C^2\otimes\mathbb C^k$ ist
$|\Omega\rangle\langle\Omega|=[E_{ij}]_{i,j=1}^2\ge0$, und $\theta_2$ bildet dies auf
$F=[E_{ji}]_{i,j=1}^2$ ab, das $e_1\otimes e_2-e_2\otimes e_1$ auf sein Negatives abbildet. $\square$

**Satz 6.2 (Choi).** Für eine lineare Abbildung $\Phi:M_k\to M_{k'}$ sind äquivalent:
(a) $\Phi$ ist CP; (b) $\Phi$ ist $k$-positiv; (c) die **Choi-Matrix**
$$
J(\Phi)=\sum_{i,j=1}^kE_{ij}\otimes\Phi(E_{ij})=[\Phi(E_{ij})]_{i,j=1}^k\in M_k(M_{k'})
$$
ist positiv; (d) es gibt $K_1,\dots,K_r\in M_{k'\times k}$ (**Kraus-Operatoren**) mit
$\Phi(A)=\sum_lK_lAK_l^*$, und man kann $r=\operatorname{rank}J(\Phi)\le kk'$ wählen.

*Beweis.* (a)$\Rightarrow$(b) ist trivial. (b)$\Rightarrow$(c): Für
$\Omega=\sum_ie_i\otimes e_i\in\mathbb C^k\otimes\mathbb C^k$ ist
$|\Omega\rangle\langle\Omega|=[E_{ij}]\ge0$ und $J(\Phi)=\Phi_k([E_{ij}])$. (c)$\Rightarrow$(d): Nach
Satz 1.3 ist $J(\Phi)=\sum_{l=1}^rv_lv_l^*$ mit $r=\operatorname{rank}J(\Phi)$ (die Eigenvektoren zu den
von null verschiedenen Eigenwerten, skaliert mit deren Quadratwurzeln). Jedes
$v\in\mathbb C^k\otimes\mathbb C^{k'}$ lässt sich als $v=\sum_ie_i\otimes Ke_i$ mit genau einem
$K\in M_{k'\times k}$ schreiben ($Ke_i$ ist der $i$-te Block von $v$), und der $(i,j)$-Block von
$vv^*$ ist $Ke_ie_j^*K^*=KE_{ij}K^*$. Schreibt man jedes $v_l$ so mit $K_l$ und vergleicht die
Blöcke, erhält man $\Phi(E_{ij})=\sum_lK_lE_{ij}K_l^*$ und daraus (d) wegen der Linearität.
(d)$\Rightarrow$(a): $\Phi_n([A_{ij}])=\sum_l(1_n\otimes K_l)[A_{ij}](1_n\otimes K_l)^*\ge0$. $\square$

**Korollar 6.3.** Ist $\Phi(A)=\sum_{l=1}^sL_lAL_l^*$, so gilt $s\ge\operatorname{rank}J(\Phi)$. Die
minimale Anzahl von Kraus-Operatoren, der **Kraus-Rang**, ist also $\operatorname{rank}J(\Phi)$.

*Beweis.* Liest man die letzte Rechnung rückwärts, so ist $J(\Phi)=\sum_{l=1}^sw_lw_l^*$ mit
$w_l=\sum_ie_i\otimes L_le_i$, eine Summe von $s$ Matrizen vom Rang $\le1$. $\square$

Eine Abbildung zwischen direkten Summen ist eine Matrix von Abbildungen zwischen Blöcken. Es setze
$\iota_c:M_{k_c}\to\mathcal M$ eine Matrix in den Block $c$, $\pi_d:\mathcal N\to M_{k'_d}$ greife den
Block $d$ heraus, und es sei $\Phi_{dc}=\pi_d\circ\Phi\circ\iota_c:M_{k_c}\to M_{k'_d}$. Wegen
$A=\sum_c\iota_c(A_c)$ ist $\Phi(A)_d=\sum_c\Phi_{dc}(A_c)$.

**Proposition 6.4.** $\Phi:\mathcal M\to\mathcal N$ ist genau dann $n$-positiv, wenn jedes
$\Phi_{dc}$ es ist. Folglich ist $\Phi$ genau dann CP, wenn es $K^{dc}_l\in M_{k'_d\times k_c}$,
$l=1,\dots,r$, gibt mit
$$
\Phi(A)_d=\sum_c\sum_{l=1}^rK^{dc}_l\,A_c\,K^{dc*}_l ,
$$
und $r=\max_{d,c}\operatorname{rank}J(\Phi_{dc})$ genügt. $\Phi$ ist genau dann TP, wenn
$\sum_{d,l}K^{dc*}_lK^{dc}_l=1_{k_c}$ für jedes $c$ gilt.

*Beweis.* Sind alle $\Phi_{dc}$ $n$-positiv und ist $[A_{ij}]\ge0$, so ist jedes
$X^c=[(A_{ij})_c]\ge0$, und der Block $d$ von $\Phi_n([A_{ij}])$ ist $\sum_c(\Phi_{dc})_n(X^c)\ge0$.
Umgekehrt ist für $[X_{ij}]\ge0$ in $M_n(M_{k_c})$ auch $[\iota_c(X_{ij})]\ge0$ in $M_n(\mathcal M)$
(sein Block $c$ ist $[X_{ij}]$, die übrigen verschwinden), und der Block $d$ seines Bildes ist
$[\Phi_{dc}(X_{ij})]$. Die Kraus-Form folgt aus Satz 6.2 für jedes Paar $(d,c)$, aufgefüllt mit
Nulloperatoren auf ein gemeinsames $r$. Schließlich ist wegen der Zyklizität
$\operatorname{Tr}\Phi(A)=\sum_c\operatorname{Tr}\big(A_c\sum_{d,l}K^{dc*}_lK^{dc}_l\big)$, und dies ist
genau dann für alle $A$ gleich $\sum_c\operatorname{Tr}A_c$, wenn jede dieser Summen gleich $1_{k_c}$
ist (man setze Matrixeinheiten ein: $\operatorname{Tr}(E_{ij}X)=X_{ji}$). $\square$

Die Bibliothek hat zwei Klassen. `Channel` speichert **sektorerhaltende** Abbildungen ($C'=C$, $k'=k$,
$\Phi_{dc}=0$ für $d\ne c$), deren Kraus-Operatoren $K_l=\bigoplus_cK^{cc}_l$ in $\mathcal M$ liegen,
als Tensor der Form `(batch, r, C, k_max, k_max)`; `InterSectorChannel` speichert die allgemeine
Abbildung aus Proposition 6.4 mit der Form `(batch, r, C', C, k'_max, k_max)`. `Channel.choi()` liefert
$J(\Phi_{cc})$ für jeden Sektor und `Channel.superoperator()` die Matrix $S_c$ von $\Phi_{cc}$ auf
zeilenweise angeordneten Vektoren, $\operatorname{vec}(X)_{(a,b)}=X_{ab}$. Beide enthalten die Zahlen
$\Phi(E_{ij})_{ab}=J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$ in unterschiedlicher Anordnung.
`Channel.from_superoperator` diagonalisiert die Choi-Matrizen wie im Beweis von Satz 6.2 und weist
Abbildungen zurück, die nicht CP sind.

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
```

> **Physikalische Bedeutung.** Ein System ist nie allein: Es kann mit einer unberührten Ancilla mit
> $n$ Niveaus korreliert sein. Ein physikalischer Prozess wirkt auf das Paar als
> $\Phi_n=\mathrm{id}_n\otimes\Phi$ und muss Zustände auf Zustände abbilden; also muss $\Phi$
> vollständig positiv sein. Die Transposition ist positiv, bildet aber den verschränkten Zustand
> $|\Omega\rangle\langle\Omega|/2$ auf einen Operator mit Eigenwert $-1/2$ ab. Die Choi-Matrix ist das
> Bild des maximal verschränkten Zustands: Ein einziger Zustand entscheidet über vollständige
> Positivität.

## 6.2 Stinespring-Dilatation und die Freiheit in den Kraus-Operatoren

Für $\Phi(A)=\sum_{l=1}^rK_lAK_l^*$ von $M_k$ nach $M_{k'}$ setze man $\Phi^*(B)=\sum_lK_l^*BK_l$; wegen
der Zyklizität ist $\operatorname{Tr}(\Phi(A)B)=\operatorname{Tr}(A\Phi^*(B))$, also ist $\Phi^*$ die
duale Abbildung zu $\Phi$ bezüglich der Spur (Abschnitt 6.3) und hängt nicht von der Wahl der $K_l$ ab.
Auf $\mathbb C^r\otimes\mathbb C^{k'}$ (die **Umgebung** $\mathbb C^r$ zuerst) sei
$\operatorname{Tr}_EX=\sum_l(e_l^*\otimes1)X(e_l\otimes1)$.

**Satz 6.5 (Stinespring, endlichdimensional).** Sei
$V=\sum_le_l\otimes K_l:\mathbb C^k\to\mathbb C^r\otimes\mathbb C^{k'}$. Dann gilt
$$
\Phi(A)=\operatorname{Tr}_E(VAV^*),\qquad\Phi^*(B)=V^*(1_r\otimes B)V,
$$
und $\Phi$ ist genau dann TP, wenn $V$ eine Isometrie ist. Umgekehrt definiert jedes $V$ durch die
erste Formel eine CP-Abbildung.

*Beweis.* $(e_l^*\otimes1)V=K_l$, also $\operatorname{Tr}_E(VAV^*)=\sum_lK_lAK_l^*$; umgekehrt ist
$\operatorname{Tr}_E(VAV^*)$ eine Kraus-Form mit $K_l=(e_l^*\otimes1)V$. Ferner ist
$V^*(1\otimes B)V=\sum_{l,l'}\langle e_l,e_{l'}\rangle K_l^*BK_{l'}=\Phi^*(B)$ und
$V^*V=\sum_lK_l^*K_l$, was genau dann $1$ ist, wenn $\Phi$ TP ist (Proposition 6.4). $\square$

Allgemein (hier ohne Beweis) ist eine CP-Abbildung $\Psi$ von einer unitalen C\*-Algebra
$\mathcal A$ in die beschränkten Operatoren auf einem Hilbertraum $H$ von der Form
$\Psi(a)=V^*\pi(a)V$ mit einer Darstellung $\pi$ von $\mathcal A$ auf einem Hilbertraum $K$ und einem
beschränkten $V:H\to K$ (W. F. Stinespring, Proc. AMS 6 (1955) 211; V. Paulsen, *Completely Bounded
Maps and Operator Algebras*, Cambridge 2002, Thm. 4.1); oben ist $\pi(B)=1_r\otimes B$. Für eine
Abbildung zwischen direkten Summen macht
$V_c=\sum_{d,l}e_l\otimes\iota'_dK^{dc}_l:\mathbb C^{k_c}\to\mathbb C^r\otimes\bigoplus_d\mathbb C^{k'_d}$
($\iota'_d$ die Inklusion des $d$-ten Summanden) $\Phi(A)_d$ zum $(d,d)$-Block von
$\sum_c\operatorname{Tr}_E(V_cA_cV_c^*)$, und $\Phi$ ist genau dann TP, wenn jedes $V_c$ eine
Isometrie ist; so ziehen `random_channel` und `random_inter_sector_channel` Kanäle (Abschnitt 6.5).

**Lemma 6.6.** Sind $X,Y\in M_{p\times n}$ mit $XX^*=YY^*$, so ist $Y=XW$ mit einem unitären
$W\in M_n$.

*Beweis.* $\|X^*z\|^2=\langle z,XX^*z\rangle=\|Y^*z\|^2$ für alle $z$, also ist $X^*z\mapsto Y^*z$ eine
wohldefinierte lineare Isometrie vom Bild von $X^*$ auf das Bild von $Y^*$. Wie im Beweis von Satz 1.8
lässt sie sich zu einem unitären $U$ mit $UX^*=Y^*$ fortsetzen; man nehme $W=U^*$. $\square$

**Satz 6.7 (unitäre Freiheit).** Seien $K_1,\dots,K_n$ und $L_1,\dots,L_n$ in $M_{k'\times k}$ (die
kürzere Familie wird mit Nullen aufgefüllt). Dann gilt $\sum_lK_lAK_l^*=\sum_jL_jAL_j^*$ für alle $A$
genau dann, wenn $L_j=\sum_lu_{jl}K_l$ mit einem unitären $u\in M_n$.

*Beweis.* Ist $L_j=\sum_lu_{jl}K_l$, so gilt
$\sum_jL_jAL_j^*=\sum_{l,l'}(u^*u)_{l'l}K_lAK_{l'}^*=\sum_lK_lAK_l^*$. Umgekehrt haben gleiche
Abbildungen gleiche Choi-Matrizen, $\sum_lv_lv_l^*=\sum_jw_jw_j^*$ mit $v_l=\sum_ie_i\otimes K_le_i$
und $w_j=\sum_ie_i\otimes L_je_i$. Lemma 6.6, angewandt auf die Matrizen mit den Spalten $v_l$ und
$w_j$, liefert $w_j=\sum_lW_{lj}v_l$; da $K\mapsto\sum_ie_i\otimes Ke_i$ linear und injektiv ist, folgt
$L_j=\sum_lW_{lj}K_l$, und $u=W^T$ ist unitär. $\square$

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

> **Physikalische Bedeutung.** Jeder Kanal ist eine isometrische Wechselwirkung des Systems mit einer
> Umgebung in einem festen Zustand, gefolgt vom Ignorieren der Umgebung. Eine Messung der Umgebung in
> der Basis $e_l$ liefert das Ergebnis $l$ mit Wahrscheinlichkeit $\operatorname{Tr}K_l\rho K_l^*$ und
> hinterlässt den Zustand $K_l\rho K_l^*/\operatorname{Tr}(K_l\rho K_l^*)$: Die Kraus-Operatoren sind
> die möglichen „Quantensprünge“. Eine andere Basis ergibt andere Sprünge, aber denselben Kanal
> (Satz 6.7); die kleinste Umgebung hat die Dimension $\operatorname{rank}J$.

## 6.3 Duale Abbildungen bezüglich einer Spur

**Definition.** Seien $\operatorname{tr}_w$, $\operatorname{tr}_v$ treue Spuren auf $\mathcal M$,
$\mathcal N$. Die **duale Abbildung** einer linearen Abbildung $\Phi:\mathcal M\to\mathcal N$
bezüglich $(w,v)$ ist die lineare Abbildung $\Phi^\sharp:\mathcal N\to\mathcal M$ mit
$$
\operatorname{tr}_v\big(\Phi(x)\,y\big)=\operatorname{tr}_w\big(x\,\Phi^\sharp(y)\big)\qquad(x\in\mathcal M,\ y\in\mathcal N).
$$
Sie existiert und ist eindeutig: Für festes $y$ ist $x\mapsto\operatorname{tr}_v(\Phi(x)y)$ ein
lineares Funktional, nach Satz 3.2 also gleich $\operatorname{tr}_w(\rho x)=\operatorname{tr}_w(x\rho)$
für genau ein $\rho=:\Phi^\sharp(y)$, und die Eindeutigkeit macht $\Phi^\sharp$ linear. Für die
ungewichteten Spuren schreiben wir $\Phi^*$. Seien $z_w=\sum_cw_cP_c$ und $z_v=\sum_dv_dQ_d$ (zentral
und invertierbar).

**Proposition 6.8.** Hat $\Phi$ die Kraus-Operatoren $K^{dc}_l$ (Proposition 6.4), so gilt
$$
\Phi^\sharp(y)_c=\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_l\,y_d\,K^{dc}_l,\qquad\text{i.e.}\qquad
\Phi^\sharp(y)=z_w^{-1}\,\Phi^*(z_v\,y).
$$
Insbesondere ist $\Phi^\sharp$ CP, mit den Kraus-Operatoren $(v_d/w_c)^{1/2}K^{dc*}_l$.

*Beweis.* Wegen der Zyklizität der Spur in jedem Block ist
$$
\operatorname{tr}_v(\Phi(x)y)=\sum_dv_d\sum_{c,l}\operatorname{Tr}\big(K^{dc}_lx_cK^{dc*}_ly_d\big)
=\sum_cw_c\operatorname{Tr}\Big(x_c\sum_d\frac{v_d}{w_c}\sum_lK^{dc*}_ly_dK^{dc}_l\Big),
$$
und die Eindeutigkeit liefert die Formel. Für $v=w=1$ lautet sie $\Phi^*(y)_c=\sum_{d,l}K^{dc*}_ly_dK^{dc}_l$;
die Faktoren $v_d$ und $1/w_c$ sind die Multiplikation mit $z_v$ und $z_w^{-1}$. $\square$

$\Phi$ heißt **spurerhaltend bezüglich $(w,v)$**, wenn $\operatorname{tr}_v\Phi(x)=\operatorname{tr}_w(x)$
für alle $x$ gilt.

**Korollar 6.9.** (a) $\Phi$ ist genau dann TP bezüglich $(w,v)$, wenn $\Phi^\sharp(1)=1$, d. h.
$\sum_{d,l}(v_d/w_c)K^{dc*}_lK^{dc}_l=1_{k_c}$ für jedes $c$. (b) $\Phi$ ist genau dann unital, wenn
$\Phi^\sharp$ TP bezüglich $(v,w)$ ist. (c) Die duale Abbildung zu $\Phi^\sharp$ bezüglich $(v,w)$ ist
$\Phi$. (d) Ist $\Phi:\mathcal M\to\mathcal M$ sektorerhaltend und $v=w$, so ist $\Phi^\sharp=\Phi^*$
für jedes $w$, und $\Phi$ ist genau dann für eine Spur TP, wenn es für alle TP ist.

*Beweis.* (a) $\operatorname{tr}_v\Phi(x)=\operatorname{tr}_w(x\Phi^\sharp(1))$, und
$\operatorname{tr}_w(x(\Phi^\sharp(1)-1))=0$ gilt genau dann für alle $x$, wenn $\Phi^\sharp(1)=1$ (man
nehme $x=(\Phi^\sharp(1)-1)^*$; $\operatorname{tr}_w$ ist treu); dann verwende man Proposition 6.8.
(b) $\operatorname{tr}_w\Phi^\sharp(y)=\operatorname{tr}_v(\Phi(1)y)$ ist genau dann für alle $y$ gleich
$\operatorname{tr}_v(y)$, wenn $\Phi(1)=1$. (c) $\operatorname{tr}_w(\Phi^\sharp(y)x)=\operatorname{tr}_v(y\Phi(x))$
wegen der Spureigenschaft; man verwende die Eindeutigkeit. (d) Es kommen nur die Paare $d=c$ vor, für
die $v_d/w_c=1$ ist. $\square$

Für Abbildungen, die den Sektor wechseln, gilt (d) nicht. Sei $\mathcal M=\mathcal N=M_1\oplus M_2$
(Sektor $0$: kein Teilchen; Sektor $1$: ein Teilchen auf zwei Plätzen) und $\Phi$ der
**Teilchenverlust**, $\Phi(x)=(x_0+\operatorname{Tr}x_1)\oplus0$, mit den Kraus-Operatoren
$K^{00}=1$, $K^{01}_l=\langle l|$. Er ist TP für die ungewichteten Spuren. Für `Tr_norm` ist
$w=v=(1,\tfrac12)$, $\Phi^\sharp(1)=1\oplus2\cdot1_2$ und
$\operatorname{tr}_v\Phi(x)=x_0+\operatorname{Tr}x_1\ne x_0+\tfrac12\operatorname{Tr}x_1=\operatorname{tr}_w(x)$.
In der Bibliothek ist `Channel.adjoint()` gleich $\Phi^*$ und `InterSectorChannel.adjoint(trace=...)`
gleich $\Phi^\sharp$ für dieselbe benannte Spur auf beiden Algebren (`'blunt'`, `'norm'`, `'tau_vN'`;
die Gewichte $1/(Ck_c)$ und $1/(C'k'_d)$ von `tau_vN` berücksichtigen auch die Anzahl der Sektoren).
Wir prüfen Proposition 6.8 für eine zufällige Abbildung $M_2\oplus M_3\to M_1\oplus M_2\oplus M_2$ und
die drei Spuren:

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

> **Physikalische Bedeutung.** Das spurunabhängige Objekt ist die Heisenberg-Abbildung
> $\Psi:\mathcal N\to\mathcal M$ auf den Observablen: Ein Zustand $\omega$ vor dem Prozess wird danach
> zu $\omega\circ\Psi$, und $\Psi$ ist CP und unital. Dichten hängen von der Spur ab (Satz 3.2); die
> Abbildung, die die Dichte von $\omega$ auf die von $\omega\circ\Psi$ abbildet, ist die duale
> Abbildung zu $\Psi$ bezüglich der gewählten Spuren, und sie ist nach Korollar 6.9(b) TP bezüglich
> dieser Spuren. Wirken die Algebren mit Vielfachheiten, so hat die physikalische Spur die Gewichte
> $w_c=m_c$ (Kapitel 2). Der Teilchenverlust oben erhält die Wahrscheinlichkeit bezüglich der
> ungewichteten Spur, die Zustände zählt, nicht aber bezüglich `Tr_norm`, das die Sektoren gleich
> gewichtet.

## 6.4 Dekohärenz

Die Standardkanäle der Bibliothek wirken innerhalb jedes Blocks ($p,\gamma\in[0,1]$):

| Kanal | $\Phi(\rho)_c$ | Kraus-Operatoren im Block $c$ |
|---|---|---|
| `dephasing_channel` | $(1-p)\rho_c+p\operatorname{diag}\rho_c$ | $\sqrt{1-p}\,1$, $\sqrt p\,E_{ii}$ |
| `depolarizing_channel` | $(1-p)\rho_c+p\operatorname{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1-p}\,1$, $\sqrt{p/k_c}\,E_{ij}$ |
| `center_expectation` | $\operatorname{Tr}(\rho_c)\,1/k_c$ | $\sqrt{1/k_c}\,E_{ij}$ |
| `amplitude_damping_channel` | Niveaus $i>0$ zerfallen nach $0$ | $E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}$, $\sqrt\gamma\,E_{0i}$ |
| `random_mixed_unitary_channel` | $\sum_jp_jU_j\rho U_j^*$ | $\sqrt{p_j}\,U_j$ |

Die Kraus-Formen prüft man durch direkte Rechnung, z. B.
$\sum_{ij}\frac1kE_{ij}\rho E_{ji}=\frac1k\sum_{ij}\rho_{jj}E_{ii}$. Alle fünf sind TP, und alle außer
der Amplitudendämpfung mit $\gamma>0$ sind unital. Die bedingte Erwartung $E$ auf das Zentrum
(Proposition 2.8) ist die Depolarisierung mit $p=1$. Anders als die Messung des Sektors,
$X\mapsto\sum_cP_cXP_c$, die auf $\mathcal M$ die Identität ist (Abschnitt 3.3), behält sie nur die
$p_c$.

**Satz 6.10 (Fixpunkte unitaler Kanäle).** Sei $\Phi(X)=\sum_lK_lXK_l^*$ mit $K_l\in\mathcal M$
spurerhaltend und unital. Dann gilt $\Phi(X)=X$ genau dann, wenn $X$ mit allen $K_l$ und $K_l^*$
kommutiert. Insbesondere bilden die Fixpunkte eine von-Neumann-Unteralgebra von $\mathcal M$.

*Beweis.* Kommutiert $X$ mit allen $K_l$, so ist $\Phi(X)=X\sum_lK_lK_l^*=X$. Umgekehrt ist
$\Phi(X^*)=\Phi(X)^*$, sodass die selbstadjungierten Teile $(X+X^*)/2$ und $(X-X^*)/2i$ eines
Fixpunkts Fixpunkte sind; sei also $X=X^*$. Mit $Y_l=K_lX-XK_l$ ist
$$
\sum_lY_lY_l^*=\sum_l\big(K_lX^2K_l^*-K_lXK_l^*X-XK_lXK_l^*+XK_lK_l^*X\big)=\Phi(X^2)-X^2,
$$
und die Spur davon ist $0$, weil $\Phi$ TP ist. Da $\operatorname{Tr}Y_lY_l^*=\|Y_l\|_2^2\ge0$, ist
jedes $Y_l=0$; durch Übergang zu den Adjungierten kommutiert $X$ auch mit $K_l^*$. Diese $X$ bilden
den Kommutanten einer selbstadjungierten Menge, geschnitten mit $\mathcal M$: eine $*$-Algebra, die
$1$ enthält, also eine von-Neumann-Algebra (Satz 2.2). $\square$

Also erhält die Dephasierung ($p>0$) genau die Diagonale jedes Blocks, und die Depolarisierung
($p>0$) und $E$ erhalten genau das Zentrum. Ohne Unitalität gilt der Satz nicht: Die
Amplitudendämpfung lässt $E_{00}$ fest, das nicht mit dem Sprung $E_{01}$ kommutiert (unten haben ihre
Fixpunkte und der Kommutant beide die Dimension $1$ pro Block, aber es sind $\mathbb CE_{00}$ und
$\mathbb C1$).

**Satz 6.11 (unitale Kanäle erhöhen die Entropie).** Sei $\Phi:\mathcal M\to\mathcal N$ positiv,
spurerhaltend und unital, und seien $\rho$ und $\Phi(\rho)$ Zustände mit den Eigenwerten $\lambda_i$
und $\mu_j$. Dann ist $\mu=\Delta\lambda$ mit einer doppelt stochastischen Matrix $\Delta$, also
$\sum_jf(\mu_j)\ge\sum_if(\lambda_i)$ für jede konkave Funktion $f$ auf $[0,1]$. Insbesondere gilt
$S(\Phi(\rho))\ge S(\rho)$ und $\operatorname{Tr}\Phi(\rho)^2\le\operatorname{Tr}\rho^2$.

*Beweis.* Man diagonalisiere jeden Block: $\rho=\sum_i\lambda_i|e_i\rangle\langle e_i|$ und
$\Phi(\rho)=\sum_j\mu_j|f_j\rangle\langle f_j|$ mit Orthonormalbasen von $\bigoplus_c\mathbb C^{k_c}$
und $\bigoplus_d\mathbb C^{k'_d}$, deren Vektoren jeweils in einem einzigen Summanden liegen, sodass
$|e_i\rangle\langle e_i|\in\mathcal M$. Dann ist $\mu_j=\langle f_j,\Phi(\rho)f_j\rangle=\sum_i\Delta_{ji}\lambda_i$
mit $\Delta_{ji}=\langle f_j,\Phi(|e_i\rangle\langle e_i|)f_j\rangle\ge0$,
$\sum_j\Delta_{ji}=\operatorname{Tr}\Phi(|e_i\rangle\langle e_i|)=1$ und
$\sum_i\Delta_{ji}=\langle f_j,\Phi(1)f_j\rangle=1$ (also ist $\Delta$ quadratisch: man summiere alle
Einträge). Nach der endlichen Jensenschen Ungleichung $f(\sum_i\Delta_{ji}\lambda_i)\ge\sum_i\Delta_{ji}f(\lambda_i)$
(durch Induktion nach der Anzahl der Terme aus der Definition der Konkavität, wie in der Lösung von
Aufgabe 2 in Kapitel 3) folgt
$\sum_jf(\mu_j)\ge\sum_j\sum_i\Delta_{ji}f(\lambda_i)=\sum_if(\lambda_i)$. Man wähle $f(t)=-t\log t$
und $f(t)=-t^2$. $\square$

Der klassische Teil einer Abbildung zwischen Sektoren ist ihre **Übergangsmatrix**
$T_{dc}(\rho)=\operatorname{Tr}\Phi_{dc}(\rho_c)/\operatorname{Tr}\rho_c$ ($\rho_c\ne0$); die Methode
`transition_matrix()` der Bibliothek ist $T=T(1)$, also für den maximal gemischten Zustand jedes
Eingangssektors.

**Proposition 6.12.** Sei $\Phi:\mathcal M\to\mathcal N$ positiv und TP, $\rho$ ein Zustand mit
$p_c=\operatorname{Tr}\rho_c>0$ und $q_d=\operatorname{Tr}\Phi(\rho)_d$. Dann ist $T(\rho)$
spaltenstochastisch und $q=T(\rho)p$. Mit den bedingten Erwartungen $E_{\mathcal M}$, $E_{\mathcal N}$
auf die Zentren gilt
$$
E_{\mathcal N}\Phi E_{\mathcal M}(\rho)=\bigoplus_d\frac{(Tp)_d}{k'_d}\,1_{k'_d},\qquad T=T(1).
$$

*Beweis.* $T_{dc}\ge0$, da $\Phi_{dc}$ positiv ist, und
$\sum_dT_{dc}=\operatorname{Tr}\Phi(\iota_c\rho_c)/\operatorname{Tr}\rho_c=1$. Ferner ist
$q_d=\sum_c\operatorname{Tr}\Phi_{dc}(\rho_c)=\sum_cT_{dc}(\rho)p_c$. Schließlich ist
$E_{\mathcal M}\rho=\bigoplus_cp_c1_{k_c}/k_c$, also hat
$\Phi(E_{\mathcal M}\rho)_d=\sum_cp_c\Phi_{dc}(1_{k_c}/k_c)$ die Spur $\sum_cT_{dc}p_c$, und
$E_{\mathcal N}$ ersetzt den Block $d$ durch seine Spur mal $1/k'_d$. $\square$

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

> **Physikalische Bedeutung.** Für ein Qubit multipliziert die Dephasierung die Kohärenz $\rho_{01}$
> mit $1-p$ und lässt die Besetzungen unverändert (reine Dephasierung, Zeit $T_\varphi$); die
> Amplitudendämpfung mit $1-\gamma=e^{-t/T_1}$ lässt das angeregte Niveau zerfallen und multipliziert
> die Kohärenz mit $\sqrt{1-\gamma}=e^{-t/2T_1}$, woraus $1/T_2=1/(2T_1)+1/T_\varphi$ folgt. Satz 6.10
> ist die algebraische Form dekohärenzfreier Information: Was mit allen Sprungoperatoren kommutiert,
> überlebt (Kapitel 10). Proposition 6.12 besagt, dass eine Abbildung zwischen Sektoren, durch das
> Zentrum betrachtet, eine klassische Markow-Matrix auf den Ladungen ist.

## 6.5 Zufällige Kanäle

Eine **Haar-zufällige Isometrie** $V:\mathbb C^k\to\mathbb C^N$ ist, wie in der Bibliothek
(`random_isometry`), der Faktor $Q$ von $G=QR$, wobei $G$ eine $N\times k$-Ginibre-Matrix ist
(unabhängige Einträge mit der Dichte $\pi^{-1}e^{-|z|^2}$), $Q^*Q=1$ und $R$ eine obere
Dreiecksmatrix mit positiver Diagonale. `random_channel(alg, r)` zieht solche
$V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}$ unabhängig in jedem Block und verwendet die
Kraus-Operatoren aus Satz 6.5.

**Lemma 6.13.** $V$ hat die Verteilung der ersten $k$ Spalten einer Haar-verteilten unitären Matrix
aus $U(N)$. Insbesondere hat $UV$ für jedes unitäre $U\in M_N$ dieselbe Verteilung wie $V$.

*Beweis.* Sei $Z$ eine $N\times N$-Ginibre-Matrix mit den ersten $k$ Spalten $G$ und
$E\in M_{N\times k}$ die ersten $k$ Spalten von $1_N$. Nach Satz 5.3 ist $Z$ fast sicher invertierbar
und $Z=Q_ZR_Z$ mit Haar-verteiltem $Q_Z$. Da $R_Z$ eine obere Dreiecksmatrix ist, ist
$G=ZE=(Q_ZE)(E^*R_ZE)$ eine Faktorisierung der verlangten Art, und zwar die einzige: Ist
$G=QR=Q'R'$, so gilt $R^*R=G^*G=R'^*R'$, also ist $R'R^{-1}$ unitär und eine obere Dreiecksmatrix mit
positiver Diagonale, folglich gleich $1$ (Beweis von Lemma 5.2). Also ist $V=Q_ZE$, und $UV=(UQ_Z)E$
hat die Verteilung von $V$, weil $UQ_Z$ die von $Q_Z$ hat. $\square$

**Lemma 6.14.** Sei $u$ ein zufälliger Einheitsvektor in $\mathbb C^N$, dessen Verteilung unter allen
unitären Matrizen invariant ist. Dann ist $\mathbb E[u_a\bar u_b]=\delta_{ab}/N$ und
$$
\mathbb E\big[u_a\bar u_b\,u_c\bar u_d\big]=\frac{\delta_{ab}\delta_{cd}+\delta_{ad}\delta_{cb}}{N(N+1)} .
$$

*Beweis.* Die Invarianz unter $\operatorname{diag}(e^{i\theta_1},\dots,e^{i\theta_N})$ multipliziert
die Momente mit $e^{i(\theta_a-\theta_b)}$ bzw. $e^{i(\theta_a-\theta_b+\theta_c-\theta_d)}$; sie
verschwinden also, außer wenn $a=b$ bzw. $\{a,c\}=\{b,d\}$ als Multimengen. Wegen der Invarianz unter
Permutationsmatrizen ist $\mathbb E|u_a|^2=1/N$, und es bleiben nur $\alpha=\mathbb E|u_1|^4$ und
$\beta=\mathbb E|u_1|^2|u_2|^2$. Erstens liefert $\sum_c|u_c|^2=1$ die Gleichung
$\alpha+(N-1)\beta=\mathbb E|u_1|^2=1/N$. Zweitens ist $(u_1+u_2)/\sqrt2$ die erste Koordinate von
$Uu$ für ein unitäres $U$, also $\mathbb E|u_1+u_2|^4=4\alpha$. Mit $s=|u_1|^2+|u_2|^2$ und
$z=u_1\bar u_2$ ist $|u_1+u_2|^4=s^2+2s(z+\bar z)+(z+\bar z)^2$; die Invarianzen $u_1\mapsto-u_1$
und $u_1\mapsto iu_1$ liefern $\mathbb E[s(z+\bar z)]=0$ und $\mathbb Ez^2=0$, also
$4\alpha=(2\alpha+2\beta)+2\beta$, d. h. $\alpha=2\beta$. Daher ist $\beta=1/(N(N+1))$ und
$\alpha=2\beta$, und das ist die Formel in allen Fällen. (Hier wurde $N\ge2$ benutzt; für $N=1$ ist
$|u_1|=1$, und beide Seiten sind gleich $1$.) $\square$

**Satz 6.15 (Ausgangsreinheit eines zufälligen Kanals).** Sei $\psi$ ein reiner Zustand in einem
Sektor der Größe $k$ und $\Phi$ von `random_channel` mit Kraus-Rang $r$ gezogen. Dann gilt
$$
\mathbb E\operatorname{Tr}\Phi(\psi)^2=\frac{k+r}{kr+1},\qquad\mathbb E\,\Phi(\psi)=\frac1k1_k .
$$

*Beweis.* Im Sektor von $\psi$ ist $\Phi(\psi)=\operatorname{Tr}_E|u\rangle\langle u|$ mit $u=V\psi$,
einem Einheitsvektor, dessen Verteilung nach Lemma 6.13 unitär invariant ist. Man indiziere
$\mathbb C^r\otimes\mathbb C^k$ durch Paare $(e,x)$. Dann ist $\Phi(\psi)_{xy}=\sum_eu_{ex}\bar u_{ey}$
mit dem Mittelwert $r\delta_{xy}/(kr)$, und
$\operatorname{Tr}\Phi(\psi)^2=\sum_{x,y,e,f}u_{ex}\bar u_{ey}u_{fy}\bar u_{fx}$. Nach Lemma 6.14 mit
$N=kr$ erzwingt der erste Term $x=y$ ($kr^2$ Indextupel) und der zweite $e=f$ ($k^2r$ Tupel), also
ist der Mittelwert $(kr^2+k^2r)/(kr(kr+1))$. $\square$

Für $r=1$ ist der Kanal unitär und die Reinheit gleich $1$; für $r\to\infty$ strebt sie gegen $1/k$,
die Reinheit des maximal gemischten Zustands. Der Beweis hat nur die unitäre Invarianz der
Verteilung von $u$ benutzt; er liefert daher auch Lubkins mittlere Reinheit des reduzierten Zustands
eines zufälligen reinen Zustands $u=g/\|g\|$ mit Gaußschem $g$, die Kapitel 7, Satz 7.11, mit
Gaußschen Integralen beweist; die Zahlen $\alpha,\beta$ aus Lemma 6.14 sind die aus Korollar 5.11
für die Einträge einer Haar-verteilten unitären Matrix.

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

> **Physikalische Bedeutung.** Ein zufälliger Kanal modelliert generisches, unstrukturiertes Rauschen
> aus einer Umgebung der Dimension $r$. Schon $r=4$ drückt die Reinheit eines reinen Qutrit-Zustands
> im Mittel auf $7/13\approx0.54$, verglichen mit $1$ ohne Rauschen und dem Minimum $1/3$. In einem
> System mit Superauswahl wirkt das Rauschen in jedem Sektor getrennt, und maßgeblich ist die
> Dimension des Sektors, nicht die des ganzen Raumes.

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

## Aufgaben

1. (Werner–Holevo) Man zeige, dass $\Phi(\rho)=\big(\operatorname{Tr}(\rho)1-\rho^T\big)/(k-1)$ auf
   $M_k$, $k\ge2$, ein unitaler Kanal mit Kraus-Rang $k(k-1)/2$ ist, obwohl er aus der Transposition
   gebildet ist. Man prüfe es mit `Channel.from_superoperator`.
2. Man zeige, dass der Dephasierungskanal auf $M_k$ für $0<p\le1$ den Kraus-Rang $k$ hat, einen
   weniger als die Anzahl der Kraus-Operatoren, die die Bibliothek speichert.
3. Sei $X$ die zyklische Verschiebung auf $\mathbb C^k$ und $Z=\operatorname{diag}(1,\omega,\dots,\omega^{k-1})$,
   $\omega=e^{2\pi i/k}$. Man zeige $\frac1{k^2}\sum_{a,b=0}^{k-1}X^aZ^b\rho Z^{-b}X^{-a}=\operatorname{Tr}(\rho)1/k$:
   Die bedingte Erwartung auf das Zentrum von $M_k$ ist ein gemischt-unitärer Kanal.
4. Man zeige $(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp$, wenn auf der mittleren Algebra
   dieselbe Spur verwendet wird, und prüfe es für die drei Spuren der Bibliothek.

<details>
<summary>Lösungen</summary>

1. Mit dem Swap-Operator $F$ ist $J(\theta)=\sum_{ij}E_{ij}\otimes E_{ji}=F$ und $J(\rho\mapsto\operatorname{Tr}(\rho)1)=1\otimes1$,
   also $J(\Phi)=(1-F)/(k-1)=2P_-/(k-1)$, wobei $P_-=(1-F)/2$ die Projektion auf die
   antisymmetrischen Vektoren ist ($F=F^*$, $F^2=1$). Also ist $J\ge0$ mit Rang $k(k-1)/2$ (Satz 6.2,
   Korollar 6.3), und $\operatorname{Tr}\Phi(\rho)=\operatorname{Tr}\rho$, $\Phi(1)=1$ sieht man direkt.

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

2. $J=(1-p)|\Omega\rangle\langle\Omega|+p\sum_iE_{ii}\otimes E_{ii}$ verschwindet außerhalb des
   Aufspanns der $e_i\otimes e_i$, und auf diesem $k$-dimensionalen Raum ist es die Matrix
   $(1-p)\mathbf 1\mathbf 1^T+p1$ ($\mathbf 1$ der Vektor aus Einsen) mit den Eigenwerten $p$ und
   $(1-p)k+p$, die alle positiv sind. Also ist $\operatorname{rank}J=k$, während die Bibliothek
   $1+k_{\max}$ Operatoren speichert.

```python
ranks = [int((torch.linalg.eigvalsh(J_c[0]) > 1e-10).sum()) for J_c in zoo['dephasing'].choi()]
print(ranks, zoo['dephasing'].kraus_rank)                 # [2, 3] versus 4 stored operators
assert ranks == alg.k_factors
```

3. $\operatorname{Tr}((X^aZ^b)^*X^{a'}Z^{b'})=k\,\delta_{aa'}\delta_{bb'}$: $Z^{-b}X^{a'-a}Z^{b'}$ hat
   verschwindende Diagonale, außer wenn $a=a'$, und $\operatorname{Tr}Z^{b'-b}=\sum_j\omega^{j(b'-b)}=k\delta_{bb'}$.
   Also bilden die $B_{ab}=X^aZ^b/\sqrt k$ eine Orthonormalbasis von $M_k$ bezüglich
   $\operatorname{Tr}(A^*B)$, und je zwei Orthonormalbasen gehen durch eine unitäre Matrix
   auseinander hervor; nach Satz 6.7 ist $\sum_mB_m\rho B_m^*$ für alle Orthonormalbasen gleich. Für
   die Basis $E_{ab}$ ist es $\sum_{ab}E_{ab}\rho E_{ba}=\operatorname{Tr}(\rho)1$. Man teile durch $k$.

```python
Xs = torch.roll(torch.eye(k, dtype=dt), 1, dims=0)
Zs = torch.diag(torch.exp(2j * math.pi * torch.arange(k, dtype=torch.float64) / k)).to(dt)
W = [torch.linalg.matrix_power(Xs, a) @ torch.linalg.matrix_power(Zs, b) / k for a in range(k) for b in range(k)]
twirl = channels.kraus_channel([M3.operator(Wab[None, None]) for Wab in W])
r3op = states.random_density_matrix(M3)
assert close(twirl(r3op), channels.center_expectation(M3)(r3op))
```

4. $\operatorname{tr}(\Psi(\Phi(x))z)=\operatorname{tr}(\Phi(x)\Psi^\sharp(z))=\operatorname{tr}(x\,\Phi^\sharp(\Psi^\sharp(z)))$
   für alle $x$ und $z$; die Eindeutigkeit der dualen Abbildung liefert die Behauptung.

```python
Psi2 = channels.random_inter_sector_channel(alg_out, alg, kraus_rank=2)
for t in trace:
    assert close((Psi2 @ Psi).adjoint(trace=t)(x), Psi.adjoint(trace=t)(Psi2.adjoint(trace=t)(x)))
```

</details>

**Weiter:** [Kapitel 7](07_tensor_products.md) behandelt Tensorprodukte und Verschränkung: die partielle
Spur – einen Kanal zwischen verschiedenen Algebren, dessen duale Abbildung die Einbettung
$A\mapsto A\otimes1$ ist –, die Page-Formel und die Fusion von Ladungen.
