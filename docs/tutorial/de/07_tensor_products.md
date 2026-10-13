[English](../07_tensor_products.md) · [Español](../es/07_tensor_products.md) · [Français](../fr/07_tensor_products.md) · **Deutsch** · [中文](../zh/07_tensor_products.md) · [日本語](../ja/07_tensor_products.md) · [Русский](../ru/07_tensor_products.md)

# 7. Tensorprodukte und Verschränkung

Zwei Quantensysteme zusammen werden durch das Tensorprodukt ihrer Algebren beschrieben. Mit
Sektoren gibt es zwei Produkte: das einfache, in dem jeder Teil seine eigene Ladung behält, und
das *fusionierte*, in dem nur die Gesamtladung erhalten ist. Dieses Kapitel beweist die Struktur
beider, definiert die partielle Spur als duale Abbildung einer Einbettung und untersucht die
Verschränkung.

**Sie lernen,** wie sich Sektoren multiplizieren und wie die Bibliothek sie anordnet, warum die
partielle Spur ein Kanal und, normiert, eine bedingte Erwartung von endlichem Index ist, wie der
reduzierte Zustand eines ladungserhaltenden Zustands aussieht, wie die Verschränkung in einen
Teilchenzahlanteil und einen Konfigurationsanteil zerfällt, und wie groß die mittlere Reinheit
und die mittlere Entropie zufälliger Zustände sind.

Im Folgenden wirke $\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)$ auf
$H_1=\bigoplus_c\mathbb C^{k_c}$ mit minimalen zentralen Projektionen $P_c$, und
$\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)$ wirke auf $H_2=\bigoplus_d\mathbb C^{l_d}$
mit minimalen zentralen Projektionen $Q_d$. $\operatorname{Tr}$ ist die Spur des Raumes, auf dem ein
Operator wirkt; auf $\mathcal M_i$ ist sie die ungewichtete Spur aus Kapitel 3. Für $X$ auf
$\mathbb C^k\otimes\mathbb C^l$ schreiben wir $X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle$
und identifizieren $\mathbb C^k\otimes\mathbb C^l$ mit $\mathbb C^{kl}$ durch $e_a\otimes e_b\mapsto e_{(a-1)l+b}$
(die Kronecker-Konvention), sodass $(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}$ gilt.

## 7.1 Das Tensorprodukt zweier Algebren

Für $*$-Algebren $\mathcal M_1$ auf $H_1$ und $\mathcal M_2$ auf $H_2$ ist das **Tensorprodukt**
$\mathcal M_1\otimes\mathcal M_2$ die lineare Hülle der Operatoren $A\otimes B$ auf
$H_1\otimes H_2$ mit $A\in\mathcal M_1$, $B\in\mathcal M_2$. Wegen
$(A\otimes B)(A'\otimes B')=AA'\otimes BB'$, $(A\otimes B)^*=A^*\otimes B^*$ und $1\otimes1=1$ ist es
eine $*$-Algebra und damit eine von-Neumann-Algebra (Satz 2.2).

**Satz 7.1.** (a) $H_1\otimes H_2=\bigoplus_{c,d}\mathbb C^{k_c}\otimes\mathbb C^{l_d}$ ist die
orthogonale Summe der Bilder der Projektionen $P_c\otimes Q_d$, und
$$
\mathcal M_1\otimes\mathcal M_2=\bigoplus_{c,d}M_{k_c}\otimes M_{l_d}\cong\bigoplus_{c,d}M_{k_cl_d}(\mathbb C):
$$
es besteht aus allen Operatoren $X=\sum_{c,d}(P_c\otimes Q_d)X(P_c\otimes Q_d)$ mit beliebigen Blöcken
auf den $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$. (b) Seine minimalen zentralen Projektionen sind die
$C_1C_2$ Projektionen $P_c\otimes Q_d$, und es gilt
$\dim(\mathcal M_1\otimes\mathcal M_2)=\dim\mathcal M_1\dim\mathcal M_2$.
(c) Für Spuren $\operatorname{tr}_w$ auf $\mathcal M_1$ und $\operatorname{tr}_v$ auf $\mathcal M_2$ (Satz 3.1) gilt
$\operatorname{tr}_w(A)\operatorname{tr}_v(B)=\operatorname{tr}_{w\otimes v}(A\otimes B)$ mit
$(w\otimes v)_{cd}=w_cv_d$. Insbesondere sind `Tr_blunt`, `Tr_norm` und `tau_vN` des Produkts die
Produkte derselben Spuren der Faktoren.

*Beweis.* (a) Die $P_c\otimes Q_d$ sind paarweise orthogonale Projektionen mit Summe $1\otimes1$,
und das Bild von $P_c\otimes Q_d$ ist $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$. Da $P_c$ und $Q_d$ zentral
sind, gilt $(P_c\otimes Q_d)(A\otimes B)=AP_c\otimes BQ_d=(A\otimes B)(P_c\otimes Q_d)$, also ist jedes
Element von $\mathcal M_1\otimes\mathcal M_2$ blockdiagonal. Umgekehrt ist mit den Matrixeinheiten
$E^c_{aa'}\in\mathcal M_1$ und $F^d_{bb'}\in\mathcal M_2$ der Blöcke $c$ und $d$
$E^c_{aa'}\otimes F^d_{bb'}=|e_a\otimes e_b\rangle\langle e_{a'}\otimes e_{b'}|$ innerhalb von
$\mathbb C^{k_c}\otimes\mathbb C^{l_d}$; dies sind die Matrixeinheiten von $M_{k_cl_d}(\mathbb C)$, also
liegt jeder im Block $(c,d)$ getragene Operator in $\mathcal M_1\otimes\mathcal M_2$. (b) Das Zentrum
einer direkten Summe voller Matrixalgebren wird von den Projektionen auf die Summanden aufgespannt
(Kapitel 2, Aufgabe 1), und $\sum_{c,d}k_c^2l_d^2=\sum_ck_c^2\sum_dl_d^2$. (c)
$\operatorname{Tr}(X\otimes Y)=\sum_{a,b}X_{aa}Y_{bb}=\operatorname{Tr}X\operatorname{Tr}Y$, also
$\operatorname{tr}_{w\otimes v}(A\otimes B)=\sum_{c,d}w_cv_d\operatorname{Tr}A_c\operatorname{Tr}B_d=\operatorname{tr}_w(A)\operatorname{tr}_v(B)$.
Die Gewichte $1$, $1/k_c$, $1/(C_1k_c)$ der drei Spuren (Abschnitt 3.1) multiplizieren sich zu $1$,
$1/(k_cl_d)$, $1/(C_1C_2k_cl_d)$, den Gewichten derselben Spuren auf dem Produkt. $\square$

`tensor_product(alg1, alg2)` konstruiert diese Algebra. Der Sektor $(c,d)$ hat den Index
$s=cC_2+d$ (im Code ab $0$ gezählt) und die Ladung $(q_c,r_d)$, und innerhalb des Sektors hat
$e_a\otimes e_b$ den Index $a\,l_d+b$. `kron(A, B, alg12)` berechnet $A\otimes B$. Die
Sektorreihenfolge unterscheidet sich von der Kronecker-Reihenfolge von $H_1\otimes H_2$ um eine
Permutation der Basis, die `kron_rows` unten explizit konstruiert; mit ihr vergleichen wir mit
Kronecker-Produkten voller Matrizen.

```python
import math
import numpy as np
import torch
from torch_vn_algebra import (TypeIAlgebra, DensityMatrix, states, composite,
                              tensor_product, fused_tensor_product, kron, partial_trace)

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def full(X):
    """An operator of (+)_c M_{k_c} as a block-diagonal matrix on (+)_c C^{k_c}, shape (batch, D, D)."""
    return torch.stack([torch.block_diag(*[x[c, :k, :k] for c, k in enumerate(X.algebra.k_factors)])
                        for x in X.matrix])

def kron_rows(alg12):
    """For every sector of a product algebra, the positions of its basis vectors in H1 (x) H2 (Kronecker order)."""
    k1, k2 = alg12.tensor_factors[0].k_factors, alg12.tensor_factors[1].k_factors
    return [torch.tensor([(sum(k1[:c]) + a) * sum(k2) + sum(k2[:d]) + b
                          for c, d, off in pairs for a in range(k1[c]) for b in range(k2[d])])
            for pairs in alg12.product_blocks]              # pair blocks (c, d, offset), by offset

def full12(X):
    """An operator of a (plain or fused) product algebra as a matrix on H1 (x) H2 in Kronecker order."""
    f1, f2 = X.algebra.tensor_factors
    N = sum(f1.k_factors) * sum(f2.k_factors)
    out = torch.zeros(X.matrix.shape[0], N, N, dtype=X.matrix.dtype)
    for s, r in enumerate(kron_rows(X.algebra)):
        out[:, r[:, None], r[None, :]] = X.matrix[:, s, :len(r), :len(r)]
    return out

def random_operator(alg):
    return alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])

alg1 = TypeIAlgebra([1, 2], [1, 2], precision='double')            # C (+) M_2
alg2 = TypeIAlgebra([2, 3], [2, 3], precision='double')            # M_2 (+) M_3
alg12 = tensor_product(alg1, alg2)
print(alg12.k_factors, alg12.charges, sum(k * k for k in alg12.k_factors))     # 4 sectors, dimension 65
assert sum(k * k for k in alg12.k_factors) == 5 * 13                           # Theorem 7.1(b)

A, A2, B, B2 = random_operator(alg1), random_operator(alg1), random_operator(alg2), random_operator(alg2)
AB = kron(A, B, alg12)
assert close(full12(AB)[0], torch.kron(full(A)[0], full(B)[0]))               # Theorem 7.1(a)
assert close(AB @ kron(A2, B2, alg12), kron(A @ A2, B @ B2, alg12))
for tr in ('Tr_blunt', 'Tr_norm', 'tau_vN'):                                   # Theorem 7.1(c)
    assert (getattr(AB, tr)() - getattr(A, tr)() * getattr(B, tr)()).abs().max() < 1e-10
```

> **Physikalische Bedeutung.** Gehorcht jedes von zwei Systemen seiner eigenen Superauswahlregel –
> zwei Kästen, die keine Teilchen austauschen können –, so hat das Gesamtsystem einen Sektor für
> jedes *Paar* von Ladungen, und keine Observable verbindet verschiedene Paare. Ein gemeinsamer
> Sektor zählt mit dem Produkt der Gewichte.

## 7.2 Die partielle Spur

**Definition.** Sei $\rho\in\mathcal M_1\otimes\mathcal M_2$. Die **partielle Spur**
$\operatorname{Tr}_2\rho\in\mathcal M_1$ ist das Element mit
$$
\operatorname{Tr}\big((\operatorname{Tr}_2\rho)\,A\big)=\operatorname{Tr}\big(\rho\,(A\otimes1)\big)\qquad\text{for all }A\in\mathcal M_1 .
$$
Sie existiert und ist eindeutig nach Satz 3.2, angewandt auf das lineare Funktional
$A\mapsto\operatorname{Tr}(\rho(A\otimes1))$ und die treue Spur $\operatorname{Tr}$ von $\mathcal M_1$: Die
partielle Spur ist die duale Abbildung der Einbettung $A\mapsto A\otimes1$ bezüglich der
ungewichteten Spuren. $\operatorname{Tr}_1$ wird ebenso definiert. Für die vollen Algebren
$\mathcal M_i=B(H_i)$ ist dies die **partielle Spur des Hilbertraums**
$\operatorname{Tr}_{H_2}:B(H_1\otimes H_2)\to B(H_1)$.

**Proposition 7.2.** (a) Für $\rho=\bigoplus_{c,d}\rho_{cd}$ gilt
$$
(\operatorname{Tr}_2\rho)_c=\sum_d\operatorname{tr}_{l_d}\rho_{cd},\qquad
(\operatorname{tr}_lX)_{aa'}=\sum_{b=1}^lX_{(a,b),(a',b)}\quad(X\text{ on }\mathbb C^k\otimes\mathbb C^l),
$$
und $\operatorname{Tr}_2\rho=\operatorname{Tr}_{H_2}\rho$, wobei $(\operatorname{Tr}_{H_2}X)_{ii'}=\sum_jX_{(i,j),(i',j)}$.
(b) $\operatorname{Tr}_2(A\otimes B)=(\operatorname{Tr}B)\,A$.
(c) $\operatorname{Tr}_2\big((A\otimes1)X(A'\otimes1)\big)=A\,(\operatorname{Tr}_2X)\,A'$ für $A,A'\in\mathcal M_1$.
(d) $\operatorname{Tr}(\operatorname{Tr}_2X)=\operatorname{Tr}X$, und $\operatorname{Tr}_2$ bildet Zustände auf Zustände ab.

*Beweis.* (a) Mit $A=\bigoplus_cA_c$ ist $A\otimes1=\bigoplus_{c,d}A_c\otimes1_{l_d}$ und
$$
\operatorname{Tr}\big(\rho(A\otimes1)\big)=\sum_{c,d}\sum_{a,a',b}(\rho_{cd})_{(a,b),(a',b)}(A_c)_{a'a}
=\sum_c\operatorname{Tr}\Big(\Big(\sum_d\operatorname{tr}_{l_d}\rho_{cd}\Big)A_c\Big),
$$
also erfüllt das Element auf der rechten Seite die Definition und ist gleich $\operatorname{Tr}_2\rho$.
Dieselbe Rechnung für $B(H_i)$ liefert die Formel für $\operatorname{Tr}_{H_2}$; für blockdiagonales
$\rho$ verschwinden die Einträge $\rho_{(i,j),(i',j)}$ mit $i,i'$ in verschiedenen Blöcken, also hat
$\operatorname{Tr}_{H_2}\rho$ die Blöcke aus (a). (b) $\operatorname{tr}_l(X\otimes Y)=(\operatorname{Tr}Y)X$. (c) Für
$B\in\mathcal M_1$ ist
$\operatorname{Tr}\big(\operatorname{Tr}_2((A\otimes1)X(A'\otimes1))B\big)=\operatorname{Tr}\big(X(A'BA\otimes1)\big)
=\operatorname{Tr}\big(A(\operatorname{Tr}_2X)A'B\big)$; man verwende die Eindeutigkeit. (d) Man nehme $A=1$. Ist
$\rho$ ein Zustand, so ist $A\mapsto\operatorname{Tr}(\rho(A\otimes1))$ positiv, da
$A^*A\otimes1=(A\otimes1)^*(A\otimes1)$, und hat bei $A=1$ den Wert $1$; also ist seine Dichte
$\operatorname{Tr}_2\rho$ ein Zustand (Satz 3.2). $\square$

Kanäle sind das Thema von Kapitel 6; hier brauchen wir nur diese Definition: Eine lineare
Abbildung $\Phi$ heißt **vollständig positiv**, wenn $\mathrm{id}_n\otimes\Phi$ für jedes $n$ auf
$M_n(\mathbb C)\otimes(\cdot)$ positiv ist, und ein **Kanal**, wenn sie außerdem spurerhaltend ist.

**Satz 7.3.** Es laufe $j$ über die Standardbasis von $H_2$, und es sei
$K_j=1\otimes\langle j|:H_1\otimes H_2\to H_1$, $K_j(x\otimes y)=\langle j,y\rangle x$. Dann gilt
$$
\operatorname{Tr}_{H_2}X=\sum_jK_jXK_j^*,\qquad\sum_jK_j^*K_j=1,\qquad\sum_jK_j^*AK_j=A\otimes1 .
$$
Folglich ist $\operatorname{Tr}_2$ ein Kanal von $\mathcal M_1\otimes\mathcal M_2$ nach $\mathcal M_1$, und
seine duale Abbildung bezüglich der ungewichteten Spuren ist $A\mapsto A\otimes1$.

*Beweis.* $K_j^*x=x\otimes e_j$, also $\langle e_a,K_jXK_j^*e_{a'}\rangle=X_{(a,j),(a',j)}$; Summation
über $j$ liefert Proposition 7.2(a). Ferner ist
$\sum_jK_j^*K_j(x\otimes y)=\sum_j\langle j,y\rangle x\otimes e_j=x\otimes y$
und $\sum_jK_j^*AK_j=\sum_jA\otimes|j\rangle\langle j|=A\otimes1$. Für $Y\ge0$ ist
$(\mathrm{id}_n\otimes\operatorname{Tr}_{H_2})(Y)=\sum_j(1_n\otimes K_j)Y(1_n\otimes K_j)^*\ge0$, denn
$\langle x,TYT^*x\rangle=\langle T^*x,YT^*x\rangle\ge0$; dies bleibt auf der Unteralgebra
$\mathcal M_1\otimes\mathcal M_2$ gültig. Schließlich ist
$\operatorname{Tr}\sum_jK_jXK_j^*A=\operatorname{Tr}X\sum_jK_j^*AK_j$, was für $A=1$ die Spurerhaltung und im
Allgemeinen die Dualität ist. $\square$

`partial_trace(X, keep=1)` formt jeden Block $\rho_{cd}$ in die Gestalt $(k_c,l_d,k_c,l_d)$ um und
bildet die Spur; `composite.partial_trace_channel(alg12, keep=1)` ist der Kanal aus Satz 7.3 mit
den Kraus-Operatoren $1_{k_c}\otimes\langle j|$ vom Sektor $(c,d)$ in den Sektor $c$. Wir prüfen beide
an der Definition, an einem unabhängigen `numpy.einsum` und an der Dualität:

```python
rho = states.random_density_matrix(alg12, batch_size=3)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
assert isinstance(rho1, DensityMatrix)                                      # states give states
R = full12(rho).numpy().reshape(3, 3, 5, 3, 5)                              # indices (batch, i, j, i', j')
assert np.abs(full(rho1).numpy() - np.einsum('xijkj->xik', R)).max() < 1e-12
assert np.abs(full(rho2).numpy() - np.einsum('xijil->xjl', R)).max() < 1e-12
one2 = alg2.identity()
assert (states.expectation(rho, kron(A, one2, alg12)) - states.expectation(rho1, A)).abs().max() < 1e-12
assert close(partial_trace(AB, keep=1), A * B.trace.item())                 # Proposition 7.2(b)

Ptr = composite.partial_trace_channel(alg12, keep=1)
print(Ptr)
assert Ptr.is_trace_preserving()
assert close(Ptr(rho), rho1)
assert close(Ptr.adjoint('blunt')(A), kron(A, one2, alg12))                 # Theorem 7.3: the dual
```

> **Physikalische Bedeutung.** $\operatorname{Tr}_2\rho$ ist der **reduzierte Zustand**: alles, was ein
> Beobachter erfährt, der nur das erste System messen kann – die Definition besagt, dass er jeden
> Erwartungswert einer Observablen $A\otimes1$ reproduziert. Gemeinsame Zustände mit gleichen
> reduzierten Zuständen unterscheiden sich in ihren Korrelationen.

**Die partielle Spur als bedingte Erwartung.** Sei $\mathcal M_2=M_m(\mathbb C)$, sodass die Sektoren
von $\mathcal M_1\otimes M_m$ die von $\mathcal M_1$ sind (Satz 7.1), sei $\iota(A)=A\otimes1_m$ und
$\mathcal E=\frac1m\operatorname{Tr}_2:\mathcal M_1\otimes M_m\to\mathcal M_1$, und seien $e_{ab}$ die
Matrixeinheiten von $M_m$.

**Proposition 7.4.** (a) $\iota\circ\mathcal E$ ist eine bedingte Erwartung von
$\mathcal M_1\otimes M_m$ auf $\mathcal M_1\otimes1$ (Abschnitt 2.6), und für jede treue Spur
$\operatorname{tr}_w$ von $\mathcal M_1\otimes M_m$ ist sie die einzige, die $\operatorname{tr}_w$ erhält. (b) Die
$m^2$ Elemente $u_{ab}=\sqrt m\,(1\otimes e_{ab})$ erfüllen
$X=\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*X))$ für alle $X$, und es gilt
$\sum_{a,b}u_{ab}u_{ab}^*=m^2\,1$. (c) Jede Familie $(v_k)$ mit
$X=\sum_kv_k\,\iota(\mathcal E(v_k^*X))$ für alle $X$ erfüllt $\sum_kv_kv_k^*=m^2\,1$.

*Beweis.* (a) $\mathcal E$ ist positiv (Satz 7.3), ebenso $\iota$
($T^*T\otimes1=(T\otimes1)^*(T\otimes1)$); nach Proposition 7.2(b) gilt $\mathcal E(A\otimes1)=A$, und die
Bimoduleigenschaft ist Proposition 7.2(c). Wir schreiben $\operatorname{tr}_w$ auch für die Spur mit
denselben Gewichten $w_c$ auf $\mathcal M_1$; dann ist $\operatorname{tr}_w(Y\otimes1)=m\operatorname{tr}_w(Y)$. Die
Definition von $\operatorname{Tr}_2$, angewandt auf $AP_c$, mit $w_c$ multipliziert und über $c$ summiert,
ergibt $\operatorname{tr}_w(X(A\otimes1))=\operatorname{tr}_w((\operatorname{Tr}_2X)A)$, und $A=1$ ergibt
$\operatorname{tr}_w(\iota\mathcal E(X))=\operatorname{tr}_w(X)$. Ist $F=\iota\circ F_0$ eine weitere bedingte Erwartung
auf $\mathcal M_1\otimes1$, die $\operatorname{tr}_w$ erhält, so gilt
$\operatorname{tr}_w(X(A\otimes1))=\operatorname{tr}_w(F(X)(A\otimes1))=m\operatorname{tr}_w(F_0(X)A)$, also
$\operatorname{tr}_w((\operatorname{Tr}_2X-mF_0(X))A)=0$ für alle $A\in\mathcal M_1$, und wegen der Treue ist
$F_0=\mathcal E$ (man nehme $A=(\operatorname{Tr}_2X-mF_0(X))^*$). (b) Man schreibe
$X=\sum_{a,b}X_{ab}\otimes e_{ab}$ mit $X_{ab}\in\mathcal M_1$ (Satz 7.1(a)). Nach Proposition 7.2(b) ist
$\mathcal E(u_{ab}^*X)=\frac{\sqrt m}m\sum_{b'}\operatorname{Tr}_2(X_{ab'}\otimes e_{bb'})=X_{ab}/\sqrt m$, also
$u_{ab}\,\iota(\mathcal E(u_{ab}^*X))=X_{ab}\otimes e_{ab}$; und
$\sum_{a,b}u_{ab}u_{ab}^*=m\sum_{a,b}1\otimes e_{aa}=m^21$. (c) Nach Proposition 7.2(a) ist
$\operatorname{Tr}_2(Y^*)=(\operatorname{Tr}_2Y)^*$, also ist die Adjungierte der Entwicklung von $Y=u_{ab}$ gleich
$u_{ab}^*=\sum_k\iota(\mathcal E(u_{ab}^*v_k))\,v_k^*$, und nach (b) gilt
$\sum_{a,b}u_{ab}u_{ab}^*=\sum_k\big(\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*v_k))\big)v_k^*=\sum_kv_kv_k^*$. $\square$

Eine Familie wie in (b) ist eine *Quasibasis* von $\mathcal E$, und $\sum_kv_kv_k^*$, das nach (c)
nicht von ihr abhängt, ist der **Watatani-Index** $\operatorname{Ind}\mathcal E=m^2$ (Y. Watatani, Mem.
Amer. Math. Soc. 83, no. 424, 1990). Für einen Faktor enthält $M_k\otimes1\subset M_k\otimes M_m$ den
Raum $\mathbb C^k$ mit Vielfachheit $m$ in $\mathbb C^{km}$: Die Inklusionsmatrix ist $\Lambda=(m)$, und
$m^2=\|\Lambda\|^2=\dim(M_k\otimes M_m)/\dim M_k$ ist ihr Jones-Index $[\mathcal M:\mathcal N]$
(V. F. R. Jones, Invent. Math. 72, 1, 1983; F. M. Goodman, P. de la Harpe und V. F. R. Jones,
*Coxeter graphs and towers of algebras*, Springer, 1989). Für allgemeines $\mathcal M_1$ ist die
Inklusion die direkte Summe der $C_1$ Inklusionen $M_{k_c}\otimes1\subset M_{k_c}\otimes M_m$,
$\Lambda=m1_{C_1}$. Das Modul `torch_vn_algebra.inclusions` implementiert Inklusionen, ihre
bedingten Erwartungen und Indizes; `from_tensor_product(alg1, m)` konstruiert diese hier:

```python
from torch_vn_algebra import inclusions
m = 3
inc = inclusions.from_tensor_product(alg1, m)                # alg1 (x) 1 c alg1 (x) M_3
X = random_operator(inc.M)
for trace in ('Tr_blunt', 'Tr_norm', 'tau_vN'):
    E = inc.conditional_expectation(trace.replace('Tr_', ''))
    assert close(E(X), partial_trace(X, keep=1) * (1 / m))   # Proposition 7.4(a): the same E for every trace
    assert abs(getattr(inc.embed(E(X)), trace)() - getattr(X, trace)()).max() < 1e-10
    assert close(E.index(), inc.M.identity() * m**2)         # Proposition 7.4(c)

def unit(a, b):
    e = torch.zeros(m, m, dtype=dt)
    e[a, b] = 1
    return math.sqrt(m) * kron(alg1.identity(), inc.M.tensor_factors[1].from_blocks([e]), inc.M)

u = [unit(a, b) for a in range(m) for b in range(m)]         # the quasi-basis of Proposition 7.4(b)
assert close(sum((v @ inc.embed(E(v.adjoint() @ X))).matrix for v in u), X)
assert close(sum((v @ v.adjoint()).matrix for v in u), inc.M.identity() * m**2)
print(inc.matrix.tolist(), inc.component_indices())           # Lambda = 3 * 1_2, index 9 per component
```

## 7.3 Fusion von Ladungen

Nun trage der Sektor $c$ von $\mathcal M_1$ eine Ladung $q_c$, der Sektor $d$ von $\mathcal M_2$
eine Ladung $r_d$, und eine **Fusionsregel** $f$ liefere die Gesamtladung $f(q_c,r_d)$ – die
Addition für Teilchenzahlen, die Addition modulo $2$ für die Fermionenparität oder irgendeine
andere Abbildung. Man setze
$$
H_q=\bigoplus_{f(q_c,r_d)=q}\mathbb C^{k_c}\otimes\mathbb C^{l_d},\qquad
\Pi_q=\sum_{f(q_c,r_d)=q}P_c\otimes Q_d,\qquad n_q=\sum_{f(q_c,r_d)=q}k_cl_d .
$$
Das **fusionierte Produkt** $\mathcal M_1\otimes_f\mathcal M_2=\{\Pi_q\}'$ besteht aus den Operatoren
auf $H_1\otimes H_2$, die die Gesamtladung erhalten.

**Proposition 7.5.** (a) $\mathcal M_1\otimes_f\mathcal M_2=\bigoplus_qB(H_q)\cong\bigoplus_qM_{n_q}(\mathbb C)$
mit den minimalen zentralen Projektionen $\Pi_q$. (b) $\mathcal M_1\otimes\mathcal M_2\subseteq\mathcal M_1\otimes_f\mathcal M_2$,
mit Gleichheit genau dann, wenn jedes $H_q$ nur ein einziges Paar $(c,d)$ enthält. (c) Sind die
Ladungen reelle Zahlen und ist $f$ die Addition, so gilt $\mathcal M_1\otimes_f\mathcal M_2=\{N\}'$ mit
$N=N_1\otimes1+1\otimes N_2$, $N_1=\sum_cq_cP_c$, $N_2=\sum_dr_dQ_d$.

*Beweis.* (a) Jedes Paar liegt in genau einem $H_q$, also sind die $\Pi_q$ paarweise orthogonale
Projektionen mit Summe $1$, und $X$ kommutiert genau dann mit allen, wenn $X=\sum_q\Pi_qX\Pi_q$
gilt, mit beliebigen Blöcken in $B(H_q)$; das Zentrum ergibt sich wie in Satz 7.1(b). (b) Jedes
Element von $\mathcal M_1\otimes\mathcal M_2$ kommutiert mit jedem $P_c\otimes Q_d$, also auch mit
deren Summen $\Pi_q$. Ist jedes $H_q$ ein einzelnes $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$, so sind
beide Algebren gleich $\bigoplus_{c,d}B(\mathbb C^{k_c}\otimes\mathbb C^{l_d})$; enthält $H_q$ zwei
Paare, so enthält $B(H_q)$ Operatoren, die den Block des einen Paares in den des anderen abbilden
und nicht blockdiagonal sind. (c) $N=\sum_{c,d}(q_c+r_d)P_c\otimes Q_d=\sum_qq\,\Pi_q$, also sind die
von null verschiedenen $\Pi_q$ die Spektralprojektionen von $N$. Sie sind Polynome in $N$ (Beweis
von Satz 1.5(a)), und $N$ ist eine Linearkombination von ihnen; daher kommutiert $X$ genau dann mit
$N$, wenn es mit jedem $\Pi_q$ kommutiert. $\square$

`fused_tensor_product(alg1, alg2, fuse)` (Standard für `fuse`: die Addition) ordnet seine Sektoren
nach der fusionierten Ladung und die Paare $(c,d)$ innerhalb eines Sektors lexikographisch, jedes
mit der Basis $e_a\otimes e_b$. Beispiel: Hard-Core-Bosonen. Die teilchenzahlerhaltenden Operatoren
zweier Plätze bilden $M_1\oplus M_2\oplus M_1$ (Ladungen $0,1,2$), die dreier Plätze
$M_1\oplus M_3\oplus M_3\oplus M_1$; das fusionierte Produkt ist das von fünf Plätzen, bis auf die
Reihenfolge der Basis: $n_q=\sum_c\binom2c\binom3{q-c}=\binom5q$.

```python
bos1 = TypeIAlgebra([1, 2, 1], [1, 2, 1], precision='double', charges=[0, 1, 2])
bos2 = TypeIAlgebra([1, 3, 3, 1], [1, 3, 3, 1], precision='double', charges=[0, 1, 2, 3])
F = fused_tensor_product(bos1, bos2)                         # U(1): total particle number
P = tensor_product(bos1, bos2)                               # particle numbers conserved separately
print(F.charges, F.k_factors, sum(k * k for k in F.k_factors), sum(k * k for k in P.k_factors))
assert F.k_factors == [math.comb(5, q) for q in range(6)]

N1, N2 = full(bos1.central([0.0, 1.0, 2.0]))[0], full(bos2.central([0.0, 1.0, 2.0, 3.0]))[0]
N = torch.kron(N1, torch.eye(8, dtype=dt)) + torch.kron(torch.eye(4, dtype=dt), N2)
X = full12(F.operator(torch.randn(1, F.C, F.k_max, F.k_max, dtype=dt)))[0]
assert close(X @ N, N @ X)                                   # Proposition 7.5(c)
A, B = random_operator(bos1), random_operator(bos2)
assert close(full12(kron(A, B, F)), full12(kron(A, B, P)))   # Proposition 7.5(b)
```

Die fusionierte Algebra hat die Dimension $252$, das einfache Produkt $120$; die Differenz bilden
die Operatoren, die Teilchen zwischen den beiden Teilen verschieben, etwa Hüpfterme.

Ein Zustand von $\mathcal M_1\otimes_f\mathcal M_2$ kann Kohärenzen zwischen den Blöcken
verschiedener Paare haben, sodass Proposition 7.2 nicht anwendbar ist. Sei $E_1(X)=\sum_cP_cXP_c$
auf $B(H_1)$ (die bedingte Erwartung aus Abschnitt 2.6, Physikalische Bedeutung); sie bildet auf
$\mathcal M_1$ ab, und es gilt $\operatorname{Tr}(E_1(X)A)=\operatorname{Tr}(XA)$ für $A\in\mathcal M_1$, da
$A=\sum_cP_cAP_c$.

**Satz 7.6.** Sei $\rho$ die Dichte eines Zustands $\omega$ auf $\mathcal M_1\otimes_f\mathcal M_2$. Die
Dichte seiner Einschränkung $A\mapsto\omega(A\otimes1)$ auf $\mathcal M_1$ ist
$$
\rho_1=E_1\big(\operatorname{Tr}_{H_2}\rho\big).
$$
Der Block $(c,c')$ von $\operatorname{Tr}_{H_2}\rho$ ist $\sum_d\operatorname{tr}_{l_d}\rho_{(c,d),(c',d)}$, summiert
über die $d$ mit $f(q_c,r_d)=f(q_{c'},r_d)$, wobei $\rho_{(c,d),(c',d)}$ der Block von $\rho$ vom
Paar $(c',d)$ zum Paar $(c,d)$ ist. Folglich wirkt $E_1$ genau dann auf allen $\operatorname{Tr}_{H_2}\rho$
trivial, wenn es keine $c\neq c'$ und $d$ mit $f(q_c,r_d)=f(q_{c'},r_d)$ gibt; zum Beispiel, wenn $f$
die Addition ist und die $q_c$ paarweise verschieden sind.

*Beweis.* $A\otimes1\in\mathcal M_1\otimes_f\mathcal M_2$ (Proposition 7.5(b)) und
$\omega(A\otimes1)=\operatorname{Tr}(\rho(A\otimes1))=\operatorname{Tr}\big((\operatorname{Tr}_{H_2}\rho)A\big)
=\operatorname{Tr}\big(E_1(\operatorname{Tr}_{H_2}\rho)A\big)$. Da $E_1(\operatorname{Tr}_{H_2}\rho)\in\mathcal M_1$ ist, ist es
nach der Eindeutigkeit in Satz 3.2 die Dichte. Die Blockformel ist Proposition 7.2(a) für
$\operatorname{Tr}_{H_2}$, mit $i$ im Block $c$, $i'$ im Block $c'$ und $j$ im Block $d$; der Block von $\rho$
zwischen $(c',d)$ und $(c,d)$ verschwindet, außer wenn beide Paare im selben $H_q$ liegen. Gibt es
$c\neq c'$ und $d$ mit $f(q_c,r_d)=f(q_{c'},r_d)=q$, so nehme man Einheitsvektoren
$x\in\mathbb C^{k_c}$, $x'\in\mathbb C^{k_{c'}}$, $y\in\mathbb C^{l_d}$ und
$\psi=(x+x')\otimes y/\sqrt2\in H_q$. Dann ist $|\psi\rangle\langle\psi|\in B(H_q)$, und
$\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|=\frac12|x+x'\rangle\langle x+x'|$ hat den von null verschiedenen
Block $\frac12|x\rangle\langle x'|$, den $E_1$ entfernt. $\square$

Beide Implementierungen der Bibliothek berechnen diese Abbildung, für jede Fusionsregel:
`partial_trace` bildet nur die Spur der diagonalen Paarblöcke $\rho_{(c,d),(c,d)}$, und
`partial_trace_channel` verwendet getrennte Kraus-Indizes für Paare eines Sektors mit demselben
$c$; ihre duale Abbildung ist wieder $A\mapsto A\otimes1$. Mit der Addition bewirkt $E_1$ in unserem
Beispiel nichts; mit der Paritätsregel $f(q,r)=q+r\bmod2$ fusionieren die Ladungen $0$ und $2$ von
`bos1` mit jedem $r$ gleich, und $E_1$ entfernt echte Kohärenzen:

```python
G = fused_tensor_product(bos1, bos2, fuse=lambda q, r: (q + r) % 2)       # only the parity is conserved
print(G.k_factors)

def block_mask(alg):
    """Ones in the diagonal blocks of (+)_c C^{k_c}: X * block_mask(alg) = sum_c P_c X P_c."""
    return full(alg.operator(torch.ones(1, alg.C, alg.k_max, alg.k_max, dtype=dt)))

for alg, coherent in ((F, False), (G, True)):
    rho = states.random_density_matrix(alg, batch_size=2)
    R = full12(rho).reshape(2, 4, 8, 4, 8)
    for keep, X_H, alg_k in ((1, torch.einsum('xijkj->xik', R), bos1), (2, torch.einsum('xijil->xjl', R), bos2)):
        E_X = X_H * block_mask(alg_k)                                     # E(Tr_H rho), Theorem 7.6
        assert close(full(partial_trace(rho, keep=keep)), E_X)
        assert close(full(composite.partial_trace_channel(alg, keep=keep)(rho)), E_X)
        lost = (X_H - E_X).abs().max().item()                             # coherences removed by E
        print(f"keep={keep}: largest coherence removed: {lost:.3e}")
        assert (lost > 1e-3) if coherent else (lost == 0)
    assert close(composite.partial_trace_channel(alg, keep=1).adjoint('blunt')(A), kron(A, bos2.identity(), alg))
```

> **Physikalische Bedeutung.** In der BCS-Beschreibung eines Supraleiters ist nur die Parität der
> Elektronenzahl erhalten, sodass der globale Zustand die Teilchenzahlen $0$ und $2$ eines Gebiets
> überlagern kann. Ein Beobachter, dessen Observablen die Teilchenzahl des Gebiets erhalten, kann
> diese Kohärenz nicht sehen: Für ihn ist der reduzierte Zustand der blockdiagonal gemachte Zustand
> aus Satz 7.6, eine Mischung von Teilchenzahlen.

## 7.4 Verschränkung reiner Zustände

**Lemma 7.7 (Schmidt-Zerlegung).** Jeder Einheitsvektor $\psi\in\mathbb C^a\otimes\mathbb C^b$ lässt
sich als $\psi=\sum_{i=1}^rs_i\,u_i\otimes v_i$ schreiben, mit $s_i>0$, $\sum_is_i^2=1$ und
orthonormalen $u_i\in\mathbb C^a$, $v_i\in\mathbb C^b$. Dann ist
$\operatorname{Tr}_2|\psi\rangle\langle\psi|=\sum_is_i^2|u_i\rangle\langle u_i|$
und $\operatorname{Tr}_1|\psi\rangle\langle\psi|=\sum_is_i^2|v_i\rangle\langle v_i|$: Die beiden reduzierten
Zustände haben dieselben von null verschiedenen Eigenwerte $s_i^2$.

*Beweis.* Sei $\Psi$ die $a\times b$-Matrix $\Psi_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,\psi\rangle$.
Nach dem Spektralsatz ist $\Psi\Psi^*=\sum_{i=1}^a\lambda_i|u_i\rangle\langle u_i|$ mit einer
Orthonormalbasis $u_i$ und $\lambda_i=\|\Psi^*u_i\|^2\ge0$, so angeordnet, dass $\lambda_i>0$ genau
für $i\le r$ gilt. Man setze $s_i=\sqrt{\lambda_i}$ und $w_i=\Psi^*u_i/s_i$ für $i\le r$; dann ist
$\langle w_i,w_j\rangle=\langle u_i,\Psi\Psi^*u_j\rangle/(s_is_j)=\delta_{ij}$, und $\Psi^*u_i=0$ für $i>r$.
Also ist $\Psi=\sum_i|u_i\rangle\langle u_i|\Psi=\sum_{i\le r}|u_i\rangle\langle\Psi^*u_i|=\sum_{i\le r}s_i|u_i\rangle\langle w_i|$,
d. h. $\Psi_{\alpha\beta}=\sum_is_i(u_i)_\alpha\overline{(w_i)_\beta}$: Das ist die Behauptung mit
$v_i=\bar w_i$, und $\sum_is_i^2=\operatorname{Tr}\Psi\Psi^*=\|\psi\|^2=1$. Schließlich ist
$|\psi\rangle\langle\psi|=\sum_{i,j}s_is_j|u_i\rangle\langle u_j|\otimes|v_i\rangle\langle v_j|$, und
Proposition 7.2(b) liefert die reduzierten Zustände. $\square$

Die **Verschränkungsentropie** von $\psi$ ist
$E(\psi)=S(\operatorname{Tr}_2|\psi\rangle\langle\psi|)=S(\operatorname{Tr}_1|\psi\rangle\langle\psi|)=-\sum_is_i^2\log s_i^2$.
Sie verschwindet genau dann, wenn $r=1$, d. h. genau dann, wenn $\psi=u\otimes v$, und es gilt
$E(\psi)\le\log\min(a,b)$, mit Gleichheit genau dann, wenn der reduzierte Zustand auf dem kleineren
Faktor maximal gemischt ist (Korollar 3.6).

Mit Ladungen liegt ein reiner Zustand in einem Sektor $H_q$ des fusionierten Produkts; wir
schreiben $\psi=\sum_{(c,d)}\psi_{cd}$ mit $\psi_{cd}\in\mathbb C^{k_c}\otimes\mathbb C^{l_d}$.

**Proposition 7.8 (symmetrieaufgelöste Verschränkung).** Sei $f$ die Addition, und seien die $q_c$
paarweise verschieden und die $r_d$ paarweise verschieden, sodass es zu jedem $c$ höchstens ein
$d=d(c)$ mit $q_c+r_d=q$ gibt. Man setze $\psi_c=\psi_{c\,d(c)}$ (bzw. $0$), $p_c=\|\psi_c\|^2$ und,
für $p_c>0$, $\hat\psi_c=\psi_c/\sqrt{p_c}$. Dann liegt $\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|$ in
$\mathcal M_1$ (es kommutiert mit $N_1$) und ist gleich
$\bigoplus_cp_c\operatorname{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|$; der reduzierte Zustand auf $\mathcal M_2$
hat bei $d(c)$ die Sektorwahrscheinlichkeit $p_c$; und
$$
E(\psi)=H(p)+\sum_cp_c\,E(\hat\psi_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*Beweis.* Nach der Blockformel aus Satz 7.6 ist der Block $(c,c')$ von
$\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|$ die Summe der $\operatorname{tr}_{l_d}|\psi_{cd}\rangle\langle\psi_{c'd}|$
über die $d$ mit $q_c+r_d=q=q_{c'}+r_d$, was $c=c'$ und $d=d(c)$ erzwingt. Der reduzierte Zustand
ist also $\bigoplus_c\operatorname{Tr}_2|\psi_c\rangle\langle\psi_c|$ mit den Blockspuren $p_c$; Vertauschen der
Rollen der Faktoren liefert die Aussage für $\mathcal M_2$. Nach Lemma 7.7 ist $E(\psi)$ die
Entropie dieses Zustands, und Satz 3.3 mit $\sigma_c=\operatorname{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|$,
$S(\sigma_c)=E(\hat\psi_c)$, liefert die Formel. $\square$

$H(p)$ ist die **Teilchenzahlverschränkung** und $\sum_cp_cE(\hat\psi_c)$ die
**Konfigurationsverschränkung**. Für einen zufälligen Zustand von $2$ Teilchen auf den fünf
Plätzen berechnen wir $E(\psi)$ aus den partiellen Spuren auf beiden Seiten und aus den
Singulärwerten der $4\times8$-Matrix $\Psi$:

```python
s = F.charges.index(2)                                      # the sector of two particles, n_q = 10
n = F.k_factors[s]
psi = torch.randn(n, dtype=dt)
psi = psi / psi.norm()
mat = torch.zeros(1, F.C, F.k_max, F.k_max, dtype=dt)
mat[0, s, :n, :n] = torch.outer(psi, psi.conj())
rho = DensityMatrix(F, matrix=mat)
r1, r2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)

psi_H = torch.zeros(4 * 8, dtype=dt)
psi_H[kron_rows(F)[s]] = psi                                # psi as a vector of H1 (x) H2
lam = torch.linalg.svdvals(psi_H.reshape(4, 8)) ** 2        # squared Schmidt coefficients
lam = lam[lam > 1e-14]
E_psi = -(lam * lam.log()).sum().item()
assert abs(r1.von_neumann_entropy().item() - E_psi) < 1e-10
assert abs(r2.von_neumann_entropy().item() - E_psi) < 1e-10  # Lemma 7.7

def entropy(sigma):
    w = torch.linalg.eigvalsh(sigma).clamp(min=1e-300)
    return -(w * w.log()).sum().item()

p = r1.sector_probabilities()[0]                            # c = 0, 1, 2 particles in the first part
assert torch.allclose(r2.sector_probabilities()[0, [2, 1, 0]], p)   # d(c) = 2 - c
H_num = -(p * p.log()).sum().item()
S_conf = sum(p[c].item() * entropy(r1.matrix[0, c, :k, :k] / p[c]) for c, k in enumerate(bos1.k_factors))
print(f"E = {E_psi:.4f} = H(p) {H_num:.4f} + configurational {S_conf:.4f}")
assert abs(E_psi - H_num - S_conf) < 1e-10                  # Proposition 7.8
```

> **Physikalische Bedeutung.** Parteien, die an die Superauswahlregel der Teilchenzahl gebunden
> sind, können nur den Konfigurationsanteil als Ressource nutzen: die *Verschränkung von Teilchen*
> nach H. M. Wiseman und J. A. Vaccaro (Phys. Rev. Lett. 91, 097902, 2003). $H(p)$ wird in kalten
> Atomen als *Teilchenzahlentropie* (number entropy) gemessen (A. Lukin et al., Science 364, 256,
> 2019).

## 7.5 Zufällige Zustände und die Pagesche Formel

Sei $g\in\mathbb C^N$ mit unabhängigen Koordinaten, deren Real- und Imaginärteile unabhängig
$N(0,\tfrac12)$-verteilt sind, sodass $\mathbb E|g_i|^2=1$. Ein **Haar-zufälliger reiner Zustand** ist
$\psi=g/\|g\|$; seine Verteilung ist unitär invariant, da die Dichte $\pi^{-N}e^{-\|g\|^2}$ und das
Lebesgue-Maß es sind (sie ist die einzige invariante Verteilung auf der Sphäre, was wir nicht
verwenden).

**Lemma 7.9.** $\|g\|$ und $\psi=g/\|g\|$ sind unabhängig.

*Beweis.* Seien $h\ge0$ und $\phi\ge0$ messbar auf $(0,\infty)$ bzw. auf der Einheitssphäre
$S^{2N-1}$. Integration in Polarkoordinaten in $\mathbb R^{2N}$, mit dem Oberflächenmaß $\sigma$,
ergibt
$$
\mathbb E\big[h(\|g\|)\phi(\psi)\big]=\pi^{-N}\int_0^\infty h(t)e^{-t^2}t^{2N-1}dt\int_{S^{2N-1}}\phi\,d\sigma=:\alpha(h)\beta(\phi).
$$
Setzt man $\phi=1$, dann $h=1$, dann beides, so erhält man $\mathbb Eh(\|g\|)=\alpha(h)\beta(1)$,
$\mathbb E\phi(\psi)=\alpha(1)\beta(\phi)$ und $\alpha(1)\beta(1)=1$, also
$\mathbb E[h(\|g\|)\phi(\psi)]=\mathbb Eh(\|g\|)\,\mathbb E\phi(\psi)$. $\square$

**Lemma 7.10.** $\mathbb E\|g\|^4=N(N+1)$. Ist $G$ eine $a\times b$-Matrix mit unabhängigen, wie die
$g_i$ verteilten Einträgen (eine *Ginibre*-Matrix), so ist $\mathbb E\operatorname{Tr}(GG^*)^2=ab(a+b)$.

*Beweis.* Mit $g_i=x+iy$ ist $\mathbb E|g_i|^4=\mathbb Ex^4+2\mathbb Ex^2\mathbb Ey^2+\mathbb Ey^4=\tfrac34+\tfrac12+\tfrac34=2$,
also $\mathbb E\|g\|^4=\sum_i\mathbb E|g_i|^4+\sum_{i\ne j}\mathbb E|g_i|^2\mathbb E|g_j|^2=2N+N(N-1)$. Weiter
ist $\operatorname{Tr}(GG^*)^2=\sum_{i,i',j,j'}G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}$. Die Einträge und ihre
Konjugierten sind komplexe Linearkombinationen unabhängiger reeller Gaußscher Variablen mit
$\mathbb E[G_{ij}G_{kl}]=0$ und $\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}$, also liefert die
Formel von Isserlis (Beweis von Proposition 1.10)
$\mathbb E[G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}]=\delta_{ii'}+\delta_{jj'}$. Summation ergibt
$ab^2+a^2b$. $\square$

**Satz 7.11 (E. Lubkin, J. Math. Phys. 19, 1028, 1978).** Für ein Haar-zufälliges
$\psi\in\mathbb C^a\otimes\mathbb C^b$ und $\rho_A=\operatorname{Tr}_2|\psi\rangle\langle\psi|$ gilt
$$
\mathbb E\operatorname{Tr}\rho_A^2=\frac{a+b}{ab+1}.
$$

*Beweis.* Sei $G_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,g\rangle$, eine Ginibre-Matrix; die Matrix
von $\psi$ ist $\Psi=G/\|g\|$, und nach dem Beweis von Lemma 7.7 ist
$\rho_A=\Psi\Psi^*=GG^*/\|g\|^2$. Folglich ist
$\operatorname{Tr}(GG^*)^2=\|g\|^4\operatorname{Tr}\rho_A^2$, wobei $\operatorname{Tr}\rho_A^2\in[0,1]$ eine Funktion von
$\psi$ ist. Nach Lemma 7.9 ist $\mathbb E\operatorname{Tr}(GG^*)^2=\mathbb E\|g\|^4\,\mathbb E\operatorname{Tr}\rho_A^2$, und nach
Lemma 7.10 ist $\mathbb E\operatorname{Tr}\rho_A^2=ab(a+b)/(ab(ab+1))$. $\square$

Kapitel 6 erhält denselben Mittelwert aus den vierten Momenten eines Haar-Vektors, als
Ausgangsreinheit eines zufälligen Kanals.

**Korollar 7.12.** $\log\frac{ab+1}{a+b}\le\mathbb E\,S(\rho_A)\le\log\min(a,b)$.

*Beweis.* Nach Korollar 3.6 ist $S(\rho_A)\le\log a$, und nach Lemma 7.7 ist
$S(\rho_A)=S(\operatorname{Tr}_1|\psi\rangle\langle\psi|)\le\log b$. Die Jensensche Ungleichung für den konkaven
$\log$ mit den Gewichten $\lambda_i>0$, den von null verschiedenen Eigenwerten von $\rho_A$, liefert
$\sum_i\lambda_i\log\lambda_i\le\log\sum_i\lambda_i^2$, d. h. $S(\rho_A)\ge-\log\operatorname{Tr}\rho_A^2$; die
Jensensche Ungleichung für den konvexen $-\log$ und Satz 7.11 liefern
$\mathbb E[-\log\operatorname{Tr}\rho_A^2]\ge-\log\mathbb E\operatorname{Tr}\rho_A^2$. $\square$

Für festes $a$ und $b\to\infty$ treffen sich die Schranken: Ein kleiner Teil eines zufälligen
reinen Zustands ist fast maximal gemischt. Der exakte Mittelwert ist die Pagesche Formel, die wir
ohne Beweis angeben.

**Satz 7.13 (Pagesche Formel).** Für $a\le b$ gilt
$$
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
$$
Sie wurde von D. N. Page vermutet (Phys. Rev. Lett. 71, 1291, 1993) und von S. K. Foong und
S. Kanno (Phys. Rev. Lett. 72, 1148, 1994), J. Sánchez-Ruiz (Phys. Rev. E 52, 5653, 1995) und
S. Sen (Phys. Rev. Lett. 77, 1, 1996) bewiesen. Für $1\ll a\le b$ ist sie $\approx\log a-\frac a{2b}$;
für zwei Qubits ist sie $\frac13$.

Die Verteilung von $\rho_A$ ist das **induzierte Maß** (K. Życzkowski und H.-J. Sommers, J. Phys. A
34, 7111, 2001). `states.random_density_matrix(alg, rank=r)` bildet $G_cG_c^*$ mit
$k_c\times r$-Ginibre-Blöcken $G_c$ und teilt durch $\sum_c\operatorname{Tr}G_cG_c^*$. Für einen einzigen
Block ist dies die obige Verteilung von $\rho_A$ mit $b=r$. Für mehrere Blöcke ist es die
Einschränkung eines Haar-zufälligen reinen Zustands von $H_1\otimes\mathbb C^r$ auf $\mathcal M_1$: Man
wende Satz 7.6 mit $\mathcal M_2=M_r$ und einer konstanten Fusionsregel an, sodass
$\mathcal M_1\otimes_f\mathcal M_2=B(H_1\otimes\mathbb C^r)$, wobei $G_c$ die Zeilen von $G$ im Block $c$
sind.

**Proposition 7.14.** Für $\rho$, gezogen mit `random_density_matrix(alg, rank=r)` auf
$\bigoplus_cM_{k_c}$, und $D=\sum_ck_c$ gilt
$$
\mathbb E\operatorname{Tr}\rho^2=\frac{\sum_ck_c(k_c+r)}{D(Dr+1)} .
$$

*Beweis.* $\operatorname{Tr}\rho^2=\sum_c\operatorname{Tr}(G_cG_c^*)^2/\|g\|^4$, wobei $g\in\mathbb C^{Dr}$ die
Gesamtheit aller Einträge ist; wie in Satz 7.11 liefern Lemma 7.9 und Lemma 7.10
$\sum_ck_cr(k_c+r)/(Dr(Dr+1))$. $\square$

Wir prüfen die Sätze 7.11 und 7.13 an $10^4$ reinen Zuständen von $\mathbb C^3\otimes\mathbb C^6$, die
mit `partial_trace` reduziert werden, und am induzierten Maß, sowie Proposition 7.14 für
$M_1\oplus M_2\oplus M_3$, $r=2$:

```python
a, b, n_samples = 3, 6, 10_000
g = torch.randn(n_samples, a * b, dtype=dt)                 # complex Gaussian: E|g_i|^2 = 1
psi = g / g.norm(dim=-1, keepdim=True)
AB_alg = tensor_product(TypeIAlgebra([a], [a], precision='double'), TypeIAlgebra([b], [b], precision='double'))
pure = DensityMatrix(AB_alg, matrix=(psi[:, :, None] * psi.conj()[:, None, :])[:, None],
                     validate=False)                        # pure by construction
rho_A = partial_trace(pure, keep=1)
purity_exact = (a + b) / (a * b + 1)
page = sum(1 / k for k in range(b + 1, a * b + 1)) - (a - 1) / (2 * b)

def mean_se(x):
    return x.mean().item(), x.std().item() / math.sqrt(x.numel())

for name, rho_s in (('reduced pure states', rho_A),
                    ('induced measure', states.random_density_matrix(rho_A.algebra, 40_000, rank=b))):
    (m_p, e_p), (m_S, e_S) = mean_se(rho_s.purity()), mean_se(rho_s.von_neumann_entropy())
    print(f"{name}: E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {purity_exact:.4f}), "
          f"E S = {m_S:.4f} +- {e_S:.4f} (Page {page:.4f})")
    assert abs(m_p - purity_exact) < 5 * e_p and abs(m_S - page) < 5 * e_S     # Theorems 7.11, 7.13
    assert math.log((a * b + 1) / (a + b)) < m_S < math.log(a)                # Corollary 7.12

alg3 = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
m_p, e_p = mean_se(states.random_density_matrix(alg3, 40_000, rank=2).purity())
print(f"E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {26 / 78:.4f})")
assert abs(m_p - 26 / 78) < 5 * e_p                                          # Proposition 7.14
```

> **Physikalische Bedeutung.** Ein typischer reiner Zustand ist zwischen jedem Teil, der höchstens
> die Hälfte des Systems umfasst, und dem Rest nahezu maximal verschränkt. Wächst der Teil, so
> steigt die Verschränkung wie $\log a$ und fällt jenseits der Mitte wieder: die *Page-Kurve*, die
> beschreibt, wie Information ein Schwarzes Loch in der Hawking-Strahlung verlässt, und warum die
> Teile eines abgeschlossenen Systems thermisch aussehen.

## 7.6 Womit man vergleichen sollte

Auf einem einzelnen Hilbertraum berechnen QuTiP (`Qobj.ptrace`, `entropy_vn`, `rand_ket`, `rand_dm`
mit der Ginibre-Verteilung und einem Rang) und `numpy` (ein `reshape` und ein `einsum`, wie oben)
partielle Spuren, Entropien und induzierte Zustände und stimmen auf einfachen Produkten mit der
Bibliothek überein. Tensornetzwerk-Bibliotheken mit erhaltenen Ladungen – ITensor
(blockdünnbesetzte Tensoren mit Quantenzahlen), TeNPy (geladene Beine, Verschränkungsspektren
nach Ladung) und, für nichtabelsche Symmetrien, TensorKit.jl (Fusionsbäume) – fusionieren
Sektoren für große Systeme auf dieselbe Weise und liefern die symmetrieaufgelöste Verschränkung.
Die Bibliothek fügt die algebraische Ebene hinzu: Produkte sind von-Neumann-Algebren mit Sektoren
und Spuren, die Fusionsregel ist eine beliebige Abbildung der Bezeichnungen, die partielle Spur ist
ein Kanal mit seiner dualen Abbildung für jede Spur und der korrekten Einschränkung für fusionierte
Produkte (Sätze 7.3 und 7.6), und alles arbeitet gestapelt (im Batch). Tensornetzwerke enthält sie
nicht, und `fused_tensor_product` nimmt eine Regel mit einem einzigen Ergebnis pro Ladungspaar,
wie bei abelschen Symmetrien; eine nichtabelsche Regel wie $j_1\otimes j_2=\bigoplus_jj$ für $SU(2)$
mit ihren Clebsch-Gordan-Koeffizienten ist hier kein fusioniertes Produkt (die Bibliothek behandelt
$SU(2)$ für Spin-$\frac12$-Ketten gesondert, Kapitel 4).

## Aufgaben

1. Man berechne die duale Abbildung von $\operatorname{Tr}_2$ bezüglich $\tau_{\mathrm{vN}}$ auf
   $\mathcal M_1\otimes\mathcal M_2$ und auf $\mathcal M_1$ und überprüfe sie mit
   `partial_trace_channel(alg12).adjoint('tau_vN')`.
2. (Subadditivität.) Sei $\rho$ ein Zustand auf $\mathcal M_1\otimes\mathcal M_2$ mit invertierbaren
   reduzierten Zuständen $\rho_1,\rho_2$. Man zeige $S(\rho)\le S(\rho_1)+S(\rho_2)$, mit Gleichheit
   genau dann, wenn $\rho=\rho_1\otimes\rho_2$.
3. Sei `site` die Algebra $\mathbb C\oplus\mathbb C$ mit den Ladungen $0,1$ (ein leerer oder besetzter
   Platz), mit sich selbst durch Addition fusioniert. Man berechne $E(\psi)$ und seine beiden Anteile
   für $\psi=(|01\rangle+|10\rangle)/\sqrt2$, ein auf zwei Plätze verteiltes Teilchen, und
   beschreibe die Einschränkung des Vektorzustands von $(|00\rangle+|11\rangle)/\sqrt2$ auf die
   fusionierte Algebra.

<details>
<summary>Lösungen</summary>

1. Die Gewichte sind $1/(C_1k_c)$ auf $\mathcal M_1$ und $1/(C_1C_2k_cl_d)$ auf dem Produkt. Für
   $A\in\mathcal M_1$ gilt nach Proposition 7.2(a)
   $\tau_{\mathrm{vN}}((\operatorname{Tr}_2X)A)=\sum_c\frac1{C_1k_c}\sum_d\operatorname{Tr}(X_{cd}(A_c\otimes1))
   =\sum_{c,d}\frac1{C_1C_2k_cl_d}\operatorname{Tr}\big(X_{cd}\,C_2l_d(A_c\otimes1)\big)$, also ist die duale
   Abbildung $A\mapsto\bigoplus_{c,d}C_2l_d\,A_c\otimes1_{l_d}=A\otimes Z$ mit dem zentralen Element
   $Z=\bigoplus_dC_2l_d1_{l_d}$ von $\mathcal M_2$. Sie ist nicht unital: Von den drei Spuren der
   Bibliothek machen nur die ungewichteten $\operatorname{Tr}_2$ zur dualen Abbildung von
   $A\mapsto A\otimes1$ (allgemein leisten das $\operatorname{tr}_w$ auf $\mathcal M_1$ und die Gewichte
   $w_{cd}=w_c$ auf dem Produkt, nach derselben Rechnung).

```python
dual_tau = composite.partial_trace_channel(alg12, keep=1).adjoint('tau_vN')
A = random_operator(alg1)
Z = alg2.central([alg2.C * l for l in alg2.k_factors])
assert close(dual_tau(A), kron(A, Z, alg12))
```

2. $\rho_1\otimes\rho_2$ ist ein Zustand (Satz 7.1(c)), invertierbar, mit Eigenvektoren
   $u_i\otimes v_j$ und Eigenwerten $\lambda_i\mu_j$, also gilt
   $\log(\rho_1\otimes\rho_2)=\log\rho_1\otimes1+1\otimes\log\rho_2$ (Satz 1.5(a)). Nach der Definition
   der partiellen Spur ist $\operatorname{Tr}\rho(\log\rho_1\otimes1)=\operatorname{Tr}\rho_1\log\rho_1$, und ebenso
   für den zweiten Faktor. Folglich ist $S(\rho\|\rho_1\otimes\rho_2)=S(\rho_1)+S(\rho_2)-S(\rho)$, und
   dies ist $\ge0$, mit Gleichheit genau dann, wenn $\rho=\rho_1\otimes\rho_2$ (Korollar 3.6). Diese
   Differenz ist die Transinformation (wechselseitige Information).

```python
rho = states.random_density_matrix(alg12, batch_size=4)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
mutual = rho1.von_neumann_entropy() + rho2.von_neumann_entropy() - rho.von_neumann_entropy()
assert (mutual > 0).all()
assert torch.allclose(mutual, states.relative_entropy(rho, kron(rho1, rho2, alg12)))
```

3. Die Sektoren sind $N=0,1,2$ mit den Größen $1,2,1$. $\psi$ liegt in $N=1$, mit
   $\psi_0=|01\rangle/\sqrt2$, $\psi_1=|10\rangle/\sqrt2$, also ist $p=(\frac12,\frac12)$, und da die
   Blöcke eindimensional sind, ist $E(\hat\psi_c)=0$: Nach Proposition 7.8 ist
   $E(\psi)=H(p)=\log2$ vollständig Teilchenzahlverschränkung, und der nutzbare
   Konfigurationsanteil ist $0$. Der Vektor $(|00\rangle+|11\rangle)/\sqrt2$ überlagert $N=0$ und
   $N=2$; auf der blockdiagonalen fusionierten Algebra ist sein Vektorzustand die Mischung
   $\frac12|00\rangle\langle00|+\frac12|11\rangle\langle11|$ von Produktzuständen (Kapitel 3,
   Aufgabe 3).

```python
site = TypeIAlgebra([1, 1], [1, 1], precision='double', charges=[0, 1])
S2 = fused_tensor_product(site, site)
assert S2.k_factors == [1, 2, 1]
mat = torch.zeros(1, 3, 2, 2, dtype=dt)
mat[0, 1] = 0.5                                       # |psi><psi| with psi = (|01> + |10>)/sqrt 2
r1 = partial_trace(DensityMatrix(S2, matrix=mat), keep=1)
assert torch.allclose(r1.sector_probabilities()[0], torch.tensor([0.5, 0.5], dtype=torch.float64))
assert abs(r1.von_neumann_entropy().item() - math.log(2)) < 1e-12
```

</details>

**Weiter:** Kapitel 8 (in Vorbereitung) behandelt die Dynamik: die Schrödinger- und die
Lindblad-Gleichung in einer Algebra mit Sektoren, Sprünge zwischen Sektoren und Postselektion.
