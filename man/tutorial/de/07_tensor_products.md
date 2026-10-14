[English](../07_tensor_products.md) · [Español](../es/07_tensor_products.md) · [Français](../fr/07_tensor_products.md) · **Deutsch** · [中文](../zh/07_tensor_products.md) · [日本語](../ja/07_tensor_products.md) · [Русский](../ru/07_tensor_products.md)

# Tutorial 7. Tensorprodukte und Verschränkung

Dieses Tutorial begleitet [Kapitel 7](../../lecture_notes/de/07_tensor_products.md) der Vorlesungsnotizen, das die Struktur des einfachen und
des fusionierten Tensorprodukts zweier Algebren mit Sektoren beweist, die partielle Spur als duale
Abbildung einer Einbettung definiert und Verschränkung und zufällige Zustände untersucht. Hier bilden wir
diese Produkte mit der Bibliothek, berechnen partielle Spuren und Verschränkungsentropien und überprüfen
die Sätze numerisch.

Im Folgenden wirke $`\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)`$ auf
$`H_1=\bigoplus_c\mathbb C^{k_c}`$ mit minimalen zentralen Projektionen $`P_c`$, und
$`\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)`$ wirke auf $`H_2=\bigoplus_d\mathbb C^{l_d}`$
mit minimalen zentralen Projektionen $`Q_d`$. $`\mathrm{Tr}`$ ist die Spur des Raumes, auf dem ein
Operator wirkt; auf $`\mathcal M_i`$ ist sie die ungewichtete Spur aus Kapitel 3. Für $`X`$ auf
$`\mathbb C^k\otimes\mathbb C^l`$ schreiben wir $`X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle`$
und identifizieren $`\mathbb C^k\otimes\mathbb C^l`$ mit $`\mathbb C^{kl}`$ durch $`e_a\otimes e_b\mapsto e_{(a-1)l+b}`$
(die Kronecker-Konvention), sodass $`(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}`$ gilt.

<a id="s-7-1"></a>

## 7.1 Das Tensorprodukt zweier Algebren

*Vorlesungsnotizen:* [Abschnitt 7.1](../../lecture_notes/de/07_tensor_products.md#s-7-1): [Satz 7.1](../../lecture_notes/de/07_tensor_products.md#r-7-1) (die Struktur von $`\mathcal M_1\otimes\mathcal M_2`$).

`tensor_product(alg1, alg2)` konstruiert diese Algebra. Der Sektor $`(c,d)`$ hat den Index
$`s=cC_2+d`$ (im Code ab $`0`$ gezählt) und die Ladung $`(q_c,r_d)`$, und innerhalb des Sektors hat
$`e_a\otimes e_b`$ den Index $`a\,l_d+b`$. `kron(A, B, alg12)` berechnet $`A\otimes B`$. Die
Sektorreihenfolge unterscheidet sich von der Kronecker-Reihenfolge von $`H_1\otimes H_2`$ um eine
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

<a id="s-7-2"></a>

## 7.2 Die partielle Spur

*Vorlesungsnotizen:* [Abschnitt 7.2](../../lecture_notes/de/07_tensor_products.md#s-7-2): die Definition von $`\mathrm{Tr}_2`$, [Proposition 7.2](../../lecture_notes/de/07_tensor_products.md#r-7-2), [Satz 7.3](../../lecture_notes/de/07_tensor_products.md#r-7-3) (die partielle Spur ist ein Kanal) und [Proposition 7.4](../../lecture_notes/de/07_tensor_products.md#r-7-4) (eine bedingte Erwartung vom Index $`m^2`$).

`partial_trace(X, keep=1)` formt jeden Block $`\rho_{cd}`$ in die Gestalt $`(k_c,l_d,k_c,l_d)`$ um und
bildet die Spur; `composite.partial_trace_channel(alg12, keep=1)` ist der Kanal aus Satz 7.3 mit
den Kraus-Operatoren $`1_{k_c}\otimes\langle j|`$ vom Sektor $`(c,d)`$ in den Sektor $`c`$. Wir prüfen beide
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

Das Modul `torch_vn_algebra.inclusions` implementiert Inklusionen, ihre bedingten Erwartungen und Indizes;
`from_tensor_product(alg1, m)` konstruiert die Inklusion $`\mathcal M_1\otimes1\subset\mathcal M_1\otimes M_m`$ aus
[Proposition 7.4](../../lecture_notes/de/07_tensor_products.md#r-7-4), mit Inklusionsmatrix $`\Lambda=m1_{C_1}`$ und Index $`m^2`$
([Vorlesungsnotizen, Abschnitt 7.2](../../lecture_notes/de/07_tensor_products.md#s-7-2)):

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

<a id="s-7-3"></a>

## 7.3 Fusion von Ladungen

*Vorlesungsnotizen:* [Abschnitt 7.3](../../lecture_notes/de/07_tensor_products.md#s-7-3): [Proposition 7.5](../../lecture_notes/de/07_tensor_products.md#r-7-5) (das fusionierte Produkt) und [Satz 7.6](../../lecture_notes/de/07_tensor_products.md#r-7-6) (der reduzierte Zustand eines fusionierten Produkts).

`fused_tensor_product(alg1, alg2, fuse)` (Standard für `fuse`: die Addition) ordnet seine Sektoren
nach der fusionierten Ladung und die Paare $`(c,d)`$ innerhalb eines Sektors lexikographisch, jedes
mit der Basis $`e_a\otimes e_b`$. Beispiel: Hard-Core-Bosonen ([Vorlesungsnotizen, Abschnitt 7.3](../../lecture_notes/de/07_tensor_products.md#s-7-3)). Die
teilchenzahlerhaltenden Operatoren bilden auf zwei Plätzen $`M_1\oplus M_2\oplus M_1`$ (Ladungen $`0,1,2`$), auf drei
Plätzen $`M_1\oplus M_3\oplus M_3\oplus M_1`$, und nach [Proposition 7.5(c)](../../lecture_notes/de/07_tensor_products.md#r-7-5) ist ihr fusioniertes Produkt $`\{N\}'`$
für die Gesamtbesetzungszahl von fünf Plätzen, mit $`n_q=\binom5q`$.

```python
bos1 = TypeIAlgebra([1, 2, 1], [1, 2, 1], precision='double', charges=[0, 1, 2])
bos2 = TypeIAlgebra([1, 3, 3, 1], [1, 3, 3, 1], precision='double', charges=[0, 1, 2, 3])
F = fused_tensor_product(bos1, bos2)                         # U(1): total particle number
P = tensor_product(bos1, bos2)                               # particle numbers conserved separately
print(F.charges, F.k_factors, sum(k * k for k in F.k_factors), sum(k * k for k in P.k_factors))
assert F.k_factors == [math.comb(5, q) for q in range(6)]
assert sum(k * k for k in F.k_factors) == 252 and sum(k * k for k in P.k_factors) == 120

N1, N2 = full(bos1.central([0.0, 1.0, 2.0]))[0], full(bos2.central([0.0, 1.0, 2.0, 3.0]))[0]
N = torch.kron(N1, torch.eye(8, dtype=dt)) + torch.kron(torch.eye(4, dtype=dt), N2)
X = full12(F.operator(torch.randn(1, F.C, F.k_max, F.k_max, dtype=dt)))[0]
assert close(X @ N, N @ X)                                   # Proposition 7.5(c)
A, B = random_operator(bos1), random_operator(bos2)
assert close(full12(kron(A, B, F)), full12(kron(A, B, P)))   # Proposition 7.5(b)
```

Die fusionierte Algebra hat die Dimension $`\sum_q\binom5q^2=\binom{10}5=252`$ (wieder Vandermonde),
das einfache Produkt $`(1+4+1)(1+9+9+1)=120`$ (Satz 7.1(b)). Nach Satz 7.1(a) und Proposition 7.5(b)
besteht das einfache Produkt aus den Elementen des fusionierten, die bezüglich der Paare $`(c,d)`$
blockdiagonal sind; die übrigen $`132`$ Dimensionen werden von den Operatoren aufgespannt, die einen
Paarblock $`(c,d)`$ in einen anderen Paarblock $`(c',d')`$ desselben $`H_q`$ abbilden, $`c'\ne c`$ (also
$`d'\ne d`$): Sie ändern die Ladung jedes Teils bei erhaltener Gesamtladung, d. h. hier verschieben
sie Teilchen zwischen den beiden Teilen.

Nach der Blockformel von Satz 7.6 ist $`(\rho_1)_c=\sum_d\mathrm{tr}_{l_d}\,\rho_{(c,d),(c,d)}`$. Beide
Implementierungen der Bibliothek sind so gebaut, dass sie diese Abbildung für jede Fusionsregel
berechnen: `partial_trace` bildet nur die Spur der diagonalen Paarblöcke $`\rho_{(c,d),(c,d)}`$, und
`partial_trace_channel` verwendet getrennte Kraus-Indizes für Paare eines Sektors mit demselben
$`c`$, sodass ihre duale Abbildung wieder $`A\mapsto A\otimes1`$ ist. Der Code prüft beide, für die
beiden Fusionsregeln unten und beide Wahlen des behaltenen Faktors, gegen das unabhängig mit
`einsum` berechnete $`E_1(\mathrm{Tr}_{H_2}\,\rho)`$, und die duale Abbildung gegen $`A\otimes1`$. Mit der Addition bewirkt $`E_1`$ in unserem
Beispiel nichts; mit der Paritätsregel $`f(q,r)=q+r\bmod2`$ fusionieren die Ladungen $`0`$ und $`2`$ von
`bos1` mit jedem $`r`$ gleich, und $`E_1`$ entfernt echte Kohärenzen:

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

<a id="s-7-4"></a>

## 7.4 Verschränkung reiner Zustände

*Vorlesungsnotizen:* [Abschnitt 7.4](../../lecture_notes/de/07_tensor_products.md#s-7-4): [Lemma 7.7](../../lecture_notes/de/07_tensor_products.md#r-7-7) (Schmidt-Zerlegung) und [Proposition 7.8](../../lecture_notes/de/07_tensor_products.md#r-7-8) (symmetrieaufgelöste Verschränkung): $`E(\psi)=H(p)+\sum_cp_cE(\hat\psi_c)`$, wobei $`p_c`$ die Wahrscheinlichkeit des Sektors $`c`$ des ersten Teils ist.

$`H(p)`$ ist die **Teilchenzahlverschränkung** und $`\sum_cp_cE(\hat\psi_c)`$ die
**Konfigurationsverschränkung**. Für einen zufälligen Zustand von $`2`$ Teilchen auf den fünf
Plätzen berechnen wir $`E(\psi)`$ aus den partiellen Spuren auf beiden Seiten und aus den
Singulärwerten der $`4\times8`$-Matrix $`\Psi`$:

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

<a id="s-7-5"></a>

## 7.5 Zufällige Zustände und die Pagesche Formel

*Vorlesungsnotizen:* [Abschnitt 7.5](../../lecture_notes/de/07_tensor_products.md#s-7-5): [Satz 7.11](../../lecture_notes/de/07_tensor_products.md#r-7-11) (Lubkin), [Korollar 7.12](../../lecture_notes/de/07_tensor_products.md#r-7-12), [Satz 7.13](../../lecture_notes/de/07_tensor_products.md#r-7-13) (Pagesche Formel) und [Proposition 7.14](../../lecture_notes/de/07_tensor_products.md#r-7-14) (das induzierte Maß).

Wir prüfen die Sätze 7.11 und 7.13 an $`10^4`$ reinen Zuständen von $`\mathbb C^3\otimes\mathbb C^6`$, die
mit `partial_trace` reduziert werden, und am induzierten Maß, sowie Proposition 7.14 für
$`M_1\oplus M_2\oplus M_3`$, $`r=2`$:

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
assert math.log(a) - (a - 1) / (a * b + 1) - (a - 1) / (2 * b) <= page <= math.log(a) - (a - 1) / (2 * b)
assert abs(sum(1 / k for k in (3, 4)) - 1 / 4 - 1 / 3) < 1e-15                # two qubits: E S = 1/3

alg3 = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')
m_p, e_p = mean_se(states.random_density_matrix(alg3, 40_000, rank=2).purity())
print(f"E Tr rho^2 = {m_p:.4f} +- {e_p:.4f} (exact {26 / 78:.4f})")
assert abs(m_p - 26 / 78) < 5 * e_p                                          # Proposition 7.14
```

<a id="s-7-6"></a>

## 7.6 Womit man vergleichen sollte

Auf einem einzelnen Hilbertraum berechnen QuTiP (`Qobj.ptrace`, `entropy_vn`, `rand_ket`, `rand_dm`
mit der Ginibre-Verteilung und einem Rang) und `numpy` (ein `reshape` und ein `einsum`, wie oben)
partielle Spuren, Entropien und induzierte Zustände für einfache Produkte (die obigen
`einsum`-Prüfungen vergleichen eine solche Rechnung mit der Bibliothek). Tensornetzwerk-Bibliotheken mit erhaltenen Ladungen – ITensor
(blockdünnbesetzte Tensoren mit Quantenzahlen), TeNPy (geladene Beine, Verschränkungsspektren
nach Ladung) und, für nichtabelsche Symmetrien, TensorKit.jl (Fusionsbäume) – fusionieren
Sektoren für große Systeme auf dieselbe Weise und liefern die symmetrieaufgelöste Verschränkung.
Die Bibliothek fügt die algebraische Ebene hinzu: Produkte sind von-Neumann-Algebren mit Sektoren
und Spuren, die Fusionsregel ist eine beliebige Abbildung der Bezeichnungen, die partielle Spur ist
ein Kanal mit seiner dualen Abbildung für jede Spur und der korrekten Einschränkung für fusionierte
Produkte (Sätze 7.3 und 7.6), und alles arbeitet gestapelt (im Batch). Tensornetzwerke enthält sie
nicht, und `fused_tensor_product` nimmt eine Regel mit einem einzigen Ergebnis pro Ladungspaar,
wie bei abelschen Symmetrien; eine nichtabelsche Regel, etwa für $`SU(2)`$, wo sich zwei Spins
$`j_1,j_2`$ zu den mehreren Gesamtspins $`|j_1-j_2|,\dots,j_1+j_2`$ verbinden (Kapitel 4, Aufgabe 2), ist
hier kein fusioniertes Produkt (die Bibliothek behandelt
$`SU(2)`$ für Spin-$`\frac12`$-Ketten gesondert, Kapitel 4).

<a id="x-7"></a>

## Aufgaben

[Aufgaben](../../lecture_notes/de/07_tensor_products.md#x-7) der Vorlesungsnotizen, numerisch:

1. Man überprüfe, dass die duale Abbildung von $`\mathrm{Tr}_2`$ bezüglich $`\tau_{\mathrm{vN}}`$,
   `partial_trace_channel(alg12).adjoint('tau_vN')`, gleich $`A\mapsto A\otimes Z`$ mit $`Z=\bigoplus_dC_2l_d1_{l_d}`$ ist.
2. (Subadditivität.) Man überprüfe an zufälligen Zuständen von `alg12`, dass die gegenseitige Information
   $`S(\rho_1)+S(\rho_2)-S(\rho)`$ positiv ist und gleich $`S(\rho\|\rho_1\otimes\rho_2)`$.
3. Für ein auf zwei Plätze verteiltes Teilchen, $`\psi=(|01\rangle+|10\rangle)/\sqrt2`$, überprüfe man
   $`E(\psi)=H(p)=\log2`$.

<details>
<summary>Lösungen</summary>

*Aufgabe 1.*

```python
dual_tau = composite.partial_trace_channel(alg12, keep=1).adjoint('tau_vN')
A = random_operator(alg1)
Z = alg2.central([alg2.C * l for l in alg2.k_factors])
assert close(dual_tau(A), kron(A, Z, alg12))
```

*Aufgabe 2.*

```python
rho = states.random_density_matrix(alg12, batch_size=4)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
mutual = rho1.von_neumann_entropy() + rho2.von_neumann_entropy() - rho.von_neumann_entropy()
assert (mutual > 0).all()
assert torch.allclose(mutual, states.relative_entropy(rho, kron(rho1, rho2, alg12)))
```

*Aufgabe 3.*

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

**Weiter:** [Tutorial 8](08_dynamics.md): die Schrödinger- und die Lindblad-Gleichung in der Bibliothek.
