[English](../07_tensor_products.md) · **Español** · [Français](../fr/07_tensor_products.md) · [Deutsch](../de/07_tensor_products.md) · [中文](../zh/07_tensor_products.md) · [日本語](../ja/07_tensor_products.md) · [Русский](../ru/07_tensor_products.md)

# 7. Productos tensoriales y entrelazamiento

Este capítulo estudia el producto tensorial de dos álgebras con sectores (físicamente, el álgebra
de dos sistemas considerados conjuntamente; véanse los recuadros *Significado físico*). Con sectores
hay dos productos: el simple, en el que cada parte conserva su propia carga, y el *fusionado*, en el
que solo se conserva la carga total. Demostramos la estructura de ambos, definimos la traza parcial
como el dual de una inclusión y estudiamos el entrelazamiento.

**Aprenderá** cómo se multiplican los sectores y cómo los ordena la biblioteca, por qué la traza
parcial es un canal y, normalizada, una esperanza condicional de índice finito, cuál es el estado
reducido de un estado que conserva la carga, cómo el entrelazamiento se divide en una parte de
número y una parte configuracional, y la pureza y la entropía medias de los estados aleatorios.

En todo el capítulo, $\mathcal M_1=\bigoplus_{c=1}^{C_1}M_{k_c}(\mathbb C)$ actúa sobre
$H_1=\bigoplus_c\mathbb C^{k_c}$ con proyecciones centrales minimales $P_c$, y
$\mathcal M_2=\bigoplus_{d=1}^{C_2}M_{l_d}(\mathbb C)$ actúa sobre $H_2=\bigoplus_d\mathbb C^{l_d}$
con proyecciones centrales minimales $Q_d$. $\operatorname{Tr}$ es la traza del espacio sobre el que
actúa un operador; sobre $\mathcal M_i$ es la traza `blunt` del capítulo 3. Para $X$ sobre
$\mathbb C^k\otimes\mathbb C^l$ escribimos $X_{(a,b),(a',b')}=\langle e_a\otimes e_b,X\,e_{a'}\otimes e_{b'}\rangle$,
e identificamos $\mathbb C^k\otimes\mathbb C^l$ con $\mathbb C^{kl}$ mediante $e_a\otimes e_b\mapsto e_{(a-1)l+b}$
(el convenio de Kronecker), de modo que $(A\otimes B)_{(a,b),(a',b')}=A_{aa'}B_{bb'}$.

## 7.1 El producto tensorial de dos álgebras

Para $*$-álgebras $\mathcal M_1$ sobre $H_1$ y $\mathcal M_2$ sobre $H_2$, el **producto tensorial**
$\mathcal M_1\otimes\mathcal M_2$ es el subespacio generado por los operadores $A\otimes B$ sobre
$H_1\otimes H_2$, con $A\in\mathcal M_1$, $B\in\mathcal M_2$. Como
$(A\otimes B)(A'\otimes B')=AA'\otimes BB'$, $(A\otimes B)^*=A^*\otimes B^*$ y $1\otimes1=1$, es
una $*$-álgebra y, por tanto, un álgebra de von Neumann (teorema 2.2).

**Teorema 7.1.** (a) $H_1\otimes H_2=\bigoplus_{c,d}\mathbb C^{k_c}\otimes\mathbb C^{l_d}$, suma
ortogonal de las imágenes de las proyecciones $P_c\otimes Q_d$, y
$$
\mathcal M_1\otimes\mathcal M_2=\bigoplus_{c,d}M_{k_c}\otimes M_{l_d}\cong\bigoplus_{c,d}M_{k_cl_d}(\mathbb C):
$$
está formada por todos los operadores $X=\sum_{c,d}(P_c\otimes Q_d)X(P_c\otimes Q_d)$, con bloques
arbitrarios sobre los $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$. (b) Sus proyecciones centrales
minimales son las $C_1C_2$ proyecciones $P_c\otimes Q_d$, y
$\dim(\mathcal M_1\otimes\mathcal M_2)=\dim\mathcal M_1\dim\mathcal M_2$.
(c) Para trazas $\operatorname{tr}_w$ sobre $\mathcal M_1$ y $\operatorname{tr}_v$ sobre $\mathcal M_2$ (teorema 3.1),
$\operatorname{tr}_w(A)\operatorname{tr}_v(B)=\operatorname{tr}_{w\otimes v}(A\otimes B)$ con
$(w\otimes v)_{cd}=w_cv_d$. En particular, `Tr_blunt`, `Tr_norm` y `tau_vN` del producto son los
productos de las mismas trazas de los factores.

*Demostración.* (a) Las $P_c\otimes Q_d$ son proyecciones mutuamente ortogonales de suma
$1\otimes1$, y la imagen de $P_c\otimes Q_d$ es $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$. Como $P_c$ y
$Q_d$ son centrales, $(P_c\otimes Q_d)(A\otimes B)=AP_c\otimes BQ_d=(A\otimes B)(P_c\otimes Q_d)$, de
modo que todo elemento de $\mathcal M_1\otimes\mathcal M_2$ es diagonal por bloques. Recíprocamente,
con las unidades matriciales $E^c_{aa'}\in\mathcal M_1$ y $F^d_{bb'}\in\mathcal M_2$ de los bloques
$c$ y $d$, se tiene $E^c_{aa'}\otimes F^d_{bb'}=|e_a\otimes e_b\rangle\langle e_{a'}\otimes e_{b'}|$ dentro de
$\mathbb C^{k_c}\otimes\mathbb C^{l_d}$; estas son las unidades matriciales de $M_{k_cl_d}(\mathbb C)$,
así que todo operador con soporte en el bloque $(c,d)$ pertenece a $\mathcal M_1\otimes\mathcal M_2$.
(b) El centro de una suma directa de álgebras de matrices completas está generado por las
proyecciones sobre los sumandos (capítulo 2, ejercicio 1), y
$\sum_{c,d}k_c^2l_d^2=\sum_ck_c^2\sum_dl_d^2$. (c)
$\operatorname{Tr}(X\otimes Y)=\sum_{a,b}X_{aa}Y_{bb}=\operatorname{Tr}X\operatorname{Tr}Y$, luego
$\operatorname{tr}_{w\otimes v}(A\otimes B)=\sum_{c,d}w_cv_d\operatorname{Tr}A_c\operatorname{Tr}B_d=\operatorname{tr}_w(A)\operatorname{tr}_v(B)$.
Los pesos $1$, $1/k_c$, $1/(C_1k_c)$ de las tres trazas (sección 3.1) se multiplican y dan $1$,
$1/(k_cl_d)$, $1/(C_1C_2k_cl_d)$, los pesos de las mismas trazas sobre el producto. $\square$

`tensor_product(alg1, alg2)` construye esta álgebra. El sector $(c,d)$ tiene índice $s=cC_2+d$
(índices desde $0$ en el código) y carga $(q_c,r_d)$, y dentro de él $e_a\otimes e_b$ tiene índice
$a\,l_d+b$. `kron(A, B, alg12)` calcula $A\otimes B$. El orden de los sectores difiere del orden de
Kronecker de $H_1\otimes H_2$ en una permutación de la base, que `kron_rows` construye
explícitamente más abajo; con ella comparamos con productos de Kronecker de matrices completas.

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

> **Significado físico** (interpretación física, no un teorema de este curso). Si cada uno de dos
> sistemas obedece su propia regla de superselección
> —dos cajas que no pueden intercambiar partículas—, el sistema conjunto tiene un sector por cada
> *par* de cargas, y ningún observable conecta pares distintos. Un sector conjunto cuenta con el
> producto de los pesos (teorema 7.1(c)).

## 7.2 La traza parcial

**Definición.** Sea $\rho\in\mathcal M_1\otimes\mathcal M_2$. La **traza parcial**
$\operatorname{Tr}_2\rho\in\mathcal M_1$ es el elemento tal que
$$
\operatorname{Tr}\big((\operatorname{Tr}_2\rho)\,A\big)=\operatorname{Tr}\big(\rho\,(A\otimes1)\big)\qquad\text{for all }A\in\mathcal M_1 .
$$
Existe y es única por el teorema 3.2, aplicado a la forma lineal
$A\mapsto\operatorname{Tr}(\rho(A\otimes1))$ y a la traza fiel $\operatorname{Tr}$ de $\mathcal M_1$: la
traza parcial es el dual, respecto de las trazas `blunt`, de la inclusión $A\mapsto A\otimes1$.
$\operatorname{Tr}_1$ se define del mismo modo; para un estado $\rho$, $\operatorname{Tr}_2\rho$ y
$\operatorname{Tr}_1\rho$ son los **estados reducidos**. Para las álgebras completas $\mathcal M_i=B(H_i)$ se
trata de la **traza parcial del espacio de Hilbert** $\operatorname{Tr}_{H_2}:B(H_1\otimes H_2)\to B(H_1)$.

**Proposición 7.2.** (a) Para $\rho=\bigoplus_{c,d}\rho_{cd}$,
$$
(\operatorname{Tr}_2\rho)_c=\sum_d\operatorname{tr}_{l_d}\rho_{cd},\qquad
(\operatorname{tr}_lX)_{aa'}=\sum_{b=1}^lX_{(a,b),(a',b)}\quad(X\text{ on }\mathbb C^k\otimes\mathbb C^l),
$$
y $\operatorname{Tr}_2\rho=\operatorname{Tr}_{H_2}\rho$, donde $(\operatorname{Tr}_{H_2}X)_{ii'}=\sum_jX_{(i,j),(i',j)}$.
(b) $\operatorname{Tr}_2(A\otimes B)=(\operatorname{Tr}B)\,A$.
(c) $\operatorname{Tr}_2\big((A\otimes1)X(A'\otimes1)\big)=A\,(\operatorname{Tr}_2X)\,A'$ para $A,A'\in\mathcal M_1$.
(d) $\operatorname{Tr}(\operatorname{Tr}_2X)=\operatorname{Tr}X$, y $\operatorname{Tr}_2$ transforma estados en estados.

*Demostración.* (a) Con $A=\bigoplus_cA_c$, se tiene $A\otimes1=\bigoplus_{c,d}A_c\otimes1_{l_d}$ y
$$
\operatorname{Tr}\big(\rho(A\otimes1)\big)=\sum_{c,d}\sum_{a,a',b}(\rho_{cd})_{(a,b),(a',b)}(A_c)_{a'a}
=\sum_c\operatorname{Tr}\Big(\Big(\sum_d\operatorname{tr}_{l_d}\rho_{cd}\Big)A_c\Big),
$$
así que el elemento de la derecha satisface la definición y es $\operatorname{Tr}_2\rho$. El mismo
cálculo para $B(H_i)$ da la fórmula de $\operatorname{Tr}_{H_2}$; para $\rho$ diagonal por bloques, las
entradas $\rho_{(i,j),(i',j)}$ con $i,i'$ en bloques distintos se anulan, de modo que
$\operatorname{Tr}_{H_2}\rho$ tiene los bloques de (a). (b) $\operatorname{tr}_l(X\otimes Y)=(\operatorname{Tr}Y)X$.
(c) Para $B\in\mathcal M_1$,
$\operatorname{Tr}\big(\operatorname{Tr}_2((A\otimes1)X(A'\otimes1))B\big)=\operatorname{Tr}\big(X(A'BA\otimes1)\big)
=\operatorname{Tr}\big(A(\operatorname{Tr}_2X)A'B\big)$; se concluye por unicidad. (d) Tomemos $A=1$. Si $\rho$ es
un estado, $A\mapsto\operatorname{Tr}(\rho(A\otimes1))$ es positiva, ya que $A^*A\otimes1=(A\otimes1)^*(A\otimes1)$,
y vale $1$ en $A=1$, así que su densidad $\operatorname{Tr}_2\rho$ es un estado (teorema 3.2). $\square$

Los canales son el tema del capítulo 6; aquí solo necesitamos esta definición: una aplicación
lineal $\Phi$ es **completamente positiva** si $\mathrm{id}_n\otimes\Phi$ es positiva sobre
$M_n(\mathbb C)\otimes(\cdot)$ para todo $n$, y es un **canal** si además conserva la traza.

**Teorema 7.3.** Sea $j$ un índice que recorre la base canónica de $H_2$ y sea
$K_j=1\otimes\langle j|:H_1\otimes H_2\to H_1$, $K_j(x\otimes y)=\langle j,y\rangle x$. Entonces
$$
\operatorname{Tr}_{H_2}X=\sum_jK_jXK_j^*,\qquad\sum_jK_j^*K_j=1,\qquad\sum_jK_j^*AK_j=A\otimes1 .
$$
Por tanto, $\operatorname{Tr}_2$ es un canal de $\mathcal M_1\otimes\mathcal M_2$ en $\mathcal M_1$, y su
dual respecto de las trazas `blunt` es $A\mapsto A\otimes1$.

*Demostración.* $K_j^*x=x\otimes e_j$, luego $\langle e_a,K_jXK_j^*e_{a'}\rangle=X_{(a,j),(a',j)}$; sumando
sobre $j$ se obtiene la proposición 7.2(a). Además,
$\sum_jK_j^*K_j(x\otimes y)=\sum_j\langle j,y\rangle x\otimes e_j=x\otimes y$
y $\sum_jK_j^*AK_j=\sum_jA\otimes|j\rangle\langle j|=A\otimes1$. Para $Y\ge0$,
$(\mathrm{id}_n\otimes\operatorname{Tr}_{H_2})(Y)=\sum_j(1_n\otimes K_j)Y(1_n\otimes K_j)^*\ge0$, porque
$\langle x,TYT^*x\rangle=\langle T^*x,YT^*x\rangle\ge0$; esto se mantiene en la subálgebra
$\mathcal M_1\otimes\mathcal M_2$. Por último,
$\operatorname{Tr}\sum_jK_jXK_j^*A=\operatorname{Tr}X\sum_jK_j^*AK_j$, lo que para $A=1$ es la conservación
de la traza y, en general, la dualidad. $\square$

`partial_trace(X, keep=1)` reorganiza cada bloque $\rho_{cd}$ con forma $(k_c,l_d,k_c,l_d)$ y toma
la traza; `composite.partial_trace_channel(alg12, keep=1)` es el canal del teorema 7.3, con
operadores de Kraus $1_{k_c}\otimes\langle j|$ del sector $(c,d)$ al sector $c$. Comprobamos ambos
frente a la definición, frente a un `numpy.einsum` independiente y frente a la dualidad:

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

> **Significado físico** (interpretación física, no un teorema de este curso). El estado reducido
> $\operatorname{Tr}_2\rho$ es todo lo que averigua un observador que solo puede medir el primer sistema; la definición dice que reproduce todos los
> valores esperados de los observables $A\otimes1$. Los estados conjuntos con estados reducidos
> iguales difieren en sus correlaciones.

**La traza parcial como esperanza condicional.** Sea $\mathcal M_2=M_m(\mathbb C)$, de modo que los
sectores de $\mathcal M_1\otimes M_m$ son los de $\mathcal M_1$ (teorema 7.1), sean
$\iota(A)=A\otimes1_m$ y $\mathcal E=\frac1m\operatorname{Tr}_2:\mathcal M_1\otimes M_m\to\mathcal M_1$,
y sean $e_{ab}$ las unidades matriciales de $M_m$.

**Proposición 7.4.** (a) $\iota\circ\mathcal E$ es una esperanza condicional de
$\mathcal M_1\otimes M_m$ sobre $\mathcal M_1\otimes1$ (sección 2.6) y, para toda traza fiel
$\operatorname{tr}_w$ de $\mathcal M_1\otimes M_m$, es la única que conserva $\operatorname{tr}_w$. (b) Los
$m^2$ elementos $u_{ab}=\sqrt m\,(1\otimes e_{ab})$ satisfacen
$X=\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*X))$ para todo $X$, y
$\sum_{a,b}u_{ab}u_{ab}^*=m^2\,1$. (c) Toda familia $(v_k)$ con
$X=\sum_kv_k\,\iota(\mathcal E(v_k^*X))$ para todo $X$ cumple $\sum_kv_kv_k^*=m^2\,1$.

*Demostración.* (a) $\mathcal E$ es positiva (teorema 7.3) y también lo es $\iota$
($T^*T\otimes1=(T\otimes1)^*(T\otimes1)$); $\mathcal E(A\otimes1)=A$ por la proposición 7.2(b), y la
propiedad de bimódulo es la proposición 7.2(c). Escribamos también $\operatorname{tr}_w$ para la traza
con los mismos pesos $w_c$ sobre $\mathcal M_1$; entonces $\operatorname{tr}_w(Y\otimes1)=m\operatorname{tr}_w(Y)$. La
definición de $\operatorname{Tr}_2$ aplicada a $AP_c$, multiplicada por $w_c$ y sumada sobre $c$, da
$\operatorname{tr}_w(X(A\otimes1))=\operatorname{tr}_w((\operatorname{Tr}_2X)A)$, y $A=1$ da
$\operatorname{tr}_w(\iota\mathcal E(X))=\operatorname{tr}_w(X)$. Si $F=\iota\circ F_0$ es otra esperanza
condicional sobre $\mathcal M_1\otimes1$ que conserva $\operatorname{tr}_w$, entonces, como
$F(X(A\otimes1))=F(X)(A\otimes1)$,
$\operatorname{tr}_w(X(A\otimes1))=\operatorname{tr}_w(F(X)(A\otimes1))=m\operatorname{tr}_w(F_0(X)A)$, luego
$\operatorname{tr}_w((\operatorname{Tr}_2X-mF_0(X))A)=0$ para todo $A\in\mathcal M_1$, y $F_0=\mathcal E$ por
fidelidad (tómese $A=(\operatorname{Tr}_2X-mF_0(X))^*$). (b) Escribamos $X=\sum_{a,b}X_{ab}\otimes e_{ab}$ con
$X_{ab}\in\mathcal M_1$ (teorema 7.1(a)). Por la proposición 7.2(b),
$\mathcal E(u_{ab}^*X)=\frac{\sqrt m}m\sum_{b'}\operatorname{Tr}_2(X_{ab'}\otimes e_{bb'})=X_{ab}/\sqrt m$, luego
$u_{ab}\,\iota(\mathcal E(u_{ab}^*X))=X_{ab}\otimes e_{ab}$; y
$\sum_{a,b}u_{ab}u_{ab}^*=m\sum_{a,b}1\otimes e_{aa}=m^21$. (c) $\operatorname{Tr}_2(Y^*)=(\operatorname{Tr}_2Y)^*$ por
la proposición 7.2(a), así que el adjunto del desarrollo de $Y=u_{ab}$ es
$u_{ab}^*=\sum_k\iota(\mathcal E(u_{ab}^*v_k))\,v_k^*$ y, por (b),
$\sum_{a,b}u_{ab}u_{ab}^*=\sum_k\big(\sum_{a,b}u_{ab}\,\iota(\mathcal E(u_{ab}^*v_k))\big)v_k^*=\sum_kv_kv_k^*$. $\square$

Una familia como la de (c) es una *cuasibase* de $\mathcal E$, y $\sum_kv_kv_k^*$, que por (c) no
depende de ella y por (b) existe, es el **índice** $\operatorname{Ind}\mathcal E=m^2$ de $\mathcal E$ (una
definición; el nombre sigue a Y. Watatani, Mem. Amer. Math. Soc. 83, no. 424, 1990). La
**matriz de inclusión** de una inclusión de $\bigoplus_cM_{k_c}$ en $\bigoplus_sM_{n_s}$ tiene las entradas
$\Lambda_{sc}$, el número de copias del bloque $c$ de la subálgebra en el bloque $s$ de la mayor
(el capítulo 9 la usa para las torres de álgebras). Para un factor, $M_k\otimes1\subset M_k\otimes M_m$:
ordenando la base $e_a\otimes e_b$ de $\mathbb C^{km}$ primero por $b$, $A\otimes1_m=\operatorname{diag}(A,\dots,A)$ con
$m$ copias, así que $\Lambda=(m)$, y $m^2=\|\Lambda\|^2=\dim(M_k\otimes M_m)/\dim M_k$. (Lectura adicional,
no usada aquí: V. F. R. Jones, Invent. Math. 72, 1, 1983, sobre el índice de los subfactores, y
F. M. Goodman, P. de la Harpe y V. F. R. Jones, *Coxeter graphs and towers of algebras*, Springer,
1989, sobre las matrices de inclusión.) Para $\mathcal M_1$ general, la inclusión es la suma directa de
las $C_1$ inclusiones $M_{k_c}\otimes1\subset M_{k_c}\otimes M_m$, $\Lambda=m1_{C_1}$. El módulo
`torch_vn_algebra.inclusions` implementa las inclusiones, sus esperanzas condicionales y sus
índices; `from_tensor_product(alg1, m)` construye esta:

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

## 7.3 Fusión de cargas

Supongamos ahora que el sector $c$ de $\mathcal M_1$ lleva una carga $q_c$, el sector $d$ de
$\mathcal M_2$ una carga $r_d$, y que una **regla de fusión** $f$ da la carga total $f(q_c,r_d)$: la
suma para números de partículas, la suma módulo $2$ para la paridad fermiónica, o cualquier otra
aplicación. Pongamos
$$
H_q=\bigoplus_{f(q_c,r_d)=q}\mathbb C^{k_c}\otimes\mathbb C^{l_d},\qquad
\Pi_q=\sum_{f(q_c,r_d)=q}P_c\otimes Q_d,\qquad n_q=\sum_{f(q_c,r_d)=q}k_cl_d .
$$
El **producto fusionado** $\mathcal M_1\otimes_f\mathcal M_2=\{\Pi_q\}'$ está formado por los operadores
sobre $H_1\otimes H_2$ que conservan la carga total.

**Proposición 7.5.** (a) $\mathcal M_1\otimes_f\mathcal M_2=\bigoplus_qB(H_q)\cong\bigoplus_qM_{n_q}(\mathbb C)$,
con proyecciones centrales minimales $\Pi_q$. (b) $\mathcal M_1\otimes\mathcal M_2\subseteq\mathcal M_1\otimes_f\mathcal M_2$,
con igualdad si y solo si cada $H_q$ contiene un único par $(c,d)$. (c) Si las cargas son números
reales y $f$ es la suma, entonces $\mathcal M_1\otimes_f\mathcal M_2=\{N\}'$ con
$N=N_1\otimes1+1\otimes N_2$, $N_1=\sum_cq_cP_c$, $N_2=\sum_dr_dQ_d$.

*Demostración.* (a) Cada par está exactamente en un $H_q$, así que las $\Pi_q$ son proyecciones
mutuamente ortogonales de suma $1$, y $X$ conmuta con todas ellas si y solo si
$X=\sum_q\Pi_qX\Pi_q$, con bloques arbitrarios en $B(H_q)$; el centro es como en el teorema 7.1(b).
(b) Todo elemento de $\mathcal M_1\otimes\mathcal M_2$ conmuta con cada $P_c\otimes Q_d$ y, por tanto,
con sus sumas $\Pi_q$. Si cada $H_q$ es un único $\mathbb C^{k_c}\otimes\mathbb C^{l_d}$, ambas álgebras
son $\bigoplus_{c,d}B(\mathbb C^{k_c}\otimes\mathbb C^{l_d})$; si $H_q$ contiene dos pares, $B(H_q)$
contiene operadores que llevan el bloque de un par al del otro, y que no son diagonales por bloques.
(c) $N=\sum_{c,d}(q_c+r_d)P_c\otimes Q_d=\sum_qq\,\Pi_q$, de modo que las $\Pi_q$ no nulas son las
proyecciones espectrales de $N$. Son polinomios en $N$ (demostración del teorema 1.5(a)) y $N$ es
una combinación de ellas, así que $X$ conmuta con $N$ si y solo si conmuta con cada $\Pi_q$. $\square$

`fused_tensor_product(alg1, alg2, fuse)` (`fuse` por defecto: la suma) ordena sus sectores según la
carga fusionada y los pares $(c,d)$ dentro de un sector lexicográficamente, cada uno con la base
$e_a\otimes e_b$. Ejemplo: bosones de núcleo duro. Sobre $L$ sitios, $(\mathbb C^2)^{\otimes L}$ con el
número de ocupación $N=\sum_in^{(i)}$, $n=\operatorname{diag}(0,1)$, el subespacio propio $N=q$ está generado
por los $\binom Lq$ vectores de la base producto con $q$ sitios ocupados, así que los operadores que
conservan el número $\{N\}'$ forman $\bigoplus_qM_{\binom Lq}$ (proposición 4.10(2)): para dos sitios
$M_1\oplus M_2\oplus M_1$ (cargas $0,1,2$), para tres sitios $M_1\oplus M_3\oplus M_3\oplus M_1$. Por la
proposición 7.5(c), el producto fusionado es $\{N\}'$ para el número de ocupación total de cinco
sitios, salvo el orden de la base, y en efecto $n_q=\sum_c\binom2c\binom3{q-c}=\binom5q$ (identidad de
Vandermonde: un subconjunto de $q$ elementos de cinco sitios consta de $c$ de los dos primeros y
$q-c$ de los tres últimos).

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

El álgebra fusionada tiene dimensión $\sum_q\binom5q^2=\binom{10}5=252$ (de nuevo Vandermonde), el
producto simple $(1+4+1)(1+9+9+1)=120$ (teorema 7.1(b)). Por el teorema 7.1(a) y la proposición 7.5(b),
el producto simple consta de los elementos del fusionado que son diagonales por bloques para los
pares $(c,d)$; las $132$ dimensiones restantes están generadas por los operadores que llevan un bloque
de par $(c,d)$ a otro bloque de par $(c',d')$ del mismo $H_q$, $c'\ne c$ (luego $d'\ne d$): cambian la
carga de cada parte conservando la total, es decir, aquí mueven partículas entre las dos partes.

Un estado de $\mathcal M_1\otimes_f\mathcal M_2$ puede tener coherencias entre bloques de pares
distintos, así que la proposición 7.2 no se aplica. Sea $E_1(X)=\sum_cP_cXP_c$ sobre $B(H_1)$ (la
esperanza condicional de la sección 2.6); su imagen es $\mathcal M_1$, y
$\operatorname{Tr}(E_1(X)A)=\operatorname{Tr}(XA)$ para $A\in\mathcal M_1$, porque $A=\sum_cP_cAP_c$.

**Teorema 7.6.** Sea $\rho$ la densidad de un estado $\omega$ sobre $\mathcal M_1\otimes_f\mathcal M_2$.
La densidad de su restricción $A\mapsto\omega(A\otimes1)$ a $\mathcal M_1$ es
$$
\rho_1=E_1\big(\operatorname{Tr}_{H_2}\rho\big).
$$
El bloque $(c,c')$ de $\operatorname{Tr}_{H_2}\rho$ es $\sum_d\operatorname{tr}_{l_d}\rho_{(c,d),(c',d)}$, sumado sobre
los $d$ con $f(q_c,r_d)=f(q_{c'},r_d)$, donde $\rho_{(c,d),(c',d)}$ es el bloque de $\rho$ del par
$(c',d)$ al par $(c,d)$. Por tanto, $E_1$ actúa trivialmente sobre todos los $\operatorname{Tr}_{H_2}\rho$
si y solo si no existen $c\neq c'$ y $d$ con $f(q_c,r_d)=f(q_{c'},r_d)$; por ejemplo, si $f$ es la
suma y los $q_c$ son distintos dos a dos.

*Demostración.* $A\otimes1\in\mathcal M_1\otimes_f\mathcal M_2$ (proposición 7.5(b)) y
$\omega(A\otimes1)=\operatorname{Tr}(\rho(A\otimes1))=\operatorname{Tr}\big((\operatorname{Tr}_{H_2}\rho)A\big)
=\operatorname{Tr}\big(E_1(\operatorname{Tr}_{H_2}\rho)A\big)$. Como $E_1(\operatorname{Tr}_{H_2}\rho)\in\mathcal M_1$, es
la densidad por la unicidad del teorema 3.2. La fórmula por bloques es la proposición 7.2(a) para
$\operatorname{Tr}_{H_2}$, con $i$ en el bloque $c$, $i'$ en el bloque $c'$ y $j$ en el bloque $d$; el bloque
de $\rho$ entre $(c',d)$ y $(c,d)$ se anula salvo que ambos pares estén en el mismo $H_q$. Si
existen $c\neq c'$ y $d$ con $f(q_c,r_d)=f(q_{c'},r_d)=q$, tomemos vectores unitarios
$x\in\mathbb C^{k_c}$, $x'\in\mathbb C^{k_{c'}}$, $y\in\mathbb C^{l_d}$ y
$\psi=(x+x')\otimes y/\sqrt2\in H_q$. Entonces $|\psi\rangle\langle\psi|\in B(H_q)$, y
$\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|=\frac12|x+x'\rangle\langle x+x'|$ tiene el bloque no nulo
$\frac12|x\rangle\langle x'|$, que $E_1$ elimina. $\square$

Por la fórmula por bloques del teorema 7.6, $(\rho_1)_c=\sum_d\operatorname{tr}_{l_d}\rho_{(c,d),(c,d)}$. Las dos
implementaciones de la biblioteca están construidas para calcular esta aplicación para cualquier
regla de fusión: `partial_trace` toma la traza solo de los bloques diagonales de pares
$\rho_{(c,d),(c,d)}$, y `partial_trace_channel` usa índices de Kraus separados para los pares de un
mismo sector con el mismo $c$, de modo que su dual es de nuevo $A\mapsto A\otimes1$. El código
comprueba ambas, para las dos reglas de fusión siguientes y las dos elecciones del factor
conservado, frente a $E_1(\operatorname{Tr}_{H_2}\rho)$ calculado independientemente con `einsum`, y el dual
frente a $A\otimes1$. Con la suma, $E_1$ no hace
nada en nuestro ejemplo; con la regla de paridad $f(q,r)=q+r\bmod2$, las cargas $0$ y $2$ de `bos1`
se fusionan igual con todo $r$, y $E_1$ elimina coherencias genuinas:

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

> **Significado físico** (interpretación física, no un teorema de este curso). En la descripción
> BCS de un superconductor solo se conserva la paridad
> del número de electrones, de modo que el estado global puede superponer los números de
> partículas $0$ y $2$ de una región. Un observador cuyos observables conservan el número de
> partículas de la región no puede ver esta coherencia: para él, el estado reducido es el estado
> diagonal por bloques del teorema 7.6 (el contenido matemático de este recuadro), una mezcla de
> números de partículas.

## 7.4 Entrelazamiento de estados puros

**Lema 7.7 (descomposición de Schmidt).** Todo vector unitario $\psi\in\mathbb C^a\otimes\mathbb C^b$
puede escribirse $\psi=\sum_{i=1}^rs_i\,u_i\otimes v_i$ con $s_i>0$, $\sum_is_i^2=1$ y
$u_i\in\mathbb C^a$, $v_i\in\mathbb C^b$ ortonormales. Entonces
$\operatorname{Tr}_2|\psi\rangle\langle\psi|=\sum_is_i^2|u_i\rangle\langle u_i|$
y $\operatorname{Tr}_1|\psi\rangle\langle\psi|=\sum_is_i^2|v_i\rangle\langle v_i|$: los dos estados reducidos
tienen los mismos valores propios no nulos $s_i^2$.

*Demostración.* Sea $\Psi$ la matriz $a\times b$ dada por $\Psi_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,\psi\rangle$.
Por el teorema espectral, $\Psi\Psi^*=\sum_{i=1}^a\lambda_i|u_i\rangle\langle u_i|$ con una base
ortonormal $u_i$ y $\lambda_i=\|\Psi^*u_i\|^2\ge0$, ordenados de modo que $\lambda_i>0$ exactamente
para $i\le r$. Pongamos $s_i=\sqrt{\lambda_i}$ y $w_i=\Psi^*u_i/s_i$ para $i\le r$; entonces
$\langle w_i,w_j\rangle=\langle u_i,\Psi\Psi^*u_j\rangle/(s_is_j)=\delta_{ij}$, y $\Psi^*u_i=0$ para $i>r$.
Por tanto $\Psi=\sum_i|u_i\rangle\langle u_i|\Psi=\sum_{i\le r}|u_i\rangle\langle\Psi^*u_i|=\sum_{i\le r}s_i|u_i\rangle\langle w_i|$,
es decir, $\Psi_{\alpha\beta}=\sum_is_i(u_i)_\alpha\overline{(w_i)_\beta}$: es la afirmación con
$v_i=\bar w_i$, y $\sum_is_i^2=\operatorname{Tr}\Psi\Psi^*=\|\psi\|^2=1$. Por último,
$|\psi\rangle\langle\psi|=\sum_{i,j}s_is_j|u_i\rangle\langle u_j|\otimes|v_i\rangle\langle v_j|$, y la
proposición 7.2(b) da los estados reducidos. $\square$

La **entropía de entrelazamiento** de $\psi$ es
$E(\psi)=S(\operatorname{Tr}_2|\psi\rangle\langle\psi|)=S(\operatorname{Tr}_1|\psi\rangle\langle\psi|)=-\sum_is_i^2\log s_i^2$.
Se anula si y solo si $r=1$ (cada término $-s_i^2\log s_i^2$ es $\ge0$, y vale $0$ solo para
$s_i^2=1$), es decir, si y solo si $\psi=u\otimes v$ (entonces $\Psi=uv^T$ tiene rango $1$, y
$r=\operatorname{rank}\Psi\Psi^*=\operatorname{rank}\Psi$), y $E(\psi)\le\log\min(a,b)$, con igualdad si y solo si
el estado reducido del factor más pequeño es máximamente mixto (corolario 3.6, aplicado a ambos
estados reducidos).

Con cargas tomamos un vector unitario $\psi$ con una carga total definida, es decir, en un sector
$H_q$ del producto fusionado (una superposición de distintos $H_q$ da sobre
$\mathcal M_1\otimes_f\mathcal M_2$ el mismo estado que una mezcla, capítulo 3, ejercicio 3); escribamos
$\psi=\sum_{(c,d)}\psi_{cd}$ con $\psi_{cd}\in\mathbb C^{k_c}\otimes\mathbb C^{l_d}$.

**Proposición 7.8 (entrelazamiento resuelto por simetría).** Sea $f$ la suma, y supongamos que los
$q_c$ son distintos dos a dos y que los $r_d$ son distintos dos a dos, de modo que para cada $c$ a
lo sumo un $d=d(c)$ cumple $q_c+r_d=q$. Pongamos $\psi_c=\psi_{c\,d(c)}$ (o $0$),
$p_c=\|\psi_c\|^2$ y, para $p_c>0$, $\hat\psi_c=\psi_c/\sqrt{p_c}$. Entonces
$\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|$ pertenece a $\mathcal M_1$ (conmuta con $N_1$) y es igual a
$\bigoplus_cp_c\operatorname{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|$; el estado reducido sobre $\mathcal M_2$
tiene la probabilidad de sector $p_c$ en $d(c)$; y
$$
E(\psi)=H(p)+\sum_cp_c\,E(\hat\psi_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*Demostración.* Por la fórmula por bloques del teorema 7.6, el bloque $(c,c')$ de
$\operatorname{Tr}_{H_2}|\psi\rangle\langle\psi|$ es la suma de $\operatorname{tr}_{l_d}|\psi_{cd}\rangle\langle\psi_{c'd}|$
sobre los $d$ con $q_c+r_d=q=q_{c'}+r_d$, lo que obliga a que $c=c'$ y $d=d(c)$. Así, el estado
reducido es $\bigoplus_c\operatorname{Tr}_2|\psi_c\rangle\langle\psi_c|$, con trazas de bloque $p_c$;
intercambiando los papeles de los factores se obtiene el enunciado para $\mathcal M_2$. Por el
lema 7.7, $E(\psi)$ es la entropía de este estado, y el teorema 3.3 con
$\sigma_c=\operatorname{Tr}_2|\hat\psi_c\rangle\langle\hat\psi_c|$, $S(\sigma_c)=E(\hat\psi_c)$, da la
fórmula. $\square$

$H(p)$ es el **entrelazamiento de número** y $\sum_cp_cE(\hat\psi_c)$ el **entrelazamiento
configuracional**. Para un estado aleatorio de $2$ partículas en los cinco sitios calculamos
$E(\psi)$ a partir de las trazas parciales de ambos lados y de los valores singulares de la matriz
$\Psi$ de tamaño $4\times8$:

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

> **Significado físico** (interpretación física, no un teorema de este curso). Las partes sujetas
> a la regla de superselección del número de partículas solo pueden usar como recurso la parte configuracional: el *entrelazamiento de partículas* de
> H. M. Wiseman y J. A. Vaccaro (Phys. Rev. Lett. 91, 097902, 2003). $H(p)$ se mide en átomos fríos
> como la *entropía de número* (A. Lukin et al., Science 364, 256, 2019).

## 7.5 Estados aleatorios y la fórmula de Page

Sea $g\in\mathbb C^N$ con coordenadas independientes cuyas partes real e imaginaria son
independientes con ley $N(0,\tfrac12)$, de modo que $\mathbb E|g_i|^2=1$. Un **estado puro aleatorio
de Haar** es $\psi=g/\|g\|$; su ley es invariante unitaria, como lo son la densidad
$\pi^{-N}e^{-\|g\|^2}$ y la medida de Lebesgue (es la única ley invariante sobre la esfera, hecho que
no utilizaremos).

**Lema 7.9.** $\|g\|$ y $\psi=g/\|g\|$ son independientes.

*Demostración.* Sean $h\ge0$ y $\phi\ge0$ medibles sobre $(0,\infty)$ y sobre la esfera unidad
$S^{2N-1}$, respectivamente. La integración en coordenadas polares en $\mathbb R^{2N}$, con $\sigma$
la medida de superficie, da
$$
\mathbb E\big[h(\|g\|)\phi(\psi)\big]=\pi^{-N}\int_0^\infty h(t)e^{-t^2}t^{2N-1}dt\int_{S^{2N-1}}\phi\,d\sigma=:\alpha(h)\beta(\phi).
$$
Tomando $\phi=1$, luego $h=1$ y luego ambas, se obtiene $\mathbb Eh(\|g\|)=\alpha(h)\beta(1)$,
$\mathbb E\phi(\psi)=\alpha(1)\beta(\phi)$ y $\alpha(1)\beta(1)=1$, luego
$\mathbb E[h(\|g\|)\phi(\psi)]=\mathbb Eh(\|g\|)\,\mathbb E\phi(\psi)$. $\square$

**Lema 7.10.** $\mathbb E\|g\|^4=N(N+1)$. Si $G$ es una matriz $a\times b$ con entradas independientes
distribuidas como los $g_i$ (una matriz *de Ginibre*), entonces $\mathbb E\operatorname{Tr}(GG^*)^2=ab(a+b)$.

*Demostración.* Con $g_i=x+iy$, $\mathbb E|g_i|^4=\mathbb Ex^4+2\mathbb Ex^2\mathbb Ey^2+\mathbb Ey^4=\tfrac34+\tfrac12+\tfrac34=2$,
luego $\mathbb E\|g\|^4=\sum_i\mathbb E|g_i|^4+\sum_{i\ne j}\mathbb E|g_i|^2\mathbb E|g_j|^2=2N+N(N-1)$. Además,
$\operatorname{Tr}(GG^*)^2=\sum_{i,i',j,j'}G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}$. Las entradas y sus
conjugadas son combinaciones lineales complejas de gaussianas reales independientes con
$\mathbb E[G_{ij}G_{kl}]=0$ y $\mathbb E[G_{ij}\bar G_{kl}]=\delta_{ik}\delta_{jl}$, así que la fórmula de
Isserlis (demostración de la proposición 1.10) da
$\mathbb E[G_{ij}\bar G_{i'j}G_{i'j'}\bar G_{ij'}]=\delta_{ii'}+\delta_{jj'}$. Sumando se obtiene
$ab^2+a^2b$. $\square$

**Teorema 7.11 (E. Lubkin, J. Math. Phys. 19, 1028, 1978).** Para $\psi\in\mathbb C^a\otimes\mathbb C^b$
aleatorio de Haar y $\rho_A=\operatorname{Tr}_2|\psi\rangle\langle\psi|$,
$$
\mathbb E\operatorname{Tr}\rho_A^2=\frac{a+b}{ab+1}.
$$

*Demostración.* Sea $G_{\alpha\beta}=\langle e_\alpha\otimes e_\beta,g\rangle$, una matriz de Ginibre; la
matriz de $\psi$ es $\Psi=G/\|g\|$ y, por la demostración del lema 7.7,
$\rho_A=\Psi\Psi^*=GG^*/\|g\|^2$. Por tanto
$\operatorname{Tr}(GG^*)^2=\|g\|^4\operatorname{Tr}\rho_A^2$, donde $\operatorname{Tr}\rho_A^2\in[0,1]$ es función de
$\psi$. Por el lema 7.9, $\mathbb E\operatorname{Tr}(GG^*)^2=\mathbb E\|g\|^4\,\mathbb E\operatorname{Tr}\rho_A^2$, y por el
lema 7.10, $\mathbb E\operatorname{Tr}\rho_A^2=ab(a+b)/(ab(ab+1))$. $\square$

El capítulo 6 obtiene la misma media a partir de los cuartos momentos de un vector de Haar, como
la pureza de salida de un canal aleatorio.

**Corolario 7.12.** $\log\frac{ab+1}{a+b}\le\mathbb E\,S(\rho_A)\le\log\min(a,b)$.

*Demostración.* Por el corolario 3.6, $S(\rho_A)\le\log a$, y por el lema 7.7,
$S(\rho_A)=S(\operatorname{Tr}_1|\psi\rangle\langle\psi|)\le\log b$. La desigualdad de Jensen para la función
cóncava $\log$ con los pesos $\lambda_i>0$, los valores propios no nulos de $\rho_A$, da
$\sum_i\lambda_i\log\lambda_i\le\log\sum_i\lambda_i^2$, es decir, $S(\rho_A)\ge-\log\operatorname{Tr}\rho_A^2$; la
desigualdad de Jensen para la función convexa $-\log$ y el teorema 7.11 dan
$\mathbb E[-\log\operatorname{Tr}\rho_A^2]\ge-\log\mathbb E\operatorname{Tr}\rho_A^2$. $\square$

Para $a$ fijo y $b\to\infty$ las cotas se juntan: una parte pequeña de un estado puro aleatorio es
casi máximamente mixta. La media exacta es la fórmula de Page, que enunciamos sin demostración.

**Teorema 7.13 (fórmula de Page).** Para $a\le b$,
$$
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
$$
Fue conjeturada por D. N. Page (Phys. Rev. Lett. 71, 1291, 1993) y demostrada por S. K. Foong y
S. Kanno (Phys. Rev. Lett. 72, 1148, 1994), J. Sánchez-Ruiz (Phys. Rev. E 52, 5653, 1995) y
S. Sen (Phys. Rev. Lett. 77, 1, 1996). Para $1\ll a\le b$ vale $\approx\log a-\frac a{2b}$; para dos
cúbits vale $\frac13$.

La ley de $\rho_A$ es la **medida inducida** (K. Życzkowski y H.-J. Sommers, J. Phys. A 34,
7111, 2001). `states.random_density_matrix(alg, rank=r)` forma $G_cG_c^*$ con bloques de Ginibre
$G_c$ de tamaño $k_c\times r$ y divide por $\sum_c\operatorname{Tr}G_cG_c^*$. Para un solo bloque, es la ley
de $\rho_A$ anterior con $b=r$. Para varios bloques, es la restricción a $\mathcal M_1$ de un estado
puro aleatorio de Haar de $H_1\otimes\mathbb C^r$: aplíquese el teorema 7.6 con $\mathcal M_2=M_r$ y una
regla de fusión constante, de modo que $\mathcal M_1\otimes_f\mathcal M_2=B(H_1\otimes\mathbb C^r)$, con
$G_c$ las filas de $G$ del bloque $c$.

**Proposición 7.14.** Para $\rho$ generada por `random_density_matrix(alg, rank=r)` sobre
$\bigoplus_cM_{k_c}$, con $D=\sum_ck_c$,
$$
\mathbb E\operatorname{Tr}\rho^2=\frac{\sum_ck_c(k_c+r)}{D(Dr+1)} .
$$

*Demostración.* $\operatorname{Tr}\rho^2=\sum_c\operatorname{Tr}(G_cG_c^*)^2/\|g\|^4$, con $g\in\mathbb C^{Dr}$ el conjunto
de todas las entradas; como en el teorema 7.11, los lemas 7.9 y 7.10 dan
$\sum_ck_cr(k_c+r)/(Dr(Dr+1))$. $\square$

Comprobamos los teoremas 7.11 y 7.13 sobre $10^4$ estados puros de $\mathbb C^3\otimes\mathbb C^6$
reducidos con `partial_trace` y sobre la medida inducida, y la proposición 7.14 para
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

> **Significado físico.** Un estado puro típico está casi máximamente entrelazado entre cualquier
> parte que no supere la mitad del sistema y el resto. Al crecer la parte, el entrelazamiento
> aumenta como $\log a$ y disminuye pasada la mitad: es la *curva de Page*, que describe cómo sale
> la información de un agujero negro en la radiación de Hawking y por qué las partes de un sistema
> cerrado parecen térmicas.

## 7.6 Con qué comparar

Sobre un único espacio de Hilbert, QuTiP (`Qobj.ptrace`, `entropy_vn`, `rand_ket`, `rand_dm` con la
distribución de Ginibre y un rango) y `numpy` (un `reshape` y un `einsum`, como antes) calculan
trazas parciales, entropías y estados inducidos, y coinciden con la biblioteca en los productos
simples. Las bibliotecas de redes tensoriales con cargas conservadas —ITensor (tensores dispersos
por bloques con números cuánticos), TeNPy (patas con carga, espectros de entrelazamiento por carga)
y, para simetrías no abelianas, TensorKit.jl (árboles de fusión)— fusionan los sectores del mismo
modo para sistemas grandes y dan el entrelazamiento resuelto por simetría. La biblioteca añade la
capa algebraica: los productos son álgebras de von Neumann con sectores y trazas, la regla de
fusión es cualquier aplicación de las etiquetas, la traza parcial es un canal con su dual para cada
traza y la restricción correcta para los productos fusionados (teoremas 7.3 y 7.6), y todo ello
funciona por lotes. No incluye redes tensoriales, y `fused_tensor_product` admite una regla con un
único resultado por par de cargas, como en las simetrías abelianas; una regla no abeliana como
$j_1\otimes j_2=\bigoplus_jj$ para $SU(2)$, con sus coeficientes de Clebsch–Gordan, no es aquí un
producto fusionado (la biblioteca trata $SU(2)$ por separado para las cadenas de espín $\frac12$,
capítulo 4).

## Ejercicios

1. Calcule el dual de $\operatorname{Tr}_2$ respecto de $\tau_{\mathrm{vN}}$ tanto sobre
   $\mathcal M_1\otimes\mathcal M_2$ como sobre $\mathcal M_1$, y compruébelo con
   `partial_trace_channel(alg12).adjoint('tau_vN')`.
2. (Subaditividad.) Sea $\rho$ un estado sobre $\mathcal M_1\otimes\mathcal M_2$ con estados reducidos
   invertibles $\rho_1,\rho_2$. Demuestre que $S(\rho)\le S(\rho_1)+S(\rho_2)$, con igualdad si y
   solo si $\rho=\rho_1\otimes\rho_2$.
3. Sea `site` el álgebra $\mathbb C\oplus\mathbb C$ con cargas $0,1$ (un sitio vacío u ocupado),
   fusionada consigo misma mediante la suma. Calcule $E(\psi)$ y sus dos partes para
   $\psi=(|01\rangle+|10\rangle)/\sqrt2$, una partícula compartida por dos sitios, y describa la
   restricción del estado vectorial de $(|00\rangle+|11\rangle)/\sqrt2$ al álgebra fusionada.

<details>
<summary>Soluciones</summary>

1. Los pesos son $1/(C_1k_c)$ sobre $\mathcal M_1$ y $1/(C_1C_2k_cl_d)$ sobre el producto. Para
   $A\in\mathcal M_1$, por la proposición 7.2(a),
   $\tau_{\mathrm{vN}}((\operatorname{Tr}_2X)A)=\sum_c\frac1{C_1k_c}\sum_d\operatorname{Tr}(X_{cd}(A_c\otimes1))
   =\sum_{c,d}\frac1{C_1C_2k_cl_d}\operatorname{Tr}\big(X_{cd}\,C_2l_d(A_c\otimes1)\big)$, así que el dual es
   $A\mapsto\bigoplus_{c,d}C_2l_d\,A_c\otimes1_{l_d}=A\otimes Z$ con el elemento central
   $Z=\bigoplus_dC_2l_d1_{l_d}$ de $\mathcal M_2$. No es unital: de las tres trazas de la
   biblioteca, solo las `blunt` hacen de $\operatorname{Tr}_2$ el dual de $A\mapsto A\otimes1$ (en general
   lo hacen $\operatorname{tr}_w$ sobre $\mathcal M_1$ y los pesos $w_{cd}=w_c$ sobre el producto, por el
   mismo cálculo).

```python
dual_tau = composite.partial_trace_channel(alg12, keep=1).adjoint('tau_vN')
A = random_operator(alg1)
Z = alg2.central([alg2.C * l for l in alg2.k_factors])
assert close(dual_tau(A), kron(A, Z, alg12))
```

2. $\rho_1\otimes\rho_2$ es un estado (teorema 7.1(c)), invertible, con vectores propios
   $u_i\otimes v_j$ y valores propios $\lambda_i\mu_j$, luego
   $\log(\rho_1\otimes\rho_2)=\log\rho_1\otimes1+1\otimes\log\rho_2$ (teorema 1.5(a)). Por la
   definición de la traza parcial,
   $\operatorname{Tr}\rho(\log\rho_1\otimes1)=\operatorname{Tr}\rho_1\log\rho_1$, y análogamente para el segundo
   factor. Por tanto $S(\rho\|\rho_1\otimes\rho_2)=S(\rho_1)+S(\rho_2)-S(\rho)$, que es $\ge0$, con
   igualdad si y solo si $\rho=\rho_1\otimes\rho_2$ (corolario 3.6). Esta diferencia es la
   información mutua.

```python
rho = states.random_density_matrix(alg12, batch_size=4)
rho1, rho2 = partial_trace(rho, keep=1), partial_trace(rho, keep=2)
mutual = rho1.von_neumann_entropy() + rho2.von_neumann_entropy() - rho.von_neumann_entropy()
assert (mutual > 0).all()
assert torch.allclose(mutual, states.relative_entropy(rho, kron(rho1, rho2, alg12)))
```

3. Los sectores son $N=0,1,2$, de tamaños $1,2,1$. $\psi$ está en $N=1$, con
   $\psi_0=|01\rangle/\sqrt2$, $\psi_1=|10\rangle/\sqrt2$, luego $p=(\frac12,\frac12)$ y, como los
   bloques son unidimensionales, $E(\hat\psi_c)=0$: por la proposición 7.8, $E(\psi)=H(p)=\log2$
   es todo entrelazamiento de número, y la parte configuracional utilizable es $0$. El vector
   $(|00\rangle+|11\rangle)/\sqrt2$ superpone $N=0$ y $N=2$; sobre el álgebra fusionada, diagonal
   por bloques, su estado vectorial es la mezcla
   $\frac12|00\rangle\langle00|+\frac12|11\rangle\langle11|$ de estados producto (capítulo 3,
   ejercicio 3).

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

**Siguiente:** el [capítulo 8](08_dynamics.md) trata la dinámica: las ecuaciones de Schrödinger y de
Lindblad en un álgebra con sectores, los saltos entre sectores y la postselección.
