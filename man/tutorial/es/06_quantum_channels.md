[English](../06_quantum_channels.md) · **Español** · [Français](../fr/06_quantum_channels.md) · [Deutsch](../de/06_quantum_channels.md) · [中文](../zh/06_quantum_channels.md) · [日本語](../ja/06_quantum_channels.md) · [Русский](../ru/06_quantum_channels.md)

# Tutorial 6. Canales cuánticos

Este tutorial acompaña al [capítulo 6](../../lecture_notes/es/06_quantum_channels.md) de los apuntes de clase, que demuestra las
descripciones de las aplicaciones completamente positivas entre sumas directas de álgebras de matrices
—matrices de Choi, operadores de Kraus, la dilatación de Stinespring— y estudia los duales, los puntos
fijos, las matrices de transición y los canales aleatorios. Aquí vemos cómo almacena la biblioteca los
canales y comprobamos numéricamente estos resultados.

En todo el capítulo, $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ y
$`\mathcal N=\bigoplus_{d=1}^{C'}M_{k'_d}(\mathbb C)`$ con proyecciones centrales minimales $`P_c`$ y
$`Q_d`$; $`A=\bigoplus_cA_c`$; $`M_{k'\times k}`$ es el espacio de las matrices complejas $`k'\times k`$;
$`E_{ij}=|i\rangle\langle j|`$ son las unidades matriciales; $`\mathrm{Tr}`$ es la traza `blunt`
$`\sum_c\mathrm{Tr}\,A_c`$, y $`\mathrm{tr}_w=\sum_cw_c\mathrm{Tr}\,A_c`$ una traza fiel
(teorema 3.1). Las matrices por bloques actúan sobre $`\mathbb C^n\otimes\mathbb C^k`$, con el índice
de bloque en primer lugar.

<a id="s-6-1"></a>

## 6.1 Aplicaciones completamente positivas y teorema de Choi

*Apuntes:* [sección 6.1](../../lecture_notes/es/06_quantum_channels.md#s-6-1): [proposición 6.1](../../lecture_notes/es/06_quantum_channels.md#r-6-1) (la trasposición), [teorema 6.2](../../lecture_notes/es/06_quantum_channels.md#r-6-2) (Choi), [corolario 6.3](../../lecture_notes/es/06_quantum_channels.md#r-6-3) (rango de Kraus) y [proposición 6.4](../../lecture_notes/es/06_quantum_channels.md#r-6-4) (aplicaciones entre sumas directas).

La biblioteca tiene dos clases. `Channel` almacena las aplicaciones **que preservan los sectores**
($`C'=C`$, $`k'=k`$, $`\Phi_{dc}=0`$ para $`d\ne c`$), cuyos operadores de Kraus $`K_l=\bigoplus_cK^{cc}_l`$
están en $`\mathcal M`$, como un tensor de forma `(batch, r, C, k_max, k_max)`; `InterSectorChannel`
almacena la aplicación general de la [proposición 6.4](../../lecture_notes/es/06_quantum_channels.md#r-6-4) con forma `(batch, r, C', C, k'_max, k_max)`.
`Channel.choi()` devuelve $`J(\Phi_{cc})`$ para cada sector y `Channel.superoperator()` la matriz $`S_c`$
de $`\Phi_{cc}`$ sobre vectores ordenados por filas, $`\mathrm{vec}(X)_{(a,b)}=X_{ab}`$. Ambas
contienen los números $`\Phi(E_{ij})_{ab}=J_{(i,a),(j,b)}=S_{(a,b),(i,j)}`$ en distinto orden.
`Channel.from_superoperator` lanza `ValueError` si una matriz de Choi no es hermítica,
$`\|J-J^*\|>2\cdot10^{-10}s`$ (en doble precisión), donde $`s=\max(1,\|J\|)`$ con la norma espectral; en
caso contrario sustituye $`J`$ por su parte hermítica $`(J+J^*)/2`$, la diagonaliza como en la
demostración del [teorema 6.2](../../lecture_notes/es/06_quantum_channels.md#r-6-2), lanza `ValueError` si un valor propio es menor que $`-10^{-10}s`$, y pone
a $`0`$ los valores propios negativos restantes. Como $`J^*=J`$ si y solo si $`\Phi(X^*)=\Phi(X)^*`$ para todo $`X`$
([apuntes de clase, sección 6.1](../../lecture_notes/es/06_quantum_channels.md#s-6-1)), salvo estas tolerancias de redondeo rechaza exactamente las
aplicaciones que no preservan los adjuntos y, entre las demás, exactamente las que no son CP; por
ejemplo, se rechazan $`X\mapsto iX`$ y la trasposición.

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
try:
    Channel.from_superoperator(alg, [1j * torch.eye(k * k, dtype=dt)[None] for k in (2, 3)])
    raise AssertionError("X -> iX was accepted")
except ValueError as e:
    print(e)                                              # superoperator does not preserve Hermiticity
```

<a id="s-6-2"></a>

## 6.2 Dilatación de Stinespring y libertad en los operadores de Kraus

*Apuntes:* [sección 6.2](../../lecture_notes/es/06_quantum_channels.md#s-6-2): [teorema 6.5](../../lecture_notes/es/06_quantum_channels.md#r-6-5) (Stinespring) y [teorema 6.7](../../lecture_notes/es/06_quantum_channels.md#r-6-7) (libertad unitaria). Para el bloque $`c=1`$ del canal de la sección 6.1, el código comprueba que $`V`$ es una isometría, que $`\Phi(\rho)=\mathrm{Tr}_E(V\rho V^*)`$ y $`\Phi^*(B)=V^*(1\otimes B)V`$, y que los operadores de Kraus $`\sum_lu_{jl}K_l`$ dan el mismo canal.

En el código, el índice del entorno de $`V`$ es el primer índice del tensor de Kraus:

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

<a id="s-6-3"></a>

## 6.3 Duales respecto de una traza

*Apuntes:* [sección 6.3](../../lecture_notes/es/06_quantum_channels.md#s-6-3): el dual $`\Phi^\sharp`$ respecto de $`(w,v)`$, [proposición 6.8](../../lecture_notes/es/06_quantum_channels.md#r-6-8) y [corolario 6.9](../../lecture_notes/es/06_quantum_channels.md#r-6-9).

En la biblioteca, `Channel.adjoint()` es $`\Phi^*`$, e `InterSectorChannel.adjoint(trace=...)` es
$`\Phi^\sharp`$ para la misma traza con nombre en ambas álgebras (`'blunt'`, `'norm'`, `'tau_vN'`; los
pesos $`1/(Ck_c)`$ y $`1/(C'k'_d)`$ de `tau_vN` también tienen en cuenta el número de sectores).
Comprobamos la [proposición 6.8](../../lecture_notes/es/06_quantum_channels.md#r-6-8) para una aplicación aleatoria $`M_2\oplus M_3\to M_1\oplus M_2\oplus M_2`$
y las tres trazas, y el [corolario 6.9(a)](../../lecture_notes/es/06_quantum_channels.md#r-6-9) para la pérdida de partículas $`\Phi(x)=(x_0+\mathrm{Tr}\,x_1)\oplus0`$
sobre $`M_1\oplus M_2`$, que es TP para las trazas blunt pero no para `Tr_norm`
([apuntes de clase, sección 6.3](../../lecture_notes/es/06_quantum_channels.md#s-6-3)):

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

<a id="s-6-4"></a>

## 6.4 Decoherencia

*Apuntes:* [sección 6.4](../../lecture_notes/es/06_quantum_channels.md#s-6-4): [teorema 6.10](../../lecture_notes/es/06_quantum_channels.md#r-6-10) (puntos fijos de los canales unitales), [teorema 6.11](../../lecture_notes/es/06_quantum_channels.md#r-6-11) (los canales unitales aumentan la entropía) y [proposición 6.12](../../lecture_notes/es/06_quantum_channels.md#r-6-12) (la matriz de transición).

Los canales estándar de la biblioteca actúan dentro de cada bloque ($`p,\gamma\in[0,1]`$):

| canal | $`\Phi(\rho)_c`$ | operadores de Kraus en el bloque $`c`$ |
|---|---|---|
| `dephasing_channel` | $`(1-p)\rho_c+p\mathrm{diag}\,\rho_c`$ | $`\sqrt{1-p}\,1`$, $`\sqrt p\,E_{ii}`$ |
| `depolarizing_channel` | $`(1-p)\rho_c+p\mathrm{Tr}(\rho_c)\,1/k_c`$ | $`\sqrt{1-p}\,1`$, $`\sqrt{p/k_c}\,E_{ij}`$ |
| `center_expectation` | $`\mathrm{Tr}(\rho_c)\,1/k_c`$ | $`\sqrt{1/k_c}\,E_{ij}`$ |
| `amplitude_damping_channel` | los niveles $`i>0`$ decaen a $`0`$ | $`E_{00}+\sqrt{1-\gamma}\sum_{i>0}E_{ii}`$, $`\sqrt\gamma\,E_{0i}`$ |
| `random_mixed_unitary_channel` | $`\sum_jp_jU_j\rho U_j^*`$ | $`\sqrt{p_j}\,U_j`$ |

Para cada canal de la tabla sobre $`M_2\oplus M_3`$, el código calcula la dimensión de sus puntos fijos (el
núcleo de $`S_c-1`$) y la del conmutante de sus operadores de Kraus: coinciden para los canales unitales
([teorema 6.10](../../lecture_notes/es/06_quantum_channels.md#r-6-10)), y para el amortiguamiento de amplitud tienen la misma dimensión pero son espacios
distintos. Después comprueba el [teorema 6.11](../../lecture_notes/es/06_quantum_channels.md#r-6-11) en $`100`$ estados aleatorios, y la [proposición 6.12](../../lecture_notes/es/06_quantum_channels.md#r-6-12) para
la aplicación `Psi` de la sección 6.3 y para la pérdida de partículas.

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
assert fixed_dims(zoo['amplitude damping']) == commutant_dims(zoo['amplitude damping']) == [1, 1]

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

Bajo el amortiguamiento de amplitud completo, la entropía del estado máximamente mixto de
$`M_2\oplus M_3`$ baja de $`\log5\approx1.609`$ a $`H(\tfrac25,\tfrac35)\approx0.673`$: un canal no unital
puede purificar.

Sobre $`M_2`$, el desfase $`D_p`$ y el amortiguamiento de amplitud $`A_\gamma`$ conmutan, y su composición
multiplica las coherencias por $`(1-p)\sqrt{1-\gamma}`$; con $`1-p=e^{-t/T_\varphi}`$ y $`1-\gamma=e^{-t/T_1}`$ esto da
$`1/T_2=1/(2T_1)+1/T_\varphi`$ ([apuntes de clase, sección 6.4](../../lecture_notes/es/06_quantum_channels.md#s-6-4), *Coherencias de un cúbit*). El código
comprueba estas identidades y las leyes de composición $`D_{p(t)}D_{p(s)}=D_{p(t+s)}`$, $`A_{\gamma(t)}A_{\gamma(s)}=A_{\gamma(t+s)}`$:

```python
q2 = TypeIAlgebra([2], [2], precision='double')
p_, g_ = 0.3, 0.4
Dp, Ag = channels.dephasing_channel(q2, p_), channels.amplitude_damping_channel(q2, g_)
r2 = states.random_density_matrix(q2)
assert close(Dp(Ag(r2)), Ag(Dp(r2)))                                         # D_p A_gamma = A_gamma D_p
coh = Dp(Ag(r2)).matrix[0, 0, 0, 1] / r2.matrix[0, 0, 0, 1]
assert abs(coh - (1 - p_) * math.sqrt(1 - g_)) < 1e-12                       # coherence factor
T1, Tphi, t1, t2 = 2.0, 3.0, 0.7, 1.1
gam = lambda t: 1 - math.exp(-t / T1)
assert close(channels.amplitude_damping_channel(q2, gam(t1))(channels.amplitude_damping_channel(q2, gam(t2))(r2)),
             channels.amplitude_damping_channel(q2, gam(t1 + t2))(r2))       # A_gamma(t) A_gamma(s) = A_gamma(t+s)
pdeph = lambda t: 1 - math.exp(-t / Tphi)
assert close(channels.dephasing_channel(q2, pdeph(t1))(channels.dephasing_channel(q2, pdeph(t2))(r2)),
             channels.dephasing_channel(q2, pdeph(t1 + t2))(r2))           # D_p(t) D_p(s) = D_p(t+s)
both = channels.dephasing_channel(q2, pdeph(t1))(channels.amplitude_damping_channel(q2, gam(t1))(r2))
coh_t = (both.matrix[0, 0, 0, 1] / r2.matrix[0, 0, 0, 1]).item()
assert abs(coh_t - math.exp(-t1 * (1 / (2 * T1) + 1 / Tphi))) < 1e-12       # 1/T_2 = 1/(2T_1) + 1/T_phi
```

<a id="s-6-5"></a>

## 6.5 Canales aleatorios

*Apuntes:* [sección 6.5](../../lecture_notes/es/06_quantum_channels.md#s-6-5): [lema 6.13](../../lecture_notes/es/06_quantum_channels.md#r-6-13) y [teorema 6.15](../../lecture_notes/es/06_quantum_channels.md#r-6-15) (pureza de salida de un canal aleatorio).

Una **isometría Haar-aleatoria** $`V:\mathbb C^k\to\mathbb C^N`$ es, como en la biblioteca
(`random_isometry`), el factor $`Q`$ de $`G=QR`$, donde $`G`$ es una matriz de Ginibre $`N\times k`$
(entradas independientes con densidad $`\pi^{-1}e^{-|z|^2}`$), $`Q^*Q=1`$, y $`R`$ es triangular superior
con diagonal positiva (en la CPU la biblioteca toma una factorización de LAPACK $`G=Q'R'`$ y multiplica la
$`i`$-ésima columna de $`Q'`$ por $`R'_{ii}/|R'_{ii}|`$, lo que da esta factorización, como en el
teorema 5.3, siempre que todos los $`R'_{ii}\ne0`$, es decir, siempre que $`G`$ tenga rango $`k`$, lo que
ocurre casi seguramente por la demostración del lema 6.13; para lotes grandes en una GPU usa
Gram–Schmidt, que la da directamente). `random_channel(alg, r)` genera tales
$`V_c:\mathbb C^{k_c}\to\mathbb C^r\otimes\mathbb C^{k_c}`$ independientemente en cada bloque y usa los
operadores de Kraus del teorema 6.5.

Por el [teorema 6.15 de los apuntes de clase](../../lecture_notes/es/06_quantum_channels.md#r-6-15), para un estado puro $`\psi`$ en un sector de tamaño $`k`$ y $`\Phi`$
generado por `random_channel` con rango de Kraus $`r`$, $`\mathbb E\mathrm{Tr}\,\Phi(\psi)^2=(k+r)/(kr+1)`$ y
$`\mathbb E\,\Phi(\psi)=1_k/k`$. Estimamos ambas cantidades a partir de $`20\,000`$ muestras, con el error estándar de la media:

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

<a id="s-6-6"></a>

## 6.6 Con qué comparar

QuTiP (`to_choi`, `to_kraus`, `to_super`, `to_stinespring`), Qiskit (`qiskit.quantum_info`: las
clases `Kraus`, `Choi`, `SuperOp`, `Stinespring` y `random_quantum_channel`; modelos de ruido como
`amplitude_damping_error` en Qiskit Aer), toqito (`toqito.channel_ops`, `toqito.channel_props`) y
QETLAB (`ChoiMatrix`, `KrausOperators`, `ApplyMap`, `RandomSuperoperator`) convierten entre las
representaciones de las secciones 6.1–6.2, comprueban la positividad completa y (la mayoría) generan
canales aleatorios, para aplicaciones $`M_n\to M_m`$ entre álgebras de matrices completas.
Sus convenios difieren (apilamiento por columnas o por filas, el orden de los factores de la matriz de
Choi), así que compare los convenios antes de comparar números. Un canal que preserva los sectores
sobre $`\mathcal M`$ puede tratarse allí como un canal sobre $`M_D`$, $`D=\sum_ck_c`$, con superoperadores
de tamaño $`D^2\times D^2`$ en lugar de $`k_c^2\times k_c^2`$ por sector. La biblioteca añade la propia
estructura de bloques: operadores de Kraus almacenados por sector y por lotes, aplicaciones entre
álgebras distintas (`InterSectorChannel`), los duales para cada una de las trazas del capítulo 3 y la
matriz de transición entre sectores.

<a id="x-6"></a>

## Ejercicios

[Ejercicios](../../lecture_notes/es/06_quantum_channels.md#x-6) de los apuntes de clase, numéricamente:

1. (Werner–Holevo) Compruebe con `Channel.from_superoperator` que $`\Phi(\rho)=\big(\mathrm{Tr}(\rho)1-\rho^T\big)/(k-1)`$
   sobre $`M_3`$ es un canal unital con rango de Kraus $`3`$.
2. Compruebe que el canal de desfase sobre $`M_2\oplus M_3`$ tiene rango de Kraus $`k_c`$ en el bloque $`c`$, mientras que la
   biblioteca almacena $`1+k_{\max}`$ operadores de Kraus.
3. Compruebe que el canal unitario mixto con los operadores de Kraus $`X^aZ^b/k`$ es `center_expectation` sobre $`M_3`$.
4. Compruebe $`(\Psi\circ\Phi)^\sharp=\Phi^\sharp\circ\Psi^\sharp`$ para las tres trazas de la biblioteca.

<details>
<summary>Soluciones</summary>

*Ejercicio 1.*

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

*Ejercicio 2.*

```python
ranks = [int((torch.linalg.eigvalsh(J_c[0]) > 1e-10).sum()) for J_c in zoo['dephasing'].choi()]
print(ranks, zoo['dephasing'].kraus_rank)                 # [2, 3] versus 4 stored operators
assert ranks == alg.k_factors and zoo['dephasing'].kraus_rank == 1 + max(alg.k_factors)
```

*Ejercicio 3.*

```python
Xs = torch.roll(torch.eye(k, dtype=dt), 1, dims=0)
Zs = torch.diag(torch.exp(2j * math.pi * torch.arange(k, dtype=torch.float64) / k)).to(dt)
W = [torch.linalg.matrix_power(Xs, a) @ torch.linalg.matrix_power(Zs, b) / k for a in range(k) for b in range(k)]
twirl = channels.kraus_channel([M3.operator(Wab[None, None]) for Wab in W])
r3op = states.random_density_matrix(M3)
assert close(twirl(r3op), channels.center_expectation(M3)(r3op))
```

*Ejercicio 4.*

```python
Psi2 = channels.random_inter_sector_channel(alg_out, alg, kraus_rank=2)
for t in trace:
    assert close((Psi2 @ Psi).adjoint(trace=t)(x), Psi.adjoint(trace=t)(Psi2.adjoint(trace=t)(x)))
```

</details>

**Siguiente:** [Tutorial 7](07_tensor_products.md): productos tensoriales, trazas parciales y entrelazamiento en la biblioteca.
