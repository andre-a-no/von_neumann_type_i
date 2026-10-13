[English](../09_spin_chains.md) · **Español** · [Français](../fr/09_spin_chains.md) · [Deutsch](../de/09_spin_chains.md) · [中文](../zh/09_spin_chains.md) · [日本語](../ja/09_spin_chains.md) · [Русский](../ru/09_spin_chains.md)

# 9. Cadenas de espines

Una cadena de $L$ espines $\frac12$ vive en $(\mathbb C^2)^{\otimes L}$, y el álgebra lineal densa en dimensión
$2^L$ se detiene en $L\approx14$. Las simetrías desplazan ese límite: la magnetización divide el espacio en
$L+1$ sectores (capítulo 4), las traslaciones dividen de nuevo cada sector y $SU(2)$ reduce los bloques una vez
más. Este capítulo construye estos sectores para la cadena XXZ, demuestra lo que dicen sobre los estados
fundamentales, contrasta la biblioteca con los fermiones libres y el ansatz de Bethe, y va más allá de las
matrices densas (métodos de Krylov) y de las cadenas limpias (desorden).

**Aprenderá** por qué el álgebra de una cadena es $\bigoplus_NM_{\binom LN}$ y por qué la torre de estas
álgebras tiene como diagrama de Bratteli el triángulo de Pascal e índice $4\cos^2\frac\pi{2L+4}$, cómo se
cuentan los sectores de momento, la regla de signos de Marshall y sus consecuencias para el momento y el espín
de los estados fundamentales (Lieb–Mattis), la solución de Jordan–Wigner de la cadena XX, las ecuaciones de
Bethe de dos magnones, por qué Lanczos converge monótonamente, y cómo medir la estadística de niveles y el
entrelazamiento de fermiones libres.

En todo el capítulo, los sitios son $0,\dots,L-1$, la base de $(\mathbb C^2)^{\otimes L}$ está formada por las
configuraciones $|b\rangle=|b_0\dots b_{L-1}\rangle$ con $b_i=1$ para el espín hacia arriba (la convención de la
biblioteca: el bit $i$ de una etiqueta es el sitio $i$), $s^a_i$ ($a=x,y,z$) son las matrices de espín
$\frac12\sigma_a$ en el sitio $i$, $s^\pm_i=s^x_i\pm is^y_i$, $n_i=s^z_i+\frac12$, $S^a=\sum_is^a_i$, y
$N=S^z+\frac L2$ es el número de espines hacia arriba. Una cadena es **abierta** (enlaces $(i,i+1)$, $i<L-1$) o
**periódica** (también el enlace $(L-1,0)$; entonces $L\ge3$ y los índices se toman mód $L$). $J>0$
(antiferromagnética), $\hbar=1$, y $\log$ es el logaritmo natural.

## 9.1 La cadena XXZ y sus sectores $U(1)$

La **cadena XXZ** con anisotropía $\Delta\in\mathbb R$ y campos $h_i\in\mathbb R$ es
$$
H=J\sum_{(i,j)}\Big[\tfrac12\big(s^+_is^-_j+s^-_is^+_j\big)+\Delta\,s^z_is^z_j\Big]+\sum_ih_is^z_i ,
$$
con la suma sobre los enlaces; $\frac12(s^+_is^-_j+s^-_is^+_j)=s^x_is^x_j+s^y_is^y_j$ desplaza un espín hacia
arriba a lo largo del enlace. Para $\Delta=1$ es la cadena de Heisenberg $J\sum\mathbf s_i\cdot\mathbf s_j$; para
$\Delta=0$, la cadena XX.

**Proposición 9.1.** $[H,S^z]=0$. Por tanto, $H\in U(1)'\cong\bigoplus_{N=0}^LM_{\binom LN}(\mathbb C)$, donde el
bloque $N$ actúa sobre el subespacio generado por las configuraciones con $N$ espines hacia arriba.

*Demostración.* Los términos diagonales conmutan con la matriz diagonal $S^z$. Para $i\neq j$, $s^+_is^-_j$ conmuta
con $s^z_l$, $l\ne i,j$, y $[s^z_i+s^z_j,s^+_is^-_j]=[s^z_i,s^+_i]s^-_j+s^+_i[s^z_j,s^-_j]=s^+_is^-_j-s^+_is^-_j=0$ por
$[s^z,s^\pm]=\pm s^\pm$ (sección 4.5); para $s^-_is^+_j$ se toman adjuntos. El resto es la proposición 4.10(2). $\square$

`SpinChain(L).xxz` construye los bloques directamente a partir de las configuraciones de cada sector
(`chain.labels[N]`). Los comparamos con la matriz completa $2^L\times2^L$ construida mediante productos de
Kronecker:

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

**La torre de cadenas.** Sea $\mathcal A_L=\bigoplus_NM_{\binom LN}$, que actúa sobre $L$ sitios. Añadir un sitio,
$x\mapsto x\otimes1_2$, aplica $\mathcal A_L$ en $\mathcal A_{L+1}$, porque $x\otimes1$ conmuta con
$S^z\otimes1+1\otimes s^z$. Esto da una torre $\mathcal A_0\subset\mathcal A_1\subset\mathcal A_2\subset\cdots$.

**Proposición 9.2.** La inclusión $\mathcal A_L\subset\mathcal A_{L+1}$ tiene la matriz de inclusión
$\Lambda_L$ de tamaño $(L+2)\times(L+1)$ con $\lambda_{N,N}=\lambda_{N,N-1}=1$ y todas las demás entradas $0$; los
tamaños de los bloques obedecen la regla de Pascal $\binom{L+1}N=\binom LN+\binom L{N-1}$, el diagrama de Bratteli
de la torre es el triángulo de Pascal y el grafo bipartito de $\Lambda_L$ es un camino con $2L+3$ vértices.

*Demostración.* Las configuraciones con $N$ espines hacia arriba en $L+1$ sitios son las que tienen $N$ espines
hacia arriba en los primeros $L$ sitios y el nuevo espín hacia abajo, y las que tienen $N-1$ y el nuevo espín
hacia arriba. Sobre estos dos subespacios, $x\otimes1$ actúa como $x_N$ y $x_{N-1}$, de modo que el bloque $N$ de
$\mathcal A_{L+1}$ contiene los bloques $N$ y $N-1$ de $\mathcal A_L$ una vez cada uno; las dimensiones suman según
la regla de Pascal. En el grafo bipartito, el sector $N$ de $\mathcal A_{L+1}$ está unido a $N-1$ y $N$ de
$\mathcal A_L$, lo que da el camino $0'-0-1'-1-\cdots-L'-L-(L+1)'$ (con primas para $\mathcal A_{L+1}$). $\square$

Un camino desde la cima del diagrama de Bratteli hasta el sector $N$ del nivel $L$ es una configuración $b$; hay
$\binom LN$ de ellos, el tamaño del bloque. Para una inclusión conexa con matriz de inclusión $\Lambda$ tomamos
como **índice** $\|\Lambda\|^2$, el índice de la traza de Markov (V. Jones, *Invent. Math.* 72 (1983);
F. Goodman, P. de la Harpe, V. Jones, *Coxeter graphs and towers of algebras*, Springer 1989, cap. 2–3); para
$M_k\otimes1\subset M_k\otimes M_m$ es el $m^2$ de la sección 7.2. `Inclusion.index` lo calcula.

**Lema 9.3.** La matriz de adyacencia $A_n$ del camino con vértices $1,\dots,n$ tiene los valores propios
$2\cos\frac{\pi k}{n+1}$, $k=1,\dots,n$, con vectores propios $v_k(j)=\sin\frac{\pi jk}{n+1}$; por tanto,
$\|A_n\|=2\cos\frac\pi{n+1}$.

*Demostración.* Pongamos $v_k(0)=v_k(n+1)=0$. Entonces $(A_nv_k)(j)=v_k(j-1)+v_k(j+1)=2\cos\frac{\pi k}{n+1}\,v_k(j)$
por $\sin(a-b)+\sin(a+b)=2\sin a\cos b$, y $v_k\ne0$ porque $v_k(1)=\sin\frac{\pi k}{n+1}\neq0$. Los $n$ valores
propios son distintos ($\cos$ es estrictamente decreciente en $(0,\pi)$), luego son todos, y la norma de la matriz
simétrica $A_n$ es el mayor módulo, alcanzado en $k=1$ y $k=n$. $\square$

**Corolario 9.4.** $[\mathcal A_{L+1}:\mathcal A_L]=\|\Lambda_L\|^2=4\cos^2\frac\pi{2L+4}$, que crece hacia $4$.

*Demostración.* Ordenando los vértices a lo largo del camino de la proposición 9.2, la matriz de adyacencia
$B=\begin{pmatrix}0&\Lambda_L\\\Lambda_L^T&0\end{pmatrix}$ se convierte en $A_{2L+3}$ (una semejanza por permutación). Por el
lema 1.1, $\|B\|^2=\|B^2\|$ y $\|\Lambda_L\Lambda_L^T\|=\|\Lambda_L^T\Lambda_L\|=\|\Lambda_L\|^2$, y $B^2=\Lambda_L\Lambda_L^T\oplus\Lambda_L^T\Lambda_L$, así que
$\|B\|=\|\Lambda_L\|$; el lema 9.3 con $n=2L+3$ da el valor. $\square$

Los números $4\cos^2\frac\pi n$, $n\ge3$, son exactamente los valores menores que $4$ que puede tomar el índice de
un subfactor (Jones, *loc. cit.*); aquí aparecen los $n$ pares por la razón elemental del lema 9.3.

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

> **Significado físico.** $S^z$ se conserva siempre que la interacción sea invariante bajo rotaciones alrededor
> del eje $z$; $N$ cuenta los espines hacia arriba, o las partículas tras la transformación de Jordan–Wigner de
> la sección 9.4. En la torre, cada sector de carga de la cadena más larga se alimenta de dos sectores de la
> más corta.

## 9.2 Traslaciones, momento y la regla de signos de Marshall

En una cadena periódica, sea $T$ la traslación $T|b_0b_1\dots b_{L-1}\rangle=|b_{L-1}b_0\dots b_{L-2}\rangle$, que lleva
el sitio $i$ al sitio $i+1$. $T$ es una matriz de permutación, luego unitaria, $T^L=1$ y $Ts^a_iT^*=s^a_{i+1}$. Por
tanto, $TS^zT^*=S^z$ y $THT^*=H$ para un campo uniforme $h_i=h$: $T$ lleva enlaces a enlaces. Sus valores propios
son $e^{ik}$ con **momentos** $k=2\pi m/L$, $m=0,\dots,L-1$.

**Proposición 9.5 (sectores de momento).** Fijemos $N$. Para una configuración $a$ con $N$ espines hacia arriba,
sea $R_a$ el tamaño de su órbita $\{T^ja\}$. (a) $R_a$ divide a $L$; para $k=2\pi m/L$ con $mR_a\equiv0\pmod L$, el
vector $|a,k\rangle=\sum_{j=0}^{R_a-1}e^{-ikj}T^j|a\rangle$ es no nulo y $T|a,k\rangle=e^{ik}|a,k\rangle$, y estos $R_a$
vectores forman una base ortogonal del subespacio generado por la órbita. (b) El valor propio $e^{2\pi im/L}$ de
$T$ sobre el sector $N$ tiene multiplicidad $d_{N,m}=\#\{\text{órbitas }:mR_a\equiv0\pmod L\}$, y todo operador que
conmute con $S^z$ y $T$ —por ejemplo, la cadena XXZ invariante por traslaciones— es diagonal por bloques, en
$\bigoplus_N\bigoplus_mM_{d_{N,m}}$.

*Demostración.* (a) El estabilizador de $a$ en el grupo cíclico $\{T^j\}\cong\mathbb Z_L$ es un subgrupo, y $R_a$ es
su índice, un divisor de $L$; los vectores $T^j|a\rangle$, $0\le j<R_a$, son ortonormales y $T^{R_a}|a\rangle=|a\rangle$. Si
$e^{ikR_a}=1$, entonces $T|a,k\rangle=e^{ik}\sum_{j=1}^{R_a}e^{-ikj}T^j|a\rangle=e^{ik}|a,k\rangle$, ya que el término $j=R_a$
es igual al término $j=0$; el coeficiente de $|a\rangle$ es $1$. Los $m\in\{0,\dots,L-1\}$ admisibles son los
múltiplos de $L/R_a$, $R_a$ en total, con valores propios distintos del unitario $T$, de modo que los vectores son
ortogonales y generan el subespacio de dimensión $R_a$ generado por la órbita. (b) El sector $N$ es la suma
ortogonal de los subespacios generados por sus órbitas. Un $X$ que conmuta con $S^z$ y $T$ conmuta con sus
proyecciones espectrales $P_N$, $Q_m$ (teorema 1.5(d)), así que $X=\sum_{N,m}P_NQ_mXP_NQ_m$. $\square$

`xxz_momentum` construye los bloques en la base normalizada $|a,k\rangle$, con un representante $a$ por órbita.
Contamos las órbitas de forma independiente y comparamos los espectros con el sector denso:

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

El sector $N=6$, de $924$ estados, se divide en bloques de tamaños $80,75,78,76,78,75,80,75,78,76,78,75$.

**Lema 9.6 (Perron–Frobenius, caso simétrico).** Sea $M$ una matriz real simétrica con $M_{ij}\ge0$ para $i\ne j$
cuyo grafo ($i\sim j$ si y solo si $M_{ij}>0$) es conexo. Entonces su mayor valor propio $\lambda$ es simple, con un
vector propio cuyas entradas son todas $>0$.

*Demostración.* $\lambda=\max_{\|x\|=1}\langle x,Mx\rangle$, y el máximo se alcanza exactamente en los vectores unitarios
del espacio propio de $\lambda$ (desarrollando en una base ortonormal de vectores propios). Para tal $x$, el vector
$|x|$ de los módulos cumple $\langle|x|,M|x|\rangle\ge\operatorname{Re}\langle x,Mx\rangle=\lambda$, porque los términos
diagonales coinciden y $M_{ij}|x_i||x_j|\ge\operatorname{Re}M_{ij}\bar x_ix_j$ para $i\ne j$; así que $|x|$ es un vector
propio. Si $|x|_i=0$, entonces $0=(\lambda-M_{ii})|x|_i=\sum_{j\ne i}M_{ij}|x|_j$ obliga a $|x|_j=0$ para todos los
vecinos $j$ de $i$, y por conexión $x=0$. Así, ningún vector propio no nulo tiene una entrada nula; pero un espacio
propio de dimensión $\ge2$ contiene un vector no nulo con $x_1=0$. Por tanto, $\lambda$ es simple y $|x|>0$ es un
vector propio. $\square$

**Teorema 9.7 (regla de signos de Marshall).** Sea la cadena abierta, o periódica con $L$ par; sean $B$ los sitios
impares y $N_B(b)$ el número de espines hacia arriba de $b$ en $B$; sean $J>0$ y $\Delta,h_i\in\mathbb R$
arbitrarios. En cada sector $N$, el menor valor propio de $H_N$ es simple, y su vector propio puede elegirse con
$(-1)^{N_B(b)}\psi(b)>0$ para toda configuración $b$ con $N$ espines hacia arriba.

*Demostración.* Sea $V|b\rangle=(-1)^{N_B(b)}|b\rangle$; $V$ es diagonal y $V=V^*=V^{-1}$. Todo enlace une un sitio par
con uno impar (en el anillo porque $L$ es par), así que $s^+_is^-_j$ cambia $N_B$ en $\pm1$ y
$Vs^+_is^-_jV=-s^+_is^-_j$. Así, $M=-VH_NV$ tiene entradas fuera de la diagonal $J/2>0$ entre configuraciones que
difieren en el intercambio de los dos espines de un enlace, y $0$ en otro caso. Su grafo es conexo: los
intercambios de vecinos $(i,i+1)$, $i<L-1$, generan todas las permutaciones de los sitios y, por tanto, conectan
dos configuraciones cualesquiera con $N$ espines hacia arriba (intercambiar espines iguales no hace nada). Por el
lema 9.6, el mayor valor propio de $M$ —el opuesto del menor de $H_N$— es simple, con un vector propio
$\varphi>0$; tómese $\psi=V\varphi$. $\square$

**Corolario 9.8 (momento del estado fundamental).** En una cadena periódica con $L$ par, $J>0$ y un campo
uniforme, el estado fundamental del sector $N$ tiene momento $k=0$ si $N$ es par y $k=\pi$ si $N$ es impar.

*Demostración.* $T$ conmuta con $H_N$ y el estado fundamental $\psi=V\varphi$ es simple, así que $T\psi=e^{ik}\psi$.
Ahora bien, $(T\psi)(b)=\psi(T^{-1}b)$, y $T^{-1}b$ tiene en el sitio $i$ el espín $b_{i+1}$; como $L$ es par, $i+1$
es par para $i$ impar, luego $N_B(T^{-1}b)=N-N_B(b)$ y
$(VT\psi)(b)=(-1)^{N_B(b)}(-1)^{N-N_B(b)}\varphi(T^{-1}b)=(-1)^N\varphi(T^{-1}b)$. Además, $VT\psi=e^{ik}\varphi$. Tanto
$\varphi$ como $\varphi\circ T^{-1}$ son positivos, así que $e^{ik}=(-1)^N$. $\square$

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

> **Significado físico.** En un antiferromagneto, la red fija los signos del estado fundamental: mover un espín a
> través de un enlace cambia el signo. Esto libra al Monte Carlo cuántico del problema del signo para
> antiferromagnetos bipartitos; falla en las cadenas frustradas, en las que un acoplamiento une sitios de la
> misma subred.

## 9.3 $SU(2)$ en $\Delta=1$

Para $\Delta=1$ y $h=0$, $H=J\sum\mathbf s_i\cdot\mathbf s_j$ conmuta con $S^x,S^y,S^z$ (proposición 4.15), así que
$H\cong\bigoplus_S1_{2S+1}\otimes H_S$ con bloques $H_S\in M_{m_S}$, $m_S=\binom L{L/2-S}-\binom L{L/2-S-1}$
(proposición 4.13), y el módulo `su2` construye $H_S$ sobre los vectores de peso máximo (proposición 4.15(b)).
Sea $E_0(S)$ el menor valor propio de $H_S$.

**Teorema 9.9 (ordenación de Lieb–Mattis para cadenas).** Sea la cadena como en el teorema 9.7, $\Delta=1$,
$h_i=0$, $J>0$. Para todo espín total $S$ de $L$ espines, el estado fundamental del sector $S^z=S$ tiene espín
total $S$. Por tanto, $E_0(S)$ es un valor propio simple de $H_S$, y $E_0(S_{\min})<E_0(S_{\min}+1)<\dots<E_0(\frac L2)$
con $S_{\min}\in\{0,\frac12\}$: el estado fundamental es un único multiplete del menor espín.

*Demostración.* Sea $\psi$ el estado fundamental del sector $S^z=S$, simple por el teorema 9.7. $\mathbf S^2$
conmuta con $H$ y preserva el sector, así que $\mathbf S^2\psi$ está en el espacio fundamental unidimensional: $\psi$
tiene un espín total $S'$. Sea $\chi$ el producto de singletes $\frac1{\sqrt2}(|{\uparrow\downarrow}\rangle-|{\downarrow\uparrow}\rangle)$
sobre los pares $(0,1),(2,3),\dots,(L-2S-2,L-2S-1)$ con los $2S$ sitios restantes hacia arriba. $S^+$ anula cada
factor y $S^z\chi=S\chi$, así que $\mathbf S^2\chi=(S^-S^++(S^z)^2+S^z)\chi=S(S+1)\chi$ (sección 4.5). En cada par, el
espín hacia arriba está en el sitio par (coeficiente $+$) o en el impar ($-$), de modo que $(-1)^{N_B(b)}\chi(b)$
tiene el mismo signo para todos los $b$ de su soporte; con $V,\varphi$ como en el teorema 9.7,
$\langle\psi,\chi\rangle=\langle\varphi,V\chi\rangle\neq0$. Los vectores propios de $\mathbf S^2$ para valores propios
distintos son ortogonales, así que $S'=S$. Por la proposición 4.15(a) y el corolario 4.12, cada multiplete de
espín $S'\ge S$ corta al sector $S^z=S$ en un vector, de modo que el espectro de $H$ allí es la unión de los
espectros de los $H_{S'}$, $S'\ge S$, con multiplicidades. Su menor valor propio es simple y su vector propio tiene
espín $S$, así que es $E_0(S)$, simple en $H_S$, y no es valor propio de $H_{S'}$, $S'>S$: $E_0(S)<E_0(S')$. $\square$

Lieb y Mattis demostraron más: para una red bipartita conexa con espines de las subredes $S_A,S_B$ y acoplamientos
antiferromagnéticos entre las subredes y ferromagnéticos dentro de ellas, $E_0(S)$ crece para $S\ge|S_A-S_B|$ y
el estado fundamental tiene espín $|S_A-S_B|$ (E. Lieb, D. Mattis, *J. Math. Phys.* 3 (1962) 749; enunciado sin
demostración). La **brecha singlete–triplete** $E_0(1)-E_0(0)$ es positiva por el teorema 9.9; en el anillo,
$L(E_0(1)-E_0(0))\to\pi^2/2$ con correcciones logarítmicas (J. des Cloizeaux, J. Pearson, *Phys. Rev.* 128 (1962)
2131; I. Affleck, D. Gepner, H. Schulz, T. Ziman, *J. Phys. A* 22 (1989) 511; enunciado sin demostración).
Comparamos los tamaños de los bloques para doce espines, comprobamos el teorema 9.9 frente a los sectores $U(1)$
y calculamos la brecha:

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

Los bloques de $SU(2)$ tienen tamaños $m_S=132,297,275,154,54,11,1$ frente a $924$ del mayor bloque de $U(1)$, y
el coste $\sum k^3$ de la diagonalización baja de $2.0\cdot10^9$ a $5.3\cdot10^7$. $L$ veces la brecha crece
lentamente, $4.108,4.181,4.232,4.270$ para $L=6,\dots,12$, hacia $\pi^2/2\approx4.935$.

**Correlaciones.** Los espines de los sitios forman un operador tensorial de rango uno y, por el **teorema de
Wigner–Eckart**, sus elementos de matriz entre multipletes son un coeficiente de Clebsch–Gordan por un
**elemento de matriz reducido** $\langle S',a\|\mathbf s_i\|S,b\rangle$ independiente de los números cuánticos
magnéticos (A. R. Edmonds, *Angular Momentum in Quantum Mechanics*, Princeton 1957, cap. 5; enunciado sin
demostración). `spin_correlation` lo usa para construir $\mathbf s_i\cdot\mathbf s_j$ sobre cada $H_S$ a partir
únicamente de los vectores de peso máximo; el método `'direct'` proyecta en cambio el operador de enlace.
Comparamos ambos con el estado fundamental escrito en el sector $U(1)$:

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

Las correlaciones $\langle\mathbf s_0\cdot\mathbf s_r\rangle=0.75,-0.4489,0.1878,-0.1658,0.1209,-0.1224,0.1070,\dots$
alternan de signo, como imponen los teoremas 9.7 y 9.9: el singlete es invariante bajo rotaciones, que permutan
los ejes, y $\langle s^+_0s^+_r\rangle=0$ en el sector $S^z=0$, así que $\langle\mathbf s_0\cdot\mathbf s_r\rangle=3\langle s^x_0s^x_r\rangle=
\frac34\langle s^+_0s^-_r+s^-_0s^+_r\rangle=\frac34\sum_{b_0\ne b_r}\psi(b)\psi(b')$, donde $b'$ se obtiene de $b$
intercambiando los espines de los sitios $0$ y $r$; como $N_B(b')-N_B(b)$ vale $\pm1$ para $r$ impar y $0$ para $r$
par, cada término tiene el signo $(-1)^r$. El decaimiento es lento, como $(-1)^r\sqrt{\log r}/r$ cuando
$L\to\infty$ (Affleck et al., *loc. cit.*; enunciado sin demostración).

## 9.4 Comprobaciones exactas: fermiones libres y el ansatz de Bethe

**Jordan–Wigner.** En una cadena abierta, pongamos $c_j=\big(\prod_{l<j}(-1)^{n_l}\big)s^-_j$.

**Teorema 9.10 (la cadena XX es libre).** (a) $\{c_i,c_j^*\}=\delta_{ij}$, $\{c_i,c_j\}=0$ y $c_j^*c_j=n_j$.
(b) La cadena XX abierta es $H=\frac J2\sum_{j<L-1}(s^+_js^-_{j+1}+s^-_js^+_{j+1})=\sum_{i,j}h_{ij}c_i^*c_j$ con
$h=\frac J2A_L$. (c) Sean $u_k$ vectores propios reales ortonormales de $h$ con valores propios
$\varepsilon_k=J\cos\frac{\pi k}{L+1}$ (lema 9.3) y $b_k=\sum_ju_k(j)c_j$. Para $K=\{k_1<\dots<k_N\}$, los vectores
$b_K^*|0\rangle=b_{k_1}^*\cdots b_{k_N}^*|0\rangle$ ($|0\rangle$ con todos los espines hacia abajo) forman una base
ortonormal de vectores propios del sector $N$, con energías $\sum_{k\in K}\varepsilon_k$.

*Demostración.* (a) Los operadores en sitios distintos conmutan, y en un mismo sitio $(-1)^n$ anticonmuta con
$s^\pm$. Para $i<j$, $c_i$ y $c_j^*$ (o $c_j$) se encuentran en el sitio $i$ como $s^-_i$ y $(-1)^{n_i}$ —una
anticonmutación— y por lo demás conmutan, así que anticonmutan. Para $i=j$, las cadenas de signos tienen cuadrado
$1$ y $c_jc_j^*+c_j^*c_j=s^-_js^+_j+s^+_js^-_j=1$, $c_j^*c_j=s^+_js^-_j=n_j$. (b)
$c_j^*c_{j+1}=s^+_j(-1)^{n_j}s^-_{j+1}=s^+_js^-_{j+1}$, porque $s^+_j\neq0$ solo sobre $n_j=0$; se toman adjuntos.
(c) $\{b_k,b_l^*\}=\sum_ju_k(j)u_l(j)=\delta_{kl}$ y $\{b_k,b_l\}=0$, y $H=\sum_k\varepsilon_kb_k^*b_k$ por la
descomposición espectral de $h$. De las CAR, $[b_l^*b_l,b_k^*]=\delta_{lk}b_k^*$, así que $[H,b_k^*]=\varepsilon_kb_k^*$;
con $H|0\rangle=0$ esto da $Hb_K^*|0\rangle=(\sum_K\varepsilon_k)b_K^*|0\rangle$. Como $b_k|0\rangle=0$ (cada $c_j$ baja un
espín), llevar los aniquiladores a la derecha con las CAR da $\langle b_{K'}^*0,b_K^*0\rangle=\delta_{KK'}$ por inducción
sobre $N$. Cada $c_j^*$ aumenta $N$ en uno, así que $b_K^*|0\rangle$ está en el sector $|K|$: $\binom LN$ vectores
ortonormales en un espacio de esa dimensión. $\square$

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

**Magnones.** En el anillo con $h=0$, el estado totalmente polarizado $|0\rangle$ tiene energía $E_F=J\Delta L/4$. Un
espín hacia arriba en $x$ rompe dos enlaces: $H|x\rangle=(E_F-J\Delta)|x\rangle+\frac J2(|x-1\rangle+|x+1\rangle)$ ($L\ge3$).

**Proposición 9.11 (un magnón).** En el sector $N=1$ del anillo, $|k\rangle=\sum_xe^{-ikx}|x\rangle$, $k=2\pi m/L$,
cumple $T|k\rangle=e^{ik}|k\rangle$ y $H|k\rangle=\big(E_F+J(\cos k-\Delta)\big)|k\rangle$.

*Demostración.* $\sum_xe^{-ikx}(|x-1\rangle+|x+1\rangle)=\sum_y(e^{-ik(y+1)}+e^{-ik(y-1)})|y\rangle=2\cos k\,|k\rangle$, usando
$e^{ikL}=1$; la misma reindexación da $T|k\rangle=e^{ik}|k\rangle$. $\square$

Dos magnones se dispersan. Escribamos $|x_1,x_2\rangle$, $x_1<x_2$, para los espines hacia arriba en $x_1,x_2$; un
par de vecinos rompe dos enlaces, un par separado cuatro.

**Proposición 9.12 (dos magnones; ecuaciones de Bethe).** Sean $k_1,k_2\in\mathbb C$, $K=k_1+k_2$, y
$(A_{12},A_{21})\neq0$ con (i) $A_{12}(1+e^{iK}-2\Delta e^{ik_2})+A_{21}(1+e^{iK}-2\Delta e^{ik_1})=0$ y
(ii) $A_{21}=e^{ik_2L}A_{12}$, $A_{12}=e^{ik_1L}A_{21}$. Sean
$\psi(x_1,x_2)=A_{12}e^{i(k_1x_1+k_2x_2)}+A_{21}e^{i(k_2x_1+k_1x_2)}$ sobre $\mathbb Z^2$ y
$\Psi=\sum_{0\le x_1<x_2<L}\psi(x_1,x_2)|x_1,x_2\rangle$. Si $\Psi\neq0$, entonces
$H\Psi=\big(E_F+J(\cos k_1+\cos k_2-2\Delta)\big)\Psi$. Si los denominadores no se anulan, (i) y (ii) dicen
$$
e^{ik_1L}=-\frac{1+e^{iK}-2\Delta e^{ik_1}}{1+e^{iK}-2\Delta e^{ik_2}},\qquad
e^{ik_2L}=-\frac{1+e^{iK}-2\Delta e^{ik_2}}{1+e^{iK}-2\Delta e^{ik_1}} .
$$

*Demostración.* Sea $E=E_F+J(\cos k_1+\cos k_2-2\Delta)$. (1) Para todo $(x_1,x_2)\in\mathbb Z^2$, cada onda plana da
valores $\psi(x_1\pm1,x_2)$ y $\psi(x_1,x_2\pm1)$ que suman $2(\cos k_1+\cos k_2)\psi(x_1,x_2)$. (2) Sustituyendo,
$\psi(x,x)+\psi(x+1,x+1)-2\Delta\psi(x,x+1)$ es $e^{iKx}$ por el miembro izquierdo de (i), así que se anula para todo
$x$. (3) Por (ii), $\psi(x_2,x_1+L)=A_{12}e^{ik_2L}e^{i(k_2x_1+k_1x_2)}+A_{21}e^{ik_1L}e^{i(k_1x_1+k_2x_2)}=\psi(x_1,x_2)$, y
aplicando esto dos veces, $\psi(x_1+L,x_2+L)=\psi(x_1,x_2)$. Los representantes de una configuración $\{y_1<y_2\}$
en $D=\{x_1<x_2<x_1+L\}$ son $(y_1+nL,y_2+nL)$ y $(y_2+nL,y_1+(n+1)L)$, así que $\psi$ es una función de las
configuraciones. Sea $(x_1,x_2)\in D$. Si $2\le x_2-x_1\le L-2$, los cuatro saltos permanecen en $D$ y el término
diagonal es $E_F-2J\Delta$, así que $(H\Psi)(x_1,x_2)=(E_F-2J\Delta)\psi+\frac J2\cdot2(\cos k_1+\cos k_2)\psi=E\psi$ por (1).
Si $x_2=x_1+1$, los saltos $(x_1+1,x_2)$, $(x_1,x_2-1)$ están bloqueados y el término diagonal es $E_F-J\Delta$; por
(1) y (2),
$(H\Psi)(x_1,x_2)=(E_F-J\Delta)\psi+\frac J2\big(2(\cos k_1+\cos k_2)\psi-\psi(x_1,x_1)-\psi(x_1+1,x_1+1)\big)=E\psi$.
Si $x_2=x_1+L-1$, se usa el representante $(x_2,x_1+L)$, que es del tipo anterior ($L\ge3$). $\square$

La proposición 9.12 solo dice que las soluciones dan vectores propios; no dice que todos los vectores propios
surjan así. Escribamos $k_{1,2}=K/2\pm q$ y $c=\cos\frac K2$. Entonces
$1+e^{iK}-2\Delta e^{ik_{1,2}}=2e^{iK/2}(c-\Delta e^{\pm iq})$, así que $A_{12}=c-\Delta e^{iq}$,
$A_{21}=-(c-\Delta e^{-iq})$ resuelven (i). Multiplicar las dos ecuaciones (ii) da $e^{iKL}=1$, $K=2\pi n/L$, y
entonces (ii) equivale a $F(q)=(-1)^ne^{iqL}(c-\Delta e^{-iq})+c-\Delta e^{iq}=0$, una ecuación en una variable.
Buscamos sus raíces $q\in(0,\pi)$, donde $e^{-iqL/2}F(q)$ es real o imaginario porque $c-\Delta e^{\pm iq}$ son
complejos conjugados (estados de dispersión), y sobre las rectas $q=i\kappa$, $\pi+i\kappa$, $\kappa>0$, donde $F$ es
real (estados ligados); para $L=9$, $\Delta=0.6$ construimos $\Psi$ para cada raíz y comparamos con el sector $N=2$
de $36$ estados:

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

Para $L=9$, los $36$ estados propios provienen de las ecuaciones de Bethe: $34$ estados de dispersión y $2$ estados
ligados. La completitud se comprueba aquí, no se demuestra, y para $L$ par falla. Si $\Delta\ne0$,
$\Phi=\sum_x(-1)^x|x,x+1\rangle$ cumple $T\Phi=-\Phi$ y $H\Phi=(E_F-J\Delta)\Phi$ (los saltos a $|x-1,x+1\rangle$ y
$|x,x+2\rangle$ se cancelan por pares). Un vector de Bethe cumple $T\Psi=e^{-iK}\Psi$, ya que
$\psi(x_1-1,x_2-1)=e^{-iK}\psi(x_1,x_2)$; si $e^{iK}=-1$, entonces
$\cos k_1+\cos k_2=2\cos\frac K2\cos\frac{k_1-k_2}2=0$ y su energía es $E_F-2J\Delta$. Así, $\Phi$ es ortogonal a todos
los vectores de Bethe (las últimas líneas del código comprueban $H\Phi$ para $L=8$). Con $\Delta=0.6$ y
$4\le L\le13$, la búsqueda anterior encuentra todos los estados para $L$ impar, y todos salvo $\Phi$ para $L$ par.
Para $N>2$ magnones, el ansatz tiene $N!$ términos y $N$ ecuaciones de Bethe acopladas (H. Bethe, *Z. Phys.* 71
(1931) 205; M. Karbach, G. Müller, *Computers in Physics* 11 (1997) 36).

> **Significado físico.** Un magnón es una onda de espín: un espín invertido deslocalizado sobre el anillo con
> momento $k$. Dos magnones se dispersan con un desfase (con $q$ real) o forman un estado ligado (con $q$ complejo).
> La intuición de Bethe fue que, para $N$ magnones, la dispersión se factoriza en fases de dos cuerpos; sus
> ecuaciones dan, por ejemplo, la energía del estado fundamental por sitio $\frac14-\log2$ de la cadena de
> Heisenberg infinita (L. Hulthén, *Ark. Mat. Astron. Fys.* 26A (1938) n.º 11; enunciado sin demostración).

## 9.5 Lanczos y métodos de Krylov

El sector $N=L/2$ tiene dimensión $12870$ para $L=16$ y $184756$ para $L=20$, pero $H$ es dispersa ($O(L)$ entradas
no nulas por fila, `xxz_sparse`), y los métodos de Krylov solo necesitan productos $Hv$. Sea $H$ hermítica de
tamaño $n$ y $v\ne0$. El **espacio de Krylov** es $\mathcal K_m=\operatorname{span}\{v,Hv,\dots,H^{m-1}v\}$.

**Proposición 9.13 (Lanczos).** Supongamos $\dim\mathcal K_m=m$. Sea $q_1=v/\|v\|$ y, para $j<m$, sea $q_{j+1}$ la
componente normalizada de $Hq_j$ ortogonal a $\mathcal K_j$. Entonces $q_1,\dots,q_m$ es una base ortonormal de
$\mathcal K_m$; con $Q_m=(q_1\cdots q_m)$, $T_m=Q_m^*HQ_m$ es real, simétrica y tridiagonal, con diagonal
$\alpha_j=\langle q_j,Hq_j\rangle$ y subdiagonal $\beta_j=\langle q_{j+1},Hq_j\rangle>0$; y $HQ_m=Q_mT_m+r_me_m^T$, donde
$r_m\perp\mathcal K_m$ se anula si y solo si $\mathcal K_{m+1}=\mathcal K_m$ (en otro caso $r_m=\beta_mq_{m+1}$,
$q_{m+1}=r_m/\|r_m\|$).

*Demostración.* Supongamos $\operatorname{span}\{q_1,\dots,q_j\}=\mathcal K_j$ (cierto para $j=1$). Entonces
$H\mathcal K_j\subseteq\mathcal K_j+\mathbb CHq_j$, pues $Hq_i\in\mathcal K_{i+1}$; así que
$\mathcal K_{j+1}=\mathcal K_j+\mathbb CH^jv\subseteq\mathcal K_j+\mathbb CHq_j\subseteq\mathcal K_{j+1}$, y
$\mathcal K_{j+1}=\mathcal K_j+\mathbb CHq_j$. Para $j<m$ la dimensión crece, de modo que la componente de $Hq_j$
ortogonal a $\mathcal K_j$ es no nula y la inducción continúa; para $j=m$ esa componente, $r_m$, se anula si y solo
si $\mathcal K_{m+1}=\mathcal K_m$. Para $i>j+1$, $q_i\perp\mathcal K_{j+1}\ni Hq_j$, así que $(T_m)_{ij}=0$, y
$(T_m)_{ji}=0$ porque $T_m$ es hermítica. $\alpha_j$ es real, y $\beta_j=\langle q_{j+1},Hq_j\rangle$ es la norma de la
componente, $>0$. Así, $Hq_j=\beta_{j-1}q_{j-1}+\alpha_jq_j+\beta_jq_{j+1}$ para $j<m$ y
$Hq_m=\beta_{m-1}q_{m-1}+\alpha_mq_m+r_m$. $\square$

Los valores propios $\theta^{(m)}_1\le\dots\le\theta^{(m)}_m$ de $T_m$ son los **valores de Ritz**.

**Teorema 9.14 (convergencia monótona).** Sean $\lambda_1\le\dots\le\lambda_n$ los valores propios de $H$ y
$\dim\mathcal K_m=m$. Entonces $\theta^{(m)}_j\ge\lambda_j$, y $\theta^{(m)}_j\ge\theta^{(m+1)}_j$ si
$\dim\mathcal K_{m+1}=m+1$ ($j\le m$). Si $v$ tiene una componente no nula en el espacio propio de $\lambda_1$ y
$\mathcal K_{m+1}=\mathcal K_m$, entonces $\theta^{(m)}_1=\lambda_1$.

*Demostración.* Por Courant–Fischer (álgebra lineal estándar),
$\theta^{(m)}_j=\min_{U}\max_{x\in U,\|x\|=1}\langle x,Hx\rangle$ sobre los $U\subseteq\mathcal K_m$ de dimensión $j$, ya que
$Q_m$ es una isometría sobre $\mathcal K_m$ y $\langle y,T_my\rangle=\langle Q_my,HQ_my\rangle$. La misma fórmula sobre
$U\subseteq\mathcal K_{m+1}$ da $\theta^{(m+1)}_j$ y sobre $U\subseteq\mathbb C^n$ da $\lambda_j$; minimizar sobre más
subespacios da un valor menor. Si $\mathcal K_{m+1}=\mathcal K_m$, entonces $\mathcal K_m$ es $H$-invariante y contiene
$p(H)v$ para todo polinomio $p$, en particular $P_1v\ne0$ ($P_1$ la proyección espectral de $\lambda_1$, un
polinomio en $H$ por interpolación de Lagrange, como en el teorema 1.5); así que
$\theta^{(m)}_1\le\langle x,Hx\rangle=\lambda_1$ en $x=P_1v/\|P_1v\|\in\mathcal K_m$. $\square$

La rapidez con que $\theta^{(m)}_1\to\lambda_1$ depende de la brecha $\lambda_2-\lambda_1$ en relación con la amplitud del
espectro (cotas de Kaniel–Paige–Saad mediante polinomios de Chebyshev; Y. Saad, *Numerical Methods for Large
Eigenvalue Problems*, SIAM 2011, cap. 6; enunciado sin demostración). En coma flotante, los $q_j$ pierden la
ortogonalidad; `krylov.lanczos` reortogonaliza por completo. Para la cadena XX se conoce la energía exacta del
estado fundamental (teorema 9.10):

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

El error del menor valor de Ritz cae monótonamente, $0.42$, $3.1\cdot10^{-3}$ y $3.1\cdot10^{-6}$ tras $5$, $10$ y
$20$ pasos, y alcanza el nivel de redondeo tras $40$, en dimensión $924$. `krylov.ground_state` reinicia Lanczos
hasta que el residuo $\|H\psi-E\psi\|$ queda por debajo de $10^{-10}$. Lo comprobamos frente a las energías del
anillo de Heisenberg tabuladas con nueve cifras en `scripts/chains.py` (los valores de referencia del artículo
complementario), de forma densa hasta $L=12$ y con Lanczos para $L=16$:

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

Para $L=16$ (dimensión $12870$), Lanczos tarda bastante menos de un segundo en un hilo de CPU. La energía por
sitio, $-0.5,-0.467129,-0.456387,-0.451545,-0.448949$ para $L=4,\dots,12$, se acerca por debajo al valor
$\frac14-\log2=-0.443147$ de la cadena infinita.

**Evolución temporal.** Para $e^{-i\tau H}v$ se proyecta sobre $\mathcal K_m$: $w_\tau=\|v\|\,Q_me^{-i\tau T_m}e_1$.

**Proposición 9.15.** Sea $\dim\mathcal K_m=m$. Para todo $\tau$, $\|w_\tau\|=\|v\|$ y
$\langle w_\tau,Hw_\tau\rangle=\langle v,Hv\rangle$. Si $\mathcal K_m$ es $H$-invariante, $w_\tau=e^{-i\tau H}v$.

*Demostración.* $Q_m$ es una isometría y $e^{-i\tau T_m}$ es unitaria ($T_m$ es real simétrica), lo que da la norma.
Como $Q_m^*HQ_m=T_m$ conmuta con $e^{-i\tau T_m}$,
$\langle w_\tau,Hw_\tau\rangle=\|v\|^2\langle e^{-i\tau T_m}e_1,T_me^{-i\tau T_m}e_1\rangle=\|v\|^2\langle e_1,T_me_1\rangle=\langle v,Hv\rangle$. Si
$\mathcal K_m$ es invariante, $HQ_m=Q_mT_m$ (proposición 9.13, $r_m=0$), así que $H^jQ_m=Q_mT_m^j$ y
$e^{-i\tau H}Q_m=Q_me^{-i\tau T_m}$; se aplica a $e_1$, con $Q_me_1=v/\|v\|$. $\square$

Así, cada paso de Krylov conserva la norma y la energía, sea preciso o no; la precisión exige un $m$ grande en
comparación con $\tau\|H\|$, y `krylov.evolve` subdivide el intervalo de tiempo y estima el error a posteriori
(M. Hochbruck, C. Lubich, *SIAM J. Numer. Anal.* 34 (1997) 1911; enunciado sin demostración). Realizamos un
*quench* desde el estado de Néel $|{\uparrow\downarrow\uparrow\downarrow}\cdots\rangle$ en el anillo de Heisenberg y
seguimos la magnetización alternada $m_s=\frac2L\sum_i(-1)^i\langle s^z_i\rangle$, comparando con el propagador denso
para $L=12$:

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

$m_s$ decae desde $1$, pasando por $0.770$ y $0.279$, hasta $-0.183$ en $t=2$, para ambos tamaños con una
diferencia de $1.7\cdot10^{-5}$; en $t=4$ difieren ($0.026$ frente a $0.035$): el tamaño finito del anillo más
pequeño se hace visible.

## 9.6 Desorden: estadística de niveles y entrelazamiento

La **cadena de Heisenberg con campo aleatorio** tiene campos independientes $h_i$ uniformes en $[-W,W]$. Conserva
$S^z$ (proposición 9.1) y, genéricamente, nada más: los campos rompen las traslaciones, las reflexiones, la
inversión de espines y $SU(2)$. Por eso la estadística de niveles se calcula por sector $N$ (teorema 5.15 y
sección 5.4: mezclar sectores simula niveles independientes). Para desorden débil, los niveles se repelen como
en el GOE, $\langle r\rangle\approx0.5307$ (sección 5.3); para desorden fuerte, los estados propios se localizan
(**localización de muchos cuerpos**) y los niveles se vuelven independientes, $\langle r\rangle=2\log2-1\approx0.3863$
(proposición 5.14). Calculamos $\langle r\rangle$ en la mitad central del sector $N=L/2$, con un error estándar sobre
las realizaciones independientes (cada una da un cociente medio):

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

Para $L=10$, el cociente medio va de $0.530$ ($W=0.5$), pasando por $0.484$ ($W=2$) y $0.410$ ($W=4$), hasta
$0.384$ ($W=16$), con errores estándar de $0.003$–$0.004$. Para $L=8$, el valor en $W=0.5$ ($0.507$) está por
debajo del de $W=1$ ($0.525$): en $W=0$ la cadena es integrable y tiene todas las simetrías anteriores, y un valor
de GOE requiere un desorden grande en comparación con el espaciamiento entre niveles. La transición es amplia y se
desplaza con $L$; el escalado de tamaño finito la sitúa cerca de $W\approx3.7$ (D. Luitz, N. Laflorencie, F. Alet,
*Phys. Rev. B* 91 (2015) 081103), y se debate si una fase localizada sobrevive cuando $L\to\infty$ (J. Šuntajs,
J. Bonča, T. Prosen, L. Vidmar, *Phys. Rev. E* 102 (2020) 062144).

**Entrelazamiento de fermiones libres.** Para los estados propios de la cadena XX abierta, la entropía de
entrelazamiento de los primeros $n$ sitios está determinada por una matriz $n\times n$.

**Teorema 9.16 (fórmula de Peschel).** Sea $\Phi=b_{k_1}^*\cdots b_{k_N}^*|0\rangle$ como en el teorema 9.10,
$C_{ij}=\langle\Phi,c_i^*c_j\Phi\rangle$ y $C_A=(C_{ij})_{i,j<n}$. Entonces $0\le C_A\le1$, y la entropía de
entrelazamiento de los sitios $A=\{0,\dots,n-1\}$ es, con $\nu$ recorriendo los valores propios de $C_A$ contados
con multiplicidad,
$$
S_A=\sum_\nu \eta(\nu),\qquad \eta(\nu)=-\nu\log\nu-(1-\nu)\log(1-\nu)\quad(0\log0=0).
$$

*Demostración.* Para un vector real $x$, pongamos $b(x)=\sum_jx_jc_j$; entonces $b(x)^*b(y)^*=-b(y)^*b(x)^*$ y
$\{b(x),b(y)^*\}=\langle x,y\rangle$. Sean $U\subseteq\mathbb R^L$ el subespacio generado por los $u_{k_a}$, y $P_A$,
$P_B=1-P_A$ las proyecciones de coordenadas sobre $\{0,\dots,n-1\}$ y su complemento. Diagonalizar la forma
$\langle P_Ax,P_Ay\rangle$ sobre $U$ da una base ortonormal $\phi_1,\dots,\phi_N$ de $U$ con
$\langle P_A\phi_a,P_A\phi_b\rangle=\nu_a\delta_{ab}$ y, por tanto, $\langle P_B\phi_a,P_B\phi_b\rangle=(1-\nu_a)\delta_{ab}$. Como
$(x_1,\dots,x_N)\mapsto b(x_1)^*\cdots b(x_N)^*|0\rangle$ es multilineal y alternada,
$\Phi=\pm\prod_a\big(\sqrt{\nu_a}f_a^*+\sqrt{1-\nu_a}\,g_a^*\big)|0\rangle$ con $f_a=b(P_A\phi_a/\sqrt{\nu_a})$,
$g_a=b(P_B\phi_a/\sqrt{1-\nu_a})$ (omitiendo los términos con coeficiente nulo). Desarrollando,
$\Phi=\sum_I\pm\prod_{a\in I}\sqrt{\nu_a}\prod_{a\notin I}\sqrt{1-\nu_a}\,F_I^*G_{I^c}^*|0\rangle$ sobre los subconjuntos
$I$, con $F_I^*$, $G_{I^c}^*$ los productos ordenados. Para $j\ge n$, la cadena de signos de $c_j^*$ contiene todo
$A$: $c_j^*=(-1)^{N_A}\otimes\tilde c_j^*$ con operadores de Jordan–Wigner $\tilde c_j$ de los sitios $\ge n$,
mientras que para $i<n$, $c_i^*=\hat c_i^*\otimes1$. Como $(-1)^{N_A}$ vale $1$ sobre $|0_A\rangle$,
$F_I^*G_{I^c}^*|0\rangle=\xi_I\otimes\chi_{I^c}$ con $\xi_I=\hat F_I^*|0_A\rangle$, $\chi_{I^c}=\tilde G_{I^c}^*|0_B\rangle$, que son
familias ortonormales por el cálculo con las CAR del teorema 9.10. Esto es una descomposición de Schmidt
(lema 7.7) con pesos $p_I=\prod_{a\in I}\nu_a\prod_{a\notin I}(1-\nu_a)$, la ley de variables de Bernoulli
independientes, así que $S_A=-\sum_Ip_I\log p_I=\sum_a\eta(\nu_a)$. Por último, $C_{ij}=\sum_a\phi_a(i)\phi_a(j)$
(desarrollando $c_i^*c_j$ en los $b(\phi_a)$), así que $C_A=\sum_a(P_A\phi_a)(P_A\phi_a)^T$ tiene los valores propios
$\nu_a$ sobre los vectores ortogonales $P_A\phi_a$ y $0$ en el resto, y $\eta(0)=0$. $\square$

`entanglement_entropy` calcula el estado reducido sector a sector (conmuta con la magnetización del bloque,
proposición 7.8) y su entropía. El estado fundamental del sector $N$ es $\Phi$ con los $N$ modos más bajos, simple
porque los $\varepsilon_k$ son distintos:

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

Las entropías $0.6931,0.5199,0.7456,0.6283,0.7719,0.6560$ para $n=1,\dots,6$ oscilan con la paridad de $n$ y
crecen como $\frac16\log L$ con $n/L$ fijo, la ley de la teoría de campos conforme para una cadena abierta
(I. Peschel, *J. Phys. A* 36 (2003) L205; P. Calabrese, J. Cardy, *J. Stat. Mech.* (2004) P06002; enunciado sin
demostración).

## 9.7 Con qué comparar

QuSpin (Python) construye hamiltonianos de cadenas de espines en sectores de magnetización, momento, paridad e
inversión de espín, y usa los resolvedores dispersos de valores propios y las exponenciales de matrices de SciPy;
los códigos dedicados de diagonalización exacta (XDiag, SPINPACK, HΦ) están optimizados para grandes problemas
dispersos, y con simetrías de red y máquinas grandes tales códigos han alcanzado unos 40–50 espines; ITensor,
TeNPy y TensorKit.jl (con árboles de fusión de $SU(2)$) tratan cadenas largas con estados de producto de
matrices. Lo que aporta la biblioteca es modesto y específico: la cadena como `TypeIAlgebra`
$\bigoplus_NM_{\binom LN}$, de modo que los estados, canales, entropías y dinámicas de los capítulos 3–8 se aplican
sector a sector; la torre de inclusiones con su índice; bloques de $SU(2)$ a partir de vectores de peso máximo
con correlaciones de Wigner–Eckart; lotes de realizaciones del desorden en una sola llamada; y estados reducidos
que nunca salen de los sectores. No tiene sectores de reflexión ni de inversión de espín, tiene sectores de
momento solo para cadenas XXZ uniformes, métodos densos hasta $L\approx16$, métodos de Krylov hasta $L\approx24$
(unos $26$ en una GPU grande), y no tiene redes tensoriales.

## Ejercicios

1. (Sectores de momento a mano) Para $L=4$ y $N=2$, halle las órbitas de las seis configuraciones bajo $T$ y las
   dimensiones $d_{2,m}$, $m=0,\dots,3$.
2. (Reglas de suma) Sea $\psi$ el estado fundamental del anillo de Heisenberg con $L$ par. Demuestre que
   $E_0=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$ y $\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle=0$.
3. (Demasiada simetría) Demuestre que en el sector $N=L/2$ el vector con todos los coeficientes iguales a $1$ es
   un múltiplo de $(S^-)^{L/2}|{\uparrow\cdots\uparrow}\rangle$ y un vector propio del anillo de Heisenberg con
   energía $JL/4$, y que todo espacio de Krylov $\mathcal K_m(H,v)$ está contenido en cada subespacio
   $H$-invariante que contenga a $v$. ¿Qué devuelve `krylov.ground_state` cuando se inicia desde este vector?

<details>
<summary>Soluciones</summary>

1. $0011\to1001\to1100\to0110\to0011$ es una órbita con $R=4$, y $0101\to1010\to0101$ una con $R=2$. Por la
   proposición 9.5, la primera contribuye a todo $m$, la segunda a los $m$ con $2m\equiv0\pmod4$, es decir,
   $m=0,2$: $d_{2,m}=2,1,2,1$.

```python
assert SpinChain(4, 'periodic').momentum_dims(2) == [2, 1, 2, 1]
```

2. El estado fundamental es un singlete simple (teorema 9.9) y $T$ conmuta con $H$, así que $T\psi=e^{ik}\psi$.
   Como $\mathbf s_i\cdot\mathbf s_{i+1}=T^i(\mathbf s_0\cdot\mathbf s_1)T^{*i}$, $\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle$ es el
   mismo para los $L$ enlaces, y $E_0=J\sum_i\langle\mathbf s_i\cdot\mathbf s_{i+1}\rangle=JL\langle\mathbf s_0\cdot\mathbf s_1\rangle$.
   Como $\mathbf S^2\psi=0$,
   $0=\langle\mathbf S^2\rangle=\sum_{i,j}\langle\mathbf s_i\cdot\mathbf s_j\rangle=L\sum_r\langle\mathbf s_0\cdot\mathbf s_r\rangle$, de nuevo
   por invariancia por traslaciones.

```python
assert abs(12 * corr_we[1].item() - gs.energy.item()) < 1e-10 and abs(corr_we.sum().item()) < 1e-10
```

3. $S^-$ baja uno de los espines hacia arriba; aplicándolo $L/2$ veces al estado totalmente polarizado se llega a
   cada configuración con $L/2$ espines hacia arriba en $(L/2)!$ órdenes, todos con coeficiente $1$.
   $|{\uparrow\cdots\uparrow}\rangle$ tiene los $L$ enlaces paralelos, energía $JL/4$, y $S^-$ conmuta con $H$, así que
   el vector es un vector propio con la misma energía. Si $W$ es invariante y $v\in W$, entonces $H^jv\in W$ para
   todo $j$. Aquí $\mathcal K_m=\mathbb Cv$ para todo $m$, y `ground_state` devuelve $E=JL/4$, el extremo superior
   del espectro (cada $\mathbf s_i\cdot\mathbf s_j\le\frac14$); un vector inicial aleatorio lo evita con probabilidad
   uno.

```python
c10 = SpinChain(10, 'periodic', complex_valued=False)
E_sym, _, _ = krylov.ground_state(c10.xxz_sparse(1.0, 1.0, sector=5), v0=torch.ones(1, 252, dtype=f64))
assert abs(E_sym.item() - 2.5) < 1e-10 and E_sym.item() > RING[10] + 7     # far from E_0
```

</details>

**Siguiente:** el [capítulo 10](10_quantum_computing.md) trata la computación cuántica con sectores: subespacios
libres de decoherencia y subsistemas sin ruido, corrección de errores en álgebras de operadores, y mezcladores
que respetan un peso de Hamming conservado: el álgebra $\bigoplus_wM_{\binom nw}$ de este capítulo en un nuevo
papel.
