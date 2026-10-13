[English](../08_dynamics.md) · **Español** · [Français](../fr/08_dynamics.md) · [Deutsch](../de/08_dynamics.md) · [中文](../zh/08_dynamics.md) · [日本語](../ja/08_dynamics.md) · [Русский](../ru/08_dynamics.md)

# 8. Dinámica

Un sistema cuántico cerrado evoluciona según un grupo uniparamétrico de unitarias; uno abierto —un
sistema acoplado a un entorno que no se observa—, según un semigrupo de canales. En un álgebra con
sectores ambos tienen un rasgo invisible en una sola álgebra de matrices: un hamiltoniano de
$\mathcal M$ nunca puede cambiar el sector, de modo que las cargas se conservan, mientras que un
sistema abierto puede saltar entre sectores, por ejemplo al perder una partícula. Este capítulo
demuestra los hechos básicos sobre ambos tipos de evolución, escribe la ecuación de Lindblad con
saltos entre sectores y muestra que "no se perdió ninguna partícula" es una medición de la carga
cuyo estado postseleccionado evoluciona según un hamiltoniano no hermítico.

**Aprenderá** por qué los pesos de los sectores son constantes del movimiento, por qué los
generadores en forma de Gorini–Kossakowski–Sudarshan–Lindblad (GKSL) producen canales, cómo entran
los saltos entre sectores en la ecuación de Lindblad y cuándo los pesos de los sectores obedecen una
ecuación de tasas clásica, qué significan la evolución sin saltos y la postselección, y cómo hallar
los estados estacionarios y las tasas de relajación.

En todo el capítulo, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ actúa sobre
$\mathcal H=\bigoplus_c\mathbb C^{k_c}$, $D=\sum_ck_c$, $P_c$ son las proyecciones centrales
minimales, los estados se dan mediante sus densidades `blunt` $\rho$ con pesos de los sectores
$p_c=\operatorname{Tr}\rho_c$ (capítulo 3), $[A,B]=AB-BA$, $\{A,B\}=AB+BA$ y $\hbar=1$. Un
**superoperador** es una aplicación lineal sobre $\mathcal M$ o sobre $M_D(\mathbb C)$; todas las
normas de superoperadores son normas de operador y, por tanto, submultiplicativas. Para una
aplicación lineal $B$ sobre un espacio de dimensión finita, $e^{tB}=\sum_{n\ge0}t^nB^n/n!$.

## 8.1 Sistemas cerrados

Un **hamiltoniano** es un $H\in\mathcal M$ autoadjunto. Un vector evoluciona según la **ecuación de
Schrödinger** $i\psi'(t)=H\psi(t)$, y una densidad según la **ecuación de von Neumann**
$\rho'(t)=-i[H,\rho(t)]$. Ambas son ecuaciones lineales con coeficientes constantes sobre un espacio
de dimensión finita, y el lema siguiente es toda la teoría de ecuaciones diferenciales que
necesitamos.

**Lema 8.1.** Sea $B$ una aplicación lineal sobre un espacio normado de dimensión finita $X$.
Entonces $\frac d{dt}e^{tB}=Be^{tB}=e^{tB}B$, $e^{(s+t)B}=e^{sB}e^{tB}$,
$\|e^{tB}\|\le e^{|t|\,\|B\|}$ y, para todo $x_0\in X$, la aplicación $x(t)=e^{tB}x_0$ es la única
$x$ derivable en un intervalo $I\ni0$ con $x'=Bx$ y $x(0)=x_0$.

*Demostración.* La serie converge absolutamente, con la cota de la norma, y como serie de potencias
en $t$ puede derivarse término a término, lo que da la primera fórmula. Entonces
$\frac d{dt}(e^{tB}e^{-tB})=Be^{tB}e^{-tB}-e^{tB}Be^{-tB}=0$, luego $e^{tB}e^{-tB}=1$. Si $x'=Bx$,
entonces $\frac d{dt}(e^{-tB}x(t))=e^{-tB}(x'-Bx)=0$, luego $e^{-tB}x(t)=x_0$ y $x(t)=e^{tB}x_0$.
Por último, $t\mapsto e^{(s+t)B}x_0$ y $t\mapsto e^{tB}e^{sB}x_0$ resuelven $y'=By$ con el mismo
valor inicial. $\square$

**Proposición 8.2.** Sea $H=H^*\in\mathcal M$ y $U(t)=e^{-itH}$, definida mediante el cálculo
funcional de la sección 1.3 con $f_t(\lambda)=e^{-it\lambda}$.

1. $U(t)=\bigoplus_ce^{-itH_c}$ es una unitaria de $\mathcal M$, $U(s+t)=U(s)U(t)$ y
   $U'(t)=-iHU(t)$; $U(t)$ coincide con la serie exponencial de $-itH$.
2. $\psi(t)=U(t)\psi_0$ es la única solución de $i\psi'=H\psi$ con $\psi(0)=\psi_0$, y
   $\rho(t)=U(t)\rho_0U(t)^*$ es la única solución en $\mathcal M$ de $\rho'=-i[H,\rho]$ con
   $\rho(0)=\rho_0$. Si $\rho_0$ es un estado, también lo es cada $\rho(t)$.
3. (*Imagen de Heisenberg*) $\operatorname{Tr}(\rho(t)A)=\operatorname{Tr}(\rho_0A(t))$ con
   $A(t)=U(t)^*AU(t)$, y $A'=i[H,A]$.

*Demostración.* (1) Diagonalicemos cada bloque, $H_c=V_c\Lambda_cV_c^*$; entonces $V=\bigoplus_cV_c$
diagonaliza $H$, y $U(t)=Ve^{-it\Lambda}V^*$ es diagonal por bloques (por el teorema 1.5(a), el
resultado no depende de la unitaria que diagonaliza). Como $|f_t|^2=1$ y $f_sf_t=f_{s+t}$, el
teorema 1.5(c) da la unitariedad y la ley de grupo. Derivando $e^{-it\Lambda}$ entrada a entrada se
obtiene $U'=V(-i\Lambda)e^{-it\Lambda}V^*=-iHU$, luego $U(t)$ es la serie $e^{-itH}$ por el lema 8.1
($X=M_D$, $B$ la multiplicación por la izquierda por $-iH$). (2) $(U\psi_0)'=-iHU\psi_0$, y
$(U\rho_0U^*)'=-iHU\rho_0U^*+U\rho_0U^*(iH)$, porque $(U')^*=iU^*H$. La unicidad es el lema 8.1 con
$B=-iH$ sobre $\mathcal H$, o bien con $B=-i[H,\cdot\,]$ sobre $\mathcal M$. Positividad:
$\langle x,U\rho_0U^*x\rangle=\langle U^*x,\rho_0U^*x\rangle\ge0$; la traza, por ciclicidad. (3)
Ciclicidad de la traza, y $(U^*AU)'=iU^*HAU-iU^*AHU=i[H,A(t)]$. $\square$

**Teorema 8.3 (conservación de la carga).** Sea $t\mapsto H(t)$, definida en un intervalo $I$, con
valores en los elementos autoadjuntos de $\mathcal M$.

1. Si $\psi:I\to\mathcal H$ es derivable con $i\psi'=H(t)\psi$, entonces cada $\|P_c\psi(t)\|$ es
   constante.
2. Si $\rho:I\to\mathcal M$ es derivable con $\rho'=-i[H(t),\rho]$, entonces
   $\operatorname{Tr}(Z\rho(t))$ es constante para todo $Z$ central; en particular, lo es cada peso
   de sector $p_c(t)=\operatorname{Tr}(P_c\rho(t))$.
3. Si $H$ es constante y $A\in\mathcal M$ conmuta con $H$, entonces $\operatorname{Tr}(A\rho(t))$ es
   constante.

*Demostración.* (1) El producto escalar es antilineal en el primer argumento, así que
$\frac d{dt}\langle\psi,P_c\psi\rangle=\langle-iH\psi,P_c\psi\rangle+\langle\psi,-iP_cH\psi\rangle =i\langle\psi,HP_c\psi\rangle-i\langle\psi,P_cH\psi\rangle=0$,
pues $P_c$ es central. (2)
$\frac d{dt}\operatorname{Tr}(Z\rho)=-i\operatorname{Tr}(Z[H,\rho])=-i\operatorname{Tr}([Z,H]\rho)=0$
por ciclicidad. (3) El mismo cálculo con $A$. $\square$

> **Significado físico.** Un hamiltoniano de $\mathcal M$ conmuta con las cargas (número de
> partículas, magnetización, un número cuántico de simetría), de modo que la probabilidad de cada
> carga es una constante del movimiento, incluso para hamiltonianos forzados y dependientes del
> tiempo. Desde el punto de vista computacional, la evolución se separa en $C$ problemas
> independientes: una diagonalización cuesta del orden de $\sum_ck_c^3$ operaciones en lugar de
> $D^3$.

La biblioteca calcula $U(t)$ con una descomposición espectral por cada tamaño de bloque
(`dynamics.propagator`), y `dynamics.von_neumann` hace evolucionar una densidad hasta muchos
instantes a la vez en la base propia, donde $\rho(t)_{ij}=e^{-it(\lambda_i-\lambda_j)}\rho_{ij}$.
Comprobamos frente a `scipy.linalg.expm`:

```python
import math
import numpy as np
import scipy.integrate
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, InterSectorChannel, channels, dynamics, states

torch.manual_seed(0)
f64, c128 = torch.float64, torch.complex128
alg = TypeIAlgebra([2, 3, 4], [2, 3, 4], precision='double')     # M_2 (+) M_3 (+) M_4, D = 9
H = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
t = 0.7
U = dynamics.propagator(H, t)                                      # exp(-itH), block by block
for c, k in enumerate(alg.k_factors):                              # independent reference
    ref = sla.expm(-1j * t * H.matrix[0, c, :k, :k].numpy())
    assert np.abs(U.matrix[0, c, :k, :k].numpy() - ref).max() < 1e-12
assert (U.adjoint() @ U - alg.identity()).matrix.abs().max() < 1e-12

rho0 = states.random_density_matrix(alg)
times = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.von_neumann(rho0, H, times)                        # exact, all times at once
for s, rho_s in zip(times, traj):
    U_s = dynamics.propagator(H, s)
    assert (rho_s.matrix - (U_s @ rho0 @ U_s.adjoint()).matrix).abs().max() < 1e-12
    assert (rho_s.sector_probabilities() - rho0.sector_probabilities()).abs().max() < 1e-12   # Theorem 8.3
A = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
A_t = U.adjoint() @ A @ U                                          # Heisenberg picture
assert abs(rho0.expectation(A_t) - dynamics.von_neumann(rho0, H, [t])[0].expectation(A)) < 1e-12
```

Sobre $\mathcal M$, un vector $\psi=\bigoplus_c\psi_c$ define el estado de la mezcla
$\bigoplus_c|\psi_c\rangle\langle\psi_c|$ (capítulo 3, ejercicio 3), así que `dynamics.schrodinger`
(vectores con la disposición por bloques `(batch, C, k_max)`) y `von_neumann` deben coincidir sobre
él. Para un $H(t)$ dependiente del tiempo, la biblioteca usa el método de Runge–Kutta clásico de
paso fijo (`dynamics.rk4`, `schrodinger_rk4`): un paso multiplica $\psi$ por un polinomio en valores
de $H(t)\in\mathcal M$, por lo que nunca traslada amplitud entre sectores, pero no es exactamente
unitario, y los pesos de los sectores derivan al nivel del error de discretización. (Para un único
sector grande, `krylov.evolve` solo necesita productos matriz–vector; véase el capítulo 9.)

```python
rows = torch.tensor([[1, 1, 0, 0], [1, 1, 1, 0], [1, 1, 1, 1]], dtype=c128)   # entries of C^2 (+) C^3 (+) C^4
psi0 = torch.randn(1, alg.C, alg.k_max, dtype=c128) * rows
psi0 = psi0 / psi0.norm()
psi = dynamics.schrodinger(psi0, H, times)                         # shape (T, batch, C, k_max)
p_vec = (psi.abs() ** 2).sum(-1)                                   # ||P_c psi(t)||^2
assert (p_vec - p_vec[0]).abs().max() < 1e-12
mix = DensityMatrix(alg, matrix=psi0[..., :, None] * psi0[..., None, :].conj())   # the state of psi0 on M
outer = psi[..., :, None] * psi[..., None, :].conj()
assert all((o - r.matrix).abs().max() < 1e-12 for o, r in zip(outer, dynamics.von_neumann(mix, H, times)))

V = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
H_of_t = lambda s: H.matrix + math.cos(3 * s) * V.matrix          # driven, in M at every time
psi_rk = dynamics.schrodinger_rk4(psi0, H_of_t, times, substeps=50, algebra=alg)
drift = ((psi_rk.abs() ** 2).sum(-1) - p_vec[0]).abs().max().item()
print(f"RK4, time-dependent H: largest change of a sector weight {drift:.1e}")
assert drift < 1e-8 and (dynamics.schrodinger_rk4(psi0, H, times, substeps=50) - psi).abs().max() < 1e-7

```

Los pesos de los sectores cambian en $1.1\cdot10^{-9}$ con RK4 y en nada más allá del redondeo con
los métodos exactos.

## 8.2 Semigrupos dinámicos cuánticos

Las aplicaciones CP y TP y los canales se definen en la sección 6.1 ([capítulo
6](06_quantum_channels.md)). Además del teorema 6.2 y de la proposición 6.4 (una aplicación
$A\mapsto\sum_lK_lAK_l^*$, también entre sectores, es CP), necesitamos tres propiedades de clausura.

**Lema 8.4.** (a) Las sumas, los múltiplos no negativos y las composiciones de aplicaciones CP son
CP. Si $\Phi$ es CP sobre $M_D$ y lleva una $*$-subálgebra $\mathcal N\subseteq M_D$ (como
$\mathcal M$, o $P_c\mathcal MP_c\cong M_{k_c}$ con unidad $P_c$) en sí misma, entonces
$\Phi|_{\mathcal N}$ es CP. (b) Si unas aplicaciones CP $\Phi_m$ convergen a $\Phi$ (es decir,
$\Phi_m(A)\to\Phi(A)$ para todo $A$), entonces $\Phi$ es CP. (c) Si $s\mapsto\Psi(s)$ es una
aplicación continua de $[a,b]$ en las aplicaciones CP, entonces $\int_a^b\Psi(s)\,ds$ es CP.

*Demostración.* (a) $\Phi\mapsto\Phi_n$ es lineal y $(\Phi\Psi)_n=\Phi_n\Psi_n$; las sumas de
matrices positivas son positivas. Para la restricción: $M_n(\mathcal N)\subseteq M_n(M_D)$, y la
positividad es la misma condición en ambos (tras una permutación de la base, un elemento de
$M_n(\mathcal N)$ es una suma directa de los bloques $[(A_{ij})_c]$ y de ceros). (b) Las matrices
positivas forman un conjunto cerrado ($\langle x,Xx\rangle\ge0$ para cada $x$ es una condición
cerrada), y $(\Phi_m)_n(X)\to\Phi_n(X)$. (c) Las sumas de Riemann son CP por (a) y convergen a la
integral, ya que $\Psi$ es uniformemente continua; se aplica (b). $\square$

**Definición.** Un **semigrupo dinámico cuántico** sobre $\mathcal M$ es una familia $(T_t)_{t\ge0}$
de canales $\mathcal M\to\mathcal M$ con $T_0=\mathrm{id}$, $T_{s+t}=T_sT_t$ y $t\mapsto T_t$
continua.

**Lema 8.5.** Todo semigrupo continuo de superoperadores ($T_0=\mathrm{id}$, $T_{s+t}=T_sT_t$) es de
la forma $T_t=e^{t\mathcal L}$ para exactamente un superoperador $\mathcal L$, su **generador**, y
$\mathcal L=\lim_{t\to0}(T_t-\mathrm{id})/t$.

*Demostración.* $V_h=\frac1h\int_0^hT_s\,ds\to\mathrm{id}$ cuando $h\to0$, luego $V_h$ es invertible
para algún $h$ pequeño fijo. Ahora bien, $T_tV_h=\frac1h\int_t^{t+h}T_s\,ds$ es derivable en $t$ con
derivada $\frac1h(T_{t+h}-T_t)=T_t\frac{T_h-\mathrm{id}}h$. Por tanto $T_t$ es derivable y
$T_t'=T_t\mathcal L$ con $\mathcal L=\frac{T_h-\mathrm{id}}hV_h^{-1}$, y $T_t=e^{t\mathcal L}$ por
el lema 8.1 (con $B$ la composición por la derecha con $\mathcal L$). Unicidad: $\mathcal L=T'(0)$.
$\square$

**Teorema 8.6 (generadores GKSL).** Sean $H=H^*$ y $L_1,\dots,L_J$ en $M_D$, y
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
$$
Entonces $(e^{t\mathcal L})_{t\ge0}$ es un semigrupo dinámico cuántico sobre $M_D$. Más
precisamente, pongamos $H_{\mathrm{eff}}=H-\frac i2\sum_jL_j^*L_j$,
$\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)$ y
$\Phi(\rho)=\sum_jL_j\rho L_j^*$, de modo que $\mathcal L=\mathcal L_0+\Phi$, y definamos
$$
S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},\qquad
S_n(t)=\int_0^tS_0(t-s)\,\Phi\,S_{n-1}(s)\,ds\quad(n\ge1).
$$
Entonces cada $S_n(t)$ es CP, y $e^{t\mathcal L}=\sum_{n\ge0}S_n(t)$ (el **desarrollo de Dyson**),
con $\|S_n(t)\|\le e^{t\|\mathcal L_0\|}(t\|\Phi\|)^n/n!$.

*Demostración.* *Paso 1.* $X(t)=e^{-itH_{\mathrm{eff}}}\rho e^{itH_{\mathrm{eff}}^*}$ cumple
$X'=\mathcal L_0(X)$ y $X(0)=\rho$, luego $S_0(t)=e^{t\mathcal L_0}$ por el lema 8.1; es CP con el
único operador de Kraus $e^{-itH_{\mathrm{eff}}}$. *Paso 2.* Escribamos $a=\|\mathcal L_0\|$,
$b=\|\Phi\|$. Por inducción, las $S_n$ son continuas y
$\|S_n(t)\|\le\int_0^te^{(t-s)a}\,b\,e^{sa}\frac{(bs)^{n-1}}{(n-1)!}ds=e^{ta}\frac{(bt)^n}{n!}$, de
modo que la serie converge uniformemente en todo $[0,T]$. *Paso 3.* Como
$S_n(t)=S_0(t)\int_0^te^{-s\mathcal L_0}\Phi S_{n-1}(s)\,ds$, la regla del producto y el teorema
fundamental del cálculo dan $S_n'=\mathcal L_0S_n+\Phi S_{n-1}$ ($n\ge1$) y $S_0'=\mathcal L_0S_0$.
Por las cotas, la serie de las derivadas converge uniformemente en $[0,T]$, luego
$Y(t)=\sum_nS_n(t)$ es derivable con $Y'=\mathcal L_0Y+\Phi Y=\mathcal LY$, $Y(0)=\mathrm{id}$, e
$Y(t)=e^{t\mathcal L}$ por el lema 8.1. *Paso 4.* $\Phi$ es CP (teorema 6.2), así que, por
inducción, el integrando $s\mapsto S_0(t-s)\Phi S_{n-1}(s)$ es una familia continua de aplicaciones
CP y $S_n(t)$ es CP (lema 8.4); también lo son las sumas parciales y su límite $e^{t\mathcal L}$.
*Paso 5.* Por ciclicidad,
$\operatorname{Tr}\mathcal L(\rho)=\sum_j\big(\operatorname{Tr}(L_j\rho L_j^*)-\operatorname{Tr}(L_j^*L_j\rho)\big)=0$,
luego
$\frac d{dt}\operatorname{Tr}e^{t\mathcal L}\rho=\operatorname{Tr}\mathcal L(e^{t\mathcal L}\rho)=0$:
cada $e^{t\mathcal L}$ es TP. La ley de semigrupo y la continuidad son el lema 8.1. $\square$

**Teorema 8.7 (Gorini–Kossakowski–Sudarshan; Lindblad; enunciado sin demostración).**
Recíprocamente, el generador de todo semigrupo dinámico cuántico sobre $M_D$ tiene la forma del
teorema 8.6. [V. Gorini, A. Kossakowski, E. C. G. Sudarshan, *J. Math. Phys.* 17, 821 (1976); G.
Lindblad, *Commun. Math. Phys.* 48, 119 (1976).] Para un semigrupo continuo en norma de aplicaciones
CP unitales sobre un álgebra de von Neumann general $\mathcal M$, Christensen y Evans [*J. London
Math. Soc.* 20, 358 (1979)] demostraron que el generador es $A\mapsto\Psi(A)+K^*A+AK$ con $\Psi$ CP
y $K\in\mathcal M$.

Sobre $\mathcal M=\bigoplus_cM_{k_c}$, apliquemos esto al semigrupo $(T_t^*)$ de la imagen de
Heisenberg, cuyas aplicaciones son CP y unitales (proposición 6.8, corolario 6.9(b)) y cuyo
generador es $\mathcal L^*$. Volviendo a los duales `blunt`,
$\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*$ con $\Phi=\Psi^*$ CP sobre $\mathcal M$ (operadores de
Kraus entre sectores, proposición 6.4) y $G=K$. La preservación de la traza,
$\operatorname{Tr}\big((\Phi^*(1)+G+G^*)\rho\big)=0$ para todo $\rho\in\mathcal M$, obliga a que
$G+G^*=-\Phi^*(1)$, es decir, $G=-iH-\frac12\Phi^*(1)$ con $H=H^*\in\mathcal M$: exactamente los
generadores de la sección siguiente.

> **Significado físico.** Un sistema que interactúa débilmente con un entorno grande cuyas
> correlaciones decaen rápidamente (la aproximación de Born–Markov) pierde la memoria de su pasado:
> su estado en el instante $s+t$ se obtiene del estado en el instante $s$ mediante un canal que solo
> depende de $t$. El teorema 8.7 dice que una dinámica así es siempre una parte hamiltoniana más
> unos "saltos" $L_j$ (emisión de fotones, pérdida de partículas, desfase) con tasas fijadas por la
> normalización de los $L_j$.

La biblioteca tiene cuatro rutinas de Lindblad: `lindblad_rhs` devuelve el segundo miembro como
función para `rk4`; `lindblad_evolve` integra con RK4; `lindblad_superoperator` devuelve, para cada
sector, la matriz $k_c^2\times k_c^2$ de $\mathcal L$ sobre $\operatorname{vec}\rho_c$ ordenado por
filas; y `lindblad_channel` calcula $e^{t\mathcal L}$ exactamente, sector a sector, como un
`Channel`. Comprobamos el teorema 8.6 para saltos arbitrarios sobre $M_3$: la matriz de Choi de
$e^{t\mathcal L}$ es positiva (teorema 6.2), la aplicación es TP, se cumple la ley de semigrupo, RK4
coincide con la exponencial, y la parte de saltos $e^{t\mathcal L}-S_0(t)=\sum_{n\ge1}S_n(t)$ es a
su vez CP:

```python
M3 = TypeIAlgebra([3], [3], precision='double')
H3 = M3.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
jumps = [M3.operator(torch.randn(1, 1, 3, 3, dtype=c128) / 2) for _ in range(2)]     # arbitrary L_j
T = lambda s: dynamics.lindblad_channel(M3, H3, jumps, t=s)        # exp(sL), exact
S = lambda s: T(s).superoperator()[0][0]                           # its 9 x 9 matrix on row-major vec
choi = lambda X: X.reshape(3, 3, 3, 3).permute(2, 0, 3, 1).reshape(9, 9)   # J_{(i,a),(j,b)} = S_{(a,b),(i,j)}
assert (choi(S(0.8)) - T(0.8).choi()[0][0]).abs().max() < 1e-12
assert torch.linalg.eigvalsh(choi(S(0.8))).min() > -1e-12         # completely positive
assert T(0.8).trace_preservation_error().max() < 1e-12             # trace preserving
assert (S(0.3) @ S(0.5) - S(0.8)).abs().max() < 1e-12              # semigroup

rho3 = states.random_density_matrix(M3)
rk = dynamics.lindblad_evolve(rho3, H3, jumps, [0.0, 0.8], substeps=40)[-1]
assert (rk.matrix - T(0.8)(rho3).matrix).abs().max() < 1e-7        # RK4 against exp(tL)
assert abs(dynamics.lindblad_rhs(H3, jumps)(0.0, rho3.matrix)[0, 0].trace()) < 1e-12   # Tr L(rho) = 0

H_eff = H3.matrix - 0.5j * sum(L.adjoint().matrix @ L.matrix for L in jumps)
S0 = channels.Channel(M3, torch.linalg.matrix_exp(-0.8j * H_eff)[:, None]).superoperator()[0][0]
jump_part = torch.linalg.eigvalsh(choi(S(0.8) - S0))
print(f"smallest eigenvalue of the Choi matrix of the jump part: {jump_part.min().item():.1e}")
assert jump_part.min() > -1e-12                                    # sum_{n >= 1} S_n(t) is CP
```

## 8.3 Ecuaciones de Lindblad en un álgebra con sectores

Sobre $\mathcal M$ hay dos tipos de saltos. Un salto **dentro de los sectores** es un operador
$L_j\in\mathcal M$. Un salto **entre sectores** es una aplicación CP $\Phi:\mathcal M\to\mathcal M$
con operadores de Kraus $K^{dc}_l:\mathbb C^{k_c}\to\mathbb C^{k_d}$ (proposición 6.4, un
`InterSectorChannel`): $\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l$. Por la proposición 6.8, su
dual `blunt` cumple $\Phi^*(1)_c=\sum_{d,l}K^{dc*}_lK^{dc}_l$, un elemento diagonal por bloques de
$\mathcal M$. Con tasas $g_j,g_m>0$, el generador de `lindblad_rhs` es
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_jg_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big)
+\sum_mg_m\Big(\Phi_m(\rho)-\tfrac12\{\Phi_m^*(1),\rho\}\Big).
$$

**Proposición 8.8.** Sea $J^{dc}_l\in M_D$ la matriz con el bloque $K^{dc}_l$ en la posición $(d,c)$
y ceros en el resto.

1. $\mathcal L$ lleva $\mathcal M$ en $\mathcal M$ y es la restricción a $\mathcal M$ del generador
   GKSL $\hat{\mathcal L}$ sobre $M_D$ con saltos $\sqrt{g_j}L_j$ y $\sqrt{g_m}J^{dc}_{m,l}$ (todos
   los $m,l,d,c$). Por tanto, $(e^{t\mathcal L})_{t\ge0}$ es un semigrupo dinámico cuántico sobre
   $\mathcal M$.
2. Sea en cambio $\tilde{\mathcal L}$ el generador GKSL sobre $M_D$ en el que los bloques de un
   mismo índice de Kraus se suman en un único salto, $L_{m,l}=\sum_{d,c}J^{dc}_{m,l}$. Entonces
   $\mathcal L(\rho)=E\big(\tilde{\mathcal L}(\rho)\big)$ para $\rho\in\mathcal M$, donde
   $E(X)=\sum_dP_dXP_d$ es la esperanza condicional de $M_D$ sobre $\mathcal M$; en general,
   $\tilde{\mathcal L}(\rho)\notin\mathcal M$.

*Demostración.* Como $\rho=\sum_cP_c\rho P_c$ y $J^{dc}P_{c'}=\delta_{cc'}J^{dc}$,
$J^{dc}\rho J^{dc*}$ es el bloque $K^{dc}\rho_cK^{dc*}$ en la posición $(d,d)$, y $J^{dc*}J^{dc}$ es
$K^{dc*}K^{dc}$ en la posición $(c,c)$. Sumando sobre $d,l$ se obtienen $\Phi(\rho)$ y $\Phi^*(1)$,
lo que prueba la identidad de (1), y todos los términos están en $\mathcal M$. Así,
$\hat{\mathcal L}$ y todas sus potencias llevan $\mathcal M$ en $\mathcal M$,
$e^{t\hat{\mathcal L}}|_{\mathcal M}=e^{t\mathcal L}$, y esto es un canal de $\mathcal M$ por el
teorema 8.6 y el lema 8.4(a). (2) $L\rho L^*=\sum_{d,d',c}J^{dc}\rho J^{d'c*}$ para
$L=\sum_{d,c}J^{dc}$, y $E$ conserva los términos con $d=d'$. Para $B\in\mathcal M$,
$P_dXBP_d=P_dXP_dB$, luego $E(XB)=E(X)B$ y $E(BX)=BE(X)$; así $E\{L^*L,\rho\}=\{E(L^*L),\rho\}$, y
$E(L^*L)=\sum_cP_cL^*LP_c=\sum_{c,d}J^{dc*}J^{dc}$ porque $J^{dc*}J^{d'c}=0$ para $d\ne d'$. Por
último, $E(-i[H,\rho])=-i[H,\rho]$. $\square$

> **Significado físico.** La parte (2) dice que la ecuación de Lindblad sobre $\mathcal M$ es lo que
> queda de una ecuación sobre $M_D$ cuando las coherencias entre sectores que crean los saltos se
> borran de inmediato: el contenido de una regla de superselección. En la pérdida de partículas no
> se borra nada: el operador de aniquilación $a_j$ reduce el número de partículas exactamente en
> uno, de modo que $a_j\rho a_j^*$ y $a_j^*a_j=n_j$ son diagonales por bloques, y la ecuación sobre
> $\mathcal M$ es la usual.

**Proposición 8.9 (pesos de los sectores).** Sea $\rho(t)=e^{t\mathcal L}\rho_0$ y $p(t)$ sus pesos
de los sectores.

1. Si no hay saltos entre sectores, cada $p_c(t)$ es constante.
2. Si para todo $m$ y todo par $d\ne c$ se tiene $\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}=r_{m,dc}1_{k_c}$
   con números $r_{m,dc}$, entonces $p'=Qp$ con $Q_{dc}=\sum_mg_mr_{m,dc}$ para $d\ne c$ y
   $Q_{cc}=-\sum_{d\ne c}Q_{dc}$. Por tanto, $p(t)=e^{tQ}p(0)$ no depende de $H$, ni de los $L_j$,
   ni de los estados dentro de los sectores.

*Demostración.* $p_d'=\operatorname{Tr}(P_d\mathcal L(\rho))$. El término hamiltoniano se anula como
en el teorema 8.3. Para $L_j\in\mathcal M$,
$\operatorname{Tr}(P_dL_j\rho L_j^*)-\operatorname{Tr}(P_dL_j^*L_j\rho) =\operatorname{Tr}(L_j^*P_dL_j\rho)-\operatorname{Tr}(L_j^*L_jP_d\rho)=0$,
pues $P_d$ es central. Para un salto $\Phi$ entre sectores,
$\operatorname{Tr}(P_d\Phi(\rho))=\sum_c\operatorname{Tr}\big(\sum_lK^{dc*}_lK^{dc}_l\rho_c\big)$ y
$\operatorname{Tr}(P_d\tfrac12\{\Phi^*(1),\rho\})=\sum_{d'}\operatorname{Tr}\big(\sum_lK^{d'd*}_lK^{d'd}_l\rho_d\big)$;
los términos $c=d$ y $d'=d$ se cancelan. Esto prueba (1), y bajo la hipótesis de (2) los términos
restantes son $\sum_{c\ne d}r_{dc}p_c-\sum_{d'\ne d}r_{d'd}p_d$. El lema 8.1 da $p(t)=e^{tQ}p(0)$.
$\square$

Bajo la hipótesis de (2), para $d\ne c$, $r_{m,dc}$ es la entrada $T_{dc}(\rho)$ de la matriz de
transición de $\Phi_m$ (definida antes de la proposición 6.12) para todo estado $\rho$. Así, para un
único salto con tasa $1$, $Q=R-\operatorname{diag}(\sum_dR_{dc})_c$ con $R$ = `transition_matrix()`
(la diagonal de $R$ se cancela).

**Ejemplo: pérdida de partículas.** Dos modos bosónicos con a lo sumo dos partículas; los sectores
son $N=0,1,2$ con bases $|00\rangle$; $|10\rangle,|01\rangle$; $|20\rangle,|11\rangle,|02\rangle$,
de modo que $\mathcal M=M_1\oplus M_2\oplus M_3$. El dímero de Bose–Hubbard
$H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)$ conmuta con $N$, y cada modo pierde partículas con
tasa $\gamma$: un salto entre sectores con operadores de Kraus $a_1,a_2$, restringidos a $N\to N-1$.
Como $a_1^*a_1+a_2^*a_2=N$, se aplica la proposición 8.9(2) con $Q_{N-1,N}=\gamma N$ y, partiendo de
$N=2$,
$$
p_2(t)=q^2,\qquad p_1(t)=2q(1-q),\qquad p_0(t)=(1-q)^2,\qquad q=e^{-\gamma t},
$$
una ley binomial: cada partícula sobrevive independientemente con probabilidad $e^{-\gamma t}$,
hagan lo que hagan el efecto túnel entre modos (hopping) y la interacción (compruébese $p'=Qp$
derivando).

```python
fock = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double', charges=[0, 1, 2])   # N = 0, 1, 2
r2 = math.sqrt(2)
a1 = {(0, 1): [[1.0, 0.0]], (1, 2): [[r2, 0.0, 0.0], [0.0, 1.0, 0.0]]}           # a_1 : N -> N - 1
a2 = {(0, 1): [[0.0, 1.0]], (1, 2): [[0.0, 1.0, 0.0], [0.0, 0.0, r2]]}           # a_2 : N -> N - 1
blocks = lambda a: {dc: torch.tensor(K, dtype=f64) for dc, K in a.items()}
loss = InterSectorChannel.from_blocks(fock, fock, {dc: [blocks(a1)[dc], blocks(a2)[dc]] for dc in a1})
J, U_int = 1.0, 0.8                                                # hopping and interaction
H_bh = fock.from_blocks([torch.zeros(1, 1, dtype=f64), torch.tensor([[0.0, J], [J, 0.0]], dtype=f64),
                         torch.tensor([[U_int, r2 * J, 0.0], [r2 * J, 0.0, r2 * J], [0.0, r2 * J, U_int]], dtype=f64)])

assert (loss.adjoint()(fock.identity()).matrix - fock.central([0.0, 1.0, 2.0]).matrix).abs().max() < 1e-12
R = loss.transition_matrix()[0]                                    # r_dc (columns: from sector c)
Q = R - torch.diag(R.sum(0))                                       # Proposition 8.9(2), for rate 1
print(Q)

gamma = 0.5
psi2 = torch.randn(3, dtype=c128)
psi2 = psi2 / psi2.norm()
m0 = torch.zeros(1, 3, 3, 3, dtype=c128)
m0[0, 2] = torch.outer(psi2, psi2.conj())                          # a pure state with N = 2
ts = [0.0, 0.5, 1.0, 2.0, 4.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m0), H_bh, [loss], ts, rates=[gamma], substeps=100)
for s, rho_s in zip(ts, traj):
    q = math.exp(-gamma * s)
    p = rho_s.sector_probabilities()[0]
    assert (p - torch.tensor([(1 - q) ** 2, 2 * q * (1 - q), q ** 2], dtype=f64)).abs().max() < 1e-8
    assert (p - torch.from_numpy(sla.expm(gamma * s * Q.numpy())[:, 2])).abs().max() < 1e-8
neg = min(w.min().item() for rho_s in traj for w in rho_s.eigenvalues())
print(f"smallest eigenvalue of rho(t) under RK4: {neg:.1e}")
assert neg > -1e-7                                                 # positive up to the RK4 error
```

RK4 no es exactamente completamente positivo: los valores propios de $\rho(t)$ pueden ser negativos
al nivel del error de discretización (aquí $\approx-2\cdot10^{-8}$, con pasos de hasta $0.02$);
`lindblad_channel` es exacto, pero solo trata saltos dentro de los sectores. Por último, comprobamos
la proposición 8.8(2) para una aplicación aleatoria entre los sectores, frente al generador GKSL
sobre $M_6$ construido con numpy:

```python
Phi = channels.random_inter_sector_channel(fock, fock, kraus_rank=2)
rho = states.random_density_matrix(fock)
on_M = dynamics.lindblad_rhs(None, [Phi], algebra=fock)(0.0, rho.matrix)[0].numpy()
ks, K = fock.k_factors, Phi.kraus[0].numpy()                       # K[l, d, c] = K_l^{dc}, zero padded
embed = lambda Y: sla.block_diag(*[Y[c, :k, :k] for c, k in enumerate(ks)])    # (+)_c Y_c in M_6
Ls = [np.block([[K[l, d, c, :kd, :kc] for c, kc in enumerate(ks)] for d, kd in enumerate(ks)])
      for l in range(2)]                                           # L_l = sum_{d,c} J_l^{dc}
X = embed(rho.matrix[0].numpy())
full = sum(L @ X @ L.conj().T - 0.5 * (L.conj().T @ L @ X + X @ L.conj().T @ L) for L in Ls)   # GKSL on M_6
E_full = full * embed(np.ones((3, 3, 3)))                          # E(X) = sum_d P_d X P_d
assert np.abs(E_full - embed(on_M)).max() < 1e-12                  # Proposition 8.8(2)
print(f"coherences between sectors created on M_6: {np.abs(full - E_full).max():.2f}")
```

## 8.4 Postselección y evolución sin saltos

Reunamos los saltos del generador de la sección 8.3 en
$\Phi_{\mathrm{tot}}(\rho)=\sum_jg_jL_j\rho L_j^*+\sum_mg_m\Phi_m(\rho)$, con
$A=\Phi_{\mathrm{tot}}^*(1)=\sum_jg_jL_j^*L_j+\sum_mg_m\Phi_m^*(1)\ge0$ en $\mathcal M$, y
$H_{\mathrm{eff}}=H-\frac i2A$. Por la proposición 8.8(1) y el teorema 8.6 (cuyos $S_n$ llevan
$\mathcal M$ en $\mathcal M$, porque $H_{\mathrm{eff}}\in\mathcal M$),
$$
\rho(t)=\sum_{n\ge0}S_n(t)\rho,\qquad S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},
\qquad S_n(t)=\int_0^tS_0(t-s)\,\Phi_{\mathrm{tot}}\,S_{n-1}(s)\,ds .
$$

**Proposición 8.10.** Para un estado $\rho$, sea $w_n(t)=\operatorname{Tr}S_n(t)\rho$. Entonces
$w_n(t)\ge0$, $\sum_nw_n(t)=1$, y $w_0$ es no creciente con
$w_0'(t)=-\operatorname{Tr}\big(A\,S_0(t)\rho\big)$. Para un estado puro
$\rho=|\psi\rangle\langle\psi|$, $w_0(t)=\|e^{-itH_{\mathrm{eff}}}\psi\|^2$.

*Demostración.* $S_n(t)$ es CP, luego positiva, y $\sum_nw_n=\operatorname{Tr}\rho(t)=1$. Con
$X=S_0(t)\rho\ge0$,
$w_0'=\operatorname{Tr}(-iH_{\mathrm{eff}}X+iXH_{\mathrm{eff}}^*)=-\operatorname{Tr}(AX)$, y
$\operatorname{Tr}(AX)=\operatorname{Tr}(X^{1/2}AX^{1/2})\ge0$. La última afirmación es
$\operatorname{Tr}|\phi\rangle\langle\phi|=\|\phi\|^2$. $\square$

> **Significado físico.** Supongamos que un detector registra cada salto: cada fotón emitido, cada
> átomo perdido. Entonces $w_n(t)$ es la probabilidad de exactamente $n$ clics hasta el instante
> $t$, $S_n(t)\rho/w_n(t)$ el estado condicionado a $n$ clics, y $S_0(t)\rho/w_0(t)$ el estado
> condicionado a **ningún clic**: un resultado nulo también cambia el estado, a través del
> $H_{\mathrm{eff}}$ no hermítico. Esta interpretación como proceso de conteo [E. B. Davies,
> *Quantum Theory of Open Systems*, 1976; M. D. Srinivas y E. B. Davies, *Optica Acta* 28, 981
> (1981)] es la base del método de las trayectorias cuánticas (función de onda de Monte Carlo) [J.
> Dalibard, Y. Castin, K. Mølmer, *Phys. Rev. Lett.* 68, 580 (1992)]; aquí solo usamos su parte
> matemática.

En el caso de las pérdidas, contar los clics es lo mismo que medir la carga:

**Proposición 8.11.** Etiquetemos los sectores por $N=0,\dots,C-1$, supongamos que no hay saltos
dentro de los sectores, que todo salto reduce la carga en uno ($K^{dc}_{m,l}=0$ salvo si $d=c-1$) y
que $\rho=P_{N_0}\rho P_{N_0}$. Entonces $S_n(t)\rho=P_{N_0-n}\,\rho(t)\,P_{N_0-n}$ para
$0\le n\le N_0$ y $S_n(t)\rho=0$ para $n>N_0$. En particular, la actualización de Lüders de
$\rho(t)$ por $P_{N_0}$ ("no se perdió ninguna partícula") es el estado sin saltos:
$$
\frac{P_{N_0}\rho(t)P_{N_0}}{p_{N_0}(t)}=\frac{e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*}}{w_0(t)},
\qquad p_{N_0}(t)=w_0(t).
$$

*Demostración.* $S_0(s)$ lleva operadores con soporte en el sector $N$ a operadores con soporte en
el sector $N$, ya que $H_{\mathrm{eff}}\in\mathcal M$; $\Phi_{\mathrm{tot}}$ lleva el sector $N\ge1$
al sector $N-1$ y anula el sector $0$. Por inducción, $S_n(s)\rho$ tiene soporte en el sector
$N_0-n$ para todo $s$, y es $0$ para $n>N_0$. Comprimir $\rho(t)=\sum_nS_n(t)\rho$ con $P_{N_0-n}$
conserva exactamente el $n$-ésimo término. $\square$

La función `states.lueders_update(rho, P)` de la biblioteca devuelve
$P\rho P/\operatorname{Tr}(P\rho P)$ y $\operatorname{Tr}(\rho P)$. Supongamos ahora que los dos
modos pierden partículas con tasas distintas $\gamma_1>\gamma_2$; entonces
$A=\gamma_1n_1+\gamma_2n_2$ no es central, la proposición 8.9(2) no se aplica, y la postselección en
$N=2$ da una evolución genuinamente no hermítica. Comprobamos la proposición 8.11 para $n=0$ y,
calculando la integral de $S_1$ con la regla de Simpson, para $n=1$:

```python
g1, g2 = 0.6, 0.2                                                  # mode 1 leaks faster
loss1 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a1).items()})
loss2 = InterSectorChannel.from_blocks(fock, fock, {dc: [K] for dc, K in blocks(a2).items()})
m1 = torch.zeros(1, 3, 3, 3, dtype=c128)
m1[0, 2, 0, 0] = 1.0                                               # |20>: both particles in mode 1
ts = [0.0, 0.5, 1.0, 2.0]
traj = dynamics.lindblad_evolve(DensityMatrix(fock, matrix=m1), H_bh, [loss1, loss2], ts, rates=[g1, g2],
                                substeps=100)

n_1 = [np.diag(x) for x in ([0.0], [1.0, 0.0], [2.0, 1.0, 0.0])]   # n_1 in each sector; n_2 = N - n_1
H_eff = [H_bh.matrix[0, N, :N + 1, :N + 1].numpy() - 0.5j * (g1 * n + g2 * (N * np.eye(N + 1) - n))
         for N, n in enumerate(n_1)]

def S0(s, N, X):                                  # e^{-isH_eff} X e^{isH_eff^*} in sector N, s a time or an array
    E = torch.linalg.matrix_exp(torch.from_numpy(-1j * np.asarray(s)[..., None, None] * H_eff[N])).numpy()
    return E @ X @ E.conj().swapaxes(-1, -2)

rho_20 = m1[0, 2].numpy()
for s, rho_s in zip(ts, traj):
    w0 = np.trace(S0(s, 2, rho_20)).real                           # probability of no loss
    post, prob = states.lueders_update(rho_s, fock.central([0.0, 0.0, 1.0]))   # post-select N = 2
    assert abs(prob.item() - w0) < 1e-8                            # Proposition 8.11, n = 0
    assert np.abs(post.matrix[0, 2].numpy() - S0(s, 2, rho_20) / w0).max() < 1e-7
closed = dynamics.von_neumann(DensityMatrix(fock, matrix=m1), H_bh, [2.0])[0]   # no losses at all
dist = 0.5 * np.abs(np.linalg.eigvalsh(S0(2.0, 2, rho_20) / w0 - closed.matrix[0, 2].numpy())).sum()
print(f"t = 2: P(no loss) = {w0:.3f} (without hopping: {math.exp(-2 * g1 * 2):.3f}); "
      f"trace distance post-selected/closed system: {dist:.3f}")

K1, K2 = blocks(a1)[(1, 2)].numpy(), blocks(a2)[(1, 2)].numpy()   # S_1(2) rho = int_0^2 S_0(2-s) Phi S_0(s) rho ds
jump = lambda X: g1 * K1 @ X @ K1.T + g2 * K2 @ X @ K2.T
s_grid = np.linspace(0.0, 2.0, 201)
S1 = scipy.integrate.simpson(S0(2.0 - s_grid, 1, jump(S0(s_grid, 2, rho_20))), x=s_grid, axis=0)
assert np.abs(S1 - traj[-1].matrix[0, 1, :2, :2].numpy()).max() < 1e-8     # Proposition 8.11, n = 1
```

En $t=2$ la probabilidad de no perder ninguna partícula es $0.230$; sin efecto túnel sería
$e^{-2\gamma_1t}=0.091$, así que aquí los pesos de los sectores sí dependen del hamiltoniano. El
estado postseleccionado está a distancia de traza $0.135$ del estado del sistema cerrado: no perder
ninguna partícula es información.

## 8.5 Estados estacionarios y relajación

Un estado $\sigma$ es **estacionario** si $e^{t\mathcal L}\sigma=\sigma$ para todo $t\ge0$, o,
equivalentemente, si $\mathcal L(\sigma)=0$ (derívese en $t=0$; recíprocamente, úsese la serie).

**Proposición 8.12.** Todo semigrupo dinámico cuántico sobre $\mathcal M$ tiene un estado
estacionario.

*Demostración.* Fijemos un estado $\rho$ y sea $\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds$. Es
un estado (un límite de combinaciones convexas de estados; los estados forman un conjunto cerrado).
Por linealidad,
$T_t\sigma_\tau-\sigma_\tau=\frac1\tau\big(\int_\tau^{\tau+t}-\int_0^t\big)T_s\rho\,ds$, cuya norma
traza es a lo sumo $2t/\tau$, ya que $\|T_s\rho\|_1=1$. Los estados forman un conjunto cerrado y
acotado en un espacio de dimensión finita, así que alguna sucesión $\sigma_{\tau_m}$,
$\tau_m\to\infty$, converge a un estado $\sigma$ (Bolzano–Weierstrass), y
$T_t\sigma-\sigma=\lim_m(T_t\sigma_{\tau_m}-\sigma_{\tau_m})=0$. $\square$

**Corolario 8.13.** Si no hay saltos entre sectores, los estados estacionarios generan un espacio de
dimensión al menos $C$: cada sector tiene su propio estado estacionario.

*Demostración.* $\mathcal L$ lleva $P_c\mathcal MP_c\cong M_{k_c}$ en sí mismo ($H$ y todos los
$L_j$ son diagonales por bloques), así que la restricción de $e^{t\mathcal L}$ es un semigrupo
dinámico cuántico sobre $M_{k_c}$ (lema 8.4(a)) con un estado estacionario $\sigma_c$ por la
proposición 8.12. Los $\sigma_c$ tienen soportes disjuntos. $\square$

**Proposición 8.14 (emisión espontánea).** Sobre $M_2$ con base $|g\rangle,|e\rangle$, sean
$H=\omega|e\rangle\langle e|$, un salto $L=|g\rangle\langle e|$ y una tasa $\gamma>0$. Entonces
$$
\rho_{ee}(t)=e^{-\gamma t}\rho_{ee}(0),\qquad\rho_{ge}(t)=e^{(i\omega-\gamma/2)t}\rho_{ge}(0),
\qquad e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}},
$$
donde $\mathcal U_t(\rho)=e^{-itH}\rho e^{itH}$ y $\mathcal A_p$ es el canal de amortiguamiento de
amplitud de parámetro $p$ (sección 6.4). El espectro de $\mathcal L$ es
$\{0,-\gamma,-\gamma/2\pm i\omega\}$, $|g\rangle\langle g|$ es el único estado estacionario, y
$\tfrac12\|\rho(t)-|g\rangle\langle g|\|_1=\big(\rho_{ee}(0)^2e^{-2\gamma t}+|\rho_{ge}(0)|^2e^{-\gamma t}\big)^{1/2}$.

*Demostración.* Aquí $L^*L=|e\rangle\langle e|$ y $L\rho L^*=\rho_{ee}|g\rangle\langle g|$, luego,
entrada a entrada, $\rho_{ee}'=-\gamma\rho_{ee}$, $\rho_{gg}'=\gamma\rho_{ee}$ y
$\rho_{ge}'=(i\omega-\gamma/2)\rho_{ge}$, lo que da la solución. En la base
$(\rho_{gg},\rho_{ge},\rho_{eg},\rho_{ee})$, la matriz de $\mathcal L$ es triangular superior con
diagonal $(0,-\gamma/2+i\omega,-\gamma/2-i\omega,-\gamma)$. $\mathcal L(\sigma)=0$ obliga a
$\sigma_{ee}=\sigma_{ge}=0$. El canal de amortiguamiento de amplitud multiplica $\rho_{ee}$ por
$1-p$ y $\rho_{ge}$ por $\sqrt{1-p}$ (sus operadores de Kraus son
$|g\rangle\langle g|+\sqrt{1-p}\,|e\rangle\langle e|$ y $\sqrt p\,|g\rangle\langle e|$), y
$\mathcal U_t$ multiplica $\rho_{ge}$ por $e^{i\omega t}$; ambos son TP, lo que fija $\rho_{gg}$.
Por último, $\rho(t)-|g\rangle\langle g|$ es autoadjunta, de traza nula y con diagonal
$(-\rho_{ee}(t),\rho_{ee}(t))$, así que sus valores propios son
$\pm(\rho_{ee}(t)^2+|\rho_{ge}(t)|^2)^{1/2}$. $\square$

La coherencia decae a la mitad de la tasa de la población, la relación $T_2=2T_1$ de la sección 6.4
sin desfase puro; la tasa no nula más lenta, la **brecha espectral** $\gamma/2$, controla la
relajación. Comprobamos el espectro con `lindblad_superoperator`, el canal y la distancia de traza
con `lindblad_channel`, y el corolario 8.13 para saltos aleatorios dentro de los sectores de
$M_2\oplus M_3$:

```python
q2 = TypeIAlgebra([2], [2], precision='double')                    # basis |g>, |e>
omega, gamma = 1.3, 0.4
H_at = q2.from_blocks([torch.tensor([[0.0, 0.0], [0.0, omega]], dtype=f64)])
sigma_minus = q2.from_blocks([torch.tensor([[0.0, 1.0], [0.0, 0.0]], dtype=f64)])   # |g><e|
gen = dynamics.lindblad_superoperator(q2, H_at, [sigma_minus], rates=[gamma])[0][0]
spec = sorted(torch.linalg.eigvals(gen).tolist(), key=lambda z: (z.real, z.imag))
print([f"{z.real:+.2f}{z.imag:+.2f}i" for z in spec])
assert np.allclose(spec, [-gamma, -gamma / 2 - 1j * omega, -gamma / 2 + 1j * omega, 0.0], atol=1e-12)

psi_at = torch.tensor([0.6, 0.8], dtype=c128)                      # rho_ee(0) = 0.64, rho_ge(0) = 0.48
rho_at = DensityMatrix(q2, matrix=torch.outer(psi_at, psi_at.conj())[None, None])
ground = DensityMatrix(q2, matrix=torch.tensor([[[[1.0, 0.0], [0.0, 0.0]]]], dtype=c128))
for s in (0.5, 2.0, 10.0):
    ch = dynamics.lindblad_channel(q2, H_at, [sigma_minus], t=s, rates=[gamma])
    UA = channels.unitary_channel(dynamics.propagator(H_at, s)).compose(
        channels.amplitude_damping_channel(q2, 1 - math.exp(-gamma * s)))
    assert (ch.superoperator()[0] - UA.superoperator()[0]).abs().max() < 1e-12
    r = ch(rho_at).matrix[0, 0]
    assert abs(r[1, 1].real - 0.64 * math.exp(-gamma * s)) < 1e-12
    assert abs(r[0, 1] - 0.48 * np.exp((1j * omega - gamma / 2) * s)) < 1e-12
    exact = math.sqrt(0.64 ** 2 * math.exp(-2 * gamma * s) + 0.48 ** 2 * math.exp(-gamma * s))
    assert abs(states.trace_distance(ch(rho_at), ground).item() - exact) < 1e-12

two = TypeIAlgebra([2, 3], [2, 3], precision='double')
H_two = two.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
L_two = [two.operator(torch.randn(1, 2, 3, 3, dtype=c128)) for _ in range(2)]
kernel = [int((torch.linalg.eigvals(G_c[0]).abs() < 1e-10).sum())
          for G_c in dynamics.lindblad_superoperator(two, H_two, L_two)]
print("stationary states per sector:", kernel)
assert all(dim >= 1 for dim in kernel)                             # Corollary 8.13
```

Para estos saltos aleatorios, el código encuentra exactamente un estado estacionario por sector. Los
criterios de unicidad —por ejemplo, un estado estacionario fiel junto con
$\{H,L_j,L_j^*\}'=\mathbb C1$ [A. Frigerio, *Commun. Math. Phys.* 63, 269 (1978)]— no se demuestran
aquí; el ejercicio 4 muestra que la unicidad puede fallar incluso en un factor.

## 8.6 Con qué comparar

QuTiP (`sesolve`, `mesolve` para la ecuación de Lindblad con integradores de EDO adaptativos,
`mcsolve` para trayectorias cuánticas, `steadystate`, `liouvillian`), dynamiqs (JAX; integradores de
Schrödinger y de Lindblad por lotes y diferenciables en GPU), QuantumOptics.jl y SciPy
(`scipy.linalg.expm`, y `scipy.sparse.linalg.expm_multiply` para $e^{t\mathcal L}\rho$ con un
generador vectorizado disperso) resuelven las mismas ecuaciones sobre un único espacio de Hilbert, y
para los problemas anteriores deben dar los mismos números. Sobre $\mathcal M$ se les pasan matrices
$D\times D$ diagonales por bloques; el generador vectorizado vive entonces en un espacio de
dimensión $D^2$ (disperso, pero con los sectores acoplados en una sola resolución), salvo que el
usuario separe el problema en sectores a mano. La biblioteca incorpora la estructura de bloques:
propagadores exactos y canales de Lindblad sector a sector (superoperadores $k_c^2\times k_c^2$),
saltos entre sectores como `InterSectorChannel` con la compresión de la proposición 8.8, lotes de
hamiltonianos y de operadores de salto (para promedios sobre el desorden) y diferenciación
automática (estas dos últimas también en dynamiqs). No tiene control adaptativo del paso (RK4 es de
paso fijo y solo aproximadamente positivo), ni integrador de trayectorias, y `lindblad_channel` no
admite saltos entre sectores; para ello, la referencia son QuTiP o dynamiqs sobre $M_D$.

## Ejercicios

1. (Imagen de Heisenberg) Demuestre que
   $\operatorname{Tr}(\mathcal L(\rho)A)=\operatorname{Tr}(\rho\,\mathcal L^*(A))$ para el generador
   de la sección 8.3, con
   $\mathcal L^*(A)=i[H,A]+\sum_jg_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})+\sum_mg_m(\Phi_m^*(A)-\frac12\{\Phi_m^*(1),A\})$,
   que $\mathcal L^*(1)=0$, y que $\mathcal L^*(P_c)=0$ si no hay saltos entre sectores.
2. (Punto excepcional) En el sector $N=1$ de la sección 8.4,
   $H_{\mathrm{eff}}=\begin{pmatrix}-i\gamma_1/2&J\\J&-i\gamma_2/2\end{pmatrix}$. Halle sus valores
   propios y demuestre que para $J=|\gamma_1-\gamma_2|/4>0$ no es diagonalizable.
3. Supongamos que $N_0$ bosones en $m$ modos se pierden cada uno con tasa $\gamma$ (saltos
   $a_1,\dots,a_m$ restringidos a $N\to N-1$). Demuestre que $p_N(t)=\binom{N_0}Nq^N(1-q)^{N_0-N}$
   con $q=e^{-\gamma t}$, para todo hamiltoniano de $\mathcal M$.
4. (Desfase puro) Sobre $M_2$, sean $H=\omega|e\rangle\langle e|$ y
   $L=\sigma_z=\operatorname{diag}(1,-1)$ con tasa $\gamma$. Demuestre que las poblaciones son
   constantes, que $\rho_{ge}(t)=e^{(i\omega-2\gamma)t}\rho_{ge}(0)$, y que los estados
   estacionarios generan el espacio bidimensional de las matrices diagonales aunque $C=1$.

<details>
<summary>Soluciones</summary>

1. Por ciclicidad, $\operatorname{Tr}(L\rho L^*A)=\operatorname{Tr}(\rho L^*AL)$,
   $\operatorname{Tr}(\{L^*L,\rho\}A)=\operatorname{Tr}(\rho\{L^*L,A\})$ y
   $\operatorname{Tr}(-i[H,\rho]A)=\operatorname{Tr}(\rho\,i[H,A])$; para los saltos entre sectores,
   $\operatorname{Tr}(\Phi(\rho)A)=\operatorname{Tr}(\rho\Phi^*(A))$ es la definición de $\Phi^*$
   (sección 6.3). Entonces
   $\mathcal L^*(1)=\sum_jg_j(L_j^*L_j-L_j^*L_j)+\sum_mg_m(\Phi_m^*(1)-\Phi_m^*(1))=0$, que es una
   vez más la preservación de la traza. Sin saltos entre sectores, $P_c$ conmuta con $H$ y con todos
   los $L_j$, luego $L_j^*P_cL_j=L_j^*L_jP_c=\frac12\{L_j^*L_j,P_c\}$ y $\mathcal L^*(P_c)=0$;
   entonces $p_c'=\operatorname{Tr}(\rho\mathcal L^*(P_c))=0$, una segunda demostración de la
   proposición 8.9(1).

2. $H_{\mathrm{eff}}=-\frac i4(\gamma_1+\gamma_2)1+\begin{pmatrix}-i\delta&J\\J&i\delta\end{pmatrix}$
   con $\delta=(\gamma_1-\gamma_2)/4$; la segunda matriz tiene traza $0$ y determinante
   $\delta^2-J^2$, así que los valores propios son
   $-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}$: dos tasas de decaimiento para $J<|\delta|$,
   una tasa y dos frecuencias para $J>|\delta|$. Para $J=|\delta|>0$, la segunda matriz $Y$ cumple
   $Y^2=(J^2-\delta^2)1=0$ pero $Y\ne0$, de modo que $H_{\mathrm{eff}}$ tiene el único valor propio
   $-\frac i4(\gamma_1+\gamma_2)$ y no es múltiplo de $1$: no es diagonalizable.

```python
lam = np.linalg.eigvals(H_eff[1])
delta = (g1 - g2) / 4
root = np.sqrt(J ** 2 - delta ** 2 + 0j)
assert np.allclose(sorted(lam, key=np.real), sorted(-0.25j * (g1 + g2) + np.array([root, -root]), key=np.real))
H_ep = np.array([[-0.5j * g1, delta], [delta, -0.5j * g2]])         # J = delta: the exceptional point
Y = H_ep + 0.25j * (g1 + g2) * np.eye(2)                           # H_ep minus its only eigenvalue
assert np.abs(Y @ Y).max() < 1e-15 and np.linalg.matrix_rank(Y) == 1   # one Jordan block of size 2
```

3. En el sector $N$, $\sum_ja_j^*a_j=N1$, y el único par de sectores es $(N-1,N)$, así que la
   proposición 8.9(2) da $p_N'=-\gamma Np_N+\gamma(N+1)p_{N+1}$. Con $q'=-\gamma q$, la derivada del
   término binomial es $-\gamma Np_N+\gamma(N_0-N)\binom{N_0}Nq^{N+1}(1-q)^{N_0-N-1}$, y
   $(N+1)\binom{N_0}{N+1}=(N_0-N)\binom{N_0}N$ muestra que el segundo término es
   $\gamma(N+1)p_{N+1}$. En $t=0$ es la masa puntual en $N_0$; la unicidad es el lema 8.1.

4. $L\rho L^*-\frac12\{L^*L,\rho\}=\sigma_z\rho\sigma_z-\rho$ multiplica las entradas fuera de la
   diagonal por $-2$ y anula la diagonal, y $-i[H,\rho]_{ge}=i\omega\rho_{ge}$. Así,
   $\mathcal L(\sigma)=0$ si y solo si $\sigma_{ge}=0$ (para $\gamma>0$), y el espectro es
   $\{0,0,-2\gamma\pm i\omega\}$. Como $H$ y $\sigma_z$ conmutan y $X\mapsto\sigma_zX\sigma_z$ es
   una involución,
   $e^{t\mathcal L}(\rho)=\frac{1+e^{-2\gamma t}}2U_t\rho U_t^* +\frac{1-e^{-2\gamma t}}2U_t\sigma_z\rho\sigma_zU_t^*$
   con $U_t=e^{-itH}$: un canal unital cuyos puntos fijos son las matrices que conmutan con sus
   operadores de Kraus (teorema 6.10), es decir, las diagonales.

```python
deph = dynamics.lindblad_superoperator(q2, H_at, [q2.from_blocks([torch.diag(torch.tensor([1.0, -1.0], dtype=f64))])],
                                       rates=[gamma])[0][0]
ev = torch.linalg.eigvals(deph)
assert int((ev.abs() < 1e-12).sum()) == 2
assert np.allclose(sorted(ev.tolist(), key=lambda z: (z.real, z.imag))[:2],
                   [-2 * gamma - 1j * omega, -2 * gamma + 1j * omega], atol=1e-12)
```

</details>

**Siguiente:** el [capítulo 9](09_spin_chains.md) trata las cadenas de espines: sectores de
magnetización y de momento, evolución temporal de Lanczos y de Krylov para sectores grandes, y
desorden.
