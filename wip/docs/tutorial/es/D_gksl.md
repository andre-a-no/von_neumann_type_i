[English](../D_gksl.md) · **Español** · [Français](../fr/D_gksl.md) · [Deutsch](../de/D_gksl.md) · [中文](../zh/D_gksl.md) · [日本語](../ja/D_gksl.md) · [Русский](../ru/D_gksl.md)

# Apéndice D. Generadores de semigrupos dinámicos cuánticos

El capítulo 8 demostró que un generador en forma de Gorini–Kossakowski–Sudarshan–Lindblad (GKSL) produce un
semigrupo dinámico cuántico (teorema 8.6) y usó el recíproco (teorema 8.7) y su versión para
$\mathcal M=\bigoplus_cM_{k_c}$. Este apéndice demuestra ambos en dimensión finita, junto con la
libertad en la forma GKSL y un criterio de unicidad del estado estacionario. La herramienta es la
matriz de Choi del capítulo 6: un generador en general no es completamente positivo, pero su matriz de Choi es
positiva sobre el complemento ortogonal del vector máximamente entrelazado. El mismo tipo de argumento —positividad
de una familia de operadores y después una descomposición— demuestra en la última sección el teorema sobre la
corrección de errores en álgebras de operadores enunciado en el capítulo 10 (teorema 10.7), con un canal de
recuperación explícito.

**Aprenderá** por qué el generador de un semigrupo completamente positivo es condicionalmente completamente
positivo, cómo leer $H$ y los operadores de salto en la matriz de Choi de un generador, en qué sentido
son únicos, por qué todo generador sobre un álgebra con sectores tiene la forma usada en la sección 8.3, cuándo
es único el estado estacionario y qué álgebras de observables puede proteger un canal de recuperación.

En todo el apéndice, $M_n=M_n(\mathbb C)$, los $E_{ij}=|i\rangle\langle j|$ son sus unidades matriciales, un
**superoperador** es una aplicación lineal sobre $M_n$ (o sobre $\mathcal M$), y
$$
\Omega=\sum_{i=1}^ne_i\otimes e_i\in\mathbb C^n\otimes\mathbb C^n,\qquad
J(\Psi)=\sum_{i,j=1}^nE_{ij}\otimes\Psi(E_{ij})=[\Psi(E_{ij})]_{i,j=1}^n
$$
son el vector máximamente entrelazado y la **matriz de Choi** de un superoperador $\Psi$ (teorema 6.2;
aquí $\Psi$ no tiene por qué ser positivo). Para vectores $v,w$, $|v\rangle\langle w|$ es el operador
$x\mapsto\langle w,x\rangle v$. Un superoperador $\Psi$ **preserva la hermiticidad** si
$\Psi(A^*)=\Psi(A)^*$ para todo $A$. Los semigrupos dinámicos cuánticos y sus generadores son los de la
sección 8.2 (lema 8.5): $T_t=e^{t\mathcal L}$ con canales $T_t$.

## D.1 Matrices de Choi de las piezas básicas

**Lema D.1.** Sea $\Psi$ un superoperador sobre $M_n$.

(a) $\Psi\mapsto J(\Psi)$ es una biyección lineal de los superoperadores sobre $M_n$ en $M_n(M_n)=M_{n^2}$.
(b) Todo $v\in\mathbb C^n\otimes\mathbb C^n$ es $v=(1\otimes K)\Omega$ para exactamente un $K\in M_n$, y
$\langle(1\otimes K)\Omega,(1\otimes L)\Omega\rangle=\operatorname{Tr}(K^*L)$. En particular,
$\langle\Omega,(1\otimes K)\Omega\rangle=\operatorname{Tr}K$.
(c) Para $K,L\in M_n$, el superoperador $A\mapsto KAL^*$ tiene la matriz de Choi
$|(1\otimes K)\Omega\rangle\langle(1\otimes L)\Omega|$. En particular, $J(\mathrm{id})=|\Omega\rangle\langle\Omega|$.
(d) $\Psi$ preserva la hermiticidad si y solo si $J(\Psi)$ es autoadjunta. Toda aplicación CP preserva la hermiticidad.

*Demostración.* (a) El bloque $(i,j)$ de $J(\Psi)$ es $\Psi(E_{ij})$ y las $E_{ij}$ forman una base de $M_n$, así que
$J$ es lineal e inyectiva, y toda matriz por bloques $[X_{ij}]$ es la matriz de Choi del superoperador
con $E_{ij}\mapsto X_{ij}$. (b) $(1\otimes K)\Omega=\sum_ie_i\otimes Ke_i$ tiene como $i$-ésimo bloque $Ke_i$, la
$i$-ésima columna de $K$, que determina $K$, y todo $v$ se obtiene de esta forma. Entonces
$\langle\sum_ie_i\otimes Ke_i,\sum_je_j\otimes Le_j\rangle=\sum_i\langle Ke_i,Le_i\rangle=\operatorname{Tr}(K^*L)$;
para la última afirmación, sustitúyase $(K,L)$ por $(1,K)$. (c) El bloque $(i,j)$ del segundo miembro es
$Ke_i(Le_j)^*=KE_{ij}L^*$. (d) El bloque $(i,j)$ de $J(\Psi)^*$ es $\Psi(E_{ji})^*$. Si $\Psi$
preserva la hermiticidad, $\Psi(E_{ji})^*=\Psi(E_{ji}^*)=\Psi(E_{ij})$. Recíprocamente, si
$\Psi(E_{ji})^*=\Psi(E_{ij})$ para todos $i,j$ y $A=\sum a_{ij}E_{ij}$, entonces
$\Psi(A)^*=\sum\bar a_{ij}\Psi(E_{ij})^*=\sum\bar a_{ij}\Psi(E_{ji})=\Psi(A^*)$. Una aplicación CP tiene una forma de Kraus
(teorema 6.2), y $(\sum_lK_lAK_l^*)^*=\sum_lK_lA^*K_l^*$. $\square$

**Lema D.2 (descomposición).** Sean $q_1,\dots,q_C$ vectores ortonormales de
$\mathbb C^n\otimes\mathbb C^n$, escritos $q_c=(1\otimes Q_c)\Omega$, sea $W=\operatorname{span}\{q_c\}$, y sea
$\Psi$ un superoperador que preserva la hermiticidad con $\langle v,J(\Psi)v\rangle\ge0$ para todo $v\perp W$.
Pongamos $J=J(\Psi)$ y $\xi_c=Jq_c-\frac12\sum_d\langle q_d,Jq_c\rangle q_d$, y sea $K_c\in M_n$ el definido por
$(1\otimes K_c)\Omega=\xi_c$. Entonces existen $L_1,\dots,L_r\in M_n$, $r=\operatorname{rank}(P_W^\perp JP_W^\perp)$
con $P_W^\perp$ la proyección sobre $W^\perp$, tales que
$$
\Psi(A)=\sum_{j=1}^rL_jAL_j^*+\sum_{c=1}^C\big(K_cAQ_c^*+Q_cAK_c^*\big)\qquad(A\in M_n),
$$
$\operatorname{Tr}(Q_c^*L_j)=0$ para todos $c,j$, y $\operatorname{Tr}(L_j^*L_k)=0$ para $j\neq k$.

*Demostración.* $J$ es autoadjunta por el lema D.1(d). Sea $P=\sum_c|q_c\rangle\langle q_c|$ la proyección sobre $W$
y $P^\perp=P_W^\perp=1-P$. Desarrollando $J=(P+P^\perp)J(P+P^\perp)$ y usando $PJP^\perp=PJ-PJP$ y
$P^\perp JP=JP-PJP$,
$$
J=P^\perp JP^\perp+\big(PJ+JP-PJP\big).
$$
El primer término es autoadjunto, y $\langle v,P^\perp JP^\perp v\rangle=\langle P^\perp v,JP^\perp v\rangle\ge0$
porque $P^\perp v\perp W$; luego $P^\perp JP^\perp\ge0$. Por el teorema 1.3,
$P^\perp JP^\perp=\sum_{j=1}^r|l_j\rangle\langle l_j|$, donde los $l_j$ son los vectores propios de los valores propios
no nulos $\mu_j$ multiplicados por $\sqrt{\mu_j}$: son ortogonales entre sí y están en la imagen de
$P^\perp JP^\perp$, que es ortogonal a $W$. Escribamos $l_j=(1\otimes L_j)\Omega$; por el lema D.1(b),
$\operatorname{Tr}(Q_c^*L_j)=\langle q_c,l_j\rangle=0$ y $\operatorname{Tr}(L_j^*L_k)=\langle l_j,l_k\rangle=0$ para $j\ne k$,
y por el lema D.1(c), $P^\perp JP^\perp$ es la matriz de Choi de $A\mapsto\sum_jL_jAL_j^*$. Para el segundo término,
$JP=\sum_c|Jq_c\rangle\langle q_c|$ y $PJP=\sum_{c,d}\langle q_d,Jq_c\rangle|q_d\rangle\langle q_c|$, así que
$$
\sum_c|\xi_c\rangle\langle q_c|=JP-\tfrac12PJP,\qquad\sum_c|q_c\rangle\langle\xi_c|=PJ-\tfrac12PJP
$$
(la segunda es la adjunta de la primera; $PJP$ es autoadjunta), y su suma es $PJ+JP-PJP$. Por el
lema D.1(c), $|\xi_c\rangle\langle q_c|$ y $|q_c\rangle\langle\xi_c|$ son las matrices de Choi de
$A\mapsto K_cAQ_c^*$ y $A\mapsto Q_cAK_c^*$. Así, los dos miembros de la identidad afirmada tienen la matriz de Choi $J$,
y el lema D.1(a) da la identidad. $\square$

## D.2 Positividad completa condicional

Un superoperador $\Psi$ sobre $M_n$ es **condicionalmente completamente positivo** (CCP) si
$\langle v,J(\Psi)v\rangle\ge0$ para todo $v\perp\Omega$. Por el lema D.1(b), estos $v$ son los vectores
$(1\otimes K)\Omega$ con $\operatorname{Tr}K=0$.

**Teorema D.3.** Para un superoperador $\mathcal L$ sobre $M_n$ son equivalentes:

1. $e^{t\mathcal L}$ es CP para todo $t\ge0$;
2. $\mathcal L$ preserva la hermiticidad y es condicionalmente completamente positivo;
3. existen una aplicación CP $\Phi$ sobre $M_n$ y $G\in M_n$ con $\mathcal L(A)=\Phi(A)+GA+AG^*$ para todo $A$.

En 3 se puede tomar $\Phi(A)=\sum_{j=1}^rL_jAL_j^*$ con $r\le n^2-1$, $\operatorname{Tr}L_j=0$ y
$\operatorname{Tr}(L_j^*L_k)=0$ para $j\neq k$, y $G$ con $\operatorname{Tr}G\in\mathbb R$.

*Demostración.* (1)$\Rightarrow$(2). Por el lema 8.1, $(e^{t\mathcal L}-\mathrm{id})/t\to\mathcal L$ cuando $t\to0$, y $J$ es
lineal sobre un espacio de dimensión finita, luego continua, así que
$J(e^{t\mathcal L})=J(\mathrm{id})+t\,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)$, donde la última matriz de Choi tiende a
$J(\mathcal L)$. Para $v\perp\Omega$, $\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0$ (lema D.1(c)), así que
para $t>0$, como $J(e^{t\mathcal L})\ge0$ (teorema 6.2),
$$
0\le\frac1t\langle v,J(e^{t\mathcal L})v\rangle=\big\langle v,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)v\big\rangle
\;\xrightarrow[t\to0]{}\;\langle v,J(\mathcal L)v\rangle .
$$
Cada $e^{t\mathcal L}$ es CP, luego preserva la hermiticidad (lema D.1(d)), así que
$\mathcal L(A^*)=\lim_t(e^{t\mathcal L}(A)^*-A^*)/t=\mathcal L(A)^*$.

(2)$\Rightarrow$(3). Apliquemos el lema D.2 con $C=1$ y $q_1=\Omega/\sqrt n$, es decir, $Q_1=1/\sqrt n$:
$\mathcal L(A)=\sum_jL_jAL_j^*+GA+AG^*$ con $G=K_1/\sqrt n$, y $\Phi=\sum_jL_j\cdot L_j^*$ es CP
(teorema 6.2). Además, $\operatorname{Tr}L_j=\sqrt n\operatorname{Tr}(Q_1^*L_j)=0$, $r\le\operatorname{rank}P_W^\perp=n^2-1$,
y por el lema D.1(b)
$\operatorname{Tr}G=\langle\Omega,\xi_1\rangle/\sqrt n=\langle q_1,\xi_1\rangle=\frac12\langle q_1,J(\mathcal L)q_1\rangle\in\mathbb R$,
pues $J(\mathcal L)$ es autoadjunta.

(3)$\Rightarrow$(1). Los pasos 1–4 de la demostración del teorema 8.6 solo usan que el generador es
$\mathcal L_0+\Phi$ con $\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)$ para alguna matriz
$H_{\mathrm{eff}}$ y $\Phi$ CP. Con $H_{\mathrm{eff}}=iG$, $\mathcal L_0(\rho)=G\rho+\rho G^*$, y esos
pasos muestran que $e^{t\mathcal L}=\sum_{m\ge0}S_m(t)$ es CP. $\square$

La positividad sola no bastaría. Sean $n\ge2$, $\theta$ la trasposición sobre $M_n$ y $\mathcal L=\theta-\mathrm{id}$. Como $\theta^2=\mathrm{id}$,
$$
e^{t\mathcal L}=e^{-t}\big(\cosh t\,\mathrm{id}+\sinh t\,\theta\big)=\tfrac{1+e^{-2t}}2\,\mathrm{id}+\tfrac{1-e^{-2t}}2\,\theta ,
$$
una combinación convexa de aplicaciones positivas que preservan la traza (proposición 6.1). Pero $J(\theta)=[E_{ji}]=\sum_{i,j}E_{ij}\otimes E_{ji}$ es el operador de intercambio $e_k\otimes e_l\mapsto e_l\otimes e_k$, así que lleva
$v=e_1\otimes e_2-e_2\otimes e_1$, que es ortogonal a $\Omega$, a $-v$ (como en la demostración de la proposición 6.1); con
$\|v\|^2=2$ y $\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0$,
$\langle v,J(\mathcal L)v\rangle=-2$: $\mathcal L$ no es CCP, y en efecto
$\langle v,J(e^{t\mathcal L})v\rangle=-(1-e^{-2t})<0$ para todo $t>0$. Un semigrupo de aplicaciones positivas que preservan
la traza no tiene por qué estar formado por canales.

Comprobamos el teorema D.3 en un generador GKSL con saltos arbitrarios sobre $M_3$, construido con `lindblad_superoperator`
(la matriz $S$ de $\mathcal L$ sobre $\operatorname{vec}\rho$ ordenado por filas, de modo que $J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$ como en la
sección 6.1), y en $\theta-\mathrm{id}$:

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

## D.3 El teorema GKSL

**Teorema D.4 (Gorini–Kossakowski–Sudarshan; Lindblad; teorema 8.7).** Un superoperador $\mathcal L$ sobre
$M_n$ es el generador de un semigrupo dinámico cuántico sobre $M_n$ si y solo si existen $H=H^*$ y
$L_1,\dots,L_r$ en $M_n$ con
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_{j=1}^r\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
$$
Se pueden elegir $r\le n^2-1$, $\operatorname{Tr}H=0$, $\operatorname{Tr}L_j=0$ y $\operatorname{Tr}(L_j^*L_k)=0$ para $j\neq k$.

*Demostración.* El «si» es el teorema 8.6. Recíprocamente, sea $(T_t)$ un semigrupo dinámico cuántico con generador
$\mathcal L$, de modo que $T_t=e^{t\mathcal L}$ (lema 8.5). Como todo $T_t$ es CP, el teorema D.3 da
$\mathcal L(\rho)=\sum_jL_j\rho L_j^*+G\rho+\rho G^*$ con $r\le n^2-1$ operadores $L_j$ de traza nula y ortogonales entre sí, y
$\operatorname{Tr}G\in\mathbb R$. Pongamos $H=\frac i2(G-G^*)$ y $R=\frac12(G+G^*)$, ambos autoadjuntos; entonces
$G=-iH+R$, $G^*=iH+R$ y $G\rho+\rho G^*=-i[H,\rho]+\{R,\rho\}$. Además,
$\operatorname{Tr}H=\frac i2(\operatorname{Tr}G-\overline{\operatorname{Tr}G})=0$. Como todo $T_t$ es TP,
$\operatorname{Tr}\mathcal L(\rho)=\frac d{dt}\operatorname{Tr}T_t(\rho)\big|_{t=0}=0$ (lema 8.1), y por ciclicidad
$$
0=\operatorname{Tr}\mathcal L(\rho)=\operatorname{Tr}\Big(\rho\,\Big(\sum_jL_j^*L_j+2R\Big)\Big)\qquad\text{for all }\rho .
$$
Con $\rho=X=\sum_jL_j^*L_j+2R$, que es autoadjunto, $\operatorname{Tr}(X^*X)=0$, luego $X=0$ y
$R=-\frac12\sum_jL_j^*L_j$, que es la forma GKSL. $\square$

**Teorema D.5 (libertad en la forma GKSL).** Sean $(H;L_1,\dots,L_r)$ y $(H';L'_1,\dots,L'_s)$ datos GKSL
($H,H'$ autoadjuntos) con generadores $\mathcal L$, $\mathcal L'$ como en el teorema D.4.

(a) $\mathcal L'=\mathcal L$ en cada uno de los casos siguientes:
(i) $L'_j=L_j$ y $H'=H+c1$, $c\in\mathbb R$;
(ii) $H'=H$ y, tras completar con ceros la familia más corta, $L'_j=\sum_ku_{jk}L_k$ con una matriz unitaria $u$;
(iii) $L'_j=L_j+a_j1$ con $a_j\in\mathbb C$, y $H'=H-\frac i2\sum_j(\bar a_jL_j-a_jL_j^*)$, que es autoadjunto.

(b) Si $\operatorname{Tr}H=\operatorname{Tr}H'=0$, todos los $L_j,L'_j$ tienen traza nula y $\mathcal L=\mathcal L'$, entonces $H=H'$ y,
tras completar con ceros, $L'_j=\sum_ku_{jk}L_k$ con una matriz unitaria $u$.

(c) Por tanto, dos datos GKSL tienen el mismo generador si y solo si uno se obtiene del otro mediante un número finito
de pasos (i)–(iii). En particular, el operador $H$ de la forma normalizada del teorema D.4 y su
parte disipativa $\Phi(\rho)=\sum_jL_j\rho L_j^*$ están determinados por $\mathcal L$.

*Demostración.* (a) (i) $[c1,\rho]=0$. (ii) Por el teorema 6.7, $\Phi'=\sum_jL'_j\cdot L_j'^*$ es igual a
$\Phi=\sum_jL_j\cdot L_j^*$; el dual no depende de los operadores de Kraus (sección 6.2), así que
$\sum_jL_j'^*L'_j=\Phi^*(1)=\sum_jL_j^*L_j$. (iii) Para un solo salto $L$ y $a\in\mathbb C$, desarrollando se obtiene
$$
(L+a)\rho(L+a)^*-\tfrac12\{(L+a)^*(L+a),\rho\}-\Big(L\rho L^*-\tfrac12\{L^*L,\rho\}\Big)=\tfrac12[B,\rho],
\qquad B=\bar aL-aL^* ,
$$
ya que los términos $\bar aL\rho+a\rho L^*$ del primer producto se combinan con
$-\frac12(aL^*\rho+a\rho L^*+\bar aL\rho+\bar a\rho L)$ y los términos con $|a|^2$ se cancelan. Como $B^*=-B$, la
matriz $\frac i2B$ es autoadjunta, así que $H'=H-\frac i2\sum_jB_j$ es autoadjunto, y
$-i[H'-H,\rho]=-\frac12\sum_j[B_j,\rho]$ compensa el cambio del disipador.

(b) Sean $\Phi=\sum_jL_j\cdot L_j^*$ y $G=-iH-\frac12\Phi^*(1)$, de modo que $\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*$.
Por el lema D.1(c),
$J(\mathcal L)=\sum_j|l_j\rangle\langle l_j|+|g\rangle\langle\Omega|+|\Omega\rangle\langle g|$ con $l_j=(1\otimes L_j)\Omega$ y
$g=(1\otimes G)\Omega$, y $l_j\perp\Omega$ porque $\operatorname{Tr}L_j=0$ (lema D.1(b)). Con $P^\perp$ la proyección
sobre $\Omega^\perp$, se tiene por tanto $P^\perp J(\mathcal L)P^\perp=\sum_j|l_j\rangle\langle l_j|=J(\Phi)$, y del mismo modo
$P^\perp J(\mathcal L')P^\perp=J(\Phi')$. Así, $J(\Phi)=J(\Phi')$, $\Phi=\Phi'$ (lema D.1(a)), y el teorema 6.7 da la
matriz unitaria $u$. Entonces $\Phi^*(1)=\Phi'^*(1)$, luego $[H-H',\rho]=0$ para todo $\rho$. Una matriz $X$ que conmuta con todas
las $E_{ij}$ es un múltiplo de $1$ (compárense las entradas de $XE_{ij}$ y de $E_{ij}X$), así que $H-H'=\lambda1$, y
$\lambda=0$ tomando trazas.

(c) Por (a), los pasos (i)–(iii) no cambian el generador, y el inverso de cada paso es un paso del
mismo tipo ((iii) con $-a_j$, ya que $\bar a_jL'_j-a_jL_j'^*=\bar a_jL_j-a_jL_j^*$). Cualquier dato se lleva a la
forma normalizada mediante (iii) con $a_j=-\operatorname{Tr}L_j/n$ seguido de (i) con $c=-\operatorname{Tr}H''/n$, siendo $H''$ el
hamiltoniano tras (iii); dos datos normalizados con el mismo generador están relacionados por (ii), por (b). $\square$

La demostración del teorema D.4 es un algoritmo: $\Phi$ a partir de la compresión de $J(\mathcal L)$ a $\Omega^\perp$, después
$G$ a partir de $\xi_1$, después $H$. Lo ejecutamos sobre el generador anterior, reconstruimos el generador con la biblioteca y
comparamos con la forma normalizada de los datos originales calculada mediante el teorema D.5(iii) y (i):

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

**Teorema D.6 (unicidad del estado estacionario).** Sea $\mathcal L$ el generador GKSL del teorema D.4
con datos $(H;L_1,\dots,L_r)$. Si $\mathcal L$ tiene un estado estacionario fiel $\sigma$ ($\sigma>0$) y las únicas
matrices que conmutan con $H$ y con todos los $L_j$, $L_j^*$ son los múltiplos de $1$, entonces $\ker\mathcal L=\mathbb C\sigma$:
$\sigma$ es el único estado estacionario. Además, $\ker\mathcal L^*=\mathbb C1$.

*Demostración.* El dual para el emparejamiento de la traza, $\operatorname{Tr}(\mathcal L(\rho)A)=\operatorname{Tr}(\rho\,\mathcal L^*(A))$, es
$\mathcal L^*(A)=i[H,A]+\sum_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})$: por ciclicidad,
$\operatorname{Tr}(-i[H,\rho]A)=\operatorname{Tr}(\rho\,i[H,A])$, $\operatorname{Tr}(L\rho L^*A)=\operatorname{Tr}(\rho L^*AL)$ y
$\operatorname{Tr}(\{L^*L,\rho\}A)=\operatorname{Tr}(\rho\{L^*L,A\})$ (es el ejercicio 1 del capítulo 8).
Un desarrollo directo, en el que los términos hamiltonianos se cancelan porque $A\mapsto i[H,A]$ es una derivación, da
para todo $X$
$$
\mathcal L^*(X^*X)-X^*\mathcal L^*(X)-\mathcal L^*(X^*)X=\sum_j[L_j,X]^*[L_j,X]:
$$
para un solo salto, de los seis términos de anticonmutador cuatro se cancelan por parejas y los otros dos,
$\frac12X^*L^*LX+\frac12X^*L^*LX$, se suman, así que el primer miembro es $L^*X^*XL-X^*L^*XL-L^*X^*LX+X^*L^*LX$, que es
$(XL-LX)^*(XL-LX)$. Sea $X\in\ker\mathcal L^*$. Entonces $\mathcal L^*(X^*)=\mathcal L^*(X)^*=0$, pues $\mathcal L^*$
preserva la hermiticidad, así que $\mathcal L^*(X^*X)=\sum_j[L_j,X]^*[L_j,X]=:Y\ge0$, y
$\operatorname{Tr}(\sigma Y)=\operatorname{Tr}(\mathcal L(\sigma)X^*X)=0$. Como $\sigma>0$,
$\operatorname{Tr}(\sigma Y)=\operatorname{Tr}(\sigma^{1/2}Y\sigma^{1/2})$ con $\sigma^{1/2}Y\sigma^{1/2}\ge0$ solo se anula si
$\sigma^{1/2}Y\sigma^{1/2}=0$, es decir, $Y=0$; tomando trazas, $\sum_j\operatorname{Tr}([L_j,X]^*[L_j,X])=0$ con términos
no negativos, luego $[L_j,X]=0$ para todo $j$. El mismo argumento para $X^*\in\ker\mathcal L^*$
da $[L_j,X^*]=0$, es decir, $[L_j^*,X]=0$. Entonces $X$ conmuta con $L_j^*L_j$, la parte disipativa de $\mathcal L^*(X)$
se anula, y $0=\mathcal L^*(X)=i[H,X]$. Por hipótesis, $X\in\mathbb C1$; recíprocamente, $\mathcal L^*(1)=0$. Así,
$\ker\mathcal L^*=\mathbb C1$. En las bases $(E_{ij})$ y $(E_{ji})$ de $M_n$, que son duales para el emparejamiento
$(\rho,A)\mapsto\operatorname{Tr}(\rho A)$, la matriz de $\mathcal L^*$ es la traspuesta de la de $\mathcal L$; así que ambas tienen el
mismo rango, y $\dim\ker\mathcal L=\dim\ker\mathcal L^*=1$. Un estado de $\mathbb C\sigma$ con traza $1$ es $\sigma$. $\square$

Este es el caso de dimensión finita de un criterio de Frigerio [*Commun. Math. Phys.* 63, 269 (1978)],
usado en la sección 8.5. Ambas hipótesis son necesarias. En el ejercicio 4 del capítulo 8, el estado $\tfrac12\,1$ es fiel
y estacionario, pero $H$ y $\sigma_z$ conmutan con todas las matrices diagonales, y los estados estacionarios generan el
espacio bidimensional de las matrices diagonales. Recíprocamente, sobre $M_3$ sean $H=0$, $L_1=|e_1\rangle\langle e_3|$,
$L_2=|e_2\rangle\langle e_3|$. La $*$-álgebra generada por $L_1,L_2$ contiene $L_1L_2^*=E_{12}$, $L_1^*L_1=E_{33}$ y, por tanto,
todas las unidades matriciales, así que solo los múltiplos de $1$ conmutan con $L_j,L_j^*$; pero todo $\rho=P\rho P$ con $P=E_{11}+E_{22}$
es estacionario, pues $L_j\rho=0=\rho L_j^*$ y $\{E_{33},\rho\}=0$, así que $\dim\ker\mathcal L\ge4$; por el teorema D.6, ningún estado
estacionario es fiel. Para el generador aleatorio anterior, el código comprueba las hipótesis y la conclusión
(la matriz de $\mathcal L^*$ sobre los vectores de $A^T$ ordenados por filas es $S^T$: $\operatorname{vec}(\mathcal L^*(A)^T)=S^T\operatorname{vec}(A^T)$), y el segundo ejemplo:

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

## D.4 Álgebras con sectores

Sea $\mathcal M=\bigoplus_{c=1}^CM_{k_c}\subseteq M_D$ el álgebra de las matrices diagonales por bloques, $D=\sum_ck_c$, con
proyecciones centrales minimales $P_c$, y sea $E(X)=\sum_cP_cXP_c$ la esperanza condicional de $M_D$ sobre
$\mathcal M$. Como en la demostración de la proposición 8.8, $E(B)=B$, $E(XB)=E(X)B$ y $E(BX)=BE(X)$ para $B\in\mathcal M$,
y $E(X^*)=E(X)^*$. Para una aplicación CP $\Phi:\mathcal M\to\mathcal M$, $\Phi^*$ es su dual `blunt` (proposición 6.8).

**Teorema D.7.** Sea $\mathcal L$ el generador de un semigrupo dinámico cuántico $(T_t)$ sobre $\mathcal M$. Existen
una aplicación CP $\Phi:\mathcal M\to\mathcal M$ y $H=H^*\in\mathcal M$ con
$$
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]\qquad(\rho\in\mathcal M).
$$
Por la proposición 6.4, $\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l$ con operadores de Kraus entre los sectores;
recíprocamente, todo $\mathcal L$ de este tipo genera un semigrupo dinámico cuántico sobre $\mathcal M$ (proposición 8.8(1)).
Así, los generadores de la sección 8.3 (donde $\Phi=\sum_jg_jL_j\cdot L_j^*+\sum_mg_m\Phi_m$) son exactamente los
generadores de los semigrupos dinámicos cuánticos sobre $\mathcal M$.

*Demostración.* *Paso 1: una aplicación CP sobre $M_D$.* Por la proposición 6.4,
$T_t(\rho)_d=\sum_{c,l}K^{dc}_l(t)\rho_cK^{dc}_l(t)^*$. Sea $N^{dc}_l(t)\in M_D$ la matriz con el bloque $K^{dc}_l(t)$ en la
posición $(d,c)$ y ceros en el resto. Para $X\in M_D$, $N^{dc}_lXN^{dc*}_l$ es el bloque
$K^{dc}_lX_{cc}K^{dc*}_l$ en la posición $(d,d)$, donde $X_{cc}$ es el $c$-ésimo bloque diagonal, es decir, el $c$-ésimo bloque
de $E(X)$. Por tanto, $\hat T_t=T_t\circ E=\sum_{d,c,l}N^{dc}_l\cdot N^{dc*}_l$ es CP sobre $M_D$ (teorema 6.2), con
$\hat T_0=E$. Por el lema 8.5, $(\hat T_t-E)/t=\big((T_t-\mathrm{id})/t\big)\circ E\to\mathcal L\circ E$.

*Paso 2: positividad fuera de los vectores de los sectores.* Sean $p_c=(1\otimes P_c)\Omega\in\mathbb C^D\otimes\mathbb C^D$ y
$W=\operatorname{span}\{p_c\}$. Por el lema D.1(c), $J(E)=\sum_c|p_c\rangle\langle p_c|$, así que para $v\perp W$ y $t>0$,
$0\le\frac1t\langle v,J(\hat T_t)v\rangle=\langle v,J((\hat T_t-E)/t)v\rangle\to\langle v,J(\mathcal L\circ E)v\rangle$.
Todo $T_t$ preserva la hermiticidad (por su forma de Kraus), luego también $\mathcal L$ (como en el teorema D.3), y también $E$
y $\mathcal L\circ E$.

*Paso 3: descomposición.* Por el lema D.1(b), $\langle p_c,p_d\rangle=\operatorname{Tr}(P_cP_d)=\delta_{cd}k_c$, así que los
$q_c=p_c/\sqrt{k_c}=(1\otimes Q_c)\Omega$ con $Q_c=P_c/\sqrt{k_c}$ son ortonormales. El lema D.2 para
$\Psi=\mathcal L\circ E$ da una aplicación CP $\hat\Phi=\sum_j\hat L_j\cdot\hat L_j^*$ sobre $M_D$ y $K_c\in M_D$ con
$$
\mathcal L(E(A))=\hat\Phi(A)+\sum_c\big(K_cAP_c+P_cAK_c^*\big)/\sqrt{k_c}\qquad(A\in M_D).
$$
Para $A\in\mathcal M$, $E(A)=A$ y $AP_c=P_cA$, así que la suma es $\tilde GA+A\tilde G^*$ con
$\tilde G=\sum_cK_cP_c/\sqrt{k_c}$.

*Paso 4: compresión.* Para $A\in\mathcal M$, $\mathcal L(A)\in\mathcal M$, así que
$\mathcal L(A)=E(\mathcal L(A))=E(\hat\Phi(A))+GA+AG^*$ con $G=E(\tilde G)\in\mathcal M$. La aplicación
$\Phi=E\circ\hat\Phi$ sobre $\mathcal M$ es
$\Phi(A)_d=\sum_{c,j}(P_d\hat L_jP_c)A_c(P_d\hat L_jP_c)^*$: es CP por la proposición 6.4, con los bloques de los $\hat L_j$ como
operadores de Kraus.

*Paso 5: preservación de la traza.* Como en el teorema D.4, $\operatorname{Tr}\mathcal L(A)=0$ para $A\in\mathcal M$, y
$\operatorname{Tr}\Phi(A)=\operatorname{Tr}(A\Phi^*(1))$, así que $\operatorname{Tr}(AX)=0$ para todo $A\in\mathcal M$, donde
$X=\Phi^*(1)+G+G^*\in\mathcal M$ es autoadjunto. Con $A=X$, $X=0$ (la traza `blunt` es fiel, teorema 3.1).
Con $H=\frac i2(G-G^*)\in\mathcal M$, $G=-iH-\frac12\Phi^*(1)$ y $G\rho+\rho G^*=-i[H,\rho]-\frac12\{\Phi^*(1),\rho\}$. $\square$

Los vectores $p_c$ sustituyen a $\Omega=\sum_cp_c$ porque $\hat T_0=E$, no $\mathrm{id}$: $J(\mathcal L\circ E)$ no tiene por qué ser
positiva sobre $\Omega^\perp$, y en el ejemplo de abajo no lo es.

**Corolario D.8 (imagen de Heisenberg).** Sea $(S_t)_{t\ge0}$ una familia continua de aplicaciones CP unitales sobre
$\mathcal M$ con $S_0=\mathrm{id}$ y $S_{s+t}=S_sS_t$. Su generador es $A\mapsto\Psi(A)+K^*A+AK$ con
$\Psi$ CP sobre $\mathcal M$, $K\in\mathcal M$ y $\Psi(1)+K+K^*=0$.

*Demostración.* Sean $T_t=S_t^*$ (duales `blunt`). Son CP (proposición 6.8) y TP (corolario 6.9(b)),
$T_0=\mathrm{id}$ y $T_{s+t}=(S_tS_s)^*=S_s^*S_t^*=T_sT_t$, porque
$\operatorname{Tr}(S_tS_s(x)y)=\operatorname{Tr}(xS_s^*S_t^*(y))$ y $S_tS_s=S_{s+t}$. Como $\Psi\mapsto\Psi^*$ es lineal sobre un
espacio de dimensión finita, $t\mapsto T_t$ es continua y $(T_t-\mathrm{id})/t\to\mathcal L_S^*$, siendo $\mathcal L_S$ el
generador de $(S_t)$. El teorema D.7 da $\mathcal L_S^*(\rho)=\Phi(\rho)+G\rho+\rho G^*$ con $G+G^*=-\Phi^*(1)$, y por
ciclicidad, $\operatorname{Tr}((G\rho+\rho G^*)A)=\operatorname{Tr}(\rho(G^*A+AG))$, así que por el corolario 6.9(c)
$\mathcal L_S(A)=\Phi^*(A)+G^*A+AG$. Tómense $\Psi=\Phi^*$ (CP por la proposición 6.8) y $K=G$. $\square$

Este es el caso de dimensión finita del teorema de Christensen y Evans citado tras el teorema 8.7. Comprobamos
el teorema D.7 sobre $M_1\oplus M_2\subseteq M_3$ para un generador de `lindblad_rhs` con un salto dentro de los
sectores y un salto aleatorio entre ellos, siguiendo los cinco pasos de la demostración:

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

## D.5 Corrección de errores en álgebras de operadores

Usamos el marco de la sección 10.3. $\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*$ es un canal sobre $M_D$, $\mathcal C=P\mathbb C^D$
un código, $B(\mathcal C)=\{X\in M_D:X=PXP\}$ los operadores sobre $\mathcal C$, $\mathfrak A\subseteq B(\mathcal C)$ una $*$-álgebra
con unidad $P$ (un álgebra de von Neumann sobre $\mathcal C$ por el teorema 2.2), y $\mathfrak A'$ su conmutante en
$B(\mathcal C)$. $\mathfrak A$ es **corregible** para $\Phi$ si existe un canal $\mathcal R$ sobre $M_D$ con
$P(\mathcal R\circ\Phi)^*(A)P=A$ para todo $A\in\mathfrak A$. Si $\mathcal R$ tiene los operadores de Kraus $R_m$, entonces
$\mathcal R\circ\Phi$ tiene los operadores de Kraus $R_mE_i$, y $\sum_mR_m^*R_m=1$.

**Teorema D.9 (Bény–Kempf–Kribs; teorema 10.7).** $\mathfrak A$ es corregible para $\Phi$ si y solo si
$PE_i^*E_jP\in\mathfrak A'$ para todos $i,j$. En ese caso, el canal $\mathcal R$ construido en la demostración corrige $\mathfrak A$.

*Demostración.* *Necesidad.* Supongamos que $\mathcal R$ corrige $\mathfrak A$ y pongamos $W_{mi}=R_mE_iP$, de modo que
$\sum_{m,i}W_{mi}^*XW_{mi}=P(\mathcal R\circ\Phi)^*(X)P$ para todo $X\in M_D$. Esto es igual a $X$ para $X\in\mathfrak A$, y a $P$ para
$X=1$, porque $(\mathcal R\circ\Phi)^*$ es unital (corolario 6.9(b)). Sea $p\in\mathfrak A$ una proyección; entonces $p\le P$
y $\sum W_{mi}^*(1-p)W_{mi}=P-p$. Para $x=px$,
$\sum_{m,i}\|(1-p)W_{mi}x\|^2=\langle x,(P-p)x\rangle=0$, y para $x=(P-p)x$,
$\sum_{m,i}\|pW_{mi}x\|^2=\langle x,px\rangle=0$. Por tanto, $(1-p)W_{mi}p=0$ y $pW_{mi}(P-p)=0$, y como
$W_{mi}=W_{mi}P$,
$$
W_{mi}p=pW_{mi}p=pW_{mi}P=pW_{mi}.
$$
Todo $A\in\mathfrak A$ autoadjunto es combinación lineal de proyecciones de $\mathfrak A$: por el teorema 1.3 sobre $\mathcal C$,
$A=\sum_k\lambda_kp_k$ con $\lambda_k$ distintos y proyecciones ortogonales $p_k$, $\sum_kp_k=P$, y
$p_k=\prod_{l\neq k}(A-\lambda_lP)/(\lambda_k-\lambda_l)\in\mathfrak A$ (aplíquense ambos miembros a los vectores propios de $\mathcal C$ y a $\mathcal C^\perp$; el producto vacío es $P$). Todo elemento
de $\mathfrak A$ es $A_1+iA_2$ con $A_1,A_2\in\mathfrak A$ autoadjuntos. Así, $W_{mi}A=AW_{mi}$ para todo $A\in\mathfrak A$ y,
tomando adjuntos, $W_{mi}^*A=AW_{mi}^*$. Como $\sum_mW_{mi}^*W_{mj}=PE_i^*(\sum_mR_m^*R_m)E_jP=PE_i^*E_jP$,
$$
A\,PE_i^*E_jP=\sum_mW_{mi}^*AW_{mj}=PE_i^*E_jP\,A\qquad(A\in\mathfrak A),
$$
es decir, $PE_i^*E_jP\in\mathfrak A'$.

*Suficiencia.* Por el teorema 2.6 para $\mathfrak A$ sobre $\mathcal C$, existe una base ortonormal $|c;a,b\rangle$
($c=1,\dots,C$, $a=1,\dots,k_c$, $b=1,\dots,m_c$) de $\mathcal C$ —las imágenes de los $e_a\otimes f_b$ del $c$-ésimo sumando por
la matriz unitaria $V$— en la que $\mathfrak A$ consta de los operadores $\sum_{c,a,a',b}(A_c)_{aa'}|c;a,b\rangle\langle c;a',b|$,
$A_c\in M_{k_c}$, y $\mathfrak A'$ de los operadores $\sum_{c,a,b,b'}(Y_c)_{bb'}|c;a,b\rangle\langle c;a,b'|$, $Y_c\in M_{m_c}$.
Por tanto, la hipótesis da matrices $\Lambda^c_{ij}\in M_{m_c}$ con
$$
\langle c;a,b|E_i^*E_j|d;a',b'\rangle=\delta_{cd}\,\delta_{aa'}\,(\Lambda^c_{ij})_{bb'} .
$$
*Paso 1.* Para cada $c$ y cada par $\alpha=(i,b)$, sean $F^c_\alpha:\mathbb C^{k_c}\to\mathbb C^D$,
$F^c_\alpha e_a=E_i|c;a,b\rangle$, y $\Lambda^c_{\alpha\beta}=(\Lambda^c_{ij})_{bb'}$ para $\beta=(j,b')$. La fórmula anterior dice que
$F^{c*}_\alpha F^d_\beta=\delta_{cd}\Lambda^c_{\alpha\beta}1_{k_c}$: una condición de Knill–Laflamme para los «errores» $F^c_\alpha$ sobre el
código $\mathbb C^{k_c}$, con, además, ortogonalidad entre los sectores.

*Paso 2.* Como en la demostración del teorema 10.5, $\Lambda^c$ es autoadjunta, ya que
$\Lambda^c_{\beta\alpha}1=F^{c*}_\beta F^c_\alpha=(F^{c*}_\alpha F^c_\beta)^*=\bar\Lambda^c_{\alpha\beta}1$, y positiva, ya que
$\sum\bar z_\alpha z_\beta\Lambda^c_{\alpha\beta}1=(\sum z_\alpha F^c_\alpha)^*(\sum z_\beta F^c_\beta)\ge0$. Escribamos
$\Lambda^c=U^c\operatorname{diag}(d^c_\kappa)U^{c*}$ con $U^c$ unitaria y $d^c_\kappa\ge0$, y
$G^c_\kappa=\sum_\alpha u^c_{\alpha\kappa}F^c_\alpha$. Entonces $G^{c*}_\kappa G^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}d^c_\kappa1$, así que
$G^c_\kappa=0$ si $d^c_\kappa=0$, y en otro caso $V^c_\kappa=G^c_\kappa/\sqrt{d^c_\kappa}$ es una isometría con
$V^{c*}_\kappa V^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}1$. Invirtiendo, $F^c_\alpha=\sum_\kappa\bar u^c_{\alpha\kappa}G^c_\kappa
=\sum_{\kappa:d^c_\kappa>0}\bar u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,V^c_\kappa$, de donde
$$
V^{c'*}_\kappa F^d_\beta=\delta_{c'd}\,\bar u^d_{\beta\kappa}\sqrt{d^d_\kappa}\,1,\qquad
F^{c*}_\alpha V^{c'}_\kappa=\delta_{cc'}\,u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,1 .
$$
*Paso 3: la recuperación.* Las $V^c_\kappa V^{c*}_\kappa$ ($d^c_\kappa>0$) son proyecciones mutuamente ortogonales; sea $Q$ su
suma. Sean $\iota_c:\mathbb C^{k_c}\to\mathbb C^D$, $\iota_ce_a=|c;a,1\rangle$, una isometría, y
$$
\mathcal R(\sigma)=\sum_{c}\sum_{\kappa:d^c_\kappa>0}\iota_cV^{c*}_\kappa\,\sigma\,V^c_\kappa\iota_c^*+(1-Q)\sigma(1-Q).
$$
Es un canal: $\sum_{c,\kappa}V^c_\kappa\iota_c^*\iota_cV^{c*}_\kappa+(1-Q)^2=Q+1-Q=1$. Mide el sector $c$ y
el síndrome $\kappa$, deshace el error sobre $\mathbb C^{k_c}$ y escribe el resultado en la copia $b=1$.

*Paso 4: corrección.* $\mathcal R^*(X)=\sum_{c,\kappa}V^c_\kappa\iota_c^*X\iota_cV^{c*}_\kappa+(1-Q)X(1-Q)$, y para
$A\in\mathfrak A$ con bloques $A_c$, $\iota_c^*A\iota_c=A_c$. Como cada $F^d_\beta$ es una combinación de los $V^d_\kappa$,
$(1-Q)F^d_\beta=0$. Con el paso 2, para $\alpha=(i,b)$, $\beta=(j,b')$,
$$
\langle c;a,b|E_i^*\mathcal R^*(A)E_j|d;a',b'\rangle=\sum_{c',\kappa}\langle e_a,F^{c*}_\alpha V^{c'}_\kappa A_{c'}V^{c'*}_\kappa F^d_\beta e_{a'}\rangle
=\delta_{cd}\sum_\kappa u^c_{\alpha\kappa}d^c_\kappa\bar u^c_{\beta\kappa}(A_c)_{aa'}=\delta_{cd}(\Lambda^c_{ij})_{bb'}(A_c)_{aa'} .
$$
Tomando $j=i$ y sumando sobre $i$: por la fórmula de la hipótesis con $\sum_iE_i^*E_i=1$, $\sum_i\Lambda^c_{ii}=1_{m_c}$, así que
$$
\langle c;a,b|(\mathcal R\circ\Phi)^*(A)|d;a',b'\rangle=\delta_{cd}\,\delta_{bb'}(A_c)_{aa'}=\langle c;a,b|A|d;a',b'\rangle .
$$
Tanto $P(\mathcal R\circ\Phi)^*(A)P$ como $A$ están en $B(\mathcal C)$ y tienen la misma matriz en la base $|c;a,b\rangle$, así que son
iguales. $\square$

**Corolario D.10.** (a) Para $\mathfrak A=B(\mathcal C)$, $\mathfrak A'=\mathbb CP$, y el teorema D.9 es el teorema de
Knill–Laflamme 10.5 ($C=1$, $m_1=1$; el paso 3 es la recuperación del teorema 10.5). (b) Para $\mathfrak A=\operatorname{span}\{P_1,\dots,P_C\}$ con
proyecciones ortogonales $P_c\ne0$, $\sum_cP_c=P$, es la equivalencia de 2 y 3 en la proposición 10.6. (c) (*Subsistemas
sin ruido.*) Con la notación de la proposición 10.1(a), sean $\mathcal C=V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})$ el $c$-ésimo
sumando para el álgebra de interacción $\mathcal A$ de $\Phi$, y $\mathfrak A=V(1_{k_c}\otimes M_{m_c})V^*$ sobre $\mathcal C$, un
álgebra de la forma $M_{m}\otimes1_{k}$ tras intercambiar los dos factores. Entonces $\mathfrak A$ es corregida por
$\mathcal R=\mathrm{id}$, y $PE_i^*E_jP=V(E_{i,c}^*E_{j,c}\otimes1_{m_c})V^*\in\mathfrak A'$ ($P=VV^*$).

*Demostración.* (a) $\mathfrak A'=\mathbb CP$ es el caso $C=1$, $m_1=1$ de la descripción de $\mathfrak A'$ de arriba, y
$PE_i^*E_jP\in\mathbb CP$ es la condición de Knill–Laflamme. La corregibilidad de $B(\mathcal C)$ es la corrección en el sentido de la
sección 10.2: si $\mathcal R(\Phi(\rho))=\rho$ para todos los estados $\rho$ con soporte en $\mathcal C$, entonces
$\operatorname{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)=\operatorname{Tr}(\rho A)$ para todos ellos, y estos $\rho$ generan $B(\mathcal C)$; el
recíproco se muestra en la sección 10.3. Con $C=1$, $m_1=1$, los operadores de Kraus $\iota_1V^{1*}_\kappa$ del paso 3 son los
del teorema 10.5, porque $\iota_1\iota_1^*=P$. (b) $T\in B(\mathcal C)$ conmuta con todas las $P_c$ si y solo si
$P_cTP_d=0$ para $c\ne d$, que es la condición 2 de la proposición 10.6. (c) Todo $A\in\mathfrak A$, extendido por $0$ a los
demás sumandos, está en $\mathcal A'=V(\bigoplus_c1_{k_c}\otimes M_{m_c})V^*$ (teorema 2.6), así que $\Phi^*(A)=A$ por la proposición
10.1(a) y $PAP=A$. La fórmula para $PE_i^*E_jP$ se sigue de $V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}$, y
$E_{i,c}^*E_{j,c}\otimes1_{m_c}$ conmuta con $1_{k_c}\otimes M_{m_c}$. $\square$

Por el teorema 2.6, toda $\mathfrak A$ es $*$-isomorfa a algún $\bigoplus_cM_{k_c}$. Cuando $C\ge2$ y algún $k_c\ge2$, una
$\mathfrak A$ corregible protege la etiqueta del sector $c$ (información clásica) junto con la información cuántica de
cada $M_{k_c}$; tales códigos se llaman códigos **híbridos**. Como ejemplo tomamos $\mathfrak A=(M_2\otimes1_2)\oplus M_1$ sobre un
código aleatorio de dimensión cinco en $\mathbb C^8$ y un ruido que actúa sobre los factores gauge $\mathbb C^{m_c}$ mediante
$\sum_iB_{c,i}^*B_{c,i}=1_{m_c}$, lleva los sectores a subespacios ortogonales, revuelve $\mathcal C^\perp$ y después se mezcla mediante
una matriz unitaria aleatoria (teorema 6.7). La condición de Knill–Laflamme falla, la condición del teorema D.9 se cumple, y la
recuperación de la demostración restablece todo $A\in\mathfrak A$, pero no un operador genérico sobre $\mathcal C$:

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

Los resultados de este apéndice se deben a V. Gorini, A. Kossakowski y E. C. G. Sudarshan
[*J. Math. Phys.* 17, 821 (1976)] y G. Lindblad [*Commun. Math. Phys.* 48, 119 (1976)] (teoremas D.3–D.5),
E. Christensen y D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)] (corolario D.8 en álgebras de von Neumann
generales), A. Frigerio (teorema D.6), y C. Bény, A. Kempf y D. W. Kribs [*Phys. Rev. Lett.* 98, 100502
(2007); *Phys. Rev. A* 76, 042303 (2007)] (teorema D.9); las demostraciones anteriores son completas para los enunciados
de dimensión finita que se usan en el curso.
