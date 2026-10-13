[English](../03_traces_and_states.md) · **Español** · [Français](../fr/03_traces_and_states.md) · [Deutsch](../de/03_traces_and_states.md) · [中文](../zh/03_traces_and_states.md) · [日本語](../ja/03_traces_and_states.md) · [Русский](../ru/03_traces_and_states.md)

# 3. Trazas, estados y entropía

Sobre $M_n(\mathbb C)$ hay esencialmente una sola traza. Sobre $\bigoplus_cM_{k_c}$ hay una
familia de $C$ parámetros, y la elección importa: las densidades, las entropías y los
duales de las aplicaciones dependen de ella. Este capítulo clasifica las trazas, representa
los estados mediante densidades, demuestra las desigualdades básicas de la entropía y
termina con una desigualdad de trazas verosímil que es falsa.

**Aprenderá** las tres trazas de la biblioteca y cómo pasar de una a otra, la parte clásica
de un estado, la descomposición de la entropía de von Neumann por sectores, la desigualdad
de Klein y sus consecuencias, y el contraste de Michelson.

En todo el capítulo, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ con todos los $k_c\ge1$,
$D=\sum_ck_c$, $P_c$ son las proyecciones centrales minimales y $A=\bigoplus_cA_c$.

## 3.1 Trazas

Una **traza** sobre $\mathcal M$ es una forma lineal $\tau$ con $\tau(AB)=\tau(BA)$ para todos
los $A,B$. Es **positiva** si $\tau(A^*A)\ge0$ para todo $A$, y **fiel** si además
$\tau(A^*A)=0$ implica $A=0$.

**Teorema 3.1.** Toda traza sobre $\mathcal M$ es de la forma
$\operatorname{tr}_w(A)=\sum_cw_c\operatorname{Tr}A_c$ con $w\in\mathbb C^C$. Es positiva si y
solo si todos los $w_c\ge0$, y fiel si y solo si todos los $w_c>0$.

*Demostración.* Sean $E^c_{ij}$ las unidades matriciales del $c$-ésimo bloque. Para una
traza $\tau$,
$\tau(E^c_{ij})=\tau(E^c_{i1}E^c_{1j})=\tau(E^c_{1j}E^c_{i1})=\delta_{ij}\tau(E^c_{11})$.
Sea $w_c=\tau(E^c_{11})$; por linealidad, $\tau(A)=\sum_c\sum_i(A_c)_{ii}w_c$. Recíprocamente,
toda $\operatorname{tr}_w$ es una traza. Positividad: $E^c_{11}=(E^c_{11})^*E^c_{11}$, luego
$w_c\ge0$ es necesario, y $\operatorname{tr}_w(A^*A)=\sum_cw_c\|A_c\|_2^2$ prueba la suficiencia
y el criterio de fidelidad. $\square$

La biblioteca implementa tres trazas fieles:

| nombre | pesos $w_c$ | valor en $1$ | significado |
|---|---|---|---|
| `Tr_blunt` | $1$ | $D$ | traza de $\bigoplus_c\mathbb C^{k_c}$ (multiplicidades $m_c=1$) |
| `Tr_norm` | $1/k_c$ | $C$ | cada bloque normalizado |
| `tau_vN` | $1/(Ck_c)$ | $1$ | el estado tracial con pesos iguales $1/C$ en los sectores |

Por el capítulo 2, la traza de un espacio sobre el que $\mathcal M$ actúa con
multiplicidades $m_c$ es $\operatorname{tr}_w$ con $w_c=m_c$. Una traza fiel es, por tanto, una
elección de "cuánto cuenta cada sector".

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, DensityMatrix, states, channels

torch.manual_seed(0)
alg = TypeIAlgebra([2, 3, 5], [2, 3, 5], precision='double')    # M_2 (+) M_3 (+) M_5, D = 10, C = 3
U = alg.random_unitary_operator(4)
A = alg.operator_from_eigenvalues(lambda d: torch.randn(4, d), batch_size=4, force_self_adjoint=True)
B = U @ A                                                       # a batch of four non-normal operators
one = alg.identity()

for f, value in ((lambda X: X.Tr_blunt(), 10), (lambda X: X.Tr_norm(), 3), (lambda X: X.tau_vN(), 1)):
    assert abs(f(one).item() - value) < 1e-12                   # normalisation
    assert (f(A @ B) - f(B @ A)).abs().max() < 1e-10            # trace property
assert (A.Tr_norm() - 3 * A.tau_vN()).abs().max() < 1e-12       # Tr_norm = C tau_vN
```

## 3.2 Estados y densidades

Un **estado** es una forma lineal $\omega$ sobre $\mathcal M$ que es positiva
($\omega(A^*A)\ge0$) y normalizada ($\omega(1)=1$).

**Teorema 3.2.** Sea $\operatorname{tr}_w$ una traza fiel. Para toda forma lineal $\omega$ existe
un único $\rho_w\in\mathcal M$ con $\omega(A)=\operatorname{tr}_w(\rho_wA)$ para todo $A$. La forma
$\omega$ es positiva si y solo si $\rho_w\ge0$, y es un estado si y solo si, además,
$\operatorname{tr}_w(\rho_w)=1$. Para dos trazas fieles, las densidades están relacionadas bloque a
bloque por
$$
w_c\,(\rho_w)_c=w'_c\,(\rho_{w'})_c .
$$

*Demostración.* $\langle X,Y\rangle_w=\operatorname{tr}_w(X^*Y)$ es un producto escalar sobre
$\mathcal M$ (teorema 3.1). Por el teorema de representación de Riesz en el espacio de
dimensión finita $\mathcal M$, existe un único $R$ con
$\omega(A)=\langle R,A\rangle_w=\operatorname{tr}_w(R^*A)$; sea $\rho_w=R^*$. Si $\rho_w\ge0$, entonces
$\omega(A^*A)=\sum_cw_c\operatorname{Tr}(A_c\,(\rho_w)_c\,A_c^*)\ge0$, ya que
$A_c(\rho_w)_cA_c^*\ge0$ (corolario 1.7 con $T=(\rho_w)_c^{1/2}A_c^*$). Recíprocamente, si
$\omega$ es positiva, tomemos $x\in\mathbb C^{k_c}$ y $A$ el operador de rango uno
$|x\rangle\langle x|$ en el bloque $c$; entonces $A=A^*A/\|x\|^2$ y
$0\le\omega(A)=w_c\langle x,(\rho_w)_cx\rangle$, luego $(\rho_w)_c\ge0$. La normalización es
$\omega(1)=\operatorname{tr}_w(\rho_w)$. Por último, tanto $\sum_cw_c\operatorname{Tr}((\rho_w)_cA_c)$ como
$\sum_cw'_c\operatorname{Tr}((\rho_{w'})_cA_c)$ son iguales a $\omega(A)$ para todo $A$, y por
unicidad $w_c(\rho_w)_c=w'_c(\rho_{w'})_c$. $\square$

La biblioteca almacena un estado mediante su densidad $\rho=\rho_{\mathrm{blunt}}$ respecto de
`Tr_blunt` (clase `DensityMatrix`); `density(trace)` y `DensityMatrix.from_density`
realizan la conversión. El valor del estado no depende del convenio:

```python
rho = states.random_density_matrix(alg, batch_size=4)          # Hilbert-Schmidt random states
for trace, tr in (('blunt', lambda X: X.Tr_blunt()), ('norm', lambda X: X.Tr_norm()),
                  ('tau_vN', lambda X: X.tau_vN())):
    d = alg.operator(rho.density(trace))
    assert (tr(d @ A) - (rho @ A).Tr_blunt()).abs().max() < 1e-10       # omega(A) for all traces
    assert (tr(d) - 1).abs().max() < 1e-10                               # tr_w(rho_w) = 1
    back = DensityMatrix.from_density(alg, rho.density(trace), trace=trace)
    assert (back.matrix - rho.matrix).abs().max() < 1e-12
```

> **Significado físico.** $\omega(A)$ es el valor esperado del observable $A$. La densidad
> $\rho$ de la traza `blunt` es la matriz densidad usual de la mecánica cuántica sobre
> $\bigoplus_c\mathbb C^{k_c}$. Las densidades respecto de otras trazas aparecen cuando la
> misma álgebra actúa con multiplicidades (capítulo 2), por ejemplo sobre la parte simétrica
> y la parte no simétrica de un sistema de espines.

## 3.3 La parte clásica de un estado

La restricción de $\omega$ al centro es una distribución de probabilidad sobre los sectores,
$p_c=\omega(P_c)=\operatorname{Tr}\rho_c$. Si $p_c>0$, el **estado condicional** en el sector $c$
tiene densidad $\sigma_c=\rho_c/p_c$ sobre $\mathbb C^{k_c}$, y
$$
\rho=\bigoplus_cp_c\,\sigma_c
$$
(si $p_c=0$, entonces $\rho_c=0$, ya que $\rho_c\ge0$ y $\operatorname{Tr}\rho_c=0$, y el término se
omite).
Un estado sobre $\mathcal M$ es una mezcla clásica de estados cuánticos en los sectores. Medir
"en qué sector" no lo perturba: el canal de medición $\rho\mapsto\sum_cP_c\rho P_c$ es la
identidad sobre $\mathcal M$, porque todo elemento de $\mathcal M$ es diagonal por bloques.
Este es el contenido operativo de una regla de superselección.

```python
p = rho.sector_probabilities()                                   # shape (batch, C)
assert torch.allclose(p.sum(-1), torch.ones(4, dtype=p.dtype))
P = [alg.central([1.0 if c == d else 0.0 for d in range(alg.C)]) for c in range(alg.C)]
measure_sector = channels.measurement_channel(P)
assert (measure_sector(rho).matrix - rho.matrix).abs().max() < 1e-12   # no disturbance
```

## 3.4 Entropía

La **entropía de von Neumann** de un estado con densidad (`blunt`) $\rho$ es
$S(\rho)=-\operatorname{Tr}\rho\log\rho=-\sum_i\lambda_i\log\lambda_i$, sumando sobre los valores
propios de $\rho$, con $0\log0=0$.

**Teorema 3.3 (entropía de una mezcla de sectores).** Con $p_c$ y $\sigma_c$ como antes
(omitiendo los sectores con $p_c=0$),
$$
S(\rho)=H(p)+\sum_cp_c\,S(\sigma_c),\qquad H(p)=-\sum_cp_c\log p_c .
$$

*Demostración.* Los valores propios no nulos de $\rho$ son los números $p_c\mu$, donde $\mu$
recorre los valores propios no nulos de las $\sigma_c$ (sectores con $p_c>0$); los valores
propios nulos contribuyen con $0$ a ambos miembros. Por tanto
$-\sum_c\sum_{\mu>0}p_c\mu\log(p_c\mu)=-\sum_cp_c\log p_c\sum_{\mu>0}\mu-\sum_cp_c\sum_{\mu>0}\mu\log\mu$,
y $\sum_{\mu>0}\mu=\operatorname{Tr}\sigma_c=1$. $\square$

> **Significado físico.** La fórmula divide la incertidumbre sobre el sistema en dos partes:
> $H(p)$, la incertidumbre clásica sobre la carga (en qué sector está el sistema), y
> $\sum_cp_cS(\sigma_c)$, la incertidumbre cuántica media dentro de los sectores. Bajo una regla de
> superselección, la primera parte es la entropía de una variable aleatoria clásica: la carga.

La entropía relativa a otra traza es $S_w(\rho)=-\operatorname{tr}_w(\rho_w\log\rho_w)$. Por el
teorema 3.2, $(\rho_w)_c=\rho_c/w_c$, luego
$$
S_w(\rho)=-\sum_cw_c\operatorname{Tr}\Big(\frac{\rho_c}{w_c}\log\frac{\rho_c}{w_c}\Big)
=S(\rho)+\sum_cp_c\log w_c .
$$
La entropía es una propiedad del estado *y* de la traza: el mismo estado tiene entropías
distintas respecto de `Tr_blunt` y de `tau_vN`.

```python
S = rho.von_neumann_entropy()
sigma = [rho.matrix[:, c, :k, :k] / p[:, c, None, None] for c, k in enumerate(alg.k_factors)]
S_sec = [-(w * torch.log(w)).sum(-1) for w in (torch.linalg.eigvalsh(s) for s in sigma)]
H = -(p * torch.log(p)).sum(-1)
assert torch.allclose(S, H + sum(p[:, c] * S_sec[c] for c in range(alg.C)))    # Theorem 3.3

w_tau = torch.tensor([1 / (alg.C * k) for k in alg.k_factors], dtype=torch.float64)
d_tau = alg.operator(rho.density('tau_vN'), is_positive=True)
S_tau = -(d_tau @ d_tau.log()).tau_vN().real
assert torch.allclose(S_tau, S + (p * torch.log(w_tau)).sum(-1))
```

## 3.5 Entropía relativa y desigualdad de Klein

Para estados con densidades $\rho,\sigma$, $\sigma>0$, la **entropía relativa** es
$S(\rho\|\sigma)=\operatorname{Tr}\rho\log\rho-\operatorname{Tr}\rho\log\sigma$, donde
$\operatorname{Tr}\rho\log\rho:=\sum_i\lambda_i\log\lambda_i$, sumando sobre los valores propios
de $\rho$, con $0\log0=0$ (y análogamente en cada bloque). No depende de la traza: por el teorema 3.2,
$(\rho_w)_c=\rho_c/w_c$ y $(\sigma_w)_c=\sigma_c/w_c$, así que los valores propios de $(\rho_w)_c$ son los
de $\rho_c$ divididos por $w_c$, y tanto $\operatorname{tr}_w(\rho_w\log\rho_w)$ como
$\operatorname{tr}_w(\rho_w\log\sigma_w)$ difieren de $\operatorname{Tr}\rho\log\rho$ y de
$\operatorname{Tr}\rho\log\sigma$ en la misma cantidad $-\sum_cp_c\log w_c$ (el logaritmo de $\rho$ solo
se necesita sobre su soporte, por la convención $0\log0=0$). La diferencia no cambia.

**Lema 3.4 (desigualdad de la tangente).** Sea $f:I\to\mathbb R$ convexa en un intervalo $I$ y
derivable en $b\in I$. Entonces $f(a)\ge f(b)+(a-b)f'(b)$ para todo $a\in I$, con
desigualdad estricta para $a\neq b$ si $f$ es estrictamente convexa.

*Demostración.* Para $t\in(0,1]$, la convexidad da $f(b+t(a-b))\le(1-t)f(b)+tf(a)$, es decir,
$\big(f(b+t(a-b))-f(b)\big)/t\le f(a)-f(b)$; haciendo $t\to0$ se obtiene la desigualdad. Si
$f$ es estrictamente convexa y $a\neq b$, sea $m=(a+b)/2$; entonces
$f(a)-f(b)>2\big(f(m)-f(b)\big)\ge2(m-b)f'(b)=(a-b)f'(b)$, el primer paso por la convexidad
estricta en $m$ y el segundo por la desigualdad que acabamos de demostrar, aplicada a $m$. $\square$

**Teorema 3.5 (desigualdad de Klein).** Sea $f:I\to\mathbb R$ convexa en un intervalo $I$, y
sean $X,Y$ autoadjuntos con espectros contenidos en $I$, con $f$ derivable en cada valor
propio de $Y$. Entonces $\operatorname{Tr}[f(X)-f(Y)-(X-Y)f'(Y)]\ge0$. Si $f$ es estrictamente
convexa, la igualdad solo se da para $X=Y$.

*Demostración.* Sean $X=\sum_ia_i|u_i\rangle\langle u_i|$ e $Y=\sum_jb_j|v_j\rangle\langle v_j|$ con
bases propias ortonormales, y $c_{ij}=|\langle u_i,v_j\rangle|^2$, de modo que
$\sum_ic_{ij}=\sum_jc_{ij}=1$. Entonces $\operatorname{Tr}f(X)=\sum_{ij}c_{ij}f(a_i)$,
$\operatorname{Tr}f(Y)=\sum_{ij}c_{ij}f(b_j)$ y
$\operatorname{Tr}Xf'(Y)=\sum_{ij}c_{ij}a_if'(b_j)$, $\operatorname{Tr}Yf'(Y)=\sum_{ij}c_{ij}b_jf'(b_j)$.
Por tanto, el primer miembro es igual a $\sum_{ij}c_{ij}[f(a_i)-f(b_j)-(a_i-b_j)f'(b_j)]$,
una suma de términos no negativos por el lema 3.4. Si $f$ es
estrictamente convexa y la suma es cero, entonces $a_i=b_j$ siempre que $c_{ij}>0$, luego
$Xv_j=\sum_ia_i\langle u_i,v_j\rangle u_i=b_jv_j$ para todo $j$, y $X=Y$. $\square$

**Corolario 3.6.** $S(\rho\|\sigma)\ge0$, con igualdad si y solo si $\rho=\sigma$. En
particular, $S(\rho)\le\log D$, con igualdad si y solo si $\rho=1/D$.

*Demostración.* Apliquemos el teorema 3.5 a $f(t)=t\log t$, $f(0)=0$, que es estrictamente
convexa en $I=[0,\infty)$ ($f''(t)=1/t>0$ en $(0,\infty)$, y $f$ es continua en $0$) y
derivable en $(0,\infty)\supseteq\operatorname{spec}\sigma$ con $f'(t)=\log t+1$, con
$X=\rho$, $Y=\sigma$: el primer miembro es
$\operatorname{Tr}\rho\log\rho-\operatorname{Tr}\sigma\log\sigma-\operatorname{Tr}(\rho-\sigma)(\log\sigma+1)
=S(\rho\|\sigma)$, ya que $\operatorname{Tr}\rho=\operatorname{Tr}\sigma=1$. Para la segunda
afirmación, tomemos $\sigma=1/D$: $0\le S(\rho\|1/D)=-S(\rho)+\log D$. $\square$

**Teorema 3.7 (regla de la cadena; monotonía bajo la esperanza condicional sobre el centro).** Sean $\rho,\sigma$ estados con
$\sigma>0$, $p_c=\operatorname{Tr}\rho_c$, $q_c=\operatorname{Tr}\sigma_c>0$, $\hat\sigma_c=\sigma_c/q_c$ y,
para $p_c>0$, $\hat\rho_c=\rho_c/p_c$. Entonces, con sumas sobre los $c$ con $p_c>0$ (los
demás sectores contribuyen con $0$ por $0\log0=0$),
$$
S(\rho\|\sigma)=\sum_cp_c\log\frac{p_c}{q_c}+\sum_cp_c\,S(\hat\rho_c\|\hat\sigma_c)
\ \ge\ \sum_cp_c\log\frac{p_c}{q_c}=S(E\rho\|E\sigma),
$$
donde $E$ es la esperanza condicional sobre el centro (proposición 2.8).

*Demostración.* Para $p_c>0$, $\operatorname{Tr}\rho_c\log\rho_c=p_c\log p_c+p_c\operatorname{Tr}\hat\rho_c\log\hat\rho_c$
(valores propios $p_c\mu$, como en el teorema 3.3) y $\log\sigma_c=(\log q_c)1+\log\hat\sigma_c$,
luego $\operatorname{Tr}\rho_c\log\sigma_c=p_c\log q_c+p_c\operatorname{Tr}\hat\rho_c\log\hat\sigma_c$ y
$\operatorname{Tr}\rho_c\log\rho_c-\operatorname{Tr}\rho_c\log\sigma_c=p_c\log(p_c/q_c)+p_cS(\hat\rho_c\|\hat\sigma_c)$.
Si $p_c=0$, entonces $\rho_c=0$ y el sector contribuye con $0$.
Sumamos sobre $c$; la segunda suma es $\ge0$ por el corolario 3.6 aplicado en cada bloque.
Por último, $E\rho=\bigoplus_c(p_c/k_c)1_{k_c}$, $E\sigma=\bigoplus_c(q_c/k_c)1_{k_c}$, y
$\operatorname{Tr}E\rho(\log E\rho-\log E\sigma)=\sum_ck_c\frac{p_c}{k_c}\log\frac{p_c}{q_c}$. $\square$

Así pues, olvidar la información cuántica dentro de los sectores solo puede disminuir la
distinguibilidad de dos estados; lo que queda es la divergencia clásica de Kullback–Leibler
entre las distribuciones de los sectores.

```python
sigma_state = states.random_density_matrix(alg, batch_size=4)
q = sigma_state.sector_probabilities()
D_rel = states.relative_entropy(rho, sigma_state)
kl = (p * torch.log(p / q)).sum(-1)
E = channels.center_expectation(alg)
assert (D_rel >= kl - 1e-12).all()                                       # Theorem 3.7
assert torch.allclose(states.relative_entropy(E(rho), E(sigma_state)), kl)
assert (rho.von_neumann_entropy() <= math.log(10) + 1e-12).all()         # Corollary 3.6
assert abs(states.maximally_mixed_state(alg).von_neumann_entropy().item() - math.log(10)) < 1e-12
```

## 3.6 El contraste de Michelson

Para $A\ge0$, $A\neq0$, el **contraste de Michelson** es
$\Delta(A)=(\lambda_{\max}-\lambda_{\min})/(\lambda_{\max}+\lambda_{\min})\in[0,1]$, donde los
valores propios extremos se toman sobre todos los sectores (y $\Delta(0)=0$). Mide cuánto
dista $A$ de un múltiplo de la identidad: $\Delta(A)=0$ si y solo si $A=\lambda1$.

**Proposición 3.8.** Si $\Phi:\mathcal M\to\mathcal M$ es lineal, positiva y unital
($\Phi(1)=1$), entonces $\Delta(\Phi(A))\le\Delta(A)$ para todo $A\ge0$.

*Demostración.* $\lambda_{\min}1\le A\le\lambda_{\max}1$, y la positividad y la unitalidad dan
$\lambda_{\min}1\le\Phi(A)\le\lambda_{\max}1$. Por tanto
$\lambda_{\min}(\Phi(A))\ge\lambda_{\min}$ y $\lambda_{\max}(\Phi(A))\le\lambda_{\max}$. Si
$\Phi(A)=0$, la afirmación es trivial. En otro caso, con $r=\lambda_{\min}/\lambda_{\max}$,
$\Delta=(1-r)/(1+r)$ es decreciente en $r$, y $r(\Phi(A))\ge r(A)$. $\square$

```python
A_pos = alg.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(4, d), batch_size=4, force_positive=True)
Phi = channels.random_mixed_unitary_channel(alg, 3, batch_size=4)       # unital
assert (Phi(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
assert (E(A_pos).michelson_contrast <= A_pos.michelson_contrast + 1e-12).all()
```

## 3.7 Una desigualdad verosímil que es falsa

Para $X,Y$ positivos en $M_k$ y una unitaria $U$, ¿se cumple $|\operatorname{Tr}(XUY)|\le\operatorname{Tr}(XY)$?
Es cierto para $U=1$ y parece un enunciado del tipo de Cauchy–Schwarz. Es falso siempre que
$X$ e $Y$ no conmutan.

**Proposición 3.9.** Sean $X,Y\ge0$ en $M_k$. Entonces $\max_U|\operatorname{Tr}(XUY)|=\|YX\|_1\ge\operatorname{Tr}(XY)$, donde el máximo se
toma sobre las unitarias $U\in M_k$, con igualdad si y solo si $XY=YX$.

*Demostración.* $\operatorname{Tr}(XUY)=\operatorname{Tr}(YX\,U)$, luego el máximo es $\|YX\|_1$ por el
teorema 1.9, y $\|YX\|_1\ge|\operatorname{Tr}(YX)|=\operatorname{Tr}(XY)$, ya que
$\operatorname{Tr}(XY)=\operatorname{Tr}(X^{1/2}YX^{1/2})\ge0$, pues
$X^{1/2}YX^{1/2}=(Y^{1/2}X^{1/2})^*(Y^{1/2}X^{1/2})\ge0$ (corolario 1.7). Si $XY=YX$, entonces
$X^{1/2}$ conmuta con $Y$ (teorema 1.5(d)), luego $YX=X^{1/2}YX^{1/2}\ge0$ y
$\|YX\|_1=\operatorname{Tr}(YX)$ (para $P\ge0$, $\|P\|_1=\operatorname{Tr}P$, sección 1.4). Recíprocamente,
supongamos que $\|YX\|_1=\operatorname{Tr}(YX)$. Si $|YX|=0$, entonces $YX=0=(YX)^*=XY$. En otro
caso, escribamos $YX=W|YX|$. Por el teorema 1.9, $\operatorname{Tr}(YX\cdot1)=\|YX\|_1$ significa
que se da la igualdad en el paso de Cauchy–Schwarz para $U=1$:
$|YX|^{1/2}$ y $W|YX|^{1/2}$ son proporcionales con un factor no negativo, que, comparando
normas, es $1$; así $W|YX|^{1/2}=|YX|^{1/2}$ e $YX=W|YX|=|YX|\ge0$. Una matriz positiva es
autoadjunta, luego $YX=(YX)^*=XY$. $\square$

La biblioteca encuentra el maximizador mediante la descomposición polar. El artículo que
acompaña a la biblioteca muestra numéricamente que el muestreo de Monte Carlo de $U$ se
queda por debajo de él, y cada vez más a medida que crece $k$ —uno de sus temas
centrales—; el código siguiente lo muestra para un par concreto. La razón es que $U=1$ ya está cerca del
óptimo, mientras que para una $U$ aleatoria de Haar en dimensión alta $|\operatorname{Tr}(XUY)|$ es
típicamente mucho menor que $\operatorname{Tr}(XY)$ (una unitaria aleatoria desordena las fases), de
modo que todas las muestras quedan muy por debajo del supremo exacto.

```python
Mk = TypeIAlgebra([6], [6], precision='double')
X = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
Y = Mk.operator_from_eigenvalues(lambda d: 0.2 + torch.rand(1, d), force_positive=True)
YX = Y @ X
W = YX @ YX.abs().inverse()                                          # polar part of YX
gap = (X @ W.adjoint() @ Y).trace.abs() - (X @ Y).trace.real
assert abs(gap.item() - (YX.trace_norm() - (X @ Y).trace.real).item()) < 1e-10
assert gap.item() > 0                                                # the inequality fails
samples = Mk.random_unitary_operator(20_000)
best = (X @ samples @ Y).trace.abs().max() - (X @ Y).trace.real
print(f"exact sup {gap.item():.4f}, best of 20000 Haar samples {best.item():.4f}")
assert best.item() < gap.item()
```

## 3.8 Con qué comparar

QuTiP (`entropy_vn`, `entropy_relative`) y los paquetes de información cuántica (por
ejemplo, `qutip`, `toqito`) calculan entropías y entropías relativas de matrices densidad
sobre un único espacio de Hilbert; sobre $\bigoplus_c\mathbb C^{k_c}$ coinciden con los
convenios de la traza `blunt` que se usan aquí. Lo que no proporcionan es la familia de
trazas del teorema 3.1 ni las conversiones del teorema 3.2, que se vuelven necesarias en
cuanto un álgebra actúa con multiplicidades o los sectores deben ponderarse de forma
distinta.

## Ejercicios

1. Demuestre que la densidad del estado tracial $\tau_{\mathrm{vN}}$ (pesos iguales en los
   sectores) respecto de $\tau_{\mathrm{vN}}$ es $1$, y calcule su densidad `blunt` y su
   entropía $S$.
2. Sea $\operatorname{tr}_w$ fiel. Demuestre que $S_w(\rho)\le\log\operatorname{tr}_w(1)$, con igualdad si
   y solo si $\rho_w=1/\operatorname{tr}_w(1)$.
3. Demuestre que el canal de medición $\rho\mapsto\sum_cP_c\rho P_c$ sobre $M_D$ (¡no sobre
   $\mathcal M$!) sí perturba los estados, y calcule $S(\rho)$ antes y después para una
   superposición de dos sectores.

<details>
<summary>Soluciones</summary>

1. $\tau_{\mathrm{vN}}(A)=\tau_{\mathrm{vN}}(1\cdot A)$, luego $\rho_{\tau}=1$ por unicidad. Por
   el teorema 3.2, la densidad `blunt` es $\bigoplus_c1_{k_c}/(Ck_c)$, con $p_c=1/C$ y
   $\sigma_c=1/k_c$, luego $S=\log C+\frac1C\sum_c\log k_c$ por el teorema 3.3.

```python
tau = states.tracial_state(alg)
assert (tau.density('tau_vN') - one.matrix).abs().max() < 1e-12
assert abs(tau.von_neumann_entropy().item() - (math.log(3) + sum(math.log(k) for k in (2, 3, 5)) / 3)) < 1e-12
```

2. Las sumas recorren los $c$ con $p_c>0$. Entonces
   $S_w(\rho)=S(\rho)+\sum_cp_c\log w_c=-\sum_cp_c\log p_c+\sum_cp_cS(\sigma_c)+\sum_cp_c\log w_c
   \le\sum_cp_c\log\frac{k_cw_c}{p_c}\le\log\sum_{p_c>0}k_cw_c\le\log\operatorname{tr}_w(1)$, usando
   $S(\sigma_c)\le\log k_c$ (corolario 3.6) y la desigualdad de Jensen finita
   $\sum_cp_c\log x_c\le\log\sum_cp_cx_c$ para la función cóncava $\log$ (que se deduce de la
   definición de concavidad por inducción sobre el número de términos), con
   $x_c=k_cw_c/p_c$; es estricta salvo que todos los $x_c$ sean iguales, ya que $\log$ es
   estrictamente cóncava. La igualdad en toda la cadena obliga a que todo $p_c>0$,
   $\sigma_c=1/k_c$ y $k_cw_c/p_c$ sea constante, es decir,
   $p_c=k_cw_c/\operatorname{tr}_w(1)$, es decir, $(\rho_w)_c=\rho_c/w_c=1_{k_c}/\operatorname{tr}_w(1)$.
   Recíprocamente, para $\rho_w=1/\operatorname{tr}_w(1)$,
   $S_w(\rho)=-\operatorname{tr}_w(1)\cdot\frac1{\operatorname{tr}_w(1)}\log\frac1{\operatorname{tr}_w(1)}=\log\operatorname{tr}_w(1)$.

3. En $M_D$, una superposición pura $\psi=(\psi_1\oplus\psi_2)/\sqrt2$ de vectores unitarios
   de dos sectores tiene $S=0$; tras la medición el estado es
   $\frac12|\psi_1\rangle\langle\psi_1|\oplus\frac12|\psi_2\rangle\langle\psi_2|$, con $S=\log2$.
   Restringido a $\mathcal M$, sin embargo, el estado vectorial de $\psi$ coincide con la
   mezcla:
   $\langle\psi,A\psi\rangle=\tfrac12\langle\psi_1,A_1\psi_1\rangle+\tfrac12\langle\psi_2,A_2\psi_2\rangle$
   para todo $A\in\mathcal M$. Ningún observable de $\mathcal M$ detecta la coherencia, y ese
   es el contenido de una regla de superselección.

</details>

**Siguiente:** el [capítulo 4](04_symmetries.md) trata las simetrías y las representaciones:
cómo la acción de un grupo sobre un sistema físico da lugar a los sectores, y cómo la
biblioteca construye las álgebras correspondientes para sistemas de espines y registros de
cúbits.
