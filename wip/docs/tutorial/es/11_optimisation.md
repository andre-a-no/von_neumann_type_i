[English](../11_optimisation.md) · **Español** · [Français](../fr/11_optimisation.md) · [Deutsch](../de/11_optimisation.md) · [中文](../zh/11_optimisation.md) · [日本語](../ja/11_optimisation.md) · [Русский](../ru/11_optimisation.md)

# 11. Optimización y desigualdades de trazas

Una desigualdad de trazas es un enunciado sobre los valores extremos de un funcional de operadores. La
sección 3.7 mostró una que es verosímil y falsa (proposición 3.9), y una ejecución de muestreo que, numéricamente,
no encontró la violación. Este capítulo explica por qué el muestreo falla en dimensión alta y lo sustituye por la optimización sobre conjuntos de
restricciones parametrizados de forma exacta; allí donde los extremos se conocen en forma cerrada, los
demostramos, de modo que el optimizador se contrasta con un teorema.

**Aprenderá** cómo la biblioteca parametriza las unitarias y los espectros con un contraste prescrito, y por qué
se alcanza todo operador admisible; la desigualdad de trazas de von Neumann y los extremos exactos que proporciona;
una desigualdad de trazas con su condición de igualdad; por qué las muestras de Haar están, en media cuadrática,
lejos de los extremos; y qué demuestra una búsqueda multiarranque.

En todo el capítulo, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ con $D=\sum_ck_c\ge2$, $\operatorname{Tr}$ es la
traza `blunt` (sección 3.1) y $\|A\|_1=\sum_c\|A_c\|_1$, la norma traza de $A$ como matriz diagonal por bloques en
$M_D(\mathbb C)$. Los valores singulares de $A\in M_k(\mathbb C)$ son $s_1(A)\ge\dots\ge s_k(A)$; $x^\downarrow$ y
$x^\uparrow$ son las reordenaciones decreciente y creciente de un vector real $x$. $\Delta$ es el contraste de
Michelson (sección 3.6), y para $\delta\in[0,1)$
$$
\ell(\delta)=\frac{1-\delta}{1+\delta},\qquad 1-\ell(\delta)=\frac{2\delta}{1+\delta};
$$
un operador positivo con $\lambda_{\max}=1$ tiene contraste $\ell(\lambda_{\min})$, luego tiene contraste $\delta$ si y
solo si $\lambda_{\min}=\ell(\delta)$, ya que $\ell$ es una involución estrictamente decreciente de $[0,1]$. Los
ejemplos son los tres funcionales del artículo que acompaña a la biblioteca (sección *From sampling to
optimisation*):

| | funcional | variables |
|---|---|---|
| Exp. 1 | $z=\lvert\operatorname{Tr}(XUY)\rvert-\operatorname{Tr}(XY)$ | $X,Y\ge0$, $U\in\mathcal M$ unitaria |
| Exp. 2 | $z=\operatorname{Tr}\lvert XY\rvert-\operatorname{Tr}(\lvert X\rvert\,\lvert Y\rvert)$ | $X=UX_0$, $Y=VY_0$, $X_0,Y_0\ge0$, $U,V$ unitarias |
| Exp. 3 | $z=\operatorname{Tr}(Y\lvert X\rvert Y)-\operatorname{Tr}\lvert YXY\rvert$ | $X=X^*$, $Y\ge0$ |

con $|X|$, $|Y|$ normalizados a $\lambda_{\max}=1$ y sus contrastes como parámetros. Las **envolventes** son
$\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}$ y $\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}$, con todas las
demás variables libres.

## 11.1 Conjuntos de restricciones y sus parametrizaciones

**Definición.** Para $\delta\in[0,1)$ y $\ell=\ell(\delta)$, la **capa de contraste** es
$$
\mathcal P_\delta=\{X\in\mathcal M:X\ge0,\ \lambda_{\max}(X)=1,\ \Delta(X)=\delta\}
=\{X=X^*:\operatorname{spec}X\subseteq[\ell,1],\ 1,\ell\in\operatorname{spec}X\}.
$$
Las **casillas** de $\mathcal M$ son los $D$ pares $(c,i)$, $i\le k_c$; para $\lambda\in\mathbb R^D$,
$\operatorname{diag}(\lambda)=\bigoplus_c\operatorname{diag}(\lambda_{c,1},\dots,\lambda_{c,k_c})$.

`optimize.PositiveParam` con contraste prescrito $\delta$ sortea, en cada arranque, dos casillas distintas
$s_\pm=(c_\pm,i_\pm)$ ($s_+$ uniformemente entre las $D$ casillas y luego $s_-$ uniformemente entre las $D-1$
restantes, de modo que todo par ordenado de casillas distintas tiene probabilidad $1/(D(D-1))$) y una unitaria
aleatoria de Haar $U_0\in\mathcal M$ (`random_unitary_operator`, lema 11.10), y envía $\theta\in\mathbb R^D$ y
$G=\bigoplus_cG_c$ a
$$
X(\theta,G)=U\operatorname{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell,\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
$$
`UnitaryParam` es la aplicación $G\mapsto U_0e^{G-G^*}$ sola. A diferencia de una sigmoide, $(1+\sin\theta)/2$
alcanza $0$ y $1$ para $\theta$ finito, lo cual es conveniente porque los espectros extremales de este capítulo
(corolario 11.7, proposición 11.14) están en $\ell$ y $1$. Un contraste libre en $[a,b]$ es un parámetro más,
$\delta=a+(b-a)(1+\sin\eta)/2$; el valor por defecto de la biblioteca es $[a,b]=[0,0.999]$, y el código de abajo lo
mantiene. Un valor alcanzado con esta restricción también se alcanza sin ella, de modo que sigue siendo una cota
inferior de un supremo sin restringir (una cota superior de un ínfimo sin restringir).

**Lema 11.1.** Para toda unitaria $U_0\in\mathcal M$, $G\mapsto U_0e^{G-G^*}$ es sobreyectiva sobre el grupo
unitario de $\mathcal M$.

*Demostración.* Sea $U$ unitaria. Por el teorema 1.3 (bloque a bloque) y el corolario 1.4(2),
$U_0^*U=\bigoplus_cW_c\operatorname{diag}(e^{i\varphi_{c,j}})_jW_c^*$ con $W_c$ unitarias, $\varphi_{c,j}\in\mathbb R$.
Entonces $S=\bigoplus_cW_c\operatorname{diag}(i\varphi_{c,j})_jW_c^*$ satisface $S^*=-S$ y $e^S=U_0^*U$ (la serie
$\sum_nS^n/n!$ se calcula término a término en la base propia), y $G=S/2$ da $G-G^*=S$. $\square$

**Proposición 11.2 (la parametrización es sobreyectiva).** Para $\delta$, $s_\pm$ y $U_0$ fijos, la imagen de
$(\theta,G)\mapsto X(\theta,G)$ es $\{X\in\mathcal P_\delta:1\in\operatorname{spec}X_{c_+},\ \ell\in\operatorname{spec}X_{c_-}\}$,
y $\mathcal P_\delta$ es la unión de estas imágenes sobre todos los pares de casillas distintas.

*Demostración.* La imagen está contenida en este conjunto: $X$ tiene valores propios $\lambda_s\in[\ell,1]$, con
$\lambda_{s_\pm}$ en los bloques $c_\pm$. Recíprocamente, diagonalícese $X_c=W_c\operatorname{diag}(\mu_c)W_c^*$ y
permútese cada $\mu_c$ (sustituyendo $W_c$ por $W_c\Pi_c$, con $\Pi_c$ una matriz de permutación) de modo que
$\mu_{s_+}=1$ y $\mu_{s_-}=\ell$. Esto es posible: si $c_+=c_-$ y $\delta>0$, entonces $1\ne\ell$ son dos valores
propios del mismo bloque y ocupan posiciones distintas; si $\delta=0$, todos los valores propios valen $1$. Si
$\ell<1$, póngase $\theta_s=\arcsin\big(2(\mu_s-\ell)/(1-\ell)-1\big)$ para las demás casillas, de modo que
$\lambda=\mu$ (si $\ell=1$, cualquier $\theta$ sirve). El lema 11.1 da $G$ con $U_0e^{G-G^*}=\bigoplus_cW_c$.
*Unión:* para $X\in\mathcal P_\delta$, elíjanse bloques $c_\pm$ que contengan los valores propios $1$, $\ell$ y
casillas distintas en ellos; si $c_+=c_-$ y $\delta>0$, el bloque tiene dos valores propios diferentes, luego
$k_{c_+}\ge2$; si $\delta=0$, $X=1$ y sirven dos casillas cualesquiera. $\square$

Como todo par de casillas se sortea con probabilidad positiva, todo punto de la capa está, con probabilidad
positiva, en la imagen de un arranque, y (en aritmética exacta) la restricción se cumple en cada iterado.
`SelfAdjointParam` usa $X=U\operatorname{diag}(\sigma\lambda)U^*$ con signos $\sigma\in\{\pm1\}^D$ sorteados en cada
arranque (en cada bloque, el número de signos menos es uniforme en $\{0,\dots,k_c\}$ y sus posiciones son
uniformemente aleatorias, de modo que todo patrón de signos tiene probabilidad positiva) y modificados solo por los
cambios de signo de la sección 11.5. Entonces $|X|=U\operatorname{diag}(\lambda)U^*$ (teorema 1.5(e), ya que
$\lambda>0$), y la unión sobre signos y casillas es $\{X=X^*:|X|\in\mathcal P_\delta\}$: para un tal $X$,
diagonalícese $X_c=W_c\operatorname{diag}(\nu_c)W_c^*$, de modo que $|X|_c=W_c\operatorname{diag}(|\nu_c|)W_c^*$ con
$|\nu_{c,i}|\ge\ell>0$, y aplíquese la demostración de la proposición 11.2 a $\mu=|\nu|$ con $\sigma_s=\nu_s/|\nu_s|$.
Solo importa $S=G-G^*$.

La carta $S\mapsto U_0e^S$ solo se comporta bien cerca de $S=0$. Sean $\mathfrak u=\{S\in\mathcal M:S^*=-S\}$ y
$\mathcal M_h=\{H\in\mathcal M:H^*=H\}$, subespacios reales con $\mathcal M=\mathfrak u\oplus\mathcal M_h$
($X=\frac{X-X^*}2+\frac{X+X^*}2$). *Cerca de $0$, $S\mapsto U_0e^S$ es un homeomorfismo de una bola abierta $N\ni0$
de $\mathfrak u$ sobre un entorno de $U_0$ en el grupo unitario.* En efecto, $\Psi(S,H)=U_0e^Se^H$ es una aplicación
$C^1$ $\mathfrak u\times\mathcal M_h\to\mathcal M$ (la serie $\sum_{m,n\ge0}U_0S^mH^n/(m!\,n!)$ y la serie de las
derivadas parciales de sus términos respecto de las coordenadas reales convergen uniformemente en los conjuntos
acotados), y $e^{tS}e^{tH}=1+t(S+H)+O(t^2)$ muestra que su diferencial en $0$ es $(S,H)\mapsto U_0(S+H)$, una
biyección sobre $\mathcal M$. Por el teorema de la función inversa existen bolas abiertas $N\ni0$ en $\mathfrak u$ y
$N_h\ni0$ en $\mathcal M_h$ tales que $\Psi$ es un homeomorfismo de $N\times N_h$ sobre un abierto
$O\subseteq\mathcal M$. Si $\Psi(S,H)$ es unitaria, también lo es $e^H=e^{-S}U_0^*\Psi(S,H)$ (en una base propia de la
matriz normal $S$, cuyos valores propios son imaginarios, $e^{\pm S}$ son diagonales con entradas de módulo $1$ e
inversas una de la otra); pero $e^H$ tiene los valores propios positivos $e^{h_i}$, siendo $h_i$ los valores propios
de $H$, luego $e^{h_i}=1$ (corolario 1.4(2)) y $H=0$. Por tanto $O\cap U(\mathcal M)=U_0e^N$, donde $U(\mathcal M)$ es
el grupo unitario de $\mathcal M$: un entorno de $U_0$ en $U(\mathcal M)$, y $S\mapsto U_0e^S$ es una biyección
continua $N\to U_0e^N$ cuya inversa, la primera componente de $\Psi^{-1}$, es continua. Lejos de $0$ la carta no es
inyectiva: $e^{2\pi iP}=1=e^0$ para toda proyección $P\in\mathcal M$. Por eso `extremize` **recentra** tras cada
ronda ($U_0\leftarrow U_0e^{G-G^*}$, $G\leftarrow0$). Otras dos cartas tienen la misma diferencial $S\mapsto U_0S$ en
$0$. La carta de Cayley $S\mapsto U_0C(-S/2)$ (ejercicio 1):
$(1+tS/2)(1-tS/2)^{-1}=(1+tS/2)(1+tS/2+O(t^2))=1+tS+O(t^2)$. La retracción polar, $S\mapsto$ el factor unitario de
$A=U_0(1+S)$ (teorema 1.8): $A^*A=(1-S)(1+S)=1+S^*S\ge1$ es invertible, luego $|A|$ es invertible y el factor
unitario es $A|A|^{-1}$; con $tS$ en lugar de $S$, $|A|^{-1}=(1+t^2S^*S)^{-1/2}=1+O(t^2)$ (en una base propia de
$S^*S$), luego el factor es $U_0(1+tS)+O(t^2)$.

El código comprueba la restricción, y la proposición 11.2 de forma constructiva: los parámetros de la demostración
reproducen un objetivo.

```python
import math
import scipy.linalg
import torch
from torch_vn_algebra import TypeIAlgebra, optimize as opt

torch.manual_seed(0)
alg = TypeIAlgebra([3, 3], [3, 3], precision='double')          # M_3 (+) M_3, D = 6
ell = lambda delta: (1 - delta) / (1 + delta)                     # lambda_min for contrast delta
d = torch.tensor([0.0, 0.25, 0.5, 0.9], dtype=torch.float64)
x = opt.PositiveParam(alg, 4, delta=d).operator()                 # four starts, one contrast each
assert torch.allclose(x.michelson_contrast, d, atol=1e-12)        # the constraint holds exactly
assert torch.allclose(x.lambda_max, torch.ones(4, dtype=torch.float64), atol=1e-12)

P = opt.PositiveParam(alg, 1, delta=torch.tensor([0.5], dtype=torch.float64))
lo = ell(0.5)
mu = lo + (1 - lo) * torch.rand(1, 2, 3, dtype=torch.float64)    # target spectrum in [lo, 1] ...
(cp, ip), (cm, im) = P.spectrum.pin_max[0].tolist(), P.spectrum.pin_min[0].tolist()
mu[0, cp, ip], mu[0, cm, im] = 1.0, lo                            # ... with 1 and lo in the drawn slots
W = alg.random_unitary_operator(1).matrix                         # target eigenvectors
target = W @ torch.diag_embed(mu.to(W.dtype)) @ W.mH
with torch.no_grad():                                             # the parameters of the proof
    P.spectrum.theta.copy_(torch.asin((2 * (mu - lo) / (1 - lo) - 1).clamp(-1, 1)))
    V = P.unitary.base.mH @ W                                     # U0^* W = exp(S), G = S / 2
    S = torch.stack([torch.from_numpy(scipy.linalg.logm(V[0, c].numpy())) for c in range(2)])
    P.unitary.G.copy_(S[None] / 2)
err = (P.operator().matrix - target).abs().max().item()
print(f"slots of 1 and lo: {(cp, ip)}, {(cm, im)};  reconstruction error {err:.1e}")
assert err < 1e-10
```

> **Significado físico** (interpretación física, no un teorema de este curso). El contraste de Michelson es la
> visibilidad $(I_{\max}-I_{\min})/(I_{\max}+I_{\min})$ de las franjas de interferencia. Para $0\le X\le1$ con
> $\lambda_{\max}=1$ —un filtro, un atenuador, la eficiencia de un detector— fija el rango dinámico
> $\lambda_{\min}/\lambda_{\max}$, el inverso de un número de condición.

## 11.2 La desigualdad de trazas de von Neumann

Una matriz $(c_{ij})$ con $c_{ij}\ge0$ y todas las sumas por filas y por columnas iguales a $1$ es **doblemente
estocástica**; para una unitaria $G$, $(|G_{ij}|^2)$ y $(|G_{ji}|^2)$ son doblemente estocásticas, pues las filas y
las columnas de $G$ son vectores unitarios.

**Lema 11.3 (reordenación).** Sean $p_1\ge\dots\ge p_k$, $q_1\ge\dots\ge q_k$ reales y $(c_{ij})$ doblemente
estocástica. Entonces $\sum_ip_iq_{k+1-i}\le\sum_{i,j}p_iq_jc_{ij}\le\sum_ip_iq_i$.

*Demostración.* *Cota superior.* Para $m\le k$, sea $r_j=\sum_{i\le m}c_{ij}\in[0,1]$, de modo que $\sum_jr_j=m$.
Como $q_j\ge q_m$ para $j\le m$ y $q_j\le q_m$ para $j>m$,
$$
\sum_{j\le m}q_j-\sum_jq_jr_j=\sum_{j\le m}q_j(1-r_j)-\sum_{j>m}q_jr_j\ge q_m\Big(m-\sum_jr_j\Big)=0 .
$$
Con $x_i=\sum_jq_jc_{ij}$, esto dice que $\sum_{i\le m}x_i\le\sum_{j\le m}q_j$, con igualdad para $m=k$. Por la
suma por partes de Abel y $p_m-p_{m+1}\ge0$,
$$
\sum_ip_ix_i=p_k\sum_{i\le k}x_i+\sum_{m<k}(p_m-p_{m+1})\sum_{i\le m}x_i
\le p_k\sum_{j\le k}q_j+\sum_{m<k}(p_m-p_{m+1})\sum_{j\le m}q_j=\sum_ip_iq_i .
$$
*Cota inferior.* Aplíquese la cota superior a $p$, a la sucesión decreciente $q'_j=-q_{k+1-j}$ y a la matriz
doblemente estocástica $c'_{ij}=c_{i,k+1-j}$:
$-\sum_{ij}p_iq_jc_{ij}=\sum_{ij}p_iq'_jc'_{ij}\le\sum_ip_iq'_i=-\sum_ip_iq_{k+1-i}$. $\square$

Toda $A\in M_k(\mathbb C)$ tiene una **descomposición en valores singulares** $A=W\Sigma V^*$, con $W,V$ unitarias y
$\Sigma=\operatorname{diag}(s_i(A))$: escríbase $A=W'|A|$ (teorema 1.8), $|A|=V\Sigma V^*$ (teorema 1.3) y $W=W'V$.

**Teorema 11.4 (desigualdad de trazas de von Neumann).** Para $A,B\in M_k(\mathbb C)$,
$|\operatorname{Tr}(AB)|\le\sum_is_i(A)\,s_i(B)$.

*Demostración.* Sean $A=W_A\Sigma_AV_A^*$, $B=W_B\Sigma_BV_B^*$, $a_i=s_i(A)$, $b_j=s_j(B)$, y $G=V_A^*W_B$,
$H=V_B^*W_A$. Por ciclicidad, $\operatorname{Tr}(AB)=\operatorname{Tr}(\Sigma_AG\Sigma_BH)=\sum_{i,j}a_iG_{ij}b_jH_{ji}$,
luego, por $|g||h|\le\frac12(|g|^2+|h|^2)$ y $a_ib_j\ge0$,
$$
|\operatorname{Tr}(AB)|\le\sum_{i,j}a_ib_j\,d_{ij},\qquad d_{ij}=\tfrac12\big(|G_{ij}|^2+|H_{ji}|^2\big).
$$
$(d_{ij})$ es doblemente estocástica, y el lema 11.3 da $\sum_{ij}a_ib_jd_{ij}\le\sum_ia_ib_i$. $\square$

**Corolario 11.5.** (a) Para $P,Q\in M_k(\mathbb C)$ autoadjuntas con valores propios $p,q$:
$\sum_ip^\downarrow_iq^\uparrow_i\le\operatorname{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i$.
(b) Para $P,Q\ge0$: $\operatorname{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i$.

*Demostración.* (a) Como en la demostración del teorema 3.5, $\operatorname{Tr}(PQ)=\sum_{ij}p_iq_jc_{ij}$ con la matriz
doblemente estocástica $c_{ij}=|\langle u_i,v_j\rangle|^2$ para bases propias ordenadas de forma decreciente;
aplíquese el lema 11.3 (su cota inferior es $\sum_ip^\downarrow_iq^\downarrow_{k+1-i}=\sum_ip^\downarrow_iq^\uparrow_i$).
(b) La primera desigualdad es el teorema 1.9. Por el teorema 1.9 y el teorema 11.4,
$\|PQ\|_1=\max_U|\operatorname{Tr}(P\,QU)|\le\sum_is_i(P)s_i(QU)$, y $s_i(QU)=q^\downarrow_i$ porque
$|QU|=(U^*Q^2U)^{1/2}=U^*QU$ (teorema 1.6). $\square$

Así, para $X,Y\ge0$ en $M_k(\mathbb C)$ con espectros $x,y$, la violación $\|YX\|_1-\operatorname{Tr}(XY)$ de la
proposición 3.9 es a lo sumo $\sum_ix^\downarrow_i(y^\downarrow_i-y^\uparrow_i)$ (corolario 11.5(b) para $(Y,X)$ y la
cota inferior de (a)). Para Exp. 2, este número es la respuesta exacta:

**Proposición 11.6 (Exp. 2 con espectros fijos).** En Exp. 2, $|X|=X_0$, $|Y|=Y_0$ y
$z=\|X_0VY_0\|_1-\operatorname{Tr}(X_0Y_0)$. Si los espectros $x_c$ de $(X_0)_c$ e $y_c$ de $(Y_0)_c$ están fijos y
las bases propias y $V$ varían, entonces
$$
\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}\big(y^\downarrow_{c,i}-y^\uparrow_{c,i}\big),
$$
y ambos se alcanzan.

*Demostración.* $X^*X=X_0U^*UX_0=X_0^2$, luego $|X|=X_0$ (teorema 1.6); análogamente $|Y|=Y_0$, y $|UA|=|A|$ da
$\operatorname{Tr}|XY|=\|X_0VY_0\|_1$. Todas las variables actúan bloque a bloque y $z$ es una suma sobre los bloques,
así que podemos suponer $C=1$. Póngase $Y'=VY_0V^*\ge0$ (espectro $y$); como $|AV|=V^*|A|V$,
$\|X_0VY_0\|_1=\|X_0Y'V\|_1=\|X_0Y'\|_1$. El corolario 11.5 para los pares $(X_0,Y')$ y $(X_0,Y_0)$, junto con
$\|X_0Y'\|_1\ge\operatorname{Tr}(X_0Y')$, da
$$
\sum_ix^\downarrow_iy^\uparrow_i-\sum_ix^\downarrow_iy^\downarrow_i\le\|X_0Y'\|_1-\operatorname{Tr}(X_0Y_0)\le\sum_ix^\downarrow_iy^\downarrow_i-\sum_ix^\downarrow_iy^\uparrow_i .
$$
Sea $J$ la permutación de inversión, $Je_i=e_{k+1-i}$; entonces $J\operatorname{diag}(w)=\operatorname{diag}(w_k,\dots,w_1)J$,
y la norma traza de $\operatorname{diag}(w)J$ es $\sum_i|w_i|$. Con $X_0=\operatorname{diag}(x^\downarrow)$ y $V=J$:
para $Y_0=\operatorname{diag}(y^\uparrow)$, $X_0VY_0=\operatorname{diag}(x^\downarrow_iy^\downarrow_i)J$ y
$\operatorname{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\uparrow_i$ (la cota superior); para $Y_0=\operatorname{diag}(y^\downarrow)$,
$X_0VY_0=\operatorname{diag}(x^\downarrow_iy^\uparrow_i)J$ y $\operatorname{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\downarrow_i$
(la cota inferior). $\square$

**Corolario 11.7 (envolvente de Exp. 2).** Sean $\ell_X,\ell_Y\in(0,1]$ y supongamos que $k_c\ge2$ para algún $c$. Sobre todos los
$X_0,Y_0\ge0$ con $\lambda_{\max}=1$, $\lambda_{\min}(X_0)=\ell_X$, $\lambda_{\min}(Y_0)=\ell_Y$, y todas las unitarias,
$$
\sup z=-\inf z=\sum_c\Big\lfloor\frac{k_c}2\Big\rfloor(1-\ell_X)(1-\ell_Y).
$$

*Demostración.* En un bloque, escríbase $x_i=x^\downarrow_{c,i}$, $y_i=y^\downarrow_{c,i}$, $k=k_c$, de modo que
$y^\uparrow_{c,i}=y_{k+1-i}$. Emparejando $i$ con $k+1-i$,
$x_i(y_i-y_{k+1-i})+x_{k+1-i}(y_{k+1-i}-y_i)=(x_i-x_{k+1-i})(y_i-y_{k+1-i})$, y para $k$ impar el índice central
aporta $x_i(y_i-y_i)=0$. Por tanto
$$
\sum_ix_i(y_i-y_{k+1-i})=\sum_{i\le\lfloor k/2\rfloor}(x_i-x_{k+1-i})(y_i-y_{k+1-i})\le\Big\lfloor\frac k2\Big\rfloor(1-\ell_X)(1-\ell_Y),
$$
ya que $x_i\in[\ell_X,1]$, $y_i\in[\ell_Y,1]$ son decrecientes, de modo que ambos factores están en $[0,1-\ell_X]$ y
$[0,1-\ell_Y]$. La igualdad se da si en cada bloque los $\lfloor k_c/2\rfloor$ valores propios mayores valen $1$ y los
$\lfloor k_c/2\rfloor$ menores valen $\ell_X$ (resp. $\ell_Y$), siendo el central ($k_c$ impar) arbitrario en
$[\ell_X,1]$ (resp. $[\ell_Y,1]$); estos espectros tienen extremos globales $1$ y $\ell_X$ ($\ell_Y$) porque $k_c\ge2$ para algún $c$. La
proposición 11.6 permite concluir. $\square$

Aquí $1-\ell_X=2\delta/(1+\delta)$ para $\Delta(|X|)=\delta$; para $k_c$ impar, el valor propio central es irrelevante.
El código comprueba el teorema 11.4, el corolario 11.5 y las configuraciones extremales de la proposición 11.6 para
$k=3$, $C=2$.

```python
M5 = TypeIAlgebra([5], [5], precision='double')
A = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))         # Ginibre, no structure
B = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))
sA, sB = torch.linalg.svdvals(A.matrix[:, 0]), torch.linalg.svdvals(B.matrix[:, 0])   # decreasing
assert ((A @ B).trace.abs() <= (sA * sB).sum(-1) + 1e-12).all()             # Theorem 11.4
Ph, Qh = (M5.operator_from_eigenvalues(lambda k: torch.randn(2000, k), batch_size=2000, force_self_adjoint=True)
          for _ in range(2))
p, q = Ph.eigenvalues()[0], Qh.eigenvalues()[0]                             # increasing order
t = (Ph @ Qh).trace.real
assert (t <= (p * q).sum(-1) + 1e-12).all() and (t >= (p * q.flip(-1)).sum(-1) - 1e-12).all()   # Cor. 11.5(a)
Pp, Qp = Ph @ Ph, Qh @ Qh                                                   # positive, spectra p^2, q^2
nuc = (Pp @ Qp).trace_norm()
assert ((Pp @ Qp).trace.real <= nuc + 1e-12).all()                          # Cor. 11.5(b)
assert (nuc <= ((p ** 2).sort(-1)[0] * (q ** 2).sort(-1)[0]).sum(-1) + 1e-10).all()

alg3 = TypeIAlgebra([3, 3], [3, 3], precision='double')
x_dn = torch.tensor([1.0, 0.6, ell(0.5)], dtype=torch.float64)              # Delta = 0.5, middle value 0.6
y_dn = torch.tensor([1.0, 0.7, ell(0.8)], dtype=torch.float64)              # Delta = 0.8, middle value 0.7
J = torch.eye(3, dtype=torch.float64).flip(0)                               # reversal permutation
blocks = lambda m: alg3.from_blocks([m, m])                                 # the same block in both sectors

def z_exp2(X0, V, Y0):
    return ((X0 @ V @ Y0).trace_norm() - (X0 @ Y0).trace.real).item()

X0 = blocks(torch.diag(x_dn))
z_max = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn.flip(0))))
z_min = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn)))
exact = 2 * (x_dn * (y_dn - y_dn.flip(0))).sum().item()                     # Proposition 11.6
envelope = 2 * (1 - ell(0.5)) * (1 - ell(0.8))                              # Corollary 11.7, floor(3/2) = 1
print(f"Exp. 2, k = 3, C = 2: z_max = {z_max:.6f}, z_min = {z_min:.6f}, formula {envelope:.6f}")
assert abs(z_max - exact) < 1e-12 and abs(z_min + exact) < 1e-12 and abs(exact - envelope) < 1e-12
```

## 11.3 Una desigualdad de trazas con condición de igualdad

Para $X=X^*$, sea $X_\pm=f_\pm(X)$ con $f_+(t)=\max(t,0)$, $f_-(t)=\max(-t,0)$. Entonces $X_\pm\ge0$, $X=X_+-X_-$,
$|X|=X_++X_-$ y $X_+X_-=0$ (teorema 1.5(c), y 1.5(e) para $|X|=(X^2)^{1/2}=h(X)$ con
$h(t)=\sqrt{t^2}=|t|=f_+(t)+f_-(t)$; corolario 1.4(3) para la positividad).

**Lema 11.8.** Para $A,B\ge0$ en $M_n(\mathbb C)$, $\|A-B\|_1\le\operatorname{Tr}A+\operatorname{Tr}B$, con igualdad si y
solo si $AB=0$.

*Demostración.* Sean $H=A-B$ y $Q=g(H)$ con $g=1$ en $(0,\infty)$ y $g=0$ en el resto; $Q$ es una proyección, ya que
$g=\bar g=g^2$ (teorema 1.5(c)). Como $|t|=t(2g(t)-1)$, se tiene $|H|=H(2Q-1)$ (teorema 1.5(c),(e)), y, usando
$\operatorname{Tr}(AR)=\operatorname{Tr}(RAR)$ para las proyecciones $R=Q,1-Q$ (ciclicidad),
$$
\|H\|_1=\operatorname{Tr}\big(A(2Q-1)\big)-\operatorname{Tr}\big(B(2Q-1)\big)
=\operatorname{Tr}A+\operatorname{Tr}B-2\operatorname{Tr}\big((1-Q)A(1-Q)\big)-2\operatorname{Tr}(QBQ).
$$
Las trazas restadas son $\ge0$, y $\operatorname{Tr}(RAR)=\|A^{1/2}R\|_2^2$ se anula si y solo si $AR=0$; así, la
igualdad se da si y solo si $A(1-Q)=0$ y $BQ=0$, y entonces $AB=AQB=A(BQ)^*=0$. Recíprocamente, si $AB=0$,
entonces $BA=0$, $(A-B)^2=A^2+B^2=(A+B)^2$, y el teorema 1.6 da $|A-B|=A+B$. $\square$

**Proposición 11.9 (Exp. 3).** Para $X=X^*$ e $Y\ge0$ en $\mathcal M$,
$\operatorname{Tr}|YXY|\le\operatorname{Tr}(Y|X|Y)$, con igualdad si y solo si $YX_+Y^2X_-Y=0$, en particular si
$XY=YX$.

*Demostración.* En $M_D(\mathbb C)$, $A=YX_+Y$ y $B=YX_-Y$ son positivos (corolario 1.7 con $T=X_\pm^{1/2}Y$),
$YXY=A-B$ y $\operatorname{Tr}(Y|X|Y)=\operatorname{Tr}A+\operatorname{Tr}B$; aplíquese el lema 11.8. Si $XY=YX$,
entonces $Y$ conmuta con $X_\pm$ (teorema 1.5(d)), y $AB=Y^4X_+X_-=0$. $\square$

Así, $z\ge0$ en Exp. 3, y $\underline z=0$ (tómese $Y=1$). La condición de igualdad es más débil que la
conmutación: para $X=\operatorname{diag}(1,-1,0)$ e $Y$ con filas $(1,0,a)$, $(0,1,0)$, $(a,0,1)$, $|a|<1$, se tiene
$AB=(Ye_1)\langle Ye_1,Ye_2\rangle(Ye_2)^*=0$, mientras que $XY\ne YX$ para $a\ne0$.

```python
def z_exp3(X, Y):
    return (Y @ X.abs() @ Y).trace.real - (Y @ X @ Y).trace_norm()

Xs = alg3.operator_from_eigenvalues(lambda k: torch.randn(500, k), batch_size=500, force_self_adjoint=True)
Yp = alg3.operator_from_eigenvalues(lambda k: torch.rand(500, k), batch_size=500, force_positive=True)
assert (z_exp3(Xs, Yp) >= -1e-12).all()                                     # Proposition 11.9

M3 = TypeIAlgebra([3], [3], precision='double')
Xe = M3.operator(torch.diag(torch.tensor([1.0, -1.0, 0.0]))[None, None])
Ye = M3.operator(torch.tensor([[1.0, 0, 0.5], [0, 1, 0], [0.5, 0, 1]])[None, None])
print(f"equality without commutation: z = {z_exp3(Xe, Ye).item():.1e}, "
      f"||XY - YX|| = {(Xe @ Ye - Ye @ Xe).operator_norm().item():.2f}")
assert abs(z_exp3(Xe, Ye).item()) < 1e-12 and (Xe @ Ye - Ye @ Xe).operator_norm().item() > 0.1
```

## 11.4 Muestreo frente a optimización

¿Qué tamaño tiene típicamente $|\operatorname{Tr}(AU)|$, comparado con su máximo $\|A\|_1$? La biblioteca muestrea $U$
como el factor $Q$ de la descomposición QR $Z=QR$, $R_{ii}>0$, de una **matriz de Ginibre** $Z$ (entradas
independientes con densidad $\pi^{-1}e^{-|z|^2}$; la receta de Mezzadri), de forma independiente en cada bloque.

**Lema 11.10.** $Q$ está distribuida según la medida de Haar en $U(k)$; en particular, para toda unitaria fija
$V\in M_k(\mathbb C)$, $VQ$ tiene la ley de $Q$. La unitaria $\bigoplus_cQ_c$ con bloques independientes está
distribuida según la medida de Haar en el grupo unitario de $\mathcal M$.

*Demostración.* Es el teorema 5.3 (la invariancia por la izquierda es el paso de equivariancia de su demostración;
la ley se identifica con la medida de Haar mediante el lema de unicidad 5.1, cuyo paso de teoría de la medida es el
corolario A.2) y el corolario 5.4. $\square$

**Proposición 11.11.** Para $U\in U(k)$ de Haar y $A\in M_k(\mathbb C)$: $\mathbb E\operatorname{Tr}(AU)=0$ y
$\mathbb E|\operatorname{Tr}(AU)|^2=\|A\|_2^2/k$.

*Demostración.* $\mathbb EU_{ij}=0$ por el corolario 5.5. La proposición 5.7 (con $C=1$) para
$X=|e_j\rangle\langle e_m|$ da la fórmula de Weingarten de orden uno,
$\mathbb E[U_{ij}\bar U_{lm}]=\mathbb E(UXU^*)_{il}=\delta_{il}\delta_{jm}/k$; para $k\ge2$ también se deduce de la
fórmula de orden dos, el teorema 5.10, haciendo $i_2=i'_2$, $j_2=j'_2$ y sumando sobre ambos, ya que
$\sum_{i_2,j_2}|U_{i_2j_2}|^2=k$. Con $\operatorname{Tr}(AU)=\sum_{ij}A_{ji}U_{ij}$,
$\mathbb E|\operatorname{Tr}(AU)|^2=\sum_{ijlm}A_{ji}\bar A_{ml}\mathbb E[U_{ij}\bar U_{lm}]=\frac1k\sum_{ij}|A_{ji}|^2$. $\square$

**Corolario 11.12.** Sean $X,Y\ge0$ en $M_k(\mathbb C)$ con espectros en $[a,b]$, $a>0$, y $U$ de Haar. Entonces
$$
\frac{\big(\mathbb E|\operatorname{Tr}(XUY)|^2\big)^{1/2}}{\max_V|\operatorname{Tr}(XVY)|}=\frac{\|YX\|_2}{\sqrt k\,\|YX\|_1}\le\min\Big(\frac1{\sqrt k},\frac{b^2}{a^2k}\Big),
\qquad
\mathbb P\big(|\operatorname{Tr}(XUY)|\ge\operatorname{Tr}(XY)\big)\le\frac{b^4}{a^4k^2}.
$$

*Demostración.* La igualdad es la proposición 11.11 para $A=YX$ (pues
$\operatorname{Tr}(XUY)=\operatorname{Tr}(AU)$) junto con la proposición 3.9. Como
$\|A\|_2^2=\sum_is_i^2\le(\sum_is_i)^2=\|A\|_1^2$ (sección 1.4), el cociente es $\le1/\sqrt k$. Como
$\|Av\|=\||A|v\|$ (sección 1.4), diagonalizando $|A|$ se obtiene $s_k(A)=\min_{\|v\|=1}\|Av\|$ y
$s_1(A)=\max_{\|v\|=1}\|Av\|$; y $\|Xv\|^2=\langle v,X^2v\rangle\in[a^2,b^2]$ para $v$ unitario (corolario 1.4(3)
para $X^2-a^2$ y $b^2-X^2$), y lo mismo para $Y$. Por tanto $a^2\le s_k(A)\le s_1(A)\le b^2$, $\|A\|_2^2\le kb^4$ y
$\|A\|_1\ge ka^2$. Por último,
$\operatorname{Tr}(XY)=\operatorname{Tr}(X^{1/2}YX^{1/2})\ge a\operatorname{Tr}X\ge ka^2$, y la desigualdad de Markov para
$|\operatorname{Tr}(XUY)|^2$ acota la probabilidad por $(kb^4/k)/(ka^2)^2$. $\square$

Así, $N$ muestras de Haar encuentran una violación de $|\operatorname{Tr}(XUY)|\le\operatorname{Tr}(XY)$ con probabilidad
a lo sumo $Nb^4/(a^4k^2)$ (una unión de $N$ sucesos), aunque existe una para todo par que no conmuta
(proposición 3.9). La cota solo es informativa para $k>\sqrt N\,b^2/a^2$ (para $b/a=6$ y $N=2000$: $k>1609$);
numéricamente, las muestras pasan por alto la violación ya para $k$ mucho menores (véase abajo). El código comprueba
la proposición 11.11 por Monte Carlo (la media muestral está a menos de cinco errores estándar del valor exacto, en
esta ejecución con semilla fija) e imprime, para espectros uniformes en $[0.2,1.2]$ ($b/a=6$), el cociente del
corolario 11.12, $\sup z$ y el mejor $z$ de $2000$ muestras, ambos relativos a $\operatorname{Tr}(XY)$.

```python
def random_pair(alg_k):
    pos = lambda k: 0.2 + torch.rand(1, k)                                  # spectrum in [a, b] = [0.2, 1.2]
    return tuple(alg_k.operator_from_eigenvalues(pos, force_positive=True) for _ in range(2))

M6 = TypeIAlgebra([6], [6], precision='double')
Xf, Yf = random_pair(M6)
N = 20_000
t2 = (Yf @ Xf @ M6.random_unitary_operator(N)).trace.abs() ** 2
m2, err, exact2 = t2.mean().item(), t2.std().item() / math.sqrt(N), (Yf @ Xf).frobenius_norm().item() ** 2 / 6
print(f"E|Tr(AU)|^2 = {m2:.4f} +- {err:.4f}   (exact ||A||_2^2 / k = {exact2:.4f})")
assert abs(m2 - exact2) < 5 * err                                           # Proposition 11.11

print(" k   rms/max   sup z/TrXY   best of 2000 samples")
for k in (2, 4, 8, 16, 32):
    Mk = TypeIAlgebra([k], [k], precision='double')
    Xk, Yk = random_pair(Mk)
    trXY, nuc = (Xk @ Yk).trace.real.item(), (Yk @ Xk).trace_norm().item()
    ratio = (Yk @ Xk).frobenius_norm().item() / (math.sqrt(k) * nuc)
    best = (Xk @ Mk.random_unitary_operator(2000) @ Yk).trace.abs().max().item()
    print(f"{k:2d}   {ratio:.4f}    {(nuc - trXY) / trXY:.2e}     {(best - trXY) / trXY:+.3f}")
    assert ratio <= min(1 / math.sqrt(k), 36 / k) + 1e-12                   # Corollary 11.12
    assert nuc > trXY and best < nuc and (k < 4 or best < trXY)            # from k = 4 on: no violation found
```

*Observación numérica* (la tabla impresa; afirmada solo para esta ejecución con semilla fija): el cociente se
mantiene por debajo de la cota demostrada $\min(1/\sqrt k,b^2/(a^2k))$ y se reduce aproximadamente a la mitad cada
vez que $k$ se duplica, es decir, decrece como $1/k$, y a partir de $k=4$ la mejor de $2000$ muestras ni siquiera
alcanza $z=0$. Una búsqueda por gradiente, en cambio, se acerca a ambos extremos en la ejecución de abajo (de nuevo,
una observación numérica). Para el par `Xf, Yf` se conocen ambos extremos sobre $U$:

**Lema 11.13.** Si $k\ge2$, entonces $\min_U|\operatorname{Tr}(AU)|=0$ para toda $A\in M_k(\mathbb C)$.

*Demostración.* Sean $A=W\Sigma V^*$ y $P$ el desplazamiento cíclico, $Pe_i=e_{(i\bmod k)+1}$, cuya diagonal se anula
para $k\ge2$. Para $U=VPW^*$, $\operatorname{Tr}(AU)=\operatorname{Tr}(\Sigma P)=0$. $\square$

Así, $\inf_Uz=-\operatorname{Tr}(XY)$ y $\sup_Uz=\|YX\|_1-\operatorname{Tr}(XY)$ (proposición 3.9). Una sola llamada a
`extremize` aproxima ambos, con $12$ arranques que maximizan $z$ y $12$ que maximizan $-z$ (con una precisión relativa
de $10^{-3}$ en esta ejecución, afirmada abajo). Por la proposición 11.11, la media cuadrática de
$|\operatorname{Tr}(XUY)|$ es $\|YX\|_2/\sqrt k$; con las cotas de la demostración del corolario 11.12 (espectros en
$[a,b]$), $\mathbb E|\operatorname{Tr}(XUY)|^2\le b^4$ y $\operatorname{Tr}(XY)\ge ka^2$, de modo que la desigualdad de
Markov da $\mathbb P\big(z\ge-(1-\varepsilon)\operatorname{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)$ para
$\varepsilon>0$: para $k$ grande, $z$ está cerca de su ínfimo con probabilidad alta.
*Observación numérica* (afirmada para esta ejecución con semilla fija): el muestreo llega a menos de un $1\,\%$ del
ínfimo, pero ni siquiera alcanza $z=0$, mientras que el supremo es positivo:

```python
trXY = (Xf @ Yf).trace.real
sup_exact, inf_exact = ((Yf @ Xf).trace_norm() - trXY).item(), -trXY.item()
U6 = opt.UnitaryParam(M6, 24)                                               # 24 starts
pm = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(12)   # 12 maximise z, 12 maximise -z
best = opt.extremize(lambda: pm * ((Xf @ U6.operator() @ Yf).trace.abs() - trXY), [U6],
                     maximize=True, rounds=6, steps=30).reshape(2, 12).max(1)[0]
best_up, best_down = best[0].item(), -best[1].item()
z_mc = (Xf @ M6.random_unitary_operator(N) @ Yf).trace.abs() - trXY
print(f"sup z: exact {sup_exact:.6e}, ascent {best_up:.6e}, best sample {z_mc.max().item():+.4f}")
print(f"inf z: exact {inf_exact:.6f}, descent {best_down:.6f}, worst sample {z_mc.min().item():+.6f}")
assert best_up <= sup_exact + 1e-12 and best_down >= inf_exact - 1e-12    # attained values: one-sided
assert sup_exact - best_up < 1e-3 * sup_exact and best_down - inf_exact < 1e-3 * abs(inf_exact)
assert z_mc.min().item() - inf_exact < 1e-2 * abs(inf_exact)               # sampling finds the typical side ...
assert z_mc.max().item() < 0                                               # ... but misses the violation
```

> **Significado físico** (interpretación física, no un teorema de este curso). El corolario 11.12 es un ejemplo
> débil y demostrable de lo que la mecánica estadística llama **tipicidad**: en dimensión alta se espera que casi
> todo estado o toda evolución dé aproximadamente el mismo valor de una magnitud macroscópica. Una desigualdad
> óptima se decide en configuraciones finamente ajustadas —bases propias alineadas, espectros de dos niveles— que
> una unitaria aleatoria (un circuito aleatorio, una matriz de dispersión aleatoria) rara vez produce.

## 11.5 Búsqueda multiarranque de envolventes

`optimize.extremize(objective, params, maximize)` ejecuta Adam en `rounds` rondas de `steps` pasos (tasa de
aprendizaje `lr * decay**round`), recentra las cartas unitarias tras cada ronda y devuelve, **para cada
arranque**, el mejor valor alcanzado. Cada elemento del lote es un arranque independiente con sus propias casillas,
signos, $U_0$ y contraste. Los signos de `SelfAdjointParam` son discretos: entre rondas, `discrete_step` prueba
ocho cambios de signo individuales por arranque, uno tras otro, y conserva cada cambio que mejora el objetivo.

Como todo iterado es admisible (proposición 11.2), todo valor devuelto se alcanza: una maximización da una **cota
inferior** del supremo, una minimización una cota superior del ínfimo, y nada más: una búsqueda local puede
detenerse en un extremo local. La otra dirección requiere una demostración, como el corolario 11.7 y la
proposición 11.9, o una relajación convexa (sección 11.6). Una envolvente calculada se contrasta (i) con formas
cerradas; (ii) con configuraciones construidas a mano, como sumas directas de soluciones más pequeñas, por debajo
de las cuales no debe quedar; (iii) con muestras estructuradas (espectros de dos niveles, bases de permutación) que
se acercan a los extremos en toda dimensión, a diferencia de las muestras de Haar (`scripts/bounds_checks.py`).

**Proposición 11.14 (envolvente inferior de Exp. 1).** Supongamos que todos los $k_c\ge2$. Sobre
$X\in\mathcal P_\delta$, $Y\ge0$ con $\lambda_{\max}(Y)=1$ y unitarias $U\in\mathcal M$,
$\ \min z=-\big(D-2\delta/(1+\delta)\big)$.

*Demostración.* $Y\le1$ da $X^{1/2}YX^{1/2}\le X$, luego $\operatorname{Tr}(XY)\le\operatorname{Tr}X\le(D-1)+\ell(\delta)$,
pues un valor propio de $X$ es $\ell(\delta)$ y los demás son $\le1$; por tanto
$z\ge-\operatorname{Tr}(XY)\ge-(D-2\delta/(1+\delta))$. La igualdad se da para $X$ con el valor propio $\ell(\delta)$
una vez y $1$ en el resto, $Y=1$, y $U=\bigoplus_cU_c$ con $\operatorname{Tr}(X_cU_c)=0$ (lema 11.13). $\square$

El código comprueba la proposición 11.14 en $M_2\oplus M_2$, el corolario 11.7 en $M_3\oplus M_3$ (ambos contrastes
prescritos, $k=3$ impar) y la proposición 11.9; supremos e ínfimos comparten una llamada, como antes. La palabra
clave `is_self_adjoint` selecciona el algoritmo para la norma traza (valores singulares, o valores propios para
$YXY$) sin una comprobación en cada paso. La última comprobación es del tipo (ii): el $z$ de Exp. 3 es aditivo
respecto de las sumas directas, y el contraste de una suma directa es el mayor contraste de sus sumandos
(ejercicio 3), luego
$$
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\max_{\delta'\le\delta}\overline z_{M_2}(\delta').
$$

```python
torch.manual_seed(1)
S = 6                                                                       # starts per grid point
grid = torch.tensor([0.0, 0.5, 0.9], dtype=torch.float64)
alg22 = TypeIAlgebra([2, 2], [2, 2], precision='double')                    # D = 4

# Exp. 1, lower envelope along Delta(X); Delta(Y) free in [0, 0.999]
X1, Y1, U1 = (opt.PositiveParam(alg22, 3 * S, delta=grid.repeat_interleave(S)),
              opt.PositiveParam(alg22, 3 * S), opt.UnitaryParam(alg22, 3 * S))
def z1():
    x, y = X1.operator(), Y1.operator()
    return (x @ U1.operator() @ y).trace.abs() - (x @ y).trace.real
inf1 = opt.extremize(z1, [X1, Y1, U1], maximize=False, rounds=4, steps=30).reshape(3, S).min(1)[0]
exact1 = -(4 - 2 * grid / (1 + grid))                                       # Proposition 11.14
print("Exp. 1  inf:", [f"{v:.4f}" for v in inf1.tolist()], " exact:", [f"{v:.4f}" for v in exact1.tolist()])
assert (inf1 >= exact1 - 1e-12).all() and (inf1 - exact1).max() < 4e-3

# Exp. 2 on M_3 (+) M_3, Delta(X0) = 0.5, Delta(Y0) = 0.8: sup and inf in one call
sign = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(S)
X2, Y2, V2 = (opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.5, dtype=torch.float64)),
              opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.8, dtype=torch.float64)),
              opt.UnitaryParam(alg3, 2 * S))
def z2():
    x0, y0 = X2.operator(), Y2.operator()
    M = alg3.operator((x0 @ V2.operator() @ y0).matrix, is_self_adjoint=False)
    return sign * (M.trace_norm() - (x0 @ y0).trace.real)
best2 = opt.extremize(z2, [X2, Y2, V2], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]
sup2, inf2 = best2[0].item(), -best2[1].item()
print(f"Exp. 2  sup {sup2:.6f}, inf {inf2:.6f}, exact +-{envelope:.6f}")
assert sup2 <= envelope + 1e-12 and inf2 >= -envelope - 1e-12               # Corollary 11.7: one-sided ...
assert envelope - sup2 < 1e-3 * envelope and inf2 + envelope < 1e-3 * envelope   # ... and reached

# Exp. 3: maximise sign * z at the contrasts `deltas` of |X|; Delta(Y) free in [0, 0.999]
def exp3(alg_, deltas, sign):
    X3, Y3 = opt.SelfAdjointParam(alg_, len(deltas), delta=deltas), opt.PositiveParam(alg_, len(deltas))
    def z3():
        y = Y3.operator()
        M = alg_.operator((y @ X3.operator() @ y).matrix, is_self_adjoint=True)
        return sign * ((y @ X3.abs_operator() @ y).trace.real - M.trace_norm())
    return opt.extremize(z3, [X3, Y3], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]

best3 = exp3(alg22, torch.full((2 * S,), 0.5, dtype=torch.float64), sign)
sup_M22, inf3 = best3[0].item(), -best3[1].item()
s_M2 = exp3(TypeIAlgebra([2], [2], precision='double'),
            torch.tensor([0.0, 0.5], dtype=torch.float64).repeat_interleave(S), 1.0)
direct_sum = (s_M2[1] + s_M2.max()).item()                                 # delta = 0.5, best delta' <= 0.5
conj = [(1 - dl) * (2 * 0.999 / 1.999) ** 2 for dl in (0.0, 0.5)]          # unproven closed form, d = 0.999
print(f"Exp. 3  inf {inf3:.1e};  sup for M_2 at delta = 0, 0.5: {s_M2[0]:.5f}, {s_M2[1]:.5f}"
      f"  (closed form {conj[0]:.5f}, {conj[1]:.5f})")
print(f"        sup for M_2 (+) M_2 at 0.5: search {sup_M22:.5f}, direct sum {direct_sum:.5f}")
assert inf3 >= -1e-12 and inf3 < 1e-6                                       # Proposition 11.9
assert sup_M22 >= direct_sum - 5e-3 * direct_sum                           # (ii), up to 0.5 %
```

La búsqueda en $M_2\oplus M_2$ se queda alrededor de un $0.2\,\%$ por debajo de la suma directa: no ha convergido
del todo, y la mejor cota inferior de $\overline z_{M_2\oplus M_2}(0.5)$ es el número mayor; detectar esto es el
propósito de la comprobación (ii). Los valores para $M_2$ coinciden hasta $10^{-4}$ con
$(1-\delta)\big(2d/(1+d)\big)^2=(1-\delta)(1-\ell(d))^2$, donde $d=0.999$ es el mayor contraste admitido de $Y$:
una forma cerrada que el artículo que acompaña a la biblioteca observó numéricamente (en álgebras reales) y no pudo
demostrar; el código no la afirma. Para bloques mayores, las envolventes superiores de Exp. 1 y Exp. 3 no se
conocen en forma cerrada, hasta donde sabemos; las curvas calculadas son entonces cotas inferiores y conjeturas,
nada más.

## 11.6 Con qué comparar

**Manopt**, **pymanopt** y **geoopt** proporcionan variedades matriciales (variedades de Stiefel, grupos
ortogonales y unitarios) con retracciones como la QR, la polar o la de Cayley, y resolvedores riemannianos
(gradientes conjugados, regiones de confianza, Adam riemanniano) que suelen ser más eficientes que Adam en una carta
recentrada. **cvxpy** con un resolvedor SDP resuelve problemas convexos en variables operador hasta la
optimalidad global, con un certificado dual salvo la tolerancia del resolvedor, y las relajaciones semidefinidas
dan cotas rigurosas en la dirección que una búsqueda local no puede dar. **NCTSSOS**, **ncpol2sdpa** y
**NCSOStools** acotan valores propios y trazas de polinomios no conmutativos en todas las dimensiones a la vez
(sección 0.3). La biblioteca añade parametrizaciones exactas de las restricciones espectrales y de contraste de
$\bigoplus_cM_{k_c}$ (los extremos pueden estar en cualquier sector), unitarias que respetan los bloques, y un único
cálculo por lotes para muchos arranques y valores de la restricción, en CPU o en GPU. Produce candidatos a
extremizadores y cotas inferiores de supremos; el resto deben aportarlo demostraciones o relajaciones.

## Ejercicios

1. (**Transformada de Cayley.**) Para $S^*=-S$, demuestre que $C(S)=(1-S)(1+S)^{-1}$ es unitaria y que
   $S\mapsto C(S)$ es una biyección sobre las unitarias sin valor propio $-1$. Así, a diferencia de la aplicación
   exponencial (lema 11.1), una carta de Cayley $U_0C(S)$ no alcanza $-U_0$.
2. (**Medida de Haar real.**) Para $O$ distribuida según la medida de Haar en $O(k)$ (QR de una matriz de Ginibre
   real) y $A$ real, demuestre que $\mathbb E\operatorname{Tr}(AO)^2=\|A\|_2^2/k$, y compruébelo con un álgebra real.
3. (**Sumas directas.**) Demuestre que el $z$ de Exp. 3 es aditivo respecto de las sumas directas y que
   $\Delta(A\oplus A')=\max(\Delta(A),\Delta(A'))$ para $A,A'\ge0$ con $\lambda_{\max}=1$. Deduzca la desigualdad de
   la suma directa de la sección 11.5 (con $\Delta(Y)$ libre en un intervalo fijo $[0,b]$).

<details>
<summary>Soluciones</summary>

1. $S$ tiene valores propios imaginarios puros, luego $1\pm S$ son invertibles; conmutan, y $(1+S)^*=1-S$. Por
   tanto $C^*C=(1-S)^{-1}(1+S)(1-S)(1+S)^{-1}=1$, y $1+C=2(1+S)^{-1}$ es invertible, luego
   $-1\notin\operatorname{spec}C$. Recíprocamente, para una unitaria $U$ con $1+U$ invertible,
   $S=(1-U)(1+U)^{-1}$ es antiautoadjunta: $S^*=(1+U^*)^{-1}(1-U^*)=(U+1)^{-1}U\,U^*(U-1)=-S$. Como $C$ conmuta
   con $S$, $C(1+S)=1-S$ equivale a $S(1+C)=1-C$, es decir, $S=(1-C)(1+C)^{-1}$: las dos aplicaciones son
   inversas una de la otra.

```python
Sk = torch.randn(4, 4, dtype=torch.complex128)
Sk = Sk - Sk.mH                                                   # skew-adjoint
I4 = torch.eye(4, dtype=torch.complex128)
Cay = (I4 - Sk) @ torch.linalg.inv(I4 + Sk)
assert torch.allclose(Cay.mH @ Cay, I4, atol=1e-12)
assert torch.allclose((I4 - Cay) @ torch.linalg.inv(I4 + Cay), Sk, atol=1e-10)
assert torch.linalg.eigvals(Cay).add(1).abs().min() > 1e-6        # -1 is never an eigenvalue
```

2. El teorema 5.3 y su demostración valen para una matriz gaussiana real (densidad $\propto e^{-\|Z\|_2^2/2}$), el
   grupo ortogonal y $V$ ortogonal (la observación que sigue al teorema 5.3), luego $VO$ tiene la ley de $O$.
   La multiplicación por la izquierda por $\operatorname{diag}(\pm1)$ da $\mathbb E[O_{ij}O_{lm}]=0$ para $i\ne l$,
   las permutaciones muestran que $\mathbb E[O_{ij}O_{im}]$ no depende de $i$, y $\sum_iO_{ij}O_{im}=\delta_{jm}$;
   luego $\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k$ y $\mathbb E\operatorname{Tr}(AO)^2=\sum_{ij}A_{ji}^2/k$.

```python
Mr = TypeIAlgebra([5], [5], complex_valued=False, precision='double')
Ar = Mr.operator(torch.randn(1, 1, 5, 5, dtype=torch.float64))
tr2 = (Ar @ Mr.random_unitary_operator(20_000)).trace ** 2
assert abs(tr2.mean().item() - Ar.frobenius_norm().item() ** 2 / 5) < 5 * tr2.std().item() / math.sqrt(20_000)
```

3. Las trazas y las normas traza de operadores diagonales por bloques son sumas sobre los bloques, y
   $|X\oplus X'|=|X|\oplus|X'|$; $\lambda_{\min}(A\oplus A')=\min(\lambda_{\min}(A),\lambda_{\min}(A'))$, y $\Delta$
   decrece con $\lambda_{\min}$. Así, si $(X,Y)$ es admisible en $\delta$ y $(X',Y')$ en $\delta'\le\delta$, entonces
   $(X\oplus X',Y\oplus Y')$ es admisible en $\delta$, con $z=z(X,Y)+z(X',Y')$; tómense supremos sobre ambos pares y
   sobre $\delta'\le\delta$.

</details>

**Siguiente:** el [capítulo 12](12_equivariant_learning.md) trata el aprendizaje automático equivariante: las capas que conmutan
con una simetría son elementos de un conmutante (capítulo 2), construidos con los mismos operadores por lotes y
diferenciables.
