[English](../11_optimisation.md) · **Español** · [Français](../fr/11_optimisation.md) · [Deutsch](../de/11_optimisation.md) · [中文](../zh/11_optimisation.md) · [日本語](../ja/11_optimisation.md) · [Русский](../ru/11_optimisation.md)

# 11. Optimización y desigualdades de trazas

Una desigualdad de trazas es un enunciado sobre los valores extremos de un funcional de operadores. La
sección 3.7 mostró una que es verosímil y falsa (proposición 3.9); una ejecución de muestreo ([tutorial, sección 3.7](../../tutorial/es/03_traces_and_states.md#s-3-7)) no encontró numéricamente la violación. Este capítulo explica por qué el muestreo falla en dimensión alta y lo sustituye por la optimización sobre conjuntos de
restricciones parametrizados de forma exacta; allí donde los extremos se conocen en forma cerrada, los
demostramos, de modo que el optimizador se contrasta con un teorema.

**Aprenderá** cómo se parametrizan las unitarias y los espectros con un contraste prescrito, y por qué
se alcanza todo operador admisible; la desigualdad de trazas de von Neumann y los extremos exactos que proporciona;
una desigualdad de trazas con su condición de igualdad; por qué las muestras de Haar están, en media cuadrática,
lejos de los extremos; y qué demuestra una búsqueda multiarranque.
Las comprobaciones numéricas de este capítulo, y el optimizador de la biblioteca, están en el
[Tutorial 11](../../tutorial/es/11_optimisation.md).

En todo el capítulo, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ con $D=\sum_ck_c\ge2$, $\mathrm{Tr}$ es la
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
| Exp. 1 | $z=\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)$ | $X,Y\ge0$, $U\in\mathcal M$ unitaria |
| Exp. 2 | $z=\mathrm{Tr}\lvert XY\rvert-\mathrm{Tr}(\lvert X\rvert\,\lvert Y\rvert)$ | $X=UX_0$, $Y=VY_0$, $X_0,Y_0\ge0$, $U,V$ unitarias |
| Exp. 3 | $z=\mathrm{Tr}(Y\lvert X\rvert Y)-\mathrm{Tr}\lvert YXY\rvert$ | $X=X^*$, $Y\ge0$ |

con $|X|$, $|Y|$ normalizados a $\lambda_{\max}=1$ y sus contrastes como parámetros. Las **envolventes** son
$\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}$ y $\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}$, con todas las
demás variables libres.

<a id="s-11-1"></a>

## 11.1 Conjuntos de restricciones y sus parametrizaciones

**Definición.** Para $\delta\in[0,1)$ y $\ell=\ell(\delta)$, la **capa de contraste** es
$$
\mathcal P_\delta=\{X\in\mathcal M:X\ge0,\ \lambda_{\max}(X)=1,\ \Delta(X)=\delta\}
=\{X=X^*:\mathrm{spec}\,X\subseteq[\ell,1],\ 1,\ell\in\mathrm{spec}\,X\}.
$$
Las **casillas** de $\mathcal M$ son los $D$ pares $(c,i)$, $i\le k_c$; para $\lambda\in\mathbb R^D$,
$\mathrm{diag}(\lambda)=\bigoplus_c\mathrm{diag}(\lambda_{c,1},\dots,\lambda_{c,k_c})$.

**La parametrización.** Para un contraste prescrito $\delta$, una búsqueda sortea, en cada arranque, dos casillas distintas
$s_\pm=(c_\pm,i_\pm)$ ($s_+$ uniformemente entre las $D$ casillas y luego $s_-$ uniformemente entre las $D-1$
restantes, de modo que todo par ordenado de casillas distintas tiene probabilidad $1/(D(D-1))$) y una unitaria
aleatoria de Haar $U_0\in\mathcal M$ (lema 11.10), y envía $\theta\in\mathbb R^D$ y
$G=\bigoplus_cG_c$ a
$$
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell,\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
$$
Las unitarias solas se parametrizan mediante $G\mapsto U_0e^{G-G^*}$. A diferencia de una sigmoide, $(1+\sin\theta)/2$
alcanza $0$ y $1$ para $\theta$ finito, lo cual es conveniente porque los espectros extremales de este capítulo
(corolario 11.7, proposición 11.14) están en $\ell$ y $1$. Un contraste libre en $[a,b]$ es un parámetro más,
$\delta=a+(b-a)(1+\sin\eta)/2$, por ejemplo con $[a,b]=[0,0.999]$. Un valor alcanzado con esta restricción también se alcanza sin ella, de modo que sigue siendo una cota
inferior de un supremo sin restringir (una cota superior de un ínfimo sin restringir).

*En la biblioteca:* esta parametrización es `optimize.PositiveParam`, y $G\mapsto U_0e^{G-G^*}$ es `optimize.UnitaryParam` — [tutorial, sección 11.1](../../tutorial/es/11_optimisation.md#s-11-1).

<a id="r-11-1"></a>

**Lema 11.1.** Para toda unitaria $U_0\in\mathcal M$, $G\mapsto U_0e^{G-G^*}$ es sobreyectiva sobre el grupo
unitario de $\mathcal M$.

*Demostración.* Sea $U$ unitaria. Por el teorema 1.3 (bloque a bloque) y el corolario 1.4(2),
$U_0^*U=\bigoplus_cW_c\mathrm{diag}(e^{i\varphi_{c,j}})_jW_c^*$ con $W_c$ unitarias, $\varphi_{c,j}\in\mathbb R$.
Entonces $S=\bigoplus_cW_c\mathrm{diag}(i\varphi_{c,j})_jW_c^*$ satisface $S^*=-S$ y $e^S=U_0^*U$ (la serie
$\sum_nS^n/n!$ se calcula término a término en la base propia), y $G=S/2$ da $G-G^*=S$. $\square$

<a id="r-11-2"></a>

**Proposición 11.2 (la parametrización es sobreyectiva).** Para $\delta$, $s_\pm$ y $U_0$ fijos, la imagen de
$(\theta,G)\mapsto X(\theta,G)$ es $\{X\in\mathcal P_\delta:1\in\mathrm{spec}\,X_{c_+},\ \ell\in\mathrm{spec}\,X_{c_-}\}$,
y $\mathcal P_\delta$ es la unión de estas imágenes sobre todos los pares de casillas distintas.

*Demostración.* La imagen está contenida en este conjunto: $X$ tiene valores propios $\lambda_s\in[\ell,1]$, con
$\lambda_{s_\pm}$ en los bloques $c_\pm$. Recíprocamente, diagonalícese $X_c=W_c\mathrm{diag}(\mu_c)W_c^*$ y
permútese cada $\mu_c$ (sustituyendo $W_c$ por $W_c\Pi_c$, con $\Pi_c$ una matriz de permutación) de modo que
$\mu_{s_+}=1$ y $\mu_{s_-}=\ell$. Esto es posible: si $c_+=c_-$ y $\delta>0$, entonces $1\ne\ell$ son dos valores
propios del mismo bloque y ocupan posiciones distintas; si $\delta=0$, todos los valores propios valen $1$. Si
$\ell<1$, póngase $\theta_s=\arcsin\big(2(\mu_s-\ell)/(1-\ell)-1\big)$ para las demás casillas, de modo que
$\lambda=\mu$ (si $\ell=1$, cualquier $\theta$ sirve). El lema 11.1 da $G$ con $U_0e^{G-G^*}=\bigoplus_cW_c$.
*Unión:* para $X\in\mathcal P_\delta$, elíjanse bloques $c_\pm$ que contengan los valores propios $1$, $\ell$ y
casillas distintas en ellos; si $c_+=c_-$ y $\delta>0$, el bloque tiene dos valores propios diferentes, luego
$k_{c_+}\ge2$; si $\delta=0$, $X=1$ y sirven dos casillas cualesquiera. $\square$

*Cálculo:* la proposición 11.2 de forma constructiva: los parámetros de la demostración reproducen un operador objetivo — [tutorial, sección 11.1](../../tutorial/es/11_optimisation.md#s-11-1).

Como todo par de casillas se sortea con probabilidad positiva, todo punto de la capa está, con probabilidad
positiva, en la imagen de un arranque, y (en aritmética exacta) la restricción se cumple en cada iterado.
Para $X$ autoadjunto (Exp. 3) se usa $X=U\mathrm{diag}(\sigma\lambda)U^*$ con signos $\sigma\in\{\pm1\}^D$ sorteados en cada
arranque (en cada bloque, el número de signos menos es uniforme en $\{0,\dots,k_c\}$ y sus posiciones son
uniformemente aleatorias, de modo que todo patrón de signos tiene probabilidad positiva) y modificados solo por los
cambios de signo de la sección 11.5. Entonces $|X|=U\mathrm{diag}(\lambda)U^*$ (teorema 1.5(e), ya que
$\lambda>0$), y la unión sobre signos y casillas es $\{X=X^*:|X|\in\mathcal P_\delta\}$: para un tal $X$,
diagonalícese $X_c=W_c\mathrm{diag}(\nu_c)W_c^*$, de modo que $|X|_c=W_c\mathrm{diag}(|\nu_c|)W_c^*$ con
$|\nu_{c,i}|\ge\ell>0$, y aplíquese la demostración de la proposición 11.2 a $\mu=|\nu|$ con $\sigma_s=\nu_s/|\nu_s|$.
Solo importa $S=G-G^*$.

*En la biblioteca:* `optimize.SelfAdjointParam` — [tutorial, sección 11.1](../../tutorial/es/11_optimisation.md#s-11-1).

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
inyectiva: $e^{2\pi iP}=1=e^0$ para toda proyección $P\in\mathcal M$. Por eso una búsqueda **recentra** la carta tras cada ronda ($U_0\leftarrow U_0e^{G-G^*}$, $G\leftarrow0$). Otras dos cartas tienen la misma diferencial $S\mapsto U_0S$ en
$0$. La carta de Cayley $S\mapsto U_0C(-S/2)$ (ejercicio 1):
$(1+tS/2)(1-tS/2)^{-1}=(1+tS/2)(1+tS/2+O(t^2))=1+tS+O(t^2)$. La retracción polar, $S\mapsto$ el factor unitario de
$A=U_0(1+S)$ (teorema 1.8): $A^*A=(1-S)(1+S)=1+S^*S\ge1$ es invertible, luego $|A|$ es invertible y el factor
unitario es $A|A|^{-1}$; con $tS$ en lugar de $S$, $|A|^{-1}=(1+t^2S^*S)^{-1/2}=1+O(t^2)$ (en una base propia de
$S^*S$), luego el factor es $U_0(1+tS)+O(t^2)$.

> **Significado físico** (interpretación física, no un teorema de este curso). El contraste de Michelson es la
> visibilidad $(I_{\max}-I_{\min})/(I_{\max}+I_{\min})$ de las franjas de interferencia. Para $0\le X\le1$ con
> $\lambda_{\max}=1$ —un filtro, un atenuador, la eficiencia de un detector— fija el rango dinámico
> $\lambda_{\min}/\lambda_{\max}$, el inverso de un número de condición.

<a id="s-11-2"></a>

## 11.2 La desigualdad de trazas de von Neumann

Una matriz $(c_{ij})$ con $c_{ij}\ge0$ y todas las sumas por filas y por columnas iguales a $1$ es **doblemente
estocástica**; para una unitaria $G$, $(|G_{ij}|^2)$ y $(|G_{ji}|^2)$ son doblemente estocásticas, pues las filas y
las columnas de $G$ son vectores unitarios.

<a id="r-11-3"></a>

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
$\Sigma=\mathrm{diag}(s_i(A))$: escríbase $A=W'|A|$ (teorema 1.8), $|A|=V\Sigma V^*$ (teorema 1.3) y $W=W'V$.

<a id="r-11-4"></a>

**Teorema 11.4 (desigualdad de trazas de von Neumann).** Para $A,B\in M_k(\mathbb C)$,
$|\mathrm{Tr}(AB)|\le\sum_is_i(A)\,s_i(B)$.

*Demostración.* Sean $A=W_A\Sigma_AV_A^*$, $B=W_B\Sigma_BV_B^*$, $a_i=s_i(A)$, $b_j=s_j(B)$, y $G=V_A^*W_B$,
$H=V_B^*W_A$. Por ciclicidad, $\mathrm{Tr}(AB)=\mathrm{Tr}(\Sigma_AG\Sigma_BH)=\sum_{i,j}a_iG_{ij}b_jH_{ji}$,
luego, por $|g||h|\le\frac12(|g|^2+|h|^2)$ y $a_ib_j\ge0$,
$$
|\mathrm{Tr}(AB)|\le\sum_{i,j}a_ib_j\,d_{ij},\qquad d_{ij}=\tfrac12\big(|G_{ij}|^2+|H_{ji}|^2\big).
$$
$(d_{ij})$ es doblemente estocástica, y el lema 11.3 da $\sum_{ij}a_ib_jd_{ij}\le\sum_ia_ib_i$. $\square$

<a id="r-11-5"></a>

**Corolario 11.5.** (a) Para $P,Q\in M_k(\mathbb C)$ autoadjuntas con valores propios $p,q$:
$\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i$.
(b) Para $P,Q\ge0$: $\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i$.

*Demostración.* (a) Como en la demostración del teorema 3.5, $\mathrm{Tr}(PQ)=\sum_{ij}p_iq_jc_{ij}$ con la matriz
doblemente estocástica $c_{ij}=|\langle u_i,v_j\rangle|^2$ para bases propias ordenadas de forma decreciente;
aplíquese el lema 11.3 (su cota inferior es $\sum_ip^\downarrow_iq^\downarrow_{k+1-i}=\sum_ip^\downarrow_iq^\uparrow_i$).
(b) La primera desigualdad es el teorema 1.9. Por el teorema 1.9 y el teorema 11.4,
$\|PQ\|_1=\max_U|\mathrm{Tr}(P\,QU)|\le\sum_is_i(P)s_i(QU)$, y $s_i(QU)=q^\downarrow_i$ porque
$|QU|=(U^*Q^2U)^{1/2}=U^*QU$ (teorema 1.6). $\square$

*Cálculo:* el teorema 11.4 y el corolario 11.5 en matrices aleatorias — [tutorial, sección 11.2](../../tutorial/es/11_optimisation.md#s-11-2).

Así, para $X,Y\ge0$ en $M_k(\mathbb C)$ con espectros $x,y$, la violación $\|YX\|_1-\mathrm{Tr}(XY)$ de la
proposición 3.9 es a lo sumo $\sum_ix^\downarrow_i(y^\downarrow_i-y^\uparrow_i)$ (corolario 11.5(b) para $(Y,X)$ y la
cota inferior de (a)). Para Exp. 2, este número es la respuesta exacta:

<a id="r-11-6"></a>

**Proposición 11.6 (Exp. 2 con espectros fijos).** En Exp. 2, $|X|=X_0$, $|Y|=Y_0$ y
$z=\|X_0VY_0\|_1-\mathrm{Tr}(X_0Y_0)$. Si los espectros $x_c$ de $(X_0)_c$ e $y_c$ de $(Y_0)_c$ están fijos y
las bases propias y $V$ varían, entonces
$$
\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}\big(y^\downarrow_{c,i}-y^\uparrow_{c,i}\big),
$$
y ambos se alcanzan.

*Demostración.* $X^*X=X_0U^*UX_0=X_0^2$, luego $|X|=X_0$ (teorema 1.6); análogamente $|Y|=Y_0$, y $|UA|=|A|$ da
$\mathrm{Tr}|XY|=\|X_0VY_0\|_1$. Todas las variables actúan bloque a bloque y $z$ es una suma sobre los bloques,
así que podemos suponer $C=1$. Póngase $Y'=VY_0V^*\ge0$ (espectro $y$); como $|AV|=V^*|A|V$,
$\|X_0VY_0\|_1=\|X_0Y'V\|_1=\|X_0Y'\|_1$. El corolario 11.5 para los pares $(X_0,Y')$ y $(X_0,Y_0)$, junto con
$\|X_0Y'\|_1\ge\mathrm{Tr}(X_0Y')$, da
$$
\sum_ix^\downarrow_iy^\uparrow_i-\sum_ix^\downarrow_iy^\downarrow_i\le\|X_0Y'\|_1-\mathrm{Tr}(X_0Y_0)\le\sum_ix^\downarrow_iy^\downarrow_i-\sum_ix^\downarrow_iy^\uparrow_i .
$$
Sea $J$ la permutación de inversión, $Je_i=e_{k+1-i}$; entonces $J\mathrm{diag}(w)=\mathrm{diag}(w_k,\dots,w_1)J$,
y la norma traza de $\mathrm{diag}(w)J$ es $\sum_i|w_i|$. Con $X_0=\mathrm{diag}(x^\downarrow)$ y $V=J$:
para $Y_0=\mathrm{diag}(y^\uparrow)$, $X_0VY_0=\mathrm{diag}(x^\downarrow_iy^\downarrow_i)J$ y
$\mathrm{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\uparrow_i$ (la cota superior); para $Y_0=\mathrm{diag}(y^\downarrow)$,
$X_0VY_0=\mathrm{diag}(x^\downarrow_iy^\uparrow_i)J$ y $\mathrm{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\downarrow_i$
(la cota inferior). $\square$

<a id="r-11-7"></a>

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

*Cálculo:* las configuraciones extremales de la proposición 11.6 y la envolvente del corolario 11.7 para $k=3$, $C=2$ — [tutorial, sección 11.2](../../tutorial/es/11_optimisation.md#s-11-2); la envolvente hallada por una búsqueda multiarranque — [tutorial, sección 11.5](../../tutorial/es/11_optimisation.md#s-11-5).

<a id="s-11-3"></a>

## 11.3 Una desigualdad de trazas con condición de igualdad

Para $X=X^*$, sea $X_\pm=f_\pm(X)$ con $f_+(t)=\max(t,0)$, $f_-(t)=\max(-t,0)$. Entonces $X_\pm\ge0$, $X=X_+-X_-$,
$|X|=X_++X_-$ y $X_+X_-=0$ (teorema 1.5(c), y 1.5(e) para $|X|=(X^2)^{1/2}=h(X)$ con
$h(t)=\sqrt{t^2}=|t|=f_+(t)+f_-(t)$; corolario 1.4(3) para la positividad).

<a id="r-11-8"></a>

**Lema 11.8.** Para $A,B\ge0$ en $M_n(\mathbb C)$, $\|A-B\|_1\le\mathrm{Tr}\,A+\mathrm{Tr}\,B$, con igualdad si y
solo si $AB=0$.

*Demostración.* Sean $H=A-B$ y $Q=g(H)$ con $g=1$ en $(0,\infty)$ y $g=0$ en el resto; $Q$ es una proyección, ya que
$g=\bar g=g^2$ (teorema 1.5(c)). Como $|t|=t(2g(t)-1)$, se tiene $|H|=H(2Q-1)$ (teorema 1.5(c),(e)), y, usando
$\mathrm{Tr}(AR)=\mathrm{Tr}(RAR)$ para las proyecciones $R=Q,1-Q$ (ciclicidad),
$$
\|H\|_1=\mathrm{Tr}\big(A(2Q-1)\big)-\mathrm{Tr}\big(B(2Q-1)\big)
=\mathrm{Tr}\,A+\mathrm{Tr}\,B-2\mathrm{Tr}\big((1-Q)A(1-Q)\big)-2\mathrm{Tr}(QBQ).
$$
Las trazas restadas son $\ge0$, y $\mathrm{Tr}(RAR)=\|A^{1/2}R\|_2^2$ se anula si y solo si $AR=0$; así, la
igualdad se da si y solo si $A(1-Q)=0$ y $BQ=0$, y entonces $AB=AQB=A(BQ)^*=0$. Recíprocamente, si $AB=0$,
entonces $BA=0$, $(A-B)^2=A^2+B^2=(A+B)^2$, y el teorema 1.6 da $|A-B|=A+B$. $\square$

<a id="r-11-9"></a>

**Proposición 11.9 (Exp. 3).** Para $X=X^*$ e $Y\ge0$ en $\mathcal M$,
$\mathrm{Tr}|YXY|\le\mathrm{Tr}(Y|X|Y)$, con igualdad si y solo si $YX_+Y^2X_-Y=0$, en particular si
$XY=YX$.

*Demostración.* En $M_D(\mathbb C)$, $A=YX_+Y$ y $B=YX_-Y$ son positivos (corolario 1.7 con $T=X_\pm^{1/2}Y$),
$YXY=A-B$ y $\mathrm{Tr}(Y|X|Y)=\mathrm{Tr}\,A+\mathrm{Tr}\,B$; aplíquese el lema 11.8. Si $XY=YX$,
entonces $Y$ conmuta con $X_\pm$ (teorema 1.5(d)), y $AB=Y^4X_+X_-=0$. $\square$

Así, $z\ge0$ en Exp. 3, y $\underline z=0$ (tómese $Y=1$). La condición de igualdad es más débil que la
conmutación: para $X=\mathrm{diag}(1,-1,0)$ e $Y$ con filas $(1,0,a)$, $(0,1,0)$, $(a,0,1)$, $|a|<1$, se tiene
$AB=(Ye_1)\langle Ye_1,Ye_2\rangle(Ye_2)^*=0$, mientras que $XY\ne YX$ para $a\ne0$.

*Cálculo:* la proposición 11.9 en pares aleatorios, y el ejemplo de igualdad sin conmutación — [tutorial, sección 11.3](../../tutorial/es/11_optimisation.md#s-11-3).

<a id="s-11-4"></a>

## 11.4 Muestreo frente a optimización

¿Qué tamaño tiene típicamente $|\mathrm{Tr}(AU)|$, comparado con su máximo $\|A\|_1$? Muestreamos $U$
como el factor $Q$ de la descomposición QR $Z=QR$, $R_{ii}>0$, de una **matriz de Ginibre** $Z$ (entradas
independientes con densidad $\pi^{-1}e^{-|z|^2}$; la receta de Mezzadri), de forma independiente en cada bloque.

<a id="r-11-10"></a>

**Lema 11.10.** $Q$ está distribuida según la medida de Haar en $U(k)$; en particular, para toda unitaria fija
$V\in M_k(\mathbb C)$, $VQ$ tiene la ley de $Q$. La unitaria $\bigoplus_cQ_c$ con bloques independientes está
distribuida según la medida de Haar en el grupo unitario de $\mathcal M$.

*Demostración.* Es el teorema 5.3 (la invariancia por la izquierda es el paso de equivariancia de su demostración;
la ley se identifica con la medida de Haar mediante el lema de unicidad 5.1, cuyo paso de teoría de la medida es el
corolario A.2) y el corolario 5.4. $\square$

*En la biblioteca:* este muestreador es `TypeIAlgebra.random_unitary_operator` — [tutorial, sección 11.4](../../tutorial/es/11_optimisation.md#s-11-4).

<a id="r-11-11"></a>

**Proposición 11.11.** Para $U\in U(k)$ de Haar y $A\in M_k(\mathbb C)$: $\mathbb E\mathrm{Tr}(AU)=0$ y
$\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k$.

*Demostración.* $\mathbb EU_{ij}=0$ por el corolario 5.5. La proposición 5.7 (con $C=1$) para
$X=|e_j\rangle\langle e_m|$ da la fórmula de Weingarten de orden uno,
$\mathbb E[U_{ij}\bar U_{lm}]=\mathbb E(UXU^*)_{il}=\delta_{il}\delta_{jm}/k$; para $k\ge2$ también se deduce de la
fórmula de orden dos, el teorema 5.10, haciendo $i_2=i'_2$, $j_2=j'_2$ y sumando sobre ambos, ya que
$\sum_{i_2,j_2}|U_{i_2j_2}|^2=k$. Con $\mathrm{Tr}(AU)=\sum_{ij}A_{ji}U_{ij}$,
$\mathbb E|\mathrm{Tr}(AU)|^2=\sum_{ijlm}A_{ji}\bar A_{ml}\mathbb E[U_{ij}\bar U_{lm}]=\frac1k\sum_{ij}|A_{ji}|^2$. $\square$

<a id="r-11-12"></a>

**Corolario 11.12.** Sean $X,Y\ge0$ en $M_k(\mathbb C)$ con espectros en $[a,b]$, $a>0$, y $U$ de Haar. Entonces
$$
\frac{\big(\mathbb E|\mathrm{Tr}(XUY)|^2\big)^{1/2}}{\max_V|\mathrm{Tr}(XVY)|}=\frac{\|YX\|_2}{\sqrt k\,\|YX\|_1}\le\min\Big(\frac1{\sqrt k},\frac{b^2}{a^2k}\Big),
\qquad
\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le\frac{b^4}{a^4k^2}.
$$

*Demostración.* La igualdad es la proposición 11.11 para $A=YX$ (pues
$\mathrm{Tr}(XUY)=\mathrm{Tr}(AU)$) junto con la proposición 3.9. Como
$\|A\|_2^2=\sum_is_i^2\le(\sum_is_i)^2=\|A\|_1^2$ (sección 1.4), el cociente es $\le1/\sqrt k$. Como
$\|Av\|=\||A|v\|$ (sección 1.4), diagonalizando $|A|$ se obtiene $s_k(A)=\min_{\|v\|=1}\|Av\|$ y
$s_1(A)=\max_{\|v\|=1}\|Av\|$; y $\|Xv\|^2=\langle v,X^2v\rangle\in[a^2,b^2]$ para $v$ unitario (corolario 1.4(3)
para $X^2-a^2$ y $b^2-X^2$), y lo mismo para $Y$. Por tanto $a^2\le s_k(A)\le s_1(A)\le b^2$, $\|A\|_2^2\le kb^4$ y
$\|A\|_1\ge ka^2$. Por último,
$\mathrm{Tr}(XY)=\mathrm{Tr}(X^{1/2}YX^{1/2})\ge a\mathrm{Tr}\,X\ge ka^2$, y la desigualdad de Markov para
$|\mathrm{Tr}(XUY)|^2$ acota la probabilidad por $(kb^4/k)/(ka^2)^2$. $\square$

Así, $N$ muestras de Haar encuentran una violación de $|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)$ con probabilidad
a lo sumo $Nb^4/(a^4k^2)$ (una unión de $N$ sucesos), aunque existe una para todo par que no conmuta
(proposición 3.9). La cota solo es informativa para $k>\sqrt N\,b^2/a^2$ (para $b/a=6$ y $N=2000$: $k>1609$).

*Cálculo:* la proposición 11.11 por Monte Carlo, y el cociente del corolario 11.12 frente a la mejor de $2000$ muestras para $k\le32$ — [tutorial, sección 11.4](../../tutorial/es/11_optimisation.md#s-11-4).

<a id="r-11-13"></a>

**Lema 11.13.** Si $k\ge2$, entonces $\min_U|\mathrm{Tr}(AU)|=0$ para toda $A\in M_k(\mathbb C)$.

*Demostración.* Sean $A=W\Sigma V^*$ y $P$ el desplazamiento cíclico, $Pe_i=e_{(i\bmod k)+1}$, cuya diagonal se anula
para $k\ge2$. Para $U=VPW^*$, $\mathrm{Tr}(AU)=\mathrm{Tr}(\Sigma P)=0$. $\square$

Así, $\inf_Uz=-\mathrm{Tr}(XY)$ y $\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)$ (proposición 3.9). Por la proposición 11.11, la media cuadrática de
$|\mathrm{Tr}(XUY)|$ es $\|YX\|_2/\sqrt k$; con las cotas de la demostración del corolario 11.12 (espectros en
$[a,b]$), $\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4$ y $\mathrm{Tr}(XY)\ge ka^2$, de modo que la desigualdad de
Markov da $\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)$ para
$\varepsilon>0$: para $k$ grande, $z$ está cerca de su ínfimo con probabilidad alta.

*Cálculo:* ambos extremos sobre $U$ por ascenso y descenso de gradiente, frente a $20\,000$ muestras de Haar — [tutorial, sección 11.4](../../tutorial/es/11_optimisation.md#s-11-4).

> **Significado físico** (interpretación física, no un teorema de este curso). El corolario 11.12 es un ejemplo
> débil y demostrable de lo que la mecánica estadística llama **tipicidad**: en dimensión alta se espera que casi
> todo estado o toda evolución dé aproximadamente el mismo valor de una magnitud macroscópica. Una desigualdad
> óptima se decide en configuraciones finamente ajustadas —bases propias alineadas, espectros de dos niveles— que
> una unitaria aleatoria (un circuito aleatorio, una matriz de dispersión aleatoria) rara vez produce.

<a id="s-11-5"></a>

## 11.5 Búsqueda multiarranque de envolventes

Una **búsqueda multiarranque** ejecuta un método de gradiente, por rondas, desde muchos arranques independientes,
cada uno con sus propias casillas, signos, $U_0$ y contraste; recentra las cartas unitarias tras cada ronda y devuelve,
para cada arranque, el mejor valor del objetivo en los iterados. Los signos de la parametrización autoadjunta son
discretos; entre rondas se modifican mediante cambios de signo individuales en casillas aleatorias, y cada uno se
conserva si mejora el objetivo.

*En la biblioteca:* `optimize.extremize` — [tutorial, sección 11.5](../../tutorial/es/11_optimisation.md#s-11-5).

En aritmética exacta todo iterado es admisible (la proposición 11.2 y las observaciones que la siguen), así que todo
valor devuelto es el valor de $z$ en un punto admisible; en coma flotante es ese valor salvo errores de redondeo. Por tanto, una maximización da una **cota inferior** del supremo, una
minimización una cota superior del ínfimo, y nada más: una búsqueda local puede detenerse en un extremo local. La
otra dirección requiere una demostración, como el corolario 11.7 y la proposición 11.9, o una relajación convexa
(sección 11.6). Una envolvente calculada se contrasta (i) con formas cerradas; (ii) con configuraciones construidas a
mano, como sumas directas de soluciones más pequeñas, por debajo de las cuales no debe quedar; (iii) con muestras
estructuradas (espectros de dos niveles, bases de permutación), cuyos valores son otros valores alcanzados.

<a id="r-11-14"></a>

**Proposición 11.14 (envolvente inferior de Exp. 1).** Supongamos que todos los $k_c\ge2$. Sobre
$X\in\mathcal P_\delta$, $Y\ge0$ con $\lambda_{\max}(Y)=1$ y unitarias $U\in\mathcal M$,
$\ \min z=-\big(D-2\delta/(1+\delta)\big)$.

*Demostración.* $Y\le1$ da $X^{1/2}YX^{1/2}\le X$, luego $\mathrm{Tr}(XY)\le\mathrm{Tr}\,X\le(D-1)+\ell(\delta)$,
pues un valor propio de $X$ es $\ell(\delta)$ y los demás son $\le1$; por tanto
$z\ge-\mathrm{Tr}(XY)\ge-(D-2\delta/(1+\delta))$. La igualdad se da para $X$ con el valor propio $\ell(\delta)$
una vez y $1$ en el resto, $Y=1$, y $U=\bigoplus_cU_c$ con $\mathrm{Tr}(X_cU_c)=0$ (lema 11.13). $\square$

Las comprobaciones del tipo (ii) pueden usar sumas directas: el $z$ de Exp. 3 es aditivo respecto de las sumas directas, y el contraste de una suma directa es el mayor contraste
de sus sumandos (ejercicio 3), luego, para las envolventes de Exp. 3 con $\Delta(Y)$ libre en un intervalo fijo $[0,b]$,
$$
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
$$
en particular, la suma de dos valores alcanzados en $M_2$ con contrastes $\delta$ y $\delta'\le\delta$ es un valor
alcanzado en $M_2\oplus M_2$ con contraste $\delta$.

*Cálculo:* la proposición 11.14 en $M_2\oplus M_2$, el corolario 11.7 en $M_3\oplus M_3$ y la proposición 11.9 mediante búsqueda multiarranque, con una comprobación del tipo (ii) — [tutorial, sección 11.5](../../tutorial/es/11_optimisation.md#s-11-5).

<a id="s-11-6"></a>

## 11.6 Con qué comparar

*En la biblioteca:* la comparación con otro software, incluidos resolvedores SDP, cuyas relajaciones convexas acotan un supremo por arriba — [tutorial, sección 11.6](../../tutorial/es/11_optimisation.md#s-11-6).

<a id="x-11"></a>

## Ejercicios

1. (**Transformada de Cayley.**) Para $S^*=-S$, demuestre que $C(S)=(1-S)(1+S)^{-1}$ es unitaria y que
   $S\mapsto C(S)$ es una biyección sobre las unitarias sin valor propio $-1$. Así, a diferencia de la aplicación
   exponencial (lema 11.1), una carta de Cayley $U_0C(S)$ no alcanza $-U_0$.
2. (**Medida de Haar real.**) Para $O$ distribuida según la medida de Haar en $O(k)$ (QR de una matriz de Ginibre
   real) y $A$ real, demuestre que $\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k$, y compruébelo con un álgebra real.
3. (**Sumas directas.**) Demuestre que el $z$ de Exp. 3 es aditivo respecto de las sumas directas y que
   $\Delta(A\oplus A')=\max(\Delta(A),\Delta(A'))$ para $A,A'\ge0$ con $\lambda_{\max}=1$. Deduzca la desigualdad de
   la suma directa de la sección 11.5 (con $\Delta(Y)$ libre en un intervalo fijo $[0,b]$).

<details>
<summary>Soluciones</summary>

1. $S$ es normal e $iS$ es autoadjunta, luego $S$ tiene valores propios imaginarios puros (corolario 1.4(2)) y
   $1\pm S$ son invertibles; conmutan, y $(1+S)^*=1-S$. Por tanto $C^*C=(1-S)^{-1}(1+S)(1-S)(1+S)^{-1}=1$, y
   $1+C=2(1+S)^{-1}$ es invertible, luego $-1\notin\mathrm{spec}\,C$. Recíprocamente, para una unitaria $U$ con
   $1+U$ invertible, $S=(1-U)(1+U)^{-1}$ es antiautoadjunta:
   $S^*=(1+U^*)^{-1}(1-U^*)=(U+1)^{-1}U\,U^*(U-1)=(U+1)^{-1}(U-1)=-S$, ya que $U\pm1$ conmutan. Para $C$ y $S$ que
   conmutan, $C(1+S)=1-S$ equivale a $S(1+C)=1-C$ (ambas dicen $C+CS+S=1$). Ahora bien, $C(S)$ conmuta con $S$, y
   $S=(1-U)(1+U)^{-1}$ conmuta con $U$; así, $U=C(S)$ si y solo si $S=(1-U)(1+U)^{-1}$: las dos aplicaciones son
   inversas una de la otra.

2. Por el teorema B.17 (también el lema C.2), el factor QR con signos corregidos $O$ de una matriz gaussiana real,
   que la biblioteca muestrea para un álgebra real, está distribuido según la medida de Haar en $O(k)$, luego $VO$
   tiene la ley de $O$ para $V\in O(k)$.
   La multiplicación por la izquierda por $\mathrm{diag}(\pm1)$ da $\mathbb E[O_{ij}O_{lm}]=0$ para $i\ne l$,
   las permutaciones muestran que $\mathbb E[O_{ij}O_{im}]$ no depende de $i$, y $\sum_iO_{ij}O_{im}=\delta_{jm}$;
   luego $\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k$ y $\mathbb E\mathrm{Tr}(AO)^2=\sum_{ij}A_{ji}^2/k$.

3. Las trazas y las normas traza de operadores diagonales por bloques son sumas sobre los bloques, y
   $|X\oplus X'|=|X|\oplus|X'|$; $\lambda_{\min}(A\oplus A')=\min(\lambda_{\min}(A),\lambda_{\min}(A'))$, y $\Delta$
   decrece con $\lambda_{\min}$. Así, si $(X,Y)$ es admisible en $\delta$ y $(X',Y')$ en $\delta'\le\delta$, entonces
   $(X\oplus X',Y\oplus Y')$ es admisible en $\delta$, con $z=z(X,Y)+z(X',Y')$; tómense supremos sobre ambos pares y
   sobre $\delta'\le\delta$.

</details>

*Cálculo:* los ejercicios 1 y 2, numéricamente — [tutorial, Ejercicios](../../tutorial/es/11_optimisation.md#x-11).

**Siguiente:** el [capítulo 12](12_equivariant_learning.md) trata el aprendizaje automático equivariante: las capas que conmutan
con una simetría son elementos de un conmutante (capítulo 2), construidos con los mismos operadores por lotes y
diferenciables.
