[English](../10_quantum_computing.md) · **Español** · [Français](../fr/10_quantum_computing.md) · [Deutsch](../de/10_quantum_computing.md) · [中文](../zh/10_quantum_computing.md) · [日本語](../ja/10_quantum_computing.md) · [Русский](../ru/10_quantum_computing.md)

# 10. Computación cuántica

El ruido en un ordenador cuántico es un canal, y un canal genera un álgebra: la $*$-álgebra de sus
operadores de Kraus. Su estructura (capítulo 2) dice qué información sobrevive intacta al ruido;
una condición sobre los productos de operadores de Kraus, en la que interviene otra álgebra, dice
qué información puede restaurar una recuperación. La última parte trata algoritmos que conservan
una carga: sus estados viven en el álgebra de los sectores de carga, y el amortiguamiento de
amplitud los traslada de un sector a otro.

**Aprenderá** por qué un canal unital y su dual tienen los mismos puntos fijos (el conmutante del
ruido), el teorema de Knill–Laflamme con una demostración completa, cómo proteger solo la etiqueta
de un sector, y qué aspecto tienen el mezclador XY del ansatz de operadores alternantes cuántico, el
amortiguamiento de amplitud y la postselección en el álgebra $\bigoplus_wM_{\binom nw}$ de los
sectores de peso de Hamming.
Las comprobaciones numéricas de este capítulo están en el
[Tutorial 10](../../tutorial/es/10_quantum_computing.md).

En todo el capítulo, un **canal** sobre $M_D(\mathbb C)$ es una aplicación $\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*$
con **operadores de Kraus** que cumplen $\sum_iE_i^*E_i=1$ (preserva la traza); es **unital** si
además $\sum_iE_iE_i^*=1$. Su **dual** $\Phi^*(T)=\sum_iE_i^*TE_i$ cumple
$\mathrm{Tr}(\Phi(\rho)T)=\mathrm{Tr}(\rho\,\Phi^*(T))$ (sección 6.2). (Por el teorema 6.2,
toda aplicación completamente positiva tiene esta forma; aquí los canales vienen *dados* por
operadores de Kraus.) El **álgebra de interacción** $\mathcal A$ de $\Phi$ es la $*$-álgebra
generada por $1$ y las $E_i$, el espacio generado por todos los productos de las $E_i$ y las $E_i^*$;
una matriz conmuta con $\mathcal A$ si y solo si conmuta con todas las $E_i$ y $E_i^*$. Un estado
$\rho$ tiene **soporte** en un subespacio $K$ con proyección $P_K$ si $\rho=P_K\rho P_K$. Sobre $n$
cúbits, $D=2^n$; $X,Y,Z$ son las matrices de Pauli, $Z=\mathrm{diag}(1,-1)$, $X_j$ es $X$ sobre el cúbit $j$ (numerados $1,\dots,n$;
la sección 10.4 sigue a la biblioteca, sitios $0,\dots,n-1$), $|b\rangle$, $b\in\{0,1\}^n$, la base
computacional, $|b|$ el **peso de Hamming**, y $P_w$ la proyección sobre $H_w=\mathrm{span}\{|b\rangle:|b|=w\}$.

<a id="s-10-1"></a>

## 10.1 Álgebras de ruido, subsistemas sin ruido y puntos fijos

Por el teorema 2.6 existe una unitaria $V$ con $V^*\mathcal AV=\bigoplus_cM_{k_c}\otimes1_{m_c}$ y
$V^*\mathcal A'V=\bigoplus_c1_{k_c}\otimes M_{m_c}$: cada operador de Kraus actúa solo sobre el primer
factor tensorial de cada sumando.

<a id="r-10-1"></a>

**Proposición 10.1 (subsistemas sin ruido; subespacios libres de decoherencia).** (a) Escribamos
$V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}$, y sea $\rho=V(\sigma\otimes\tau)V^*$, situado en el sumando
$c$, para estados $\sigma$ sobre $\mathbb C^{k_c}$ y $\tau$ sobre $\mathbb C^{m_c}$. Entonces
$\Phi(\rho)=V(\Phi_c(\sigma)\otimes\tau)V^*$ con el canal $\Phi_c(\sigma)=\sum_iE_{i,c}\sigma E_{i,c}^*$ sobre
$M_{k_c}$: el estado $\tau$ queda intacto. Además, $\Phi^*(T)=T$ para todo $T\in\mathcal A'$.
(b) Si $K$ es un subespacio con $E_iP_K=e_iP_K$, $e_i\in\mathbb C$, para todo $i$, entonces $\Phi(\rho)=\rho$
para todo estado $\rho$ con soporte en $K$.

*Demostración.* (a) $\Phi(\rho)=V\sum_i(E_{i,c}\sigma E_{i,c}^*\otimes\tau)V^*$, y $\Phi_c$ preserva la
traza porque el sumando $c$ de $V^*(\sum_iE_i^*E_i)V=1$ es $(\sum_iE_{i,c}^*E_{i,c})\otimes1_{m_c}$. Para
$T\in\mathcal A'$, $\Phi^*(T)=T\sum_iE_i^*E_i=T$. (b) $P_KE_i^*=(E_iP_K)^*=\bar e_iP_K$, luego
$E_i\rho E_i^*=E_iP_K\rho P_KE_i^*=|e_i|^2\rho$, y $\sum_i|e_i|^2P_K=P_K\sum_iE_i^*E_iP_K=P_K$. $\square$

Un sumando con $k_c=1$ es el caso (b): un **subespacio libre de decoherencia** (Zanardi–Rasetti 1997,
Lidar–Chuang–Whaley 1998). Para $k_c>1$, el factor $\mathbb C^{m_c}$ es un **subsistema sin ruido**
(Knill–Laflamme–Viola 2000): el ruido actúa, pero solo sobre el otro factor. La proposición 10.1 vale
para todo canal; para los canales unitales, los observables preservados son *exactamente* los de
$\mathcal A'$:

<a id="r-10-2"></a>

**Teorema 10.2 (puntos fijos de los canales unitales).** Sea $\Phi$ unital. Entonces
$$
\mathrm{Fix}(\Phi)=\{T\in M_D(\mathbb C):\Phi(T)=T\}=\mathcal A'=\mathrm{Fix}(\Phi^*).
$$
En particular, $\mathrm{Fix}(\Phi)$ es un álgebra de von Neumann, $*$-isomorfa a $\bigoplus_cM_{m_c}$.

*Demostración.* La primera igualdad es el teorema 6.10 para el álgebra $M_D(\mathbb C)$ (un solo
bloque): $T$ es un punto fijo si y solo si conmuta con todas las $E_i$ y $E_i^*$, es decir, si y solo
si $T\in\mathcal A'$. El dual $\Phi^*$ es de nuevo un canal unital: sus operadores de Kraus $E_i^*$
cumplen $\sum_iE_iE_i^*=1$ (preserva la traza) y $\sum_iE_i^*E_i=1$ (unital), y generan la misma
$*$-álgebra $\mathcal A$. Así, el teorema 6.10 aplicado a $\Phi^*$ da $\mathrm{Fix}(\Phi^*)=\mathcal A'$.
La última afirmación es el teorema 2.6. $\square$

*Cálculo:* $\mathrm{Fix}(\Phi)=\mathcal A'$ para el desfase colectivo, con conmutantes y puntos fijos calculados como núcleos — [tutorial, sección 10.1](../../tutorial/es/10_quantum_computing.md#s-10-1).

Para todo canal, $\mathcal A'\subseteq\mathrm{Fix}(\Phi^*)$ (proposición 10.1(a)), pero sin
unitalidad el teorema falla: para el amortiguamiento de amplitud, $\mathrm{Fix}(\Phi)=\mathbb C|0\rangle\langle0|$,
mientras que $\mathcal A'=\mathrm{Fix}(\Phi^*)=\mathbb C1$ (sección 6.4 y ejercicio 1). (Históricamente: Kribs 2003,
Arias–Gheondea–Gudder 2002.)

<a id="r-10-3"></a>

**Corolario 10.3 (desfase colectivo).** Sean $U=e^{-i\theta Z_{\mathrm{tot}}}$,
$Z_{\mathrm{tot}}=\sum_jZ_j$, con $e^{2i\theta d}\neq1$ para $d=1,\dots,n$, y
$\Phi(\rho)=\tfrac12(U\rho U^*+U^*\rho U)$. Entonces $\mathcal A=\mathrm{span}\{P_0,\dots,P_n\}$, cada $H_w$
es un subespacio libre de decoherencia, y $\mathrm{Fix}(\Phi)=\bigoplus_wB(H_w)\cong\bigoplus_wM_{\binom nw}$,
de dimensión $\sum_w\binom nw^2$.

*Demostración.* $Z_{\mathrm{tot}}|b\rangle=(n-2|b|)|b\rangle$, luego $U=\sum_we^{-i\theta(n-2w)}P_w$, con
valores propios distintos (sus cocientes son $e^{2i\theta(w-w')}$, $0<|w-w'|\le n$). Los polinomios de
Lagrange de la demostración del teorema 1.5(a) dan $P_w=p_w(U)\in\mathcal A$; recíprocamente, $U,U^*$
están en el álgebra $\mathrm{span}\{P_w\}$. Así, $\mathcal A=\mathrm{span}\{P_w\}$ ($k_c=1$,
$m_c=\binom nw$), y la proposición 10.1(b) se aplica a cada $H_w$. Si $T$ conmuta con todas las $P_w$,
entonces $P_vTP_{v'}=TP_vP_{v'}=0$ para $v\neq v'$, es decir, $T$ es diagonal por bloques, y
recíprocamente; el teorema 10.2 da el resto. $\square$

*Cálculo:* las dimensiones $5$, $70$, $70$ para $n=4$ — [tutorial, sección 10.1](../../tutorial/es/10_quantum_computing.md#s-10-1). *En la biblioteca:* $\bigoplus_wM_{\binom nw}$ es `SpinChain(n).algebra`, sobre la cual este canal es la identidad — [tutorial, sección 10.1](../../tutorial/es/10_quantum_computing.md#s-10-1).

> **Significado físico** (interpretación física, no un teorema de este curso). Un campo fluctuante que es el mismo en todos los cúbits desplaza la fase
> de $|b\rangle$ en proporción a $n-2|b|$. Los estados de un mismo peso de Hamming adquieren una fase
> común, de modo que $|01\rangle,|10\rangle$ codifican un cúbit inmune a este ruido (iones atrapados:
> Kielpinski et al., *Science* 291, 2001), mientras que las coherencias entre pesos decaen: el ruido
> impone una regla de superselección (sección 3.3).

**Decoherencia colectiva.** Llamamos **colectivo** a un canal si sus operadores de Kraus son múltiplos
de rotaciones $e^{-i\theta\,\mathbf n\cdot\mathbf S}$ ($\theta\in\mathbb R$, $\mathbf n\in\mathbb R^3$ un vector unitario),
$S_\alpha=\tfrac12\sum_j\sigma_{\alpha,j}$, $\sigma\in\{X,Y,Z\}$; esto modela un campo que fluctúa también en dirección
(el modelo, no un teorema). Como $\mathbf n\cdot\mathbf S$ es autoadjunto, $e^{-i\theta\,\mathbf n\cdot\mathbf S}$ es un
polinomio en él (demostración del teorema 1.5(a)), así que estos operadores de Kraus están en la $*$-álgebra $\mathcal B$
generada por $1,S_x,S_y,S_z$.

<a id="r-10-4"></a>

**Proposición 10.4 (subespacio singlete).** Sea $K=\ker S_x\cap\ker S_y\cap\ker S_z$. Todo canal
cuyos operadores de Kraus están en $\mathcal B$ fija todo estado con soporte en $K$.

*Demostración.* $\mathcal B$ está generada linealmente por las palabras en $S_x,S_y,S_z$ (y $1$); las
palabras de longitud $\ge1$ anulan $K$, así que cada elemento de $\mathcal B$ actúa sobre $K$ como un
escalar. Aplíquese la proposición 10.1(b). $\square$

Para $n=4$, la proposición 4.13 ([capítulo 4](04_symmetries.md), caso $L=4$) da
$\mathcal B\cong M_5\oplus M_3\oplus M_1$ con multiplicidades $1,3,2$ (espín total $2,1,0$), luego $\dim\mathcal B=35$ y
$\mathcal B'\cong M_1\oplus M_3\oplus M_2$ tiene dimensión $14$. Además, $\dim K=m_0=2$: los $S_a$ actúan sumando a sumando,
así que $K$ es la suma de sus intersecciones con los sumandos; sobre el sumando de espín $S$ actúan como $S^{(S)}_a\otimes1$
(corolario 4.12), y
$\sum_bu_b\otimes f_b$ ($f_b$ una base de $\mathbb C^{m_S}$) es anulado por todos ellos si y solo si cada $u_b$ es anulado por
todos los $S^{(S)}_a$. Esos $u$ forman un subespacio invariante del espín $S$ irreducible (teorema 4.11), que no es todo el espacio
para $S>0$ (pues $S^{(S)}_z\neq0$), luego es $0$, mientras que para $S=0$ es $\mathbb C^1$. La
proposición 10.1(a) y su demostración valen literalmente para cualquier $*$-álgebra de la forma
$V(\bigoplus_cM_{k_c}\otimes1_{m_c})V^*$ que contenga los operadores de Kraus, aquí $\mathcal B$: cuatro cúbits contienen un
cúbit libre de decoherencia ($k=1$, $m=2$) y un cútrit sin ruido ($k=3$, $m=3$) para todo canal cuyos operadores de Kraus
están en $\mathcal B$.

*Cálculo:* $\dim\mathcal B'=14$, $\dim\mathcal B=35$, $\dim K=2$, y la proposición 10.4 para rotaciones colectivas aleatorias — [tutorial, sección 10.1](../../tutorial/es/10_quantum_computing.md#s-10-1).

<a id="s-10-2"></a>

## 10.2 Corrección de errores: el teorema de Knill–Laflamme

Un **código** es un subespacio no nulo $\mathcal C\subseteq\mathbb C^D$ con proyección $P$. **Corrige** $\Phi$
si existe un canal $\mathcal R$ (la **recuperación**) con $\mathcal R(\Phi(\rho))=\rho$ para todo estado
$\rho$ con soporte en $\mathcal C$.

<a id="r-10-5"></a>

**Teorema 10.5 (Knill–Laflamme).** $\mathcal C$ corrige $\Phi$ si y solo si existen números
$\lambda_{ij}$ con
$$
PE_i^*E_jP=\lambda_{ij}P\qquad\text{for all }i,j.
$$

*Demostración.* *Necesidad.* Sean $R_m$ los operadores de Kraus de $\mathcal R$; entonces
$\mathcal R\circ\Phi$ tiene los operadores de Kraus $R_mE_i$. Para un vector unitario $\psi\in\mathcal C$,
$\sum_{m,i}R_mE_i|\psi\rangle\langle\psi|E_i^*R_m^*=|\psi\rangle\langle\psi|$, y para $\phi\perp\psi$ el valor
esperado en $\phi$ da $\sum_{m,i}|\langle\phi,R_mE_i\psi\rangle|^2=0$, luego $R_mE_i\psi\in\mathbb C\psi$.
Por tanto $R_mE_iP=c_{mi}P$: una aplicación lineal $T$ con $Tx\in\mathbb Cx$ para todo $x$ es un
múltiplo de $1$ (si $Tx=ax$, $Ty=by$, $T(x+y)=c(x+y)$, con $x,y$ independientes, entonces $a=c=b$).
Como $\sum_mR_m^*R_m=1$, $PE_i^*E_jP=\sum_m(R_mE_iP)^*(R_mE_jP)=\sum_m\bar c_{mi}c_{mj}P$.

*Suficiencia.* $\Lambda=(\lambda_{ij})$ es hermítica, pues $\lambda_{ji}P=(PE_i^*E_jP)^*=\bar\lambda_{ij}P$, y
positiva, pues $\sum_{ij}\bar c_ic_j\lambda_{ij}P=PF^*FP\ge0$ para $F=\sum_jc_jE_j$. Escribamos
$\Lambda=UDU^*$, con $U$ unitaria, $D=\mathrm{diag}(d_k)$, $d_k\ge0$ (teorema 1.3), y
$F_k=\sum_iu_{ik}E_i$. Entonces $\sum_kF_k\rho F_k^*=\sum_{ij}(UU^*)_{ij}E_i\rho E_j^*=\Phi(\rho)$ y
$PF_k^*F_lP=(U^*\Lambda U)_{kl}P=d_k\delta_{kl}P$. Si $d_k=0$, entonces $\|F_kPx\|^2=0$ para todo $x$,
luego $F_kP=0$. Si $d_k>0$, sea $V_k=F_kP/\sqrt{d_k}$; entonces $V_k^*V_l=\delta_{kl}P$, luego
$V_kP=V_k$, las $Q_k=V_kV_k^*$ son proyecciones mutuamente ortogonales y $Q=\sum_kQ_k$ es una
proyección. La recuperación
$$
\mathcal R(\sigma)=\sum_{k:d_k>0}V_k^*\sigma V_k+(1-Q)\sigma(1-Q)
$$
es un canal, ya que $\sum_kV_kV_k^*+(1-Q)^2=1$. Para $\rho=P\rho P$, $\Phi(\rho)=\sum_{d_k>0}d_kV_k\rho V_k^*$;
como $V_l^*V_k=\delta_{kl}P$ y $(1-Q)V_k=(1-Q)Q_kV_k=0$, $\mathcal R(\Phi(\rho))=\sum_kd_k\rho$, y
$\sum_kd_k=\mathrm{Tr}\,\Lambda=1$ porque $\sum_i\lambda_{ii}P=P\sum_iE_i^*E_iP=P$. $\square$

(Knill–Laflamme 1997; Nielsen–Chuang, teorema 10.1.) $\mathcal R$ mide el **síndrome** —cuál de los
subespacios ortogonales $Q_k\mathbb C^D$ contiene el estado— y devuelve el estado al código mediante
una rotación.

**Ejemplo: el código de inversión de bit.** $\mathcal C=\mathrm{span}\{|000\rangle,|111\rangle\}$ y
"a lo sumo una inversión de bit", $E_0=\sqrt{p_0}\,1$, $E_j=\sqrt{p_j}\,X_j$. Para $i\neq j$, $X_iX_j$
(con $X_0=1$) invierte uno o dos bits y envía $\mathcal C$ a $\mathcal C^\perp$; así,
$\Lambda=\mathrm{diag}(p_0,\dots,p_3)$, $V_k=X_kP$, $Q=1$ y $\mathcal R(\sigma)=\sum_{k=0}^3PX_k\sigma X_kP$.
Si cada cúbit se invierte de forma independiente con probabilidad $p$, la condición falla (los operadores
de Kraus incluyen múltiplos de $X_1X_2$ y de $X_3$, y $PX_1X_2X_3P=|000\rangle\langle111|+|111\rangle\langle000|$ no es un
múltiplo de $P$), y $\mathcal R$ envía $|0_L\rangle=|000\rangle$ a $|1_L\rangle$ siempre que se producen $\ge2$ inversiones y
lo deja en su sitio en otro caso: la fidelidad de $|0_L\rangle$ tras la recuperación es la
probabilidad de a lo sumo una inversión, $(1-p)^3+3p(1-p)^2=1-3p^2+2p^3=0.972$ para $p=0.1$.

*Cálculo:* el teorema 10.5 en este ejemplo: la matriz $\Lambda$, la recuperación $\mathcal R$ y la fidelidad $0.972$ — [tutorial, sección 10.2](../../tutorial/es/10_quantum_computing.md#s-10-2). *En la biblioteca:* los canales dados por operadores de Kraus, `channels.kraus_channel`, y su composición `R @ Phi` — [tutorial, sección 10.2](../../tutorial/es/10_quantum_computing.md#s-10-2).

<a id="s-10-3"></a>

## 10.3 Corrección de errores en álgebras de operadores

La imagen de Heisenberg permite proteger menos que el código entero. Sea $\mathfrak A$ un álgebra
de von Neumann sobre $\mathcal C=P\mathbb C^D$ (una $*$-álgebra de operadores sobre $\mathcal C$ que
contiene a $P$), los observables que se quieren proteger. $\mathfrak A$ es **corregible** para $\Phi$ si
existe un canal $\mathcal R$ con
$$
P\,(\mathcal R\circ\Phi)^*(A)\,P=A\qquad\text{for all }A\in\mathfrak A,
$$
es decir, $\mathrm{Tr}(\mathcal R(\Phi(\rho))A)=\mathrm{Tr}(\rho A)$ para todos los estados $\rho$
con soporte en $\mathcal C$ y todo $A\in\mathfrak A$ (el lado izquierdo es
$\mathrm{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(A)P)$, y un operador $Y\in B(\mathcal C)$ con $\mathrm{Tr}(\rho Y)=0$ para todos
estos $\rho$ es nulo: $\langle x,Yx\rangle=\langle Px,YPx\rangle=0$ para todo $x$, y polarización como en la demostración del
corolario 1.4). Para $\mathfrak A=B(\mathcal C)$ esto es la corrección de la sección 10.2: si
$\sigma=\mathcal R(\Phi(\rho))$ cumple $P\sigma P=\rho$ (tómense todos los $A\in B(\mathcal C)$), entonces, como
$\mathcal R\circ\Phi$ preserva la traza,
$\mathrm{Tr}(1-P)\sigma(1-P)=\mathrm{Tr}\,\sigma-\mathrm{Tr}\,P\sigma P=1-1=0$, luego $\sigma^{1/2}(1-P)=0$ y
$\sigma=\rho$; recíprocamente, $\mathcal R(\Phi(\rho))=\rho$ da la igualdad de trazas. El extremo opuesto es un
álgebra conmutativa: solo se protege una etiqueta clásica.

<a id="r-10-6"></a>

**Proposición 10.6 (corrección de la etiqueta de un sector).** Sean $P=\sum_{c=1}^CP_c$ con
proyecciones $P_c\neq0$ mutuamente ortogonales, y $\mathfrak A=\mathrm{span}\{P_1,\dots,P_C\}$. Son
equivalentes:
1. existen $N_1,\dots,N_C\ge0$, $\sum_cN_c=1$, con $\mathrm{Tr}(N_c\Phi(\rho))=1$ para todo $c$ y
   todo estado $\rho$ con soporte en $P_c\mathbb C^D$;
2. $P_cE_i^*E_jP_d=0$ para todos $i,j$ y $c\neq d$, es decir, $PE_i^*E_jP\in\mathfrak A'$, el conmutante
   de $\mathfrak A$ en $B(\mathcal C)$;
3. $\mathfrak A$ es corregible para $\Phi$.

*Demostración.* $T\in B(\mathcal C)$ conmuta con todas las $P_c$ si y solo si $P_cTP_d=0$ para
$c\neq d$ (como en el corolario 10.3).

(2)$\Rightarrow$(1). Los subespacios $S_c=\mathrm{span}\{E_iP_cx\}$ son mutuamente ortogonales,
pues $\langle E_iP_cx,E_jP_dy\rangle=\langle x,P_cE_i^*E_jP_dy\rangle=0$. Sea $N_c$ la proyección sobre
$S_c$ para $c<C$ y $N_C=1-\sum_{c<C}N_c$, también una proyección. Para $\rho$ con soporte en
$P_c\mathbb C^D$, $\Phi(\rho)=\sum_iE_iP_c\rho P_cE_i^*$ tiene su imagen en $S_c$, y $N_c$ es la
identidad sobre $S_c$ (para $c=C$ porque $S_C\perp S_{c'}$, $c'<C$), luego
$\mathrm{Tr}(N_c\Phi(\rho))=\mathrm{Tr}\,\Phi(\rho)=1$.
(1)$\Rightarrow$(2). Para un vector unitario $x\in P_c\mathbb C^D$ y $d\neq c$, los números
$\mathrm{Tr}(N_{d'}\Phi(|x\rangle\langle x|))\ge0$ suman $1$ y el correspondiente a $d'=c$ vale $1$,
luego $0=\mathrm{Tr}(N_d\Phi(|x\rangle\langle x|))=\sum_i\|N_d^{1/2}E_ix\|^2$; así $N_dE_iP_c=0$ y,
sumando sobre $d$, $N_cE_iP_c=E_iP_c$. Por tanto $P_cE_i^*E_jP_d=P_cE_i^*N_cE_jP_d=0$ para $c\neq d$.
(1),(2)$\Rightarrow$(3). Sean $\omega_c=P_c/\mathrm{Tr}\,P_c$ y $\mathcal R(\sigma)=\sum_c\mathrm{Tr}(N_c\sigma)\,\omega_c$
("medir y luego preparar"), un canal con los operadores de Kraus $\sqrt\mu\,|u\rangle\langle f|N_c^{1/2}$
para los pares propios $(\mu,u)$ de $\omega_c$ y los vectores $f$ de una base ortonormal. Como
$\mathrm{Tr}(\mathcal R(\sigma)P_d)=\mathrm{Tr}(N_d\sigma)$, $\mathcal R^*(P_d)=N_d$. Usando
$N_dE_iP_c=\delta_{cd}E_iP_c$ (del paso anterior) y después $\sum_iE_i^*E_i=1$,
$$
P(\mathcal R\circ\Phi)^*(P_d)P=P\,\Phi^*(N_d)\,P=\sum_i\sum_{c,c'}P_{c'}E_i^*N_dE_iP_c=\sum_i\sum_{c'}P_{c'}E_i^*E_iP_d=P_d\Big(\sum_iE_i^*E_i\Big)P_d=P_d .
$$
(3)$\Rightarrow$(1). Sean $N_c=\mathcal R^*(P_c)$ para $c<C$ y $N_C=\mathcal R^*(1-\sum_{c<C}P_c)$; son
positivos y suman $\mathcal R^*(1)=1$. Para $\rho$ con soporte en $P_c\mathbb C^D$ y $d<C$,
$\mathrm{Tr}(N_d\Phi(\rho))=\mathrm{Tr}(\rho\,P(\mathcal R\circ\Phi)^*(P_d)P)=\mathrm{Tr}(\rho P_d)=\delta_{cd}$,
y por tanto también $\mathrm{Tr}(N_C\Phi(\rho))=\delta_{cC}$. $\square$

<a id="r-10-7"></a>

**Teorema 10.7 (Bény–Kempf–Kribs).** Un álgebra de von Neumann
$\mathfrak A$ sobre $\mathcal C$ es corregible para $\Phi$ si y solo si $PE_i^*E_jP\in\mathfrak A'$ para
todos $i,j$, donde $\mathfrak A'$ es el conmutante de $\mathfrak A$ en $B(\mathcal C)$.

*Demostración.* Es el teorema D.9, demostrado en el [apéndice D](D_gksl.md), cuya demostración construye además
un canal de recuperación: mide el sector de $\mathfrak A\cong\bigoplus_cM_{k_c}$ y un síndrome, y deshace el error en
cada factor $\mathbb C^{k_c}$ como en la demostración del teorema 10.5. $\square$

Para $\mathfrak A=B(\mathcal C)$, $\mathfrak A'=\mathbb CP$ y esto es el teorema 10.5; la proposición
10.6 es el caso del $\mathfrak A=\mathrm{span}\{P_c\}$ conmutativo (corolario D.10(a),(b)). Para
$\mathfrak A=M_k\otimes1_m$ sobre $\mathcal C=\mathbb C^k\otimes\mathbb C^m$ se obtienen los códigos de
operadores (de subsistemas); los subsistemas sin ruido de la sección 10.1 son códigos de este tipo, corregidos
por $\mathcal R=\mathrm{id}$ (corolario D.10(c), con los dos factores tensoriales intercambiados). Por el
teorema 2.6, toda $\mathfrak A$ es $*$-isomorfa a algún $\bigoplus_cM_{k_c}$; cuando $C\ge2$ y algún $k_c\ge2$, una
$\mathfrak A$ corregible protege la etiqueta del sector $c$ junto con la información cuántica de cada $M_{k_c}$: un
código **híbrido** (observación tras el corolario D.10, con un ejemplo en el que la condición de Knill–Laflamme
falla pero se aplica el teorema 10.7). Históricamente: C. Bény, A. Kempf, D. W. Kribs, *Phys. Rev. Lett.* 98,
100502 (2007) y *Phys. Rev. A* 76, 042303 (2007); códigos de operadores: Kribs–Laflamme–Poulin,
*Phys. Rev. Lett.* 94, 180501 (2005).

**Ejemplo.** El código de repetición bajo inversiones de fase independientes (probabilidad $q$ por
cúbit) seguidas de a lo sumo una inversión de bit. Los operadores de Kraus son múltiplos de
$X_kZ^s$, $s\in\{0,1\}^3$, y cada $E_i^*E_j$ invierte a lo sumo dos bits, de modo que la condición 2 se
cumple para $P_0=|000\rangle\langle000|$, $P_1=|111\rangle\langle111|$: el bit clásico sobrevive, y se lee
por voto mayoritario, $N_0=\sum_{|b|\le1}|b\rangle\langle b|$, $N_1=1-N_0$. El cúbit no sobrevive:
$PZ_1P$ no es un múltiplo de $P$, y tras la recuperación $\mathcal R$ de la sección 10.2 el estado
$|{+_L}\rangle=(|000\rangle+|111\rangle)/\sqrt2$ sufre la inversión de fase lógica $Z_1Z_2Z_3$ exactamente
cuando se ha producido un número impar de inversiones de fase; su fidelidad es
$(1+(1-2q)^3)/2=0.8645$ para $q=0.05$. (El canal es $\Phi_{\mathrm{one}}\circ\mathcal Z$, con $\Phi_{\mathrm{one}}$ el
canal "a lo sumo una inversión de bit" del ejemplo de inversión de bit y $\mathcal Z$ las inversiones de fase
independientes; sobre $\mathcal C$ cada $Z^s$ actúa como $Z_1^{|s|}$, así que
$\mathcal Z(\rho)=\sum_sq^{|s|}(1-q)^{3-|s|}Z_1^{|s|}\rho Z_1^{|s|}$ tiene de nuevo soporte en $\mathcal C$, y
$\mathcal R\circ\Phi_{\mathrm{one}}$ es la identidad sobre tales estados por el teorema 10.5.
$Z_1$ envía $|{+_L}\rangle$ a $|{-_L}\rangle\perp|{+_L}\rangle$, y $|s|$ es par con probabilidad
$\frac12\big((1-q+q)^3+(1-q-q)^3\big)$.)

*Cálculo:* la proposición 10.6 en este ejemplo: la condición 2, el voto mayoritario $N_0,N_1$ y la fidelidad $0.8645$ de $|{+_L}\rangle$ — [tutorial, sección 10.3](../../tutorial/es/10_quantum_computing.md#s-10-3).

> **Significado físico** (interpretación física, no un teorema de este curso). Proteger un álgebra
> conmutativa significa proteger información clásica (un bit, la etiqueta de un sector de superselección),
> lo que es mucho más barato que proteger un cúbit: el ruido que conmuta con las proyecciones de los
> sectores es inofensivo (entonces $P_cE_i^*E_jP_d=E_i^*E_jP_cP_d=0$ para $c\neq d$, condición 2 de la
> proposición 10.6).

<a id="s-10-4"></a>

## 10.4 Mezcladores con restricciones

El ansatz de operadores alternantes cuántico (Hadfield et al., *Algorithms* 12, 34, 2019) prepara
$|\beta,\gamma\rangle=e^{-i\beta_pH_M}e^{-i\gamma_pH_C}\cdots e^{-i\beta_1H_M}e^{-i\gamma_1H_C}|b_0\rangle$, con $H_C$ un
hamiltoniano de coste diagonal y $H_M$ un **mezclador**. Para una restricción $|b|=w$ (bisección de
grafos, carteras de $w$ activos) se toma $|b_0|=w$ y el **mezclador XY** en un anillo,
$H_M=\tfrac14\sum_{(i,j)}(X_iX_j+Y_iY_j)$.

<a id="r-10-8"></a>

**Proposición 10.8.** $H_M$, $H_C$, $e^{-i\beta H_M}$ y $e^{-i\gamma H_C}$ están en
$\mathcal M=\{P_0,\dots,P_n\}'=\bigoplus_wB(H_w)\cong\bigoplus_wM_{\binom nw}$, luego $|\beta,\gamma\rangle\in H_w$:
el ansatz nunca sale del subespacio factible.

*Demostración.* Sobre los cúbits $i,j$, $X\otimes X+Y\otimes Y=2(|01\rangle\langle10|+|10\rangle\langle01|)$ (úsese
$Y|0\rangle=i|1\rangle$, $Y|1\rangle=-i|0\rangle$): $X_iX_j+Y_iY_j$ intercambia los bits $i,j$ cuando son
distintos y anula $|b\rangle$ en otro caso, de modo que envía cada $H_w$ en sí mismo. Un $H_C$ diagonal
conmuta con las $P_w$ diagonales, y lo mismo hacen las exponenciales (teorema 1.5(d)). $\square$

*En la biblioteca:* $\mathcal M$ es `SpinChain(n).algebra` y el mezclador XY es `SpinChain.xxz` con $\Delta=0$ — [tutorial, sección 10.4](../../tutorial/es/10_quantum_computing.md#s-10-4). *Cálculo:* la proposición 10.8 para la bisección de grafos sobre seis cúbits, con la biblioteca y con matrices densas — [tutorial, sección 10.4](../../tutorial/es/10_quantum_computing.md#s-10-4).

El **amortiguamiento de amplitud** del cúbit $i$ ($|1\rangle\to|0\rangle$ con probabilidad $g$) tiene los
operadores de Kraus $K_0^{(i)}=|0\rangle\langle0|_i+\sqrt{1-g}\,|1\rangle\langle1|_i$ y
$K_1^{(i)}=\sqrt g\,|0\rangle\langle1|_i$, que reduce el peso: una aplicación entre sectores.

*En la biblioteca:* `SpinChain.site_amplitude_damping`, un `InterSectorChannel` — [tutorial, sección 10.4](../../tutorial/es/10_quantum_computing.md#s-10-4).

<a id="r-10-9"></a>

**Proposición 10.9 (fuga y postselección).** Sean $\rho_0$ con soporte en $H_w$, $U_l$ unitarias
en $\mathcal M$, $\mathcal D$ el amortiguamiento de todos los cúbits, y
$\rho_l=\mathcal D(U_l\rho_{l-1}U_l^*)$. Entonces todo $\rho_l$ es diagonal por bloques con pesos
$\le w$, y
$$
P_w\rho_pP_w=(1-g)^{wp}\;U\rho_0U^*,\qquad U=U_p\cdots U_1 :
$$
el peso $w$ se encuentra con probabilidad $(1-g)^{wp}$, y el estado postseleccionado es el estado sin
ruido.

*Demostración.* Los operadores de Kraus de $\mathcal D$ son $K_s=\prod_iK^{(i)}_{s_i}$, $s\in\{0,1\}^n$;
$K_s$ envía $H_v$ a $H_{v-|s|}$, es decir, $P_vK_s=K_sP_{v+|s|}$ ($P_v=0$ para $v>n$). Si $\sigma$ es
diagonal por bloques con pesos $\le w$, entonces $P_vK_s\sigma K_s^*P_{v'}=K_sP_{v+|s|}\sigma P_{v'+|s|}K_s^*$
se anula para $v\neq v'$ y para $v>w-|s|$: $\mathcal D(\sigma)$ es diagonal por bloques con pesos
$\le w$, y su bloque de peso $w$ proviene solo de $s=0$. Como $K_0|b\rangle=(1-g)^{|b|/2}|b\rangle$,
$P_w\mathcal D(\sigma)P_w=(1-g)^wP_w\sigma P_w$. Las $U_l$ conmutan con las $P_v$; se concluye por
inducción sobre $l$. $\square$

*Cálculo:* la proposición 10.9 en el circuito de bisección, con postselección mediante `DensityMatrix.condition_on` — [tutorial, sección 10.4](../../tutorial/es/10_quantum_computing.md#s-10-4).

**Tiempo continuo.** Lo mismo vale para la ecuación de Lindblad (teorema 8.6, [capítulo 8](08_dynamics.md))
$$
\dot\rho=-i[H,\rho]+\gamma_1\sum_i\Big(S^-_i\rho S^+_i-\tfrac12\{S^+_iS^-_i,\rho\}\Big)+\gamma_2\sum_k\Big(L_k\rho L_k^*-\tfrac12\{L_k^*L_k,\rho\}\Big)
$$
con $H\in\mathcal M$ y saltos $L_k\in\mathcal M$ dentro de los sectores, p. ej. el desfase puro
$L_k=S^z_k$. Aquí $\sum_iS^+_iS^-_i=\sum_i|1\rangle\langle1|_i=N$, el **operador número**,
$NP_v=vP_v$: es central.

<a id="r-10-10"></a>

**Proposición 10.10.** Si $\rho(0)$ tiene soporte en $H_w$, entonces $P_w\rho(t)P_w=e^{-\gamma_1wt}\sigma(t)$,
donde $\sigma$ resuelve la misma ecuación con $\gamma_1=0$ y $\sigma(0)=\rho(0)$. Por tanto, la tasa de
aceptación es $e^{-\gamma_1wt}$, y el estado postseleccionado $\sigma(t)$ no depende de $\gamma_1$.

*Demostración.* Sea $B_{vv'}=P_v\rho P_{v'}$. Como $H,L_k,N$ conmutan con las $P_v$ y
$P_vS^-_i=S^-_iP_{v+1}$,
$\dot B_{vv'}=-i[H,B_{vv'}]+\gamma_1\sum_iS^-_iB_{v+1,v'+1}S^+_i-\tfrac{\gamma_1}2(v+v')B_{vv'}+\gamma_2\mathcal L_2(B_{vv'})$,
con $\mathcal L_2$ la parte de desfase. Los bloques con $\max(v,v')>w$ satisfacen un sistema lineal
cerrado con valores iniciales nulos, luego se anulan (unicidad en el lema 8.1). Entonces
$\dot B_{ww}=-i[H,B_{ww}]-\gamma_1wB_{ww}+\gamma_2\mathcal L_2(B_{ww})$, de modo que
$e^{\gamma_1wt}B_{ww}$ resuelve la ecuación de $\sigma$ con el mismo valor inicial, y es igual a
$\sigma$ (lema 8.1). Esa ecuación vuelve a tener forma GKSL, luego $\mathrm{Tr}\,\sigma(t)=1$
(teorema 8.6). $\square$

Compárese con la proposición 8.11, que admite un $\sum_jL_j^*L_j$ no central pero ningún salto
dentro de los sectores: aquí el término de pérdida $\gamma_1N$ es central, de modo que el factor
sin saltos es el escalar $e^{-\gamma_1wt}$ y se admite el desfase dentro del sector.

Para $\gamma_2=0$, $\sigma(t)=e^{-iHt}\rho(0)e^{iHt}$: la postselección es perfecta.

*Cálculo:* la proposición 10.10 para el mezclador XY sobre seis cúbits, con `dynamics.lindblad_evolve` — [tutorial, sección 10.4](../../tutorial/es/10_quantum_computing.md#s-10-4).

> **Significado físico** (interpretación física, no un teorema de este curso). La relajación de energía ($T_1$), uno de los errores principales de los
> cúbits superconductores, es un salto a un sector inferior. Medir el peso de Hamming —una proyección
> central, inofensiva dentro de un sector (sección 3.3)— lo detecta, y la postselección lo elimina
> al precio de la tasa de aceptación. El desfase conmuta con las proyecciones de los sectores y no
> se detecta.

<a id="s-10-5"></a>

## 10.5 Con qué comparar

*En la biblioteca:* la comparación con otro software (Qiskit, Cirq, PennyLane, QuTiP, stim) — [tutorial, sección 10.5](../../tutorial/es/10_quantum_computing.md#s-10-5).

<a id="x-10"></a>

## Ejercicios

1. (El teorema 10.2 necesita la unitalidad; compárese con la sección 6.4.) Para el amortiguamiento
   de amplitud de un cúbit, $E_0=|0\rangle\langle0|+\sqrt{1-g}\,|1\rangle\langle1|$, $E_1=\sqrt g\,|0\rangle\langle1|$,
   $0<g<1$, demuestre que $\mathcal A'=\mathbb C1$, $\mathrm{Fix}(\Phi)=\mathbb C|0\rangle\langle0|$ y
   $\mathrm{Fix}(\Phi^*)=\mathbb C1$.
2. Demuestre que $\sum_w\binom nw^2=\binom{2n}n$. Para el canal del corolario 10.3 con $\theta=\pi/2$
   (que viola su hipótesis), ¿cuáles son $\mathcal A$ y
   $\mathrm{Fix}(\Phi)$? Compruebe $\dim\mathrm{Fix}(\Phi)$ para $n=4$.
3. (Linealidad de la corrección de errores.) Sea $PE_a^*E_bP=\lambda_{ab}P$ para un número finito de
   operadores $E_a$ (no necesariamente operadores de Kraus de un canal), y constrúyase $\mathcal R$ a
   partir de $\Lambda$ como en la demostración del teorema 10.5. Demuestre que $\mathcal R$ corrige
   todo canal cuyos operadores de Kraus están en $\mathrm{span}\{E_a\}$. Deduzca que el código
   de inversión de bit corrige el error coherente $e^{-i\phi X_1}$.

<details>
<summary>Soluciones</summary>

1. $E_1$, $E_1^*$ y sus productos dan todas las unidades matriciales, así que $\mathcal A=M_2$ y
   $\mathcal A'=\mathbb C1$. Para $T=(x_{kl})$, $\Phi(T)=\begin{pmatrix}x_{00}+gx_{11}&\sqrt{1-g}\,x_{01}\\\sqrt{1-g}\,x_{10}&(1-g)x_{11}\end{pmatrix}=T$
   si y solo si $x_{11}=x_{01}=x_{10}=0$, y
   $\Phi^*(T)=\begin{pmatrix}x_{00}&\sqrt{1-g}\,x_{01}\\\sqrt{1-g}\,x_{10}&(1-g)x_{11}+gx_{00}\end{pmatrix}$,
   que es igual a $T$ si y solo si $T=x_{00}1$. ($\Phi$ no es unital, y $\Phi^*$ es unital pero no
   preserva la traza.)

2. Cuéntense los subconjuntos de $n$ elementos de $\{1,\dots,2n\}$ según el número $w$ de elementos
   en la primera mitad: $\binom{2n}n=\sum_w\binom nw\binom n{n-w}=\sum_w\binom nw^2$. Para
   $\theta=\pi/2$, $U=(-i)^n\sum_w(-1)^wP_w$ solo ve la paridad de $w$, así que
   $\mathcal A=\mathrm{span}\{P_{\mathrm{even}},P_{\mathrm{odd}}\}$ con $P_{\mathrm{even}}=\sum_{w\text{ even}}P_w$,
   y por el teorema 10.2, $\mathrm{Fix}(\Phi)$ consta de los operadores diagonales por bloques
   respecto de la paridad: dimensión $2\cdot(2^{n-1})^2=128$ para $n=4$.

3. Con $F_k=\sum_au_{ak}E_a$, $\mathrm{span}\{E_a\}=\mathrm{span}\{F_k\}$, así que cada operador de
   Kraus del canal es $G_m=\sum_kh_{mk}F_k$, y $G_mP=\sum_{d_k>0}h_{mk}\sqrt{d_k}V_k$. Como
   $V_l^*V_k=\delta_{lk}P$ y $(1-Q)V_k=0$, $\mathcal R(G_m\rho G_m^*)=\sum_l|h_{ml}|^2d_l\,\rho$ para
   $\rho=P\rho P$, luego $\mathcal R(\Phi(\rho))=c\rho$; tomando trazas, $c=1$. Para
   $E_a\in\{1,X_1,X_2,X_3\}$, $\Lambda=1$, $V_k=X_kP$, y $\mathcal R$ es la recuperación de la sección
   10.2; $e^{-i\phi X_1}=\cos\phi\,1-i\sin\phi\,X_1$ está en el espacio generado.

</details>

*Cálculo:* los ejercicios 1–3, numéricamente — [tutorial, Ejercicios](../../tutorial/es/10_quantum_computing.md#x-10).

**Siguiente:** el [capítulo 11](11_optimisation.md): optimización sobre los grupos unitarios de un
álgebra, y desigualdades de trazas.
