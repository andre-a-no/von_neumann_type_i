[English](../02_von_neumann_algebras.md) · **Español** · [Français](../fr/02_von_neumann_algebras.md) · [Deutsch](../de/02_von_neumann_algebras.md) · [中文](../zh/02_von_neumann_algebras.md) · [日本語](../ja/02_von_neumann_algebras.md) · [Русский](../ru/02_von_neumann_algebras.md)

# 2. Álgebras de von Neumann de dimensión finita

Un álgebra de von Neumann es un álgebra de operadores que es cerrada para la operación de
"tomar dos veces el conmutante". En dimensión finita se trata de una condición puramente
algebraica, y el teorema principal de este capítulo afirma que toda álgebra de este tipo es,
salvo un cambio de base unitario, una suma directa de álgebras de matrices completas, cada
una repetida con cierta multiplicidad. Esta es la estructura sobre la que está construida la
biblioteca.

**Aprenderá** el teorema del doble conmutante, el teorema de estructura con demostración
completa, cómo un grupo de simetría da lugar a un álgebra de este tipo, qué es el centro y
cómo representa todo esto la biblioteca.

En todo el capítulo, *$*$-álgebra sobre $\mathbb C^n$* significa un subespacio lineal complejo
$\mathcal M\subseteq M_n(\mathbb C)$ cerrado para productos y adjuntos y que contiene la
identidad $1$.

## 2.1 Conmutantes y el teorema del doble conmutante

Para un conjunto $S\subseteq M_n(\mathbb C)$, el **conmutante** es
$S'=\{X\in M_n(\mathbb C):XA=AX\text{ for all }A\in S\}$. Siempre es un álgebra unital;
si $S$ es cerrado para adjuntos ($S^*=S$), entonces $S'$ es una $*$-álgebra, ya que
$XA=AX$ para todo $A\in S$ implica $X^*A^*=A^*X^*$ para todo $A\in S$. Todo elemento de $S$
conmuta con todo elemento de $S'$, luego $S\subseteq S''$; aplicado a $S'$, esto da
$S'\subseteq S'''$, y tomando conmutantes en $S\subseteq S''$ (lo cual invierte las
inclusiones) se obtiene $S'''\subseteq S'$. Por tanto $S'=S'''$.

Cuando un álgebra actúa sobre un subespacio $E$ de $\mathbb C^n$, identificamos $E$ con
$\mathbb C^{\dim E}$ mediante una base ortonormal, de modo que todos los resultados
enunciados sobre $\mathbb C^n$ se le aplican.

**Definición.** Un **álgebra de von Neumann** sobre $\mathbb C^n$ es una $*$-álgebra
$\mathcal M$ con $\mathcal M=\mathcal M''$. Su **centro** es
$Z(\mathcal M)=\mathcal M\cap\mathcal M'$. Un **factor** es un álgebra de von Neumann con
$Z(\mathcal M)=\mathbb C1$.

En dimensión infinita, $\mathcal M=\mathcal M''$ es una condición topológica (clausura en la
topología débil de operadores). En dimensión finita se cumple automáticamente:

**Lema 2.1.** Sea $\mathcal M$ una $*$-álgebra sobre $\mathbb C^n$ y $K\subseteq\mathbb C^n$
un subespacio con $\mathcal MK\subseteq K$. Entonces la proyección ortogonal $P$ sobre $K$
pertenece a $\mathcal M'$.

*Demostración.* Para $A\in\mathcal M$, $AK\subseteq K$ y $A^*K\subseteq K$ (pues $A^*\in\mathcal M$).
La segunda inclusión da $AK^\perp\subseteq K^\perp$: para $y\in K^\perp$, $x\in K$,
$\langle x,Ay\rangle=\langle A^*x,y\rangle=0$. Por tanto $A$ aplica $K$ en $K$ y $K^\perp$
en $K^\perp$, lo que significa que $AP=PAP=PA$. $\square$

**Teorema 2.2 (teorema del doble conmutante).** Toda $*$-álgebra $\mathcal M$ sobre
$\mathbb C^n$ satisface $\mathcal M=\mathcal M''$.

*Demostración.* Se tiene $\mathcal M\subseteq\mathcal M''$; sea $T\in\mathcal M''$.

*Paso 1: un vector.* Para $x\in\mathbb C^n$, el subespacio $K=\mathcal Mx$ es invariante
por $\mathcal M$ y contiene a $x=1x$. Por el lema 2.1, su proyección $P$ está en $\mathcal M'$,
luego $TP=PT$ y $Tx=TPx=PTx\in K$: existe $A\in\mathcal M$ con $Tx=Ax$.

*Paso 2: todos los vectores a la vez.* En $\mathbb C^n\otimes\mathbb C^n=(\mathbb C^n)^{\oplus n}$
consideramos la $*$-álgebra $\tilde{\mathcal M}=\{A\oplus\dots\oplus A:A\in\mathcal M\}$. Una
matriz por bloques $[X_{ij}]_{i,j=1}^n$ conmuta con todo $A\oplus\dots\oplus A$ si y solo si
$X_{ij}A=AX_{ij}$ para todos los $i,j$, luego $\tilde{\mathcal M}'=\{[X_{ij}]:X_{ij}\in\mathcal M'\}$.
El operador $\tilde T=T\oplus\dots\oplus T$ conmuta con toda $[X_{ij}]$ de este tipo, porque
$TX_{ij}=X_{ij}T$. Por tanto $\tilde T\in\tilde{\mathcal M}''$. Apliquemos el paso 1 a
$\tilde{\mathcal M}$, $\tilde T$ y el vector $\xi=(e_1,\dots,e_n)$: existe $A\in\mathcal M$ con
$(Te_1,\dots,Te_n)=(Ae_1,\dots,Ae_n)$, es decir, $T=A\in\mathcal M$. $\square$

Así pues, en dimensión finita "álgebra de von Neumann" y "$*$-álgebra unital de matrices"
son lo mismo. El teorema del doble conmutante ofrece, no obstante, el punto de vista
adecuado: en física, un álgebra de observables viene *dada* muy a menudo como un
conmutante: los operadores que conmutan con una simetría (sección 2.4).

## 2.2 El centro

**Lema 2.3.** Sea $\mathcal Z$ una $*$-álgebra conmutativa sobre $\mathbb C^n$. Existen
proyecciones $P_1,\dots,P_C\in\mathcal Z$ no nulas, ortogonales dos a dos, con
$\sum_cP_c=1$, tales que $\mathcal Z=\operatorname{span}\{P_1,\dots,P_C\}$. Son las
proyecciones minimales no nulas de $\mathcal Z$ y son únicas salvo el orden.

*Demostración.* Todo $z\in\mathcal Z$ es normal ($z^*\in\mathcal Z$ conmuta con $z$). Sean
$z_1,\dots,z_r$ generadores de $\mathcal Z$ como espacio vectorial. Descomponemos
$\mathbb C^n$ en subespacios propios comunes. Los subespacios propios de $z_1$ son
ortogonales entre sí y generan $\mathbb C^n$ (teorema 1.3); son invariantes por
$z_2,\dots,z_r$ y por sus adjuntos (que conmutan con $z_1$), de modo que cada $z_j$
restringido a un subespacio propio de $z_1$ es de nuevo normal, y podemos descomponer cada
subespacio propio de $z_1$ en subespacios propios de $z_2$, y así sucesivamente. El
resultado es una descomposición ortogonal $\mathbb C^n=E_1\oplus\dots\oplus E_C$ en
subespacios no nulos sobre los que cada $z\in\mathcal Z$ actúa como un escalar $\chi_c(z)$,
donde hemos reunido los subespacios con la misma función $\chi_c:\mathcal Z\to\mathbb C$, de
modo que $\chi_c\neq\chi_d$ para $c\neq d$. Sea $P_c$ la proyección sobre $E_c$. Entonces
$z=\sum_c\chi_c(z)P_c$ para todo $z\in\mathcal Z$.

Falta ver que $P_c\in\mathcal Z$. Para $d\neq c$ elegimos $z_{cd}\in\mathcal Z$ con
$\chi_c(z_{cd})\neq\chi_d(z_{cd})$. El producto
$\prod_{d\neq c}\big(z_{cd}-\chi_d(z_{cd})1\big)/\big(\chi_c(z_{cd})-\chi_d(z_{cd})\big)$
pertenece a $\mathcal Z$, actúa como $1$ sobre $E_c$ y como $0$ sobre cada $E_d$, $d\neq c$;
es $P_c$. Minimalidad y unicidad: una proyección $\sum_c a_cP_c$ tiene $a_c\in\{0,1\}$, luego
las proyecciones no nulas de $\mathcal Z$ son las sumas de conjuntos no vacíos de las $P_c$. $\square$

Aplicadas a $\mathcal Z=Z(\mathcal M)$, las proyecciones $P_c$ son las **proyecciones
centrales minimales** de $\mathcal M$. Todo $A\in\mathcal M$ se descompone como
$A=\sum_cAP_c$ con $AP_c=P_cAP_c$, luego
$\mathcal M=\bigoplus_c\mathcal MP_c$ (una suma directa ortogonal de álgebras que actúan sobre
$E_c=P_c\mathbb C^n$).

**Lema 2.4.** Cada $\mathcal M_c=\mathcal MP_c$, como álgebra sobre $E_c$ con unidad $P_c$, es
un factor.

*Demostración.* Sea $y\in\mathcal M_c$ que conmuta con $\mathcal M_c$. Para $A\in\mathcal M$,
$Ay=AP_cy=(AP_c)y=y(AP_c)=yP_cA=yA$, usando que $P_c$ es central y que $y=yP_c=P_cy$.
Así $y\in Z(\mathcal M)$, $y=\sum_d\chi_d(y)P_d$, e $y=yP_c=\chi_c(y)P_c$. $\square$

## 2.3 El teorema de estructura

**Lema 2.5 (factores).** Sea $\mathcal N$ un factor sobre un espacio $E$. Existen enteros
$k,m\ge1$ y una unitaria $V:\mathbb C^k\otimes\mathbb C^m\to E$ con
$V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$. Además, $V^*\mathcal N'V=1_k\otimes M_m(\mathbb C)$.

*Demostración.* *Paso 1: una proyección minimal.* $\mathcal N$ contiene proyecciones no nulas
(por ejemplo, $1$). Sea $e\in\mathcal N$ una proyección no nula de rango mínimo. Afirmamos que
$e\mathcal Ne=\mathbb Ce$. Sea $x=x^*\in e\mathcal Ne$. Cada proyección espectral $q$ de $x$
asociada a un valor propio no nulo $\lambda$ es $q=p(x)$, donde $p$ es el polinomio de
Lagrange con $p(\lambda)=1$ y $p(\mu)=0$ para $\mu\in(\operatorname{spec}x\cup\{0\})\setminus\{\lambda\}$
(demostración del teorema 1.5(a)). Como $p(0)=0$ y $x=ex=xe$, se tiene $q\in\mathcal N$ y
$q=eqe$, de modo que la imagen de $q$ está contenida en la imagen de $e$. Por minimalidad,
$q=e$ o $q=0$; por tanto $x$ tiene a lo sumo un valor propio no nulo y su subespacio propio
es la imagen de $e$, luego $x=\lambda e$. Un elemento general de $e\mathcal Ne$ es
$x_1+ix_2$ con $x_{1,2}=x_{1,2}^*\in e\mathcal Ne$.

*Paso 2: $e$ genera $E$.* Sea $F=\operatorname{span}\mathcal Ne E$. Es invariante por
$\mathcal N$, luego su proyección $p_F$ pertenece a $\mathcal N'$ (lema 2.1). También es
invariante por la $*$-álgebra $\mathcal N'$: para $b\in\mathcal N'$, $bAe\xi=Aeb\xi\in F$; así,
por el lema 2.1 aplicado a $\mathcal N'$, $p_F\in\mathcal N''=\mathcal N$ (teorema 2.2). Por
tanto $p_F\in Z(\mathcal N)=\mathbb C1$, y $p_F\neq0$ (pues $eE\subseteq F$), luego $F=E$.

*Paso 3: unidades matriciales.* Llamamos *$e$-isometría* a un $v\in\mathcal N$ con $v^*v=e$;
entonces $vv^*$ es una proyección de $\mathcal N$. (En primer lugar, $v=ve$: para todo $x$,
$\|v(1-e)x\|^2=\langle(1-e)x,v^*v(1-e)x\rangle=\langle(1-e)x,e(1-e)x\rangle=0$. Por tanto
$vv^*v=ve=v$ y $(vv^*)^2=vev^*=vv^*$. Para las $v_i$ de más abajo esto da $e_iv_i=v_i$ y
$v_i^*e_i=v_i^*$.) Elegimos $e$-isometrías $v_1=e,v_2,\dots,v_k$ con imágenes
$e_i=v_iv_i^*$ ortogonales entre sí, con $k$ maximal (las imágenes son subespacios
ortogonales de la misma dimensión $\operatorname{rank}e$, luego
$k\le\dim E/\operatorname{rank}e$). Afirmamos que $\sum_ie_i=1$. En caso contrario,
$g=1-\sum_ie_i\neq0$ es una proyección de $\mathcal N$. Por el paso 2, $gAe\neq0$ para algún
$A\in\mathcal N$ (si $gAe=0$ para todo $A$, entonces $gAe\xi=0$ para todos los $A,\xi$, luego
$g$ se anula sobre $F=E$). Sea $y=gAe$. Entonces $y^*y\in e\mathcal Ne$ es positivo y no nulo,
luego $y^*y=\lambda e$ con $\lambda>0$, y $v=y/\sqrt\lambda$ es una $e$-isometría con imagen
contenida en $gE$, ortogonal a todas las $e_i$, lo que contradice la maximalidad.

Para $i\neq j$, $v_i^*v_j=v_i^*e_ie_jv_j=0$. Sea $E_{ij}=v_iv_j^*$. Entonces
$E_{ij}E_{kl}=v_iv_j^*v_kv_l^*=\delta_{jk}v_iev_l^*=\delta_{jk}E_{il}$,
$E_{ij}^*=E_{ji}$ y $\sum_iE_{ii}=1$. Para $A\in\mathcal N$,
$A=\sum_{ij}e_iAe_j=\sum_{ij}v_i(v_i^*Av_j)v_j^*$ y $v_i^*Av_j\in e\mathcal Ne=\mathbb Ce$,
luego $A=\sum_{ij}a_{ij}E_{ij}$ con $a_{ij}\in\mathbb C$. Así, $\mathcal N=\operatorname{span}\{E_{ij}\}$.

*Paso 4: la unitaria.* Sea $m=\operatorname{rank}e$ y $f_1,\dots,f_m$ una base ortonormal
de $eE$. Los vectores $v_if_r$ son ortonormales:
$\langle v_if_r,v_jf_s\rangle=\langle f_r,v_i^*v_jf_s\rangle=\delta_{ij}\langle f_r,ef_s\rangle=\delta_{ij}\delta_{rs}$.
Generan $E$, ya que $E=\bigoplus_ie_iE$ y $e_iE=v_ieE$. Definimos $V(|i\rangle\otimes f_r)=v_if_r$,
donde $|1\rangle,\dots,|k\rangle$ es la base canónica de $\mathbb C^k$ e identificamos
$\mathbb C^m$ con $eE$ mediante los $f_r$. Entonces $E_{ij}v_lf_r=\delta_{jl}v_if_r$, es decir,
$V^*E_{ij}V=|i\rangle\langle j|\otimes1_m$, y $V^*\mathcal NV=M_k(\mathbb C)\otimes1_m$.

*Paso 5: el conmutante.* Escribamos $X\in M_k\otimes M_m$ como matriz por bloques
$X=\sum_{ij}|i\rangle\langle j|\otimes X_{ij}$. Para $a,b$ fijos,
$(|a\rangle\langle b|\otimes1)X=\sum_j|a\rangle\langle j|\otimes X_{bj}$ y
$X(|a\rangle\langle b|\otimes1)=\sum_i|i\rangle\langle b|\otimes X_{ia}$. Si son iguales, al
comparar los bloques en la posición $(a,j)$ con $j\neq b$ se obtiene $X_{bj}=0$, y en la
posición $(a,b)$, $X_{bb}=X_{aa}$. Como $a,b$ son arbitrarios, $X=1_k\otimes Y$ con
$Y=X_{11}$; recíprocamente, todo $1_k\otimes Y$ conmuta con $M_k\otimes1_m$. $\square$

**Teorema 2.6 (estructura de las álgebras de von Neumann de dimensión finita).** Sea
$\mathcal M$ un álgebra de von Neumann sobre $\mathbb C^n$. Existen $C\ge1$, enteros
$k_c,m_c\ge1$ con $\sum_ck_cm_c=n$, y una unitaria
$V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n$ tales que
$$
V^*\mathcal MV=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},\qquad
V^*\mathcal M'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),
$$
y $V^*Z(\mathcal M)V=\operatorname{span}\{P_1,\dots,P_C\}$, donde $P_c$ es la proyección
sobre el $c$-ésimo sumando. Los pares $(k_c,m_c)$ son únicos salvo el orden: las $P_c$ son
únicas (lema 2.3), $k_c^2=\dim\mathcal MP_c$ y $m_c=\operatorname{rank}P_c/k_c$. En particular,
$\mathcal M\cong\bigoplus_cM_{k_c}(\mathbb C)$ como $*$-álgebras, y
$\dim\mathcal M=\sum_ck_c^2$, $\dim\mathcal M'=\sum_cm_c^2$.

*Demostración.* Por el lema 2.3 y la observación que lo sigue, $\mathcal M=\bigoplus_c\mathcal M_c$
con $\mathcal M_c$ actuando sobre $E_c$, y cada $\mathcal M_c$ es un factor (lema 2.4). El
lema 2.5 proporciona unitarias $V_c$ con $V_c^*\mathcal M_cV_c=M_{k_c}\otimes1_{m_c}$; sea
$V=\bigoplus_cV_c$. Para el conmutante: $X\in\mathcal M'$ conmuta con cada $P_c$, luego
$X=\sum_cP_cXP_c$, y $P_cXP_c$ conmuta con $\mathcal M_c$; recíprocamente, toda suma de este
tipo conmuta con $\mathcal M$. El lema 2.5 identifica $\mathcal M_c'$ sobre $E_c$. $\square$

Los números $k_c$ describen el álgebra en sí. Las **multiplicidades** $m_c$ describen cómo
está situada dentro de $M_n(\mathbb C)$: son invisibles para la estructura algebraica, pero
determinan la traza de $\mathbb C^n$ restringida a $\mathcal M$,
$$
\operatorname{Tr}_{\mathbb C^n}(A)=\sum_cm_c\operatorname{Tr}(A_c),\qquad A=V\Big(\bigoplus_cA_c\otimes1_{m_c}\Big)V^*.
$$
Por eso la biblioteca trabaja con una familia de trazas (capítulo 3).

> **Significado físico.** Las $P_c$ son **sectores de superselección** o **sectores de
> simetría**: ningún observable de $\mathcal M$ conecta sectores distintos, de modo que la
> medición de "en qué sector" es compatible con cualquier otra medición. El centro es la
> parte *clásica* del sistema. La multiplicidad $m_c$ es una degeneración: cada valor propio
> del bloque $A_c$ aparece en $\mathbb C^n$ con $m_c$ veces su multiplicidad en $A_c$.

## 2.4 Ejemplo: un grupo de simetría y su conmutante

Sea un grupo finito $G$ que actúa sobre $\mathbb C^n$ mediante una representación unitaria
$\pi$. Los operadores que conmutan con la simetría forman $\pi(G)'$, un álgebra de von
Neumann; su conmutante $\pi(G)''$ es el subespacio generado por los $\pi(g)$ (teorema 2.2
aplicado a dicho subespacio, que es una $*$-álgebra: contiene $\pi(e)=1$, los productos
$\pi(g)\pi(h)=\pi(gh)$ y los adjuntos $\pi(g)^*=\pi(g^{-1})$).

**Proposición 2.8.** Apliquemos el teorema 2.6 a $\mathcal M=\pi(G)''$. Entonces
$V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$, donde cada $\rho_c$ es una representación
unitaria irreducible de $G$ sobre $\mathbb C^{k_c}$, y $\rho_c,\rho_d$ no son equivalentes
para $c\neq d$. En consecuencia,
$$
\pi\cong\bigoplus_c\rho_c\otimes1_{m_c},\qquad
\pi(G)''\cong\bigoplus_cM_{\dim\rho_c},\qquad
\pi(G)'\cong\bigoplus_cM_{m_c}:
$$
los bloques de los operadores simétricos $\pi(G)'$ tienen como tamaños las multiplicidades
de las representaciones irreducibles, y los papeles de $k_c$ y $m_c$ se intercambian entre
un álgebra y su conmutante.

*Demostración.* $\pi(g)\in\mathcal M$, luego $V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}$
para ciertas $\rho_c(g)\in M_{k_c}$, y $g\mapsto\rho_c(g)$ es una representación unitaria
porque $\pi$ lo es. Como los $\pi(g)$ generan $\mathcal M$, los $\rho_c(g)$ generan $M_{k_c}$
y, más aún, la aplicación $\sum_ga_g\pi(g)\mapsto\big(\sum_ga_g\rho_c(g)\big)_c$ es
sobreyectiva sobre $\bigoplus_cM_{k_c}$. *Irreducibilidad:* un subespacio de
$\mathbb C^{k_c}$ invariante por todas las $\rho_c(g)$ es invariante por el subespacio que
generan, $M_{k_c}$, luego es $0$ o $\mathbb C^{k_c}$. *No equivalencia:* si
$\rho_d(g)=T\rho_c(g)T^{-1}$ para todo $g$ y alguna $T$ invertible, con $c\neq d$, entonces
todo elemento $(x_{c'})_{c'}$ de la imagen satisface $x_d=Tx_cT^{-1}$, lo que excluye el
elemento con $x_c=1$, $x_d=0$, en contradicción con la sobreyectividad. Los dos últimos
isomorfismos son el teorema 2.6. $\square$

Que la descomposición de una representación en irreducibles sea única salvo
equivalencia —de modo que las $m_c$ sean *las* multiplicidades en cualquier
descomposición— se demuestra en el capítulo 4.

Tomamos $G=S_3$ actuando por permutación de los tres factores de $(\mathbb C^2)^{\otimes3}$,
es decir, tres cúbits o tres espines $\tfrac12$. La teoría de representaciones (un cálculo
de caracteres, capítulo 4) predice: el subespacio simétrico (espín total $\tfrac32$,
dimensión $4$) porta la representación trivial con multiplicidad $4$; su complemento
(dimensión $4$) porta la representación estándar de dimensión dos con multiplicidad $2$; la
representación signo no aparece. Así, $\pi(S_3)''\cong M_1\oplus M_2$ (dimensión $5$) y
$\pi(S_3)'\cong M_4\oplus M_2$ (dimensión $20$). Comprobamos las dimensiones mediante
álgebra lineal: el conmutante de un conjunto $S$ es el núcleo de la aplicación lineal
$X\mapsto(AX-XA)_{A\in S}$.

```python
import itertools
import torch

torch.manual_seed(0)
dt = torch.complex128

def permutation_operator(perm, n_qubits=3):
    """Unitary on (C^2)^{(x)3} that permutes the tensor factors."""
    dim = 2 ** n_qubits
    P = torch.zeros(dim, dim, dtype=dt)
    for bits in itertools.product([0, 1], repeat=n_qubits):
        src = int(''.join(map(str, bits)), 2)
        dst = int(''.join(str(bits[perm[i]]) for i in range(n_qubits)), 2)
        P[dst, src] = 1
    return P

group = [permutation_operator(p) for p in itertools.permutations(range(3))]

def commutant_basis(ops, tol=1e-10):
    """Orthonormal basis (Hilbert-Schmidt) of {X : X A = A X for all A in ops}."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    # vec(AX - XA) = (A (x) 1 - 1 (x) A^T) vec(X) for row-major vec
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    _, s, Vh = torch.linalg.svd(L)
    null = Vh[(s > tol).sum():].conj().resolve_conj()   # rows span the null space
    return [v.reshape(d, d) for v in null]

def span_dim(ops, tol=1e-10):
    return int((torch.linalg.svdvals(torch.stack([A.reshape(-1) for A in ops])) > tol).sum())

comm = commutant_basis(group)
print(len(comm), span_dim(group), len(commutant_basis(comm)))
assert len(comm) == 20                                   # pi(S3)'  = M_4 (+) M_2
assert span_dim(group) == 5                              # pi(S3)'' = M_1 (+) M_2
assert len(commutant_basis(comm)) == 5                   # double commutant theorem
```

Las multiplicidades se observan en el espectro. El promedio sobre el grupo,
$X\mapsto\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$, lleva cualquier operador al conmutante
(ejercicio 3). Por tanto, un elemento autoadjunto aleatorio de $\pi(S_3)'$ tiene (casi
seguramente) $4$ valores propios simples en el subespacio simétrico y $2$ doblemente
degenerados en su complemento; el código comprueba los pares degenerados: el factor $1_2$
de $M_2\otimes1_2$.

```python
def group_average(X):
    return sum(g @ X @ g.conj().T for g in group) / len(group)

H0 = torch.randn(8, 8, dtype=dt)
Hsym = group_average(H0 + H0.conj().T)                  # a random symmetric Hamiltonian
assert all(torch.allclose(g @ Hsym, Hsym @ g) for g in group)

P_sym = sum(group) / len(group)                         # projection onto the symmetric subspace
w, Q = torch.linalg.eigh(P_sym)
Q_rest = Q[:, w.real < 0.5]                             # orthonormal basis of the complement (dim 4)
e_rest = torch.linalg.eigvalsh(Q_rest.conj().T @ Hsym @ Q_rest)
print(e_rest)
assert Q_rest.shape[1] == 4
assert torch.allclose(e_rest[0], e_rest[1]) and torch.allclose(e_rest[2], e_rest[3])
```

## 2.5 El álgebra en la biblioteca

La biblioteca almacena el álgebra abstracta $\bigoplus_cM_{k_c}(\mathbb C)$ del teorema 2.6
sobre $\bigoplus_c\mathbb C^{k_c}$, sin las multiplicidades: un operador es la lista de sus
bloques $A_c$, rellenados hasta un tamaño común. El conmutante del ejemplo es
`TypeIAlgebra([4, 2], [4, 2])`.

```python
from torch_vn_algebra import TypeIAlgebra, channels

alg = TypeIAlgebra([4, 2], [4, 2], precision='double')    # M_4 (+) M_2
print(alg.k_factors, alg.C)

A = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
print(A.matrix.shape)              # (1, 2, 4, 4): block 1 occupies the top-left 2 x 2 corner
assert A.matrix[0, 1, 2:, :].abs().max() == 0 and A.matrix[0, 1, :, 2:].abs().max() == 0
```

El centro está formado por los elementos $\sum_cz_cP_c$; `alg.central` los construye, y
conmutan con todo:

```python
def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

P0, P1 = alg.central([1.0, 0.0]), alg.central([0.0, 1.0])   # minimal central projections
Z = alg.central([2.0, -1.0])
assert P0.is_projection and P1.is_projection
assert close(P0 + P1, alg.identity())
assert close(Z @ A, A @ Z)
assert close(A @ P0 + A @ P1, A)                         # A = sum_c A P_c
```

## 2.6 Esperanzas condicionales

Una **esperanza condicional** de $\mathcal M$ sobre una subálgebra de von Neumann
$\mathcal N\subseteq\mathcal M$ es una aplicación lineal $E:\mathcal M\to\mathcal N$ que es
positiva ($A\ge0\Rightarrow E(A)\ge0$), satisface $E(B)=B$ para $B\in\mathcal N$ y es un
morfismo de $\mathcal N$-bimódulos: $E(B_1AB_2)=B_1E(A)B_2$ para $B_1,B_2\in\mathcal N$.

**Proposición 2.7.** Para $\mathcal M=\bigoplus_cM_{k_c}$ y $\mathcal N=Z(\mathcal M)$,
$$
E(A)=\sum_c\frac{\operatorname{Tr}A_c}{k_c}\,P_c
$$
es una esperanza condicional. Conserva toda traza de la forma
$\sum_cw_c\operatorname{Tr}A_c$.

*Demostración.* La linealidad es evidente; $A\ge0$ implica $\operatorname{Tr}A_c\ge0$, luego
$E(A)\ge0$. Para $B=\sum_cb_cP_c$, $E(B)=\sum_cb_ck_cP_c/k_c=B$. Para $B_i=\sum_cb^{(i)}_cP_c$,
$\operatorname{Tr}(B_1AB_2)_c=b^{(1)}_cb^{(2)}_c\operatorname{Tr}A_c$, lo que da la propiedad
de bimódulo. Por último, $\operatorname{Tr}(E(A)_c)=\operatorname{Tr}A_c$ para todo $c$. $\square$

En la biblioteca, esta aplicación es el canal `channels.center_expectation(alg)`. Al
aplicarla a operadores se obtiene su parte clásica (promediada en cada sector):

```python
E = channels.center_expectation(alg)
EA = E(A)
assert close(E(EA), EA)                                  # idempotent
assert close(E(Z @ A @ P1), Z @ EA @ P1)                 # bimodule property
A2 = alg.from_blocks([torch.randn(4, 4, dtype=dt), torch.randn(2, 2, dtype=dt)])
assert close(EA @ A2, A2 @ EA)                           # E(A) commutes with a fresh element
for c in range(alg.C):                                   # Tr E(A)_c = Tr A_c
    assert abs(EA.matrix[0, c].diagonal().sum() - A.matrix[0, c].diagonal().sum()) < 1e-10
```

El promedio sobre el grupo de la sección 2.4 es una esperanza condicional sobre el
conmutante $\pi(G)'$ dentro de $M_n(\mathbb C)$, y conserva $\operatorname{Tr}$ (ejercicio 3).
(Las esperanzas condicionales no son únicas en general: sobre $\mathbb C1$, cada estado
$\varphi$ da una, $A\mapsto\varphi(A)1$.) Ambas son casos particulares de la misma idea:
proyectar sobre los operadores que "no ven" cierta estructura.

## 2.7 Con qué comparar

Los ejemplos de la sección 2.4 pueden comprobarse con cualquier paquete de álgebra lineal;
el software de teoría de representaciones (GAP, SageMath) da la descomposición
$\pi\cong\bigoplus\rho_\lambda\otimes1_{m_\lambda}$ de forma simbólica. Los códigos de
diagonalización exacta en física (QuSpin, los módulos de simetría de QuTiP) usan la misma
estructura de bloques para diagonalizar hamiltonianos sector a sector. La biblioteca toma
la estructura como punto de partida: una vez conocidos los bloques, todos los cálculos
posteriores —funciones, trazas, estados, canales, dinámica— permanecen dentro del álgebra.

## Ejercicios

1. Demuestre que el centro de $M_k(\mathbb C)$ es $\mathbb C1$ (use las unidades matriciales
   $|i\rangle\langle j|$) y que $Z(\bigoplus_cM_{k_c})$ tiene dimensión $C$.
2. Sea $S$ el desplazamiento cíclico en $\mathbb C^3$, $Se_i=e_{i+1\bmod3}$. Determine
   $\{S\}'$ y $\{S\}''$ y compruebe numéricamente las dimensiones.
3. Para una representación unitaria $\pi$ de un grupo finito, demuestre que
   $E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*$ es una esperanza condicional de
   $M_n(\mathbb C)$ sobre $\pi(G)'$ que conserva $\operatorname{Tr}$.

<details>
<summary>Soluciones</summary>

1. Si $X$ conmuta con todas las $|i\rangle\langle j|$, entonces
   $X|i\rangle\langle j|=|i\rangle\langle j|X$; aplicando ambos miembros a $e_j$ y comparando
   se obtiene $Xe_i=X_{jj}e_i$ para todos los $i,j$, luego $X=X_{11}1$. En una suma directa,
   un elemento $\bigoplus X_c$ es central si y solo si lo es cada $X_c$, luego el centro es
   $\{\bigoplus_cz_c1_{k_c}\}$, de dimensión $C$.

2. $S$ es unitaria con valores propios $1,\omega,\omega^2$, $\omega=e^{2\pi i/3}$, todos
   simples. Un operador que conmuta con $S$ conserva sus subespacios propios, de modo que
   en la base propia es diagonal: $\{S\}'\cong\mathbb C^3$, las matrices circulantes.
   Entonces $\{S\}''$ también es $\cong\mathbb C^3$, porque en la base propia de $S$ una
   matriz que conmuta con todas las matrices diagonales es diagonal (compárense las
   entradas de $XD$ y $DX$ para una $D$ con entradas diagonales distintas); está generado
   por $1,S,S^2$ (Vandermonde: los valores propios son distintos).

```python
S = torch.roll(torch.eye(3, dtype=dt), 1, dims=0)
assert len(commutant_basis([S])) == 3
assert len(commutant_basis(commutant_basis([S]))) == 3
assert span_dim([torch.eye(3, dtype=dt), S, S @ S]) == 3
```

3. $E$ es lineal y positiva (una suma de aplicaciones $X\mapsto UXU^*$). Para $h\in G$,
   $\pi(h)E(X)\pi(h)^*=\frac1{|G|}\sum_g\pi(hg)X\pi(hg)^*=E(X)$, luego $E(X)\in\pi(G)'$. Si
   $X\in\pi(G)'$, cada término es igual a $X$, luego $E(X)=X$. Para $B_1,B_2\in\pi(G)'$,
   $\pi(g)B_1XB_2\pi(g)^*=B_1\pi(g)X\pi(g)^*B_2$, lo que da la propiedad de bimódulo.
   Por último, $\operatorname{Tr}(\pi(g)X\pi(g)^*)=\operatorname{Tr}X$ para todo $g$, luego
   $\operatorname{Tr}E(X)=\operatorname{Tr}X$.

</details>

**Siguiente:** el [capítulo 3](03_traces_and_states.md) clasifica las trazas sobre
$\bigoplus_cM_{k_c}$, introduce los estados y las densidades, y demuestra las desigualdades
de entropía y de trazas que se usan más adelante.
