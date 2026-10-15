[English](../02_von_neumann_algebras.md) · **Español** · [Français](../fr/02_von_neumann_algebras.md) · [Deutsch](../de/02_von_neumann_algebras.md) · [中文](../zh/02_von_neumann_algebras.md) · [日本語](../ja/02_von_neumann_algebras.md) · [Русский](../ru/02_von_neumann_algebras.md)

# 2. Álgebras de von Neumann de dimensión finita

Un álgebra de von Neumann es un álgebra de operadores que es cerrada para la operación de
"tomar dos veces el conmutante". En dimensión finita se trata de una condición puramente
algebraica, y el teorema principal de este capítulo afirma que toda álgebra de este tipo es,
salvo un cambio de base unitario, una suma directa de álgebras de matrices completas, cada
una repetida con cierta multiplicidad. Esta es la estructura sobre la que está construida la
biblioteca.

*En la biblioteca:* el álgebra $`\bigoplus_cM_{k_c}(\mathbb C)`$ es la clase `TypeIAlgebra` — [tutorial, sección 2.5](../../tutorial/es/02_von_neumann_algebras.md#s-2-5).

**Aprenderá** el teorema del doble conmutante, el teorema de estructura con demostración
completa, cómo un grupo de simetría da lugar a un álgebra de este tipo, qué es el centro y
qué es una esperanza condicional. Las comprobaciones numéricas de este capítulo están en el
[Tutorial 2](../../tutorial/es/02_von_neumann_algebras.md).

En todo el capítulo, *$`\ast`$-álgebra sobre* $`\mathbb C^n`$ significa un subespacio lineal complejo
$`\mathcal M\subseteq M_n(\mathbb C)`$ cerrado para productos y adjuntos y que contiene la
identidad $`1`$.

<a id="s-2-1"></a>

## 2.1 Conmutantes y el teorema del doble conmutante

Para un conjunto $`S\subseteq M_n(\mathbb C)`$, el **conmutante** es
$`S'=\{X\in M_n(\mathbb C):XA=AX\text{ para todo }A\in S\}`$. Siempre es un álgebra unital;
si $`S`$ es cerrado para adjuntos ($`S^*=S`$), entonces $`S'`$ es una $`*`$-álgebra, ya que
$`XA=AX`$ para todo $`A\in S`$ implica $`X^*A^*=A^*X^*`$ para todo $`A\in S`$. Todo elemento de $`S`$
conmuta con todo elemento de $`S'`$, luego $`S\subseteq S''`$; aplicado a $`S'`$, esto da
$`S'\subseteq S'''`$, y tomando conmutantes en $`S\subseteq S''`$ (lo cual invierte las
inclusiones) se obtiene $`S'''\subseteq S'`$. Por tanto $`S'=S'''`$.

*Cálculo:* los conmutantes como núcleos, numéricamente — [tutorial, sección 2.4](../../tutorial/es/02_von_neumann_algebras.md#s-2-4). *En la biblioteca:* `decompose.commutant_dimension` — [tutorial, sección 4.2](../../tutorial/es/04_symmetries.md#s-4-2); `decompose.commutant_basis` — [tutorial, sección 4.3](../../tutorial/es/04_symmetries.md#s-4-3); qué tolerancia del núcleo es correcta para matrices con ruido, teorema 4.25 — [tutorial, sección 4.9](../../tutorial/es/04_symmetries.md#s-4-9).

Cuando un álgebra actúa sobre un subespacio $`E`$ de $`\mathbb C^n`$, identificamos $`E`$ con
$`\mathbb C^{\dim E}`$ mediante una base ortonormal, de modo que todos los resultados
enunciados sobre $`\mathbb C^n`$ se le aplican.

**Definición.** Un **álgebra de von Neumann** sobre $`\mathbb C^n`$ es una $`*`$-álgebra
$`\mathcal M`$ con $`\mathcal M=\mathcal M''`$. Su **centro** es
$`Z(\mathcal M)=\mathcal M\cap\mathcal M'`$. Un **factor** es un álgebra de von Neumann con
$`Z(\mathcal M)=\mathbb C1`$.

La definición está pensada para operadores sobre espacios de Hilbert de dimensión infinita, que
este curso no trata. En dimensión finita la condición $`\mathcal M=\mathcal M''`$ se cumple automáticamente:

<a id="r-2-1"></a>

**Lema 2.1.** Sea $`\mathcal M`$ una $`*`$-álgebra sobre $`\mathbb C^n`$ y $`K\subseteq\mathbb C^n`$
un subespacio con $`\mathcal MK\subseteq K`$. Entonces la proyección ortogonal $`P`$ sobre $`K`$
pertenece a $`\mathcal M'`$.

*Demostración.* Para $`A\in\mathcal M`$, $`AK\subseteq K`$ y $`A^*K\subseteq K`$ (pues $`A^*\in\mathcal M`$).
La segunda inclusión da $`AK^\perp\subseteq K^\perp`$: para $`y\in K^\perp`$, $`x\in K`$,
$`\langle x,Ay\rangle=\langle A^*x,y\rangle=0`$. Por tanto $`A`$ aplica $`K`$ en $`K`$ y $`K^\perp`$
en $`K^\perp`$, lo que significa que $`AP=PAP=PA`$. $`\square`$

<a id="r-2-2"></a>

**Teorema 2.2 (teorema del doble conmutante).** Toda $`*`$-álgebra $`\mathcal M`$ sobre
$`\mathbb C^n`$ satisface $`\mathcal M=\mathcal M''`$.

*Demostración.* Se tiene $`\mathcal M\subseteq\mathcal M''`$; sea $`T\in\mathcal M''`$.

*Paso 1: un vector.* Para $`x\in\mathbb C^n`$, el subespacio $`K=\mathcal Mx`$ es invariante
por $`\mathcal M`$ y contiene a $`x=1x`$. Por el lema 2.1, su proyección $`P`$ está en $`\mathcal M'`$,
luego $`TP=PT`$ y $`Tx=TPx=PTx\in K`$: existe $`A\in\mathcal M`$ con $`Tx=Ax`$.

*Paso 2: todos los vectores a la vez.* En $`\mathbb C^n\otimes\mathbb C^n=(\mathbb C^n)^{\oplus n}`$
consideramos la $`*`$-álgebra $`\tilde{\mathcal M}=\{A\oplus\dots\oplus A:A\in\mathcal M\}`$. Una
matriz por bloques $`[X_{ij}]_{i,j=1}^n`$ conmuta con todo $`A\oplus\dots\oplus A`$ si y solo si
$`X_{ij}A=AX_{ij}`$ para todos los $`i,j`$, luego $`\tilde{\mathcal M}'=\{[X_{ij}]:X_{ij}\in\mathcal M'\}`$.
El operador $`\tilde T=T\oplus\dots\oplus T`$ conmuta con toda $`[X_{ij}]`$ de este tipo, porque
$`TX_{ij}=X_{ij}T`$. Por tanto $`\tilde T\in\tilde{\mathcal M}''`$. Apliquemos el paso 1 a
$`\tilde{\mathcal M}`$, $`\tilde T`$ y el vector $`\xi=(e_1,\dots,e_n)`$: existe $`A\in\mathcal M`$ con
$`(Te_1,\dots,Te_n)=(Ae_1,\dots,Ae_n)`$, es decir, $`T=A\in\mathcal M`$. $`\square`$

*Cálculo:* el teorema 2.2 en un ejemplo — [tutorial, sección 2.4](../../tutorial/es/02_von_neumann_algebras.md#s-2-4).

Así pues, en dimensión finita "álgebra de von Neumann" y "$`*`$-álgebra unital de matrices"
son lo mismo. El teorema del doble conmutante ofrece, no obstante, el punto de vista
adecuado: en física, un álgebra de observables viene *dada* muy a menudo como un
conmutante: los operadores que conmutan con una simetría (sección 2.4).

<a id="s-2-2"></a>

## 2.2 El centro

<a id="r-2-3"></a>

**Lema 2.3.** Sea $`\mathcal Z`$ una $`*`$-álgebra conmutativa sobre $`\mathbb C^n`$. Existen
proyecciones $`P_1,\dots,P_C\in\mathcal Z`$ no nulas, ortogonales dos a dos, con
$`\sum_cP_c=1`$, tales que $`\mathcal Z=\mathrm{span}\{P_1,\dots,P_C\}`$. Son las
proyecciones minimales no nulas de $`\mathcal Z`$ y son únicas salvo el orden.

*Demostración.* Todo $`z\in\mathcal Z`$ es normal ($`z^*\in\mathcal Z`$ conmuta con $`z`$). Sean
$`z_1,\dots,z_r`$ generadores de $`\mathcal Z`$ como espacio vectorial. Descomponemos
$`\mathbb C^n`$ en subespacios propios comunes. Los subespacios propios de $`z_1`$ son
ortogonales entre sí y generan $`\mathbb C^n`$ (teorema 1.3); son invariantes por
$`z_2,\dots,z_r`$ y por sus adjuntos (que conmutan con $`z_1`$), de modo que cada $`z_j`$
restringido a un subespacio propio de $`z_1`$ es de nuevo normal, y podemos descomponer cada
subespacio propio de $`z_1`$ en subespacios propios de $`z_2`$, y así sucesivamente. El
resultado es una descomposición ortogonal $`\mathbb C^n=E_1\oplus\dots\oplus E_C`$ en
subespacios no nulos sobre los que cada $`z\in\mathcal Z`$ actúa como un escalar $`\chi_c(z)`$,
donde hemos reunido los subespacios con la misma función $`\chi_c:\mathcal Z\to\mathbb C`$, de
modo que $`\chi_c\neq\chi_d`$ para $`c\neq d`$. Sea $`P_c`$ la proyección sobre $`E_c`$. Entonces
$`z=\sum_c\chi_c(z)P_c`$ para todo $`z\in\mathcal Z`$.

Falta ver que $`P_c\in\mathcal Z`$. Para $`d\neq c`$ elegimos $`z_{cd}\in\mathcal Z`$ con
$`\chi_c(z_{cd})\neq\chi_d(z_{cd})`$. El producto
$`\prod_{d\neq c}\big(z_{cd}-\chi_d(z_{cd})1\big)/\big(\chi_c(z_{cd})-\chi_d(z_{cd})\big)`$
pertenece a $`\mathcal Z`$, actúa como $`1`$ sobre $`E_c`$ y como $`0`$ sobre cada $`E_d`$, $`d\neq c`$;
es $`P_c`$ (para $`C=1`$, $`P_1=1`$: el producto vacío es $`1`$). Minimalidad y unicidad: una proyección $`\sum_c a_cP_c`$ tiene $`a_c\in\{0,1\}`$, luego
las proyecciones no nulas de $`\mathcal Z`$ son las sumas de conjuntos no vacíos de las $`P_c`$. $`\square`$

Aplicadas a $`\mathcal Z=Z(\mathcal M)`$, las proyecciones $`P_c`$ son las **proyecciones
centrales minimales** de $`\mathcal M`$. Todo $`A\in\mathcal M`$ se descompone como
$`A=\sum_cAP_c`$ con $`AP_c=P_cAP_c`$, luego
$`\mathcal M=\bigoplus_c\mathcal MP_c`$ (una suma directa ortogonal de álgebras que actúan sobre
$`E_c=P_c\mathbb C^n`$).

*En la biblioteca:* los elementos centrales $`\sum_cz_cP_c`$, entre ellos las proyecciones centrales minimales, se construyen con `TypeIAlgebra.central` — [tutorial, sección 2.5](../../tutorial/es/02_von_neumann_algebras.md#s-2-5).

<a id="r-2-4"></a>

**Lema 2.4.** Cada $`\mathcal M_c=\mathcal MP_c`$, como álgebra sobre $`E_c`$ con unidad $`P_c`$, es
un factor.

*Demostración.* Sea $`y\in\mathcal M_c`$ que conmuta con $`\mathcal M_c`$. Para $`A\in\mathcal M`$,
$`Ay=AP_cy=(AP_c)y=y(AP_c)=yP_cA=yA`$, usando que $`P_c`$ es central y que $`y=yP_c=P_cy`$.
Así $`y\in Z(\mathcal M)`$, $`y=\sum_d\chi_d(y)P_d`$, e $`y=yP_c=\chi_c(y)P_c`$. $`\square`$

<a id="s-2-3"></a>

## 2.3 El teorema de estructura

<a id="r-2-5"></a>

**Lema 2.5 (factores).** Sea $`\mathcal N`$ un factor sobre un espacio $`E`$. Existen enteros
$`k,m\ge1`$ y una unitaria $`V:\mathbb C^k\otimes\mathbb C^m\to E`$ con
$`V^*\mathcal NV=M_k(\mathbb C)\otimes1_m`$. Además, $`V^*\mathcal N'V=1_k\otimes M_m(\mathbb C)`$.

*Demostración.* *Paso 1: una proyección minimal.* $`\mathcal N`$ contiene proyecciones no nulas
(por ejemplo, $`1`$). Sea $`e\in\mathcal N`$ una proyección no nula de rango mínimo. Afirmamos que
$`e\mathcal Ne=\mathbb Ce`$. Sea $`x=x^*\in e\mathcal Ne`$. Cada proyección espectral $`q`$ de $`x`$
asociada a un valor propio no nulo $`\lambda`$ es $`q=p(x)`$, donde $`p`$ es el polinomio de
Lagrange con $`p(\lambda)=1`$ y $`p(\mu)=0`$ para $`\mu\in(\mathrm{spec}\,x\cup\{0\})\setminus\{\lambda\}`$
(teorema 1.5(a),(b): $`p`$ coincide sobre $`\mathrm{spec}\,x`$ con la función indicatriz de $`\{\lambda\}`$).
Como $`p(0)=0`$, $`p(t)=t\,r(t)`$ para un polinomio $`r`$, luego $`q=x\,r(x)=r(x)\,x\in\mathcal N`$, y $`x=ex=xe`$ da
$`q=eq=qe`$; así, la imagen de $`q`$ está contenida en la imagen de $`e`$. Si $`q\neq0`$, la minimalidad da
$`\mathrm{rank}\,q\ge\mathrm{rank}\,e`$, de modo que ambas imágenes coinciden y $`q=e`$. Como las
proyecciones espectrales de valores propios distintos son ortogonales, $`x`$ tiene a lo sumo un valor propio
no nulo $`\lambda`$, con proyección espectral $`e`$, luego $`x=\lambda e`$ (y $`x=0=0\cdot e`$ si no hay ninguno). Un elemento general de $`e\mathcal Ne`$ es
$`x_1+ix_2`$ con $`x_{1,2}=x_{1,2}^*\in e\mathcal Ne`$.

*Paso 2: $`e`$ genera $`E`$.* Sea $`F=\mathrm{span}\,\mathcal Ne E`$. Es invariante por
$`\mathcal N`$, luego su proyección $`p_F`$ pertenece a $`\mathcal N'`$ (lema 2.1). También es
invariante por la $`*`$-álgebra $`\mathcal N'`$: para $`b\in\mathcal N'`$, $`bAe\xi=Aeb\xi\in F`$; así,
por el lema 2.1 aplicado a $`\mathcal N'`$, $`p_F\in\mathcal N''=\mathcal N`$ (teorema 2.2). Por
tanto $`p_F\in Z(\mathcal N)=\mathbb C1`$, y $`p_F\neq0`$ (pues $`eE\subseteq F`$), luego $`F=E`$.

*Paso 3: unidades matriciales.* Llamamos *$`e`$-isometría* a un $`v\in\mathcal N`$ con $`v^*v=e`$;
entonces $`vv^*`$ es una proyección de $`\mathcal N`$. (En primer lugar, $`v=ve`$: para todo $`x`$,
$`\|v(1-e)x\|^2=\langle(1-e)x,v^*v(1-e)x\rangle=\langle(1-e)x,e(1-e)x\rangle=0`$. Por tanto
$`vv^*v=ve=v`$ y $`(vv^*)^2=vev^*=vv^*`$. Para las $`v_i`$ de más abajo esto da $`e_iv_i=v_i`$ y
$`v_i^*e_i=v_i^*`$.) Elegimos $`e`$-isometrías $`v_1=e,v_2,\dots,v_k`$ con imágenes
$`e_i=v_iv_i^*`$ ortogonales entre sí, con $`k`$ maximal (las imágenes son subespacios
ortogonales de la misma dimensión $`\mathrm{rank}\,e`$, luego
$`k\le\dim E/\mathrm{rank}\,e`$). Afirmamos que $`\sum_ie_i=1`$. En caso contrario,
$`g=1-\sum_ie_i\neq0`$ es una proyección de $`\mathcal N`$. Por el paso 2, $`gAe\neq0`$ para algún
$`A\in\mathcal N`$ (si $`gAe=0`$ para todo $`A`$, entonces $`gAe\xi=0`$ para todos los $`A,\xi`$, luego
$`g`$ se anula sobre $`F=E`$). Sea $`y=gAe`$. Entonces $`y^*y\in e\mathcal Ne`$ es positivo y no nulo,
luego $`y^*y=\lambda e`$ con $`\lambda>0`$, y $`v=y/\sqrt\lambda`$ es una $`e`$-isometría con imagen
contenida en $`gE`$, ortogonal a todas las $`e_i`$, lo que contradice la maximalidad.

Para $`i\neq j`$, $`v_i^*v_j=v_i^*e_ie_jv_j=0`$. Sea $`E_{ij}=v_iv_j^*`$. Entonces
$`E_{ij}E_{kl}=v_iv_j^*v_kv_l^*=\delta_{jk}v_iev_l^*=\delta_{jk}E_{il}`$,
$`E_{ij}^*=E_{ji}`$ y $`\sum_iE_{ii}=1`$. Para $`A\in\mathcal N`$,
$`A=\sum_{ij}e_iAe_j=\sum_{ij}v_i(v_i^*Av_j)v_j^*`$ y $`v_i^*Av_j\in e\mathcal Ne=\mathbb Ce`$,
luego $`A=\sum_{ij}a_{ij}E_{ij}`$ con $`a_{ij}\in\mathbb C`$. Así, $`\mathcal N=\mathrm{span}\{E_{ij}\}`$.

*Paso 4: la unitaria.* Sea $`m=\mathrm{rank}\,e`$ y $`f_1,\dots,f_m`$ una base ortonormal
de $`eE`$. Los vectores $`v_if_r`$ son ortonormales:
$`\langle v_if_r,v_jf_s\rangle=\langle f_r,v_i^*v_jf_s\rangle=\delta_{ij}\langle f_r,ef_s\rangle=\delta_{ij}\delta_{rs}`$.
Generan $`E`$, ya que $`E=\bigoplus_ie_iE`$ y $`e_iE=v_ieE`$. Definimos $`V(|i\rangle\otimes f_r)=v_if_r`$,
donde $`|1\rangle,\dots,|k\rangle`$ es la base canónica de $`\mathbb C^k`$ e identificamos
$`\mathbb C^m`$ con $`eE`$ mediante los $`f_r`$. Entonces $`E_{ij}v_lf_r=\delta_{jl}v_if_r`$, es decir,
$`V^*E_{ij}V=|i\rangle\langle j|\otimes1_m`$, y $`V^*\mathcal NV=M_k(\mathbb C)\otimes1_m`$.

*Paso 5: el conmutante.* Escribamos $`X\in M_k\otimes M_m`$ como matriz por bloques
$`X=\sum_{ij}|i\rangle\langle j|\otimes X_{ij}`$. Para $`a,b`$ fijos,
$`(|a\rangle\langle b|\otimes1)X=\sum_j|a\rangle\langle j|\otimes X_{bj}`$ y
$`X(|a\rangle\langle b|\otimes1)=\sum_i|i\rangle\langle b|\otimes X_{ia}`$. Si son iguales, al
comparar los bloques en la posición $`(a,j)`$ con $`j\neq b`$ se obtiene $`X_{bj}=0`$, y en la
posición $`(a,b)`$, $`X_{bb}=X_{aa}`$. Como $`a,b`$ son arbitrarios, $`X=1_k\otimes Y`$ con
$`Y=X_{11}`$; recíprocamente, todo $`1_k\otimes Y`$ conmuta con $`M_k\otimes1_m`$. $`\square`$

<a id="r-2-6"></a>

**Teorema 2.6 (estructura de las álgebras de von Neumann de dimensión finita).** Sea
$`\mathcal M`$ un álgebra de von Neumann sobre $`\mathbb C^n`$. Existen $`C\ge1`$, enteros
$`k_c,m_c\ge1`$ con $`\sum_ck_cm_c=n`$, y una unitaria
$`V:\bigoplus_c\mathbb C^{k_c}\otimes\mathbb C^{m_c}\to\mathbb C^n`$ tales que
```math
V^*\mathcal MV=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)\otimes1_{m_c},\qquad
V^*\mathcal M'V=\bigoplus_{c=1}^C1_{k_c}\otimes M_{m_c}(\mathbb C),
```
y $`V^*Z(\mathcal M)V=\mathrm{span}\{V^*P_1V,\dots,V^*P_CV\}`$, donde $`P_1,\dots,P_C`$ son
las proyecciones centrales minimales de $`\mathcal M`$ (lema 2.3) y $`V^*P_cV`$ es la proyección
sobre el $`c`$-ésimo sumando. Los pares $`(k_c,m_c)`$ son únicos salvo el orden: para *cualquier*
unitaria $`V`$ con la primera propiedad, $`V(\cdot)V^*`$ lleva las proyecciones sobre los sumandos a las $`P_c`$
(véase la demostración), de modo que $`k_c^2=\dim\mathcal MP_c`$ y $`m_c=\mathrm{rank}\,P_c/k_c`$ están
determinados por $`\mathcal M`$. En particular,
$`\mathcal M\cong\bigoplus_cM_{k_c}(\mathbb C)`$ como $`*`$-álgebras, y
$`\dim\mathcal M=\sum_ck_c^2`$, $`\dim\mathcal M'=\sum_cm_c^2`$.

*Demostración.* Por el lema 2.3 y la observación que lo sigue, $`\mathcal M=\bigoplus_c\mathcal M_c`$
con $`\mathcal M_c`$ actuando sobre $`E_c`$, y cada $`\mathcal M_c`$ es un factor (lema 2.4). El
lema 2.5 proporciona unitarias $`V_c`$ con $`V_c^*\mathcal M_cV_c=M_{k_c}\otimes1_{m_c}`$; sea
$`V=\bigoplus_cV_c`$. Para el conmutante: $`X\in\mathcal M'`$ conmuta con cada $`P_c`$, luego
$`X=\sum_cP_cXP_c`$, y $`P_cXP_c`$ conmuta con $`\mathcal M_c`$; recíprocamente, toda suma de este
tipo conmuta con $`\mathcal M`$. El lema 2.5 identifica $`\mathcal M_c'`$ sobre $`E_c`$, y la conjugación por $`V`$
da la fórmula para $`V^*\mathcal M'V`$.

*Unicidad.* Sea $`V`$ cualquier unitaria con $`V^*\mathcal MV=\mathcal B:=\bigoplus_{c=1}^{C'}M_{k_c}\otimes1_{m_c}`$ y
sea $`Q_c`$ la proyección sobre el $`c`$-ésimo sumando. Un elemento $`\bigoplus_cX_c`$ de $`\mathcal B`$ ($`X_c\in M_{k_c}\otimes1_{m_c}`$)
es central si y solo si cada $`X_c`$ conmuta con $`M_{k_c}\otimes1_{m_c}`$, es decir (lema 2.5, paso 5), $`X_c\in(M_{k_c}\otimes1_{m_c})\cap(1_{k_c}\otimes M_{m_c})=\mathbb C1`$
(un elemento $`a\otimes1=1\otimes b`$ tiene todos los bloques diagonales iguales a $`b`$ y todos los bloques no diagonales nulos, luego $`a`$ es un escalar).
Por tanto $`Z(\mathcal B)=\mathrm{span}\{Q_1,\dots,Q_{C'}\}`$, y $`Z(\mathcal M)=VZ(\mathcal B)V^*`$ está generado por las proyecciones
$`VQ_cV^*`$, no nulas, ortogonales dos a dos y de suma $`1`$. Cada una de ellas es una proyección no nula de $`Z(\mathcal M)`$, luego (lema 2.3)
una suma de un conjunto no vacío $`S_c`$ de las $`P_d`$; los $`S_c`$ son disjuntos (ortogonalidad) y cubren todos los $`d`$ (suma $`1`$), y como
$`C'=\dim Z(\mathcal M)=C`$, cada $`S_c`$ consta de un solo elemento. Así, tras renumerar, $`VQ_cV^*=P_c`$, de donde
$`\dim\mathcal MP_c=\dim(M_{k_c}\otimes1_{m_c})=k_c^2`$ y $`\mathrm{rank}\,P_c=k_cm_c`$. $`\square`$

*En la biblioteca:* el álgebra $`\bigoplus_cM_{k_c}(\mathbb C)`$ es un `TypeIAlgebra` — [tutorial, sección 2.5](../../tutorial/es/02_von_neumann_algebras.md#s-2-5); la unitaria $`V`$ y los pares $`(k_c,m_c)`$ los calcula `decompose.wedderburn` (con las interfaces `from_group`, `from_generators`) — [tutorial, sección 4.3](../../tutorial/es/04_symmetries.md#s-4-3).

Los números $`k_c`$ describen el álgebra en sí. Las **multiplicidades** $`m_c`$ describen cómo
está situada dentro de $`M_n(\mathbb C)`$: son invisibles para la estructura algebraica, pero
determinan la traza de $`\mathbb C^n`$ restringida a $`\mathcal M`$,
```math
\mathrm{Tr}_{\mathbb C^n}(A)=\sum_cm_c\mathrm{Tr}(A_c),\qquad A=V\Big(\bigoplus_cA_c\otimes1_{m_c}\Big)V^*.
```
Por eso la biblioteca trabaja con una familia de trazas (capítulo 3). Del mismo modo, para $`A`$
autoadjunto, la restricción de $`A`$ a $`V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})`$ tiene los valores propios del
bloque $`A_c`$, cada uno con $`m_c`$ veces su multiplicidad en $`A_c`$, ya que $`A_c\otimes1_{m_c}`$ consiste en $`m_c`$
copias de $`A_c`$.

*En la biblioteca:* las trazas `Tr_blunt`, `Tr_norm`, `tau_vN` — [tutorial, sección 3.1](../../tutorial/es/03_traces_and_states.md#s-3-1).

> **Significado físico** (interpretación física, no un teorema de este curso). Las $`P_c`$ son **sectores de superselección** o **sectores de
> simetría**: ningún observable de $`\mathcal M`$ conecta sectores distintos, de modo que la
> medición de "en qué sector" es compatible con cualquier otra medición. El centro es la
> parte *clásica* del sistema. La multiplicidad $`m_c`$ es una degeneración de todo observable de $`\mathcal M`$.

<a id="s-2-4"></a>

## 2.4 Ejemplo: un grupo de simetría y su conmutante

Sea un grupo finito $`G`$ que actúa sobre $`\mathbb C^n`$ mediante una representación unitaria
$`\pi`$. Los operadores que conmutan con la simetría forman $`\pi(G)'`$, un álgebra de von
Neumann; su conmutante $`\pi(G)''`$ es el subespacio generado por los $`\pi(g)`$ (teorema 2.2
aplicado a dicho subespacio, que es una $`*`$-álgebra: contiene $`\pi(e)=1`$, los productos
$`\pi(g)\pi(h)=\pi(gh)`$ y los adjuntos $`\pi(g)^*=\pi(g^{-1})`$).

<a id="r-2-7"></a>

**Proposición 2.7.** Apliquemos el teorema 2.6 a $`\mathcal M=\pi(G)''`$. Entonces
$`V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}`$, donde cada $`\rho_c`$ es una representación
unitaria irreducible de $`G`$ sobre $`\mathbb C^{k_c}`$, y $`\rho_c,\rho_d`$ no son equivalentes
para $`c\neq d`$. En consecuencia,
```math
\pi\cong\bigoplus_c\rho_c\otimes1_{m_c},\qquad
\pi(G)''\cong\bigoplus_cM_{\dim\rho_c},\qquad
\pi(G)'\cong\bigoplus_cM_{m_c}:
```
los bloques de los operadores simétricos $`\pi(G)'`$ tienen como tamaños las multiplicidades
de las representaciones irreducibles, y los papeles de $`k_c`$ y $`m_c`$ se intercambian entre
un álgebra y su conmutante.

*Demostración.* $`\pi(g)\in\mathcal M`$, luego $`V^*\pi(g)V=\bigoplus_c\rho_c(g)\otimes1_{m_c}`$
para ciertas $`\rho_c(g)\in M_{k_c}`$, y $`g\mapsto\rho_c(g)`$ es una representación unitaria
porque $`\pi`$ lo es. Como los $`\pi(g)`$ generan $`\mathcal M`$, los $`\rho_c(g)`$ generan $`M_{k_c}`$
y, más aún, la aplicación $`\sum_ga_g\pi(g)\mapsto\big(\sum_ga_g\rho_c(g)\big)_c`$ es
sobreyectiva sobre $`\bigoplus_cM_{k_c}`$. *Irreducibilidad:* un subespacio de
$`\mathbb C^{k_c}`$ invariante por todas las $`\rho_c(g)`$ es invariante por el subespacio que
generan, $`M_{k_c}`$, luego es $`0`$ o $`\mathbb C^{k_c}`$. *No equivalencia:* si
$`\rho_d(g)=T\rho_c(g)T^{-1}`$ para todo $`g`$ y alguna $`T`$ invertible, con $`c\neq d`$, entonces
todo elemento $`(x_{c'})_{c'}`$ de la imagen satisface $`x_d=Tx_cT^{-1}`$, lo que excluye el
elemento con $`x_c=1`$, $`x_d=0`$, en contradicción con la sobreyectividad. Los dos últimos
isomorfismos son el teorema 2.6. $`\square`$

*En la biblioteca:* la descomposición de una representación de un grupo finito, `decompose.from_group` — [tutorial, sección 4.3](../../tutorial/es/04_symmetries.md#s-4-3).

La descomposición de una representación en irreducibles es única salvo equivalencia, de modo
que las $`m_c`$ son *las* multiplicidades en cualquier descomposición. *Demostración.* Sean
$`\rho,\rho'`$ irreducibles y $`T`$ un operador de entrelazamiento, $`T\rho(g)=\rho'(g)T`$ para todo $`g`$.
El núcleo y la imagen de $`T`$ son invariantes, luego $`T=0`$ o $`T`$ es invertible (lema de Schur); si
$`\rho'=\rho`$, lo mismo vale para $`T-\lambda1`$ con $`\lambda`$ un valor propio de $`T`$, luego
$`T=\lambda1`$. Por tanto, los operadores de entrelazamiento de $`\rho`$ en $`\rho'`$ forman un espacio de
dimensión $`1`$ si $`\rho\cong\rho'`$ (para un operador de entrelazamiento invertible $`S`$, $`T\mapsto S^{-1}T`$
los identifica con los operadores de entrelazamiento de $`\rho`$ en sí misma, que son $`\mathbb C1`$) y $`0`$ en otro
caso (uno no nulo sería invertible, una equivalencia). Si $`\pi\cong\bigoplus_j\sigma_j`$ con $`\sigma_j`$
irreducibles (se admiten repeticiones), un operador de entrelazamiento de una $`\rho`$ irreducible en
$`\bigoplus_j\sigma_j`$ es una familia $`(T_j)_j`$ de operadores de entrelazamiento de $`\rho`$ en $`\sigma_j`$, y la
composición con la equivalencia no cambia las dimensiones; así, los operadores de entrelazamiento de $`\rho`$
en $`\pi`$ forman un espacio de dimensión $`\#\{j:\sigma_j\cong\rho\}`$, un número que solo
depende de $`\pi`$ y de $`\rho`$. Para la descomposición de la proposición 2.7, cuyas $`\rho_c`$ son dos a
dos no equivalentes, este número es $`m_c`$ si $`\rho\cong\rho_c`$ y $`0`$ si $`\rho`$ no es equivalente a
ninguna $`\rho_c`$. $`\square`$

Tomamos $`G=S_3`$ actuando por permutación de los tres factores de $`(\mathbb C^2)^{\otimes3}`$,
es decir, tres cúbits o tres espines $`\tfrac12`$: $`\pi(\sigma)|x\rangle=|\sigma\cdot x\rangle`$ para $`x\in\{0,1\}^3`$, donde
$`(\sigma\cdot x)_{\sigma(i)}=x_i`$. Sea $`F=\{v:\pi(g)v=v\text{ para todo }g\}`$ el **subespacio simétrico**. La
descomposición es: $`F`$ (dimensión $`4`$) porta la representación trivial con multiplicidad $`4`$; su
complemento ortogonal (dimensión $`4`$) porta la **representación estándar** de dimensión dos (la
representación por permutaciones de $`S_3`$ sobre los vectores de $`\mathbb C^3`$ de suma de coordenadas $`0`$) con
multiplicidad $`2`$; la representación signo no aparece. Así, $`\pi(S_3)''\cong M_1\oplus M_2`$ (dimensión $`5`$) y
$`\pi(S_3)'\cong M_4\oplus M_2`$ (dimensión $`20`$).

*Demostración.* Escribamos $`V^*\pi(S_3)''V=\bigoplus_cM_{k_c}\otimes1_{m_c}`$ (teorema 2.6), de modo que
$`\sum_ck_cm_c=8`$ y $`\sum_cm_c^2=\dim\pi(S_3)'`$.

*Paso 1: $`\dim\pi(S_3)'=20`$.* Como $`\pi(g)^*|g\cdot x\rangle=|x\rangle`$, la matriz
$`\pi(g)X\pi(g)^*`$ tiene la entrada $`X_{xy}`$ en la posición $`(g\cdot x,g\cdot y)`$. Por tanto $`X\in\pi(S_3)'`$, es decir,
$`\pi(g)X\pi(g)^*=X`$ para todo $`g`$, si y solo si $`(x,y)\mapsto X_{xy}`$ es constante en las órbitas de $`S_3`$ actuando
sobre los pares $`(x,y)`$ mediante $`g\cdot(x,y)=(g\cdot x,g\cdot y)`$; las matrices indicatrices de las órbitas forman una base de
$`\pi(S_3)'`$, cuya dimensión es, por tanto, el número de órbitas. Un par $`(x,y)`$ es una terna
$`\big((x_1,y_1),(x_2,y_2),(x_3,y_3)\big)`$ de elementos del conjunto de cuatro elementos $`\{0,1\}^2`$, y $`S_3`$
permuta las tres componentes, de modo que una órbita es un multiconjunto de tamaño $`3`$ formado con cuatro elementos: $`4`$ de la forma
$`\{a,a,a\}`$, $`4\cdot3=12`$ de la forma $`\{a,a,b\}`$ con $`b\neq a`$, y $`4`$ de la forma $`\{a,b,c\}`$ con
elementos distintos, $`20`$ en total.

*Paso 2: $`k=(1,2)`$, $`m=(4,2)`$.* Como $`k_c\ge1`$, $`\sum_cm_c\le\sum_ck_cm_c=8`$. Sea $`M=\max_cm_c`$.
Entonces $`20=\sum_cm_c^2\le M\sum_cm_c\le8M`$ y $`M^2\le20`$, luego $`M\in\{3,4\}`$. Si $`M=4`$, los demás
$`m_c`$ tienen suma de cuadrados $`4`$, luego son $`(2)`$ o $`(1,1,1,1)`$. Si $`M=3`$, los demás $`m_c`$ tienen suma
de cuadrados $`11`$ y suma a lo sumo $`5`$; si todos fueran $`\le2`$, su suma de cuadrados sería a lo sumo
$`2\cdot5=10`$, así que uno de ellos es $`3`$ y los restantes tienen suma de cuadrados $`2`$: son $`(1,1)`$.
Por tanto $`m`$ es, salvo el orden, $`(4,2)`$, $`(4,1,1,1,1)`$ o $`(3,3,1,1)`$. En los dos últimos casos,
$`\sum_cm_c=8=\sum_ck_cm_c`$ obliga a que todos los $`k_c=1`$, y entonces $`\pi(S_3)''\cong\bigoplus_cM_1`$ sería
conmutativa, mientras que las transposiciones $`(12)`$ y $`(23)`$ no conmutan:
$`\pi((12))\pi((23))|001\rangle=|100\rangle`$, pero $`\pi((23))\pi((12))|001\rangle=|010\rangle`$. Así,
$`m=(4,2)`$ (numerando los sumandos de modo que $`m_1=4`$), y $`4k_1+2k_2=8`$ con $`k_1,k_2\ge1`$ da $`k=(1,2)`$. En particular,
$`\dim\pi(S_3)''=\sum_ck_c^2=5`$.

*Paso 3: la representación trivial y la representación signo.* $`\rho_2`$ no fija ningún vector no nulo $`u`$: la recta
$`\mathbb Cu`$ sería un subespacio invariante de $`\mathbb C^2`$, en contra de la irreducibilidad
(proposición 2.7). Escribamos $`V^*v=v_1\oplus v_2`$ con $`v_1\in\mathbb C^1\otimes\mathbb C^4`$ y
$`v_2=\sum_{r=1}^2u_r\otimes f_r\in\mathbb C^2\otimes\mathbb C^2`$ ($`f_1,f_2`$ la base canónica). Entonces
$`v\in F`$ si y solo si $`\rho_1(g)v_1=v_1`$ y $`\rho_2(g)u_r=u_r`$ para todos los $`g`$ y $`r`$, luego todo $`v\in F`$ tiene
$`v_2=0`$. El vector no nulo $`|000\rangle`$ está en $`F`$, así que el escalar $`\rho_1(g)`$ fija un vector no nulo:
$`\rho_1`$ es trivial, $`F=V(\mathbb C^1\otimes\mathbb C^4)`$ tiene dimensión $`4`$, y su complemento
ortogonal es $`V(\mathbb C^2\otimes\mathbb C^2)`$. La representación signo no es equivalente a
$`\rho_1`$ (pues $`\mathrm{sgn}(12)=-1`$) ni a $`\rho_2`$ (las dimensiones difieren), luego, por la unicidad de
las multiplicidades demostrada más arriba, aparece $`0`$ veces.

*Paso 4: $`\rho_2`$ es la representación estándar.* Los vectores $`e_1=|100\rangle`$, $`e_2=|010\rangle`$,
$`e_3=|001\rangle`$ son permutados por los $`\pi(\sigma)`$ igual que $`\sigma`$ permuta las posiciones, luego
$`W=\{\sum_ix_ie_i:\sum_ix_i=0\}`$ es invariante y $`\pi`$ actúa sobre él mediante la representación estándar.
$`W\perp F`$: para $`v\in F`$ y $`w=\sum_ix_ie_i\in W`$, $`\langle v,w\rangle=\langle\pi(g)v,w\rangle=\langle v,\pi(g)^{-1}w\rangle`$
para todo $`g`$, luego $`\langle v,w\rangle=\langle v,\frac16\sum_g\pi(g)w\rangle`$, y
$`\frac16\sum_g\pi(g)w=\frac13\big(\sum_ix_i\big)(e_1+e_2+e_3)=0`$ (cada $`e_i`$ se envía a cada $`e_j`$ mediante
exactamente dos de las seis permutaciones). Por tanto $`V^*W\subseteq\mathbb C^2\otimes\mathbb C^2`$. Las aplicaciones
$`T_r:\mathbb C^2\otimes\mathbb C^2\to\mathbb C^2`$, $`u\otimes f\mapsto\langle f_r,f\rangle u`$, satisfacen
$`T_r(\rho_2(g)\otimes1)=\rho_2(g)T_r`$. Elegimos $`r`$ con $`T_rV^*|_W\neq0`$ (si ambas se anularan, $`V^*W`$ sería
$`0`$). Entonces $`T_rV^*|_W`$ es un operador de entrelazamiento no nulo de $`\pi|_W`$ en $`\rho_2`$; su imagen es un subespacio
invariante no nulo de $`\mathbb C^2`$, luego es $`\mathbb C^2`$, y como $`\dim W=2`$, es biyectivo. Así, $`\rho_2`$ es
equivalente a la representación estándar (que, de paso, resulta ser irreducible). $`\square`$

*Cálculo:* las dimensiones $`20`$ y $`5`$, y el teorema 2.2 en este ejemplo — [tutorial, sección 2.4](../../tutorial/es/02_von_neumann_algebras.md#s-2-4). *En la biblioteca:* `decompose.permutation_representation`, `decompose.commutant_dimension` — [tutorial, sección 4.2](../../tutorial/es/04_symmetries.md#s-4-2).

Las multiplicidades se observan en el espectro. Todo $`K\in\pi(S_3)'`$ autoadjunto es de la forma
$`V^*KV=(1_1\otimes\kappa_1)\oplus(1_2\otimes\kappa_2)`$ con $`\kappa_1\in M_4`$, $`\kappa_2\in M_2`$ autoadjuntas
(teorema 2.6). Sobre $`F`$ sus valores propios son los de $`\kappa_1`$; sobre el complemento de $`F`$ son los
de $`\kappa_2`$, cada uno repetido dos veces, ya que $`1_2\otimes\kappa_2`$ consiste en dos copias de $`\kappa_2`$. Así, *todo* valor propio de
$`K`$ sobre el complemento de $`F`$ tiene multiplicidad par. Para un $`K`$ aleatorio se cumple algo más. El promedio sobre el
grupo, $`E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*`$, aplica $`M_8(\mathbb C)`$ sobre el conmutante y es la
identidad sobre él (ejercicio 3), de modo que, para una matriz de Ginibre $`Z`$ (entradas gaussianas complejas
estándar independientes), $`K=E(Z+Z^*)`$ tiene densidad sobre el espacio vectorial real $`\pi(S_3)'_h`$ de los elementos autoadjuntos
(teorema A.15 y la observación que lo sigue). Por la proposición 4.7(a) (una referencia anticipada; su demostración usa
el teorema 2.6 y el apéndice A, no este párrafo), con probabilidad uno los seis
valores propios de $`\kappa_1`$ y $`\kappa_2`$ son distintos dos a dos: entonces $`K`$ tiene $`4`$ valores propios simples sobre $`F`$ y
$`2`$ doblemente degenerados sobre su complemento. Por último, $`P_F=\frac1{|G|}\sum_g\pi(g)`$ es la proyección ortogonal
sobre $`F`$: $`\pi(h)P_F=P_F`$ para todo $`h`$, luego su imagen está contenida en $`F`$; $`P_Fv=v`$ para $`v\in F`$; y
$`P_F^*=P_F`$, mediante el cambio $`g\to g^{-1}`$.

*Cálculo:* los valores propios doblemente degenerados de un $`K`$ aleatorio — [tutorial, sección 2.4](../../tutorial/es/02_von_neumann_algebras.md#s-2-4).

<a id="s-2-5"></a>

## 2.5 El álgebra en la biblioteca

*En la biblioteca:* cómo `TypeIAlgebra` almacena $`\bigoplus_cM_{k_c}(\mathbb C)`$, sus operadores y su centro — [tutorial, sección 2.5](../../tutorial/es/02_von_neumann_algebras.md#s-2-5).

<a id="s-2-6"></a>

## 2.6 Esperanzas condicionales

Una **esperanza condicional** de $`\mathcal M`$ sobre una subálgebra de von Neumann
$`\mathcal N\subseteq\mathcal M`$ (una $`*`$-álgebra contenida en $`\mathcal M`$, con la misma unidad) es una
aplicación lineal $`E:\mathcal M\to\mathcal N`$ que es positiva ($`A\ge0\Rightarrow E(A)\ge0`$), satisface $`E(B)=B`$ para $`B\in\mathcal N`$ y es un
morfismo de $`\mathcal N`$-bimódulos: $`E(B_1AB_2)=B_1E(A)B_2`$ para $`B_1,B_2\in\mathcal N`$.

<a id="r-2-8"></a>

**Proposición 2.8.** Para $`\mathcal M=\bigoplus_cM_{k_c}`$ sobre $`\bigoplus_c\mathbb C^{k_c}`$, con $`P_c`$ la proyección sobre
el $`c`$-ésimo bloque, y $`\mathcal N=Z(\mathcal M)=\mathrm{span}\{P_1,\dots,P_C\}`$ (ejercicio 1),
```math
E(A)=\sum_c\frac{\mathrm{Tr}\,A_c}{k_c}\,P_c
```
es una esperanza condicional. Conserva toda traza de la forma
$`\sum_cw_c\mathrm{Tr}\,A_c`$.

*Demostración.* La linealidad es evidente; $`A\ge0`$ implica $`\mathrm{Tr}\,A_c=\sum_i\langle e_i,Ae_i\rangle\ge0`$
(suma sobre los vectores $`e_i`$ de la base canónica del $`c`$-ésimo bloque), luego $`E(A)\ge0`$. Para $`B=\sum_cb_cP_c`$, $`E(B)=\sum_cb_ck_cP_c/k_c=B`$. Para $`B_i=\sum_cb^{(i)}_cP_c`$,
$`\mathrm{Tr}(B_1AB_2)_c=b^{(1)}_cb^{(2)}_c\mathrm{Tr}\,A_c`$, lo que da la propiedad
de bimódulo. Por último, $`\mathrm{Tr}(E(A)_c)=\mathrm{Tr}\,A_c`$ para todo $`c`$. $`\square`$

*En la biblioteca:* $`E`$ es `channels.center_expectation` — [tutorial, sección 2.6](../../tutorial/es/02_von_neumann_algebras.md#s-2-6).

El promedio sobre el grupo de la sección 2.4 es una esperanza condicional sobre el
conmutante $`\pi(G)'`$ dentro de $`M_n(\mathbb C)`$, y conserva $`\mathrm{Tr}`$ (ejercicio 3). Análogamente,
en el álgebra de matrices completa $`M_D(\mathbb C)`$, $`D=\sum_ck_c`$, que contiene a $`\mathcal M=\bigoplus_cM_{k_c}`$ como
sus matrices diagonales por bloques, la aplicación $`E_{\mathcal M}(X)=\sum_cP_cXP_c`$ es una esperanza condicional sobre
$`\mathcal M`$ que conserva $`\mathrm{Tr}`$. *Demostración.* $`P_cXP_c`$ se anula fuera del $`c`$-ésimo bloque diagonal,
luego $`E_{\mathcal M}(X)\in\mathcal M`$. Para $`X\ge0`$, $`\langle v,P_cXP_cv\rangle=\langle P_cv,XP_cv\rangle\ge0`$. Un elemento
$`B\in\mathcal M`$ es diagonal por bloques, luego $`BP_c=P_cB`$ y $`B=\sum_cP_cBP_c`$; por tanto $`E_{\mathcal M}(B)=B`$ y
$`E_{\mathcal M}(B_1XB_2)=\sum_cB_1P_cXP_cB_2=B_1E_{\mathcal M}(X)B_2`$ para $`B_1,B_2\in\mathcal M`$. Por último,
$`\sum_c\mathrm{Tr}(P_cXP_c)=\sum_c\mathrm{Tr}(P_cX)=\mathrm{Tr}\,X`$. $`\square`$

*En la biblioteca:* el promedio sobre el grupo es `decompose.twirl` — [tutorial, sección 4.3](../../tutorial/es/04_symmetries.md#s-4-3).

Para un **estado** $`\omega`$ sobre $`\mathcal M`$ (un funcional lineal con $`\omega(A^*A)\ge0`$ para todo $`A`$ y $`\omega(1)=1`$;
los estados se estudian en el capítulo 3), la proposición 2.8 y la linealidad dan
$`\omega(E(A))=\sum_c\omega(P_c)\mathrm{Tr}(A_c)/k_c`$, y $`\omega\circ E=\omega`$ sobre el centro, ya que allí
$`E(z)=z`$. (Las esperanzas condicionales no son únicas en general: sobre $`\mathbb C1`$, cada estado
$`\varphi`$ da una, $`A\mapsto\varphi(A)1`$; en $`M_2(\mathbb C)`$, los estados $`A\mapsto A_{11}`$ y $`A\mapsto A_{22}`$
dan dos distintas.) Todas son casos particulares de la misma idea:
proyectar sobre los operadores que "no ven" cierta estructura.

> **Significado físico** (interpretación física, no un teorema de este curso). La aplicación $`E_{\mathcal M}`$ borra las
> coherencias entre sectores de superselección, que es lo que ocurre cuando se mide la carga (el
> sector) sin registrar el resultado. La esperanza condicional $`E`$ sobre el centro va más allá:
> solo conserva la información clásica, el peso $`\omega(P_c)`$ de cada sector, y reemplaza el estado
> dentro de cada sector por el maximalmente mixto (la traza normalizada $`\mathrm{Tr}(A_c)/k_c`$).

<a id="s-2-7"></a>

## 2.7 Con qué comparar

*Cálculo:* la comparación con otro software (GAP, SageMath, QuSpin, QuTiP) — [tutorial, sección 2.7](../../tutorial/es/02_von_neumann_algebras.md#s-2-7).

<a id="x-2"></a>

## Ejercicios

1. Demuestre que el centro de $`M_k(\mathbb C)`$ es $`\mathbb C1`$ (use las unidades matriciales
   $`|i\rangle\langle j|`$) y que $`Z(\bigoplus_cM_{k_c})`$ tiene dimensión $`C`$.
2. Sea $`S`$ el desplazamiento cíclico en $`\mathbb C^3`$, $`Se_i=e_{i+1\bmod3}`$. Determine
   $`\{S\}'`$ y $`\{S\}''`$.
3. Para una representación unitaria $`\pi`$ de un grupo finito, demuestre que
   $`E(X)=\frac1{|G|}\sum_g\pi(g)X\pi(g)^*`$ es una esperanza condicional de
   $`M_n(\mathbb C)`$ sobre $`\pi(G)'`$ que conserva $`\mathrm{Tr}`$.

<details>
<summary>Soluciones</summary>

1. Si $`X`$ conmuta con todas las $`|i\rangle\langle j|`$, entonces
   $`X|i\rangle\langle j|=|i\rangle\langle j|X`$; aplicando ambos miembros a $`e_j`$ y comparando
   se obtiene $`Xe_i=X_{jj}e_i`$ para todos los $`i,j`$, luego $`X=X_{11}1`$. En una suma directa,
   un elemento $`\bigoplus X_c`$ es central si y solo si lo es cada $`X_c`$, luego el centro es
   $`\{\bigoplus_cz_c1_{k_c}\}`$, de dimensión $`C`$.

2. $`S`$ es unitaria con valores propios $`1,\omega,\omega^2`$, $`\omega=e^{2\pi i/3}`$, todos
   simples: los vectores $`u_j=\sum_i\omega^{-ij}e_i`$ ($`j=0,1,2`$) satisfacen $`Su_j=\omega^ju_j`$ y son ortogonales.
   Un operador que conmuta con $`S`$ conserva sus subespacios propios, de modo que
   en la base propia es diagonal, y toda matriz diagonal conmuta con $`S`$: $`\{S\}'\cong\mathbb C^3`$. Contiene
   $`1,S,S^2`$, que en la base propia son $`\mathrm{diag}(\lambda^0)`$, $`\mathrm{diag}(\lambda)`$,
   $`\mathrm{diag}(\lambda^2)`$ para $`\lambda=(1,\omega,\omega^2)`$, linealmente independientes (Vandermonde: los valores
   propios son distintos); así, $`\{S\}'=\mathrm{span}\{1,S,S^2\}`$, las matrices circulantes. Entonces $`\{S\}''`$ también es
   $`\{S\}'`$, porque en la base propia de $`S`$ una matriz que conmuta con todas las matrices
   diagonales es diagonal (compárense las entradas de $`XD`$ y $`DX`$ para una $`D`$ con entradas
   diagonales distintas).

3. $`E`$ es lineal y positiva (una suma de aplicaciones $`X\mapsto UXU^*`$). Para $`h\in G`$,
   $`\pi(h)E(X)\pi(h)^*=\frac1{|G|}\sum_g\pi(hg)X\pi(hg)^*=E(X)`$, luego $`E(X)\in\pi(G)'`$. Si
   $`X\in\pi(G)'`$, cada término es igual a $`X`$, luego $`E(X)=X`$. Para $`B_1,B_2\in\pi(G)'`$,
   $`\pi(g)B_1XB_2\pi(g)^*=B_1\pi(g)X\pi(g)^*B_2`$, lo que da la propiedad de bimódulo.
   Por último, $`\mathrm{Tr}(\pi(g)X\pi(g)^*)=\mathrm{Tr}\,X`$ para todo $`g`$, luego
   $`\mathrm{Tr}\,E(X)=\mathrm{Tr}\,X`$.

*Cálculo:* el ejercicio 2, numéricamente — [tutorial, Ejercicios](../../tutorial/es/02_von_neumann_algebras.md#x-2).

</details>

**Siguiente:** el [capítulo 3](03_traces_and_states.md) clasifica las trazas sobre
$`\bigoplus_cM_{k_c}`$, introduce los estados y las densidades, y demuestra las desigualdades
de entropía y de trazas que se usan más adelante.
