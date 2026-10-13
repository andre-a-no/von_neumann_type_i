[English](../05_random_matrices.md) · **Español** · [Français](../fr/05_random_matrices.md) · [Deutsch](../de/05_random_matrices.md) · [中文](../zh/05_random_matrices.md) · [日本語](../ja/05_random_matrices.md) · [Русский](../ru/05_random_matrices.md)

# 5. Matrices aleatorias en un álgebra

La teoría de matrices aleatorias se pregunta cómo es un operador "típico". En un álgebra
$\bigoplus_cM_{k_c}$ la respuesta debe respetar los sectores: las unitarias aleatorias naturales
son las del álgebra, no las de la gran álgebra de matrices $M_D$, y los hamiltonianos aleatorios
naturales tienen bloques independientes. Este capítulo construye la medida de Haar en el grupo
unitario del álgebra, calcula exactamente sus momentos de primer y segundo orden y muestra después
qué falla en la estadística espectral cuando se ignoran los sectores.

**Aprenderá** por qué la descomposición QR de una matriz gaussiana da una unitaria de Haar (y por
qué hay que fijar las fases), que promediar sobre el grupo unitario del álgebra es la esperanza
condicional sobre el centro, la fórmula de Weingarten de orden dos, la conjetura de Wigner y el
cociente de espaciamientos, y por qué la estadística de niveles calculada sobre un espectro con
sectores mezclados se parece a la de niveles independientes.

En todo el capítulo, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ con $D=\sum_ck_c$, $P_c$ son las
proyecciones centrales minimales y $E$ es la esperanza condicional sobre el centro de la
proposición 2.8. $U(k)$ es el grupo de las matrices unitarias $k\times k$ y
$U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}$ es el grupo unitario de $\mathcal M$. Las variables
aleatorias viven en un espacio de probabilidad con esperanza $\mathbb E$; la esperanza de una matriz
aleatoria se toma entrada a entrada. Una variable *gaussiana compleja estándar* tiene partes real e
imaginaria independientes $\sim N(0,\tfrac12)$, de modo que $\mathbb E|z|^2=1$ (esto es lo que
genera `torch.randn` para un dtype complejo); una *matriz de Ginibre* tiene entradas independientes
gaussianas complejas estándar.

## 5.1 Medida de Haar

Sea $G$ un grupo compacto de matrices (un subgrupo cerrado de algún $U(n)$). Una medida de
probabilidad de Borel $\mu$ sobre $G$ es **invariante por la izquierda** si $\mu(gB)=\mu(B)$ para
todo $g\in G$ y todo boreliano $B$, o, equivalentemente, si $\int f(gh)\,d\mu(h)=\int f\,d\mu$ para
toda función de Borel acotada $f$; la **invariancia por la derecha** se define con $Bg$ y $f(hg)$.

**Lema 5.1 (unicidad).** Si $\mu$ es una medida de probabilidad de Borel invariante por la
izquierda y $\nu$ una invariante por la derecha sobre $G$, entonces $\mu=\nu$. Por tanto, $G$
admite a lo sumo una medida de probabilidad de Borel invariante por la izquierda, y si existe es
también invariante por la derecha e invariante bajo $g\mapsto g^{-1}$.

*Demostración.* Sea $f:G\to\mathbb R$ continua (y, por tanto, acotada). Por el teorema de Fubini,
$\int\!\big(\int f(gh)\,d\mu(h)\big)d\nu(g)=\int\!\big(\int f(gh)\,d\nu(g)\big)d\mu(h)$. La
integral interior de la izquierda es $\int f\,d\mu$ por la invariancia por la izquierda, y la de la
derecha es $\int f\,d\nu$ por la invariancia por la derecha; así, $\int f\,d\mu=\int f\,d\nu$. Dos
medidas de probabilidad de Borel sobre un espacio métrico con las mismas integrales de funciones
continuas son iguales: las funciones continuas $\max(0,1-n\operatorname{dist}(\cdot,K))$ decrecen
hacia la indicatriz de un cerrado $K$, de modo que las medidas coinciden sobre los cerrados
(convergencia dominada), que generan los borelianos y son estables por intersecciones (teorema
$\pi$–$\lambda$ de Dynkin). Si $\mu_1,\mu_2$ son invariantes por la izquierda, entonces
$\nu(B)=\mu_2(B^{-1})$ es invariante por la derecha, ya que $(Bg)^{-1}=g^{-1}B^{-1}$; luego
$\mu_1=\nu=\mu_2$, y $\mu_1$ es invariante por la derecha y por inversión. $\square$

**Definición.** La **medida de Haar** de $G$ es su medida de probabilidad de Borel invariante por
la izquierda. Un elemento aleatorio con esta ley está **distribuido según la medida de Haar**.

Todo grupo compacto tiene una medida de Haar (véase G. B. Folland, *A Course in Abstract Harmonic
Analysis*, 2.ª ed., CRC Press, 2016, sección 2.2; se enuncia aquí sin demostración). Para $U(k)$ no
necesitamos este teorema general: el teorema 5.3 construye la medida.

**Lema 5.2 (descomposición QR).** Toda $Z\in M_k(\mathbb C)$ invertible puede escribirse $Z=QR$
con $Q$ unitaria y $R$ triangular superior con diagonal positiva. $Q$ y $R$ son únicas y dependen
de forma continua de $Z$.

*Demostración.* Apliquemos Gram–Schmidt a las columnas $z_1,\dots,z_k$ de $Z$, que son linealmente
independientes: $q_j=w_j/\|w_j\|$ con $w_j=z_j-\sum_{i<j}\langle q_i,z_j\rangle q_i$; aquí
$w_j\neq0$ porque
$z_j\notin\operatorname{span}\{z_1,\dots,z_{j-1}\}=\operatorname{span}\{q_1,\dots,q_{j-1}\}$. Entonces
$z_j=\sum_{i<j}\langle q_i,z_j\rangle q_i+\|w_j\|q_j$, es decir, $Z=QR$ con $R_{jj}=\|w_j\|>0$, y
cada paso es continuo sobre las matrices invertibles. Si $Q_1R_1=Q_2R_2$, entonces
$T=Q_2^*Q_1=R_2R_1^{-1}$ es unitaria y triangular superior con diagonal positiva. Su inversa $T^*$
es triangular inferior, y también triangular superior (la inversa de una matriz triangular superior
es triangular superior), luego $T$ es diagonal; una unitaria diagonal con diagonal positiva es
$1$. $\square$

**Teorema 5.3 (unitarias de Haar a partir de matrices gaussianas).** Sea $Z$ una matriz de
Ginibre $k\times k$. Entonces $Z$ es invertible con probabilidad uno, y el factor $Q$ del lema 5.2
está distribuido según la medida de Haar en $U(k)$. Además, si $Z=Q'R'$ es cualquier
descomposición QR con $Q'$ unitaria y $R'$ triangular superior (por ejemplo, la que devuelve
LAPACK, cuya diagonal es real de cualquier signo), entonces $Q=Q'\Lambda$ con
$\Lambda=\operatorname{diag}(R'_{ii}/|R'_{ii}|)$ está distribuida según la medida de Haar.

*Demostración.* *Invertibilidad.* $\det Z$ es un polinomio no nulo en las $2k^2$ coordenadas
reales de $Z$, y el conjunto de ceros de un polinomio real no nulo tiene medida de Lebesgue nula
(inducción sobre el número de variables: escribamos $p=\sum_jp_j(x')x_n^j$ con algún $p_j\neq0$;
para $x'$ fuera del conjunto nulo de ceros de $p_j$ hay un número finito de raíces $x_n$; se
concluye por Fubini). $Z$ tiene densidad, luego $P(\det Z=0)=0$.

*Invariancia de la ley gaussiana.* La densidad de $Z$ respecto de la medida de Lebesgue sobre
$M_k(\mathbb C)\cong\mathbb R^{2k^2}$ es $\pi^{-k^2}e^{-\|Z\|_2^2}$. Para $V\in U(k)$, la aplicación
$Z\mapsto VZ$ es lineal real y conserva $\|Z\|_2$ (sección 1.4), así que es una aplicación
ortogonal de $\mathbb R^{2k^2}$ con jacobiano $1$, y conserva la densidad: $VZ$ tiene la misma ley
que $Z$.

*Equivariancia.* Si $Z=QR$, entonces $VZ=(VQ)R$ es la descomposición del lema 5.2 de $VZ$, luego
$Q(VZ)=VQ(Z)$. Por tanto, para todo boreliano $B\subseteq U(k)$,
$P(VQ(Z)\in B)=P(Q(VZ)\in B)=P(Q(Z)\in B)$: la ley de $Q(Z)$ es invariante por la izquierda y,
por el lema 5.1, es la medida de Haar (esto prueba también que la medida de Haar de $U(k)$
existe). Por último, $Q'R'=(Q'\Lambda)(\Lambda^*R')$, y $\Lambda^*R'$ es triangular superior con
diagonal $|R'_{ii}|>0$, luego $Q'\Lambda=Q(Z)$ por unicidad. $\square$

La misma demostración con matrices gaussianas reales da la medida de Haar del grupo ortogonal
$O(k)$. Sin la corrección $\Lambda$ el resultado *no* está distribuido según la medida de Haar: el
convenio de LAPACK para los signos de los $R'_{ii}$ introduce un sesgo, que el código siguiente
pone de manifiesto (F. Mezzadri, *How to generate random matrices from the classical compact
groups*, Notices AMS 54 (2007), 592–604).

**Corolario 5.4 (medida de Haar del álgebra).** $U(\mathcal M)=\{\bigoplus_cU_c:U_c\in U(k_c)\}\cong\prod_cU(k_c)$,
y $U=\bigoplus_cU_c$ está distribuida según la medida de Haar en $U(\mathcal M)$ si y solo si los
bloques $U_c$ son independientes y cada $U_c$ está distribuido según la medida de Haar en
$U(k_c)$.

*Demostración.* $U^*U=\bigoplus_cU_c^*U_c$, luego $U$ es unitaria si y solo si lo es cada bloque.
Sean los $U_c$ independientes y distribuidos según la medida de Haar. Para
$V=\bigoplus_cV_c\in U(\mathcal M)$, los bloques $V_cU_c$ de $VU$ son independientes (funciones de
variables independientes) y $V_cU_c$ tiene la ley de $U_c$, así que $VU$ tiene la ley de $U$. Por
tanto, la ley de $U$ es invariante por la izquierda y, por el lema 5.1, es la medida de Haar del
grupo compacto $U(\mathcal M)$. Recíprocamente, por unicidad, toda $U$ distribuida según la medida
de Haar tiene esta ley. $\square$

El grupo $U(\mathcal M)$ es mucho más pequeño que $U(D)$: para $C\ge2$, una unitaria de Haar de
$U(D)$ pertenece a $\mathcal M$ con probabilidad cero (por el teorema 5.3 y la construcción de
Gram–Schmidt del lema 5.2, su primera columna es $z/\|z\|$ para un vector gaussiano $z$, así que
todas sus entradas son no nulas casi seguramente). La medida de Haar de $SU(k)$ se obtiene con el
mismo espíritu: si $U$ es de Haar en $U(k)$ y $\theta=\arg\det U$, entonces
$V=e^{-i\theta/k}U\in SU(k)$, y para $W\in SU(k)$, $\det(WU)=\det U$ da $V(WU)=WV(U)$; como $WU$
tiene la ley de $U$, la ley de $V$ es invariante por la izquierda en $SU(k)$ y, por tanto, es la
medida de Haar de $SU(k)$ por el lema 5.1.

**Corolario 5.5 (primeros momentos de las entradas).** Si $U$ está distribuida según la medida de
Haar en $U(k)$, entonces $\mathbb EU_{ij}=0$ y $\mathbb E|U_{ij}|^2=1/k$.

*Demostración.* $-1\in U(k)$, luego $-U$ tiene la ley de $U$ y $\mathbb EU=-\mathbb EU$. Para una
matriz de permutación $P$, $PU$ tiene la ley de $U$, así que $\mathbb E|U_{ij}|^2$ no depende de
$i$; como $\sum_i|U_{ij}|^2=1$, vale $1/k$. $\square$

La biblioteca genera unitarias de Haar de $U(\mathcal M)$ con `random_unitary_operator`, mediante
una única descomposición QR por lotes por cada tamaño de bloque, para todos los sectores y todas
las muestras a la vez (en una GPU, mediante Gram–Schmidt, que da la misma $Q$ por la unicidad del
lema 5.2). La medida `'haar'` es el teorema 5.3 en cada bloque; para un álgebra compleja,
`'haar_su'` multiplica por $e^{-i\theta/k}$ como antes; `'coe'` y `'cse'` son los conjuntos
circulares de Dyson (el primero es el ejercicio 3). Comprobamos los corolarios 5.4 y 5.5, y que el
muestreador de la biblioteca es exactamente la QR corregida del teorema 5.3:

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinChain, channels

torch.manual_seed(0)
dt = torch.complex128

def close(X, Y, tol=1e-10):
    X = X.matrix if hasattr(X, 'matrix') else X
    Y = Y.matrix if hasattr(Y, 'matrix') else Y
    return (X - Y).abs().max().item() < tol

def mc(samples, exact, name):
    """Monte Carlo mean of real samples with its standard error; must agree within 5 standard errors."""
    mean, err = samples.mean().item(), samples.std().item() / math.sqrt(samples.numel())
    print(f"{name}: {mean:.4f} +- {err:.4f}   (theory {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12
    return mean

alg = TypeIAlgebra([1, 2, 3], [1, 2, 3], precision='double')        # M_1 (+) M_2 (+) M_3, D = 6
B = 100_000
U = alg.random_unitary_operator(B)                                   # Haar on U(1) x U(2) x U(3)
assert close(U.adjoint() @ U, alg.identity())                        # unitary elements of the algebra
u = [U.matrix[:, c, :k, :k] for c, k in enumerate(alg.k_factors)]   # the blocks U_c
for c, k in enumerate(alg.k_factors):                                # Corollary 5.5
    mc(u[c][:, 0, 0].real, 0.0, f"E Re (U_{c})_11")
    mc(u[c][:, 0, 0].abs() ** 2, 1 / k, f"E |(U_{c})_11|^2")
mc(u[1][:, 0, 0].abs() ** 2 * u[2][:, 0, 0].abs() ** 2, 1 / 6, "E |(U_1)_11|^2 |(U_2)_11|^2")   # independent blocks

SU = alg.random_unitary_operator(4, measure='haar_su').matrix
assert torch.allclose(torch.linalg.det(SU[:, 2, :3, :3]), torch.ones(4, dtype=dt))

torch.manual_seed(1)
Z = torch.randn(B, 2, 2, dtype=dt)                                   # Ginibre matrices
torch.manual_seed(1)
V = alg.random_unitary(2, batch_size=B)                              # the library, same random numbers
Q, R = torch.linalg.qr(Z)                                            # LAPACK: diag(R) real, either sign
d = R.diagonal(dim1=-2, dim2=-1)
assert close(V, Q * (d / d.abs()).unsqueeze(-2))                     # Theorem 5.3: V = Q Lambda
raw = Q.diagonal(dim1=-2, dim2=-1).sum(-1).real                      # Re Tr Q without the correction
print(f"E Re Tr Q without the phase correction: {raw.mean():.4f} +- {raw.std() / math.sqrt(B):.4f}")
assert abs(raw.mean()) > 20 * raw.std() / math.sqrt(B)               # far from 0: not Haar
```

Sin la corrección, $\mathbb E\operatorname{Re}\operatorname{Tr}Q\approx-0.84$ para $k=2$, mientras
que una unitaria de Haar tiene $\mathbb EU=0$.

## 5.2 Momentos y el cálculo de Weingarten

**Lema 5.6.** Si $T\in M_k(\mathbb C)$ conmuta con todas las unitarias diagonales y todas las
matrices de permutación, entonces $T\in\mathbb C1$. En particular, el conmutante de $U(k)$ es
$\mathbb C1$.

*Demostración.* Para $\Delta=\operatorname{diag}(e^{i\theta_1},\dots,e^{i\theta_k})$ con los
$e^{i\theta_j}$ distintos, $\Delta T=T\Delta$ da $(e^{i\theta_i}-e^{i\theta_j})T_{ij}=0$, luego
$T$ es diagonal. Para la matriz de permutación $P_\sigma e_j=e_{\sigma(j)}$, $P_\sigma TP_\sigma^{-1}=T$
da $T_{\sigma(j)\sigma(j)}=T_{jj}$, así que todas las entradas diagonales son iguales. $\square$

**Proposición 5.7 (el twirl de Haar es la esperanza condicional sobre el centro).** Sea $U$
distribuida según la medida de Haar en $U(\mathcal M)$. Para todo $X\in\mathcal M$,
$$
\mathbb E\,[UXU^*]=E(X)=\bigoplus_c\frac{\operatorname{Tr}X_c}{k_c}\,1_{k_c}.
$$

*Demostración.* $(UXU^*)_c=U_cX_cU_c^*$, y $U_c$ es de Haar en $U(k_c)$ (corolario 5.4). Sea
$T_c=\mathbb E[U_cX_cU_c^*]$ (las entradas están acotadas). Para $V\in U(k_c)$, $VU_c$ tiene la ley
de $U_c$, luego $VT_cV^*=\mathbb E[(VU_c)X_c(VU_c)^*]=T_c$. Por el lema 5.6, $T_c=\lambda1$, y
$k_c\lambda=\operatorname{Tr}T_c=\mathbb E\operatorname{Tr}(U_cX_cU_c^*)=\operatorname{Tr}X_c$. $\square$

La normalización es la de la proposición 2.8 y la de `channels.center_expectation`. Si
$X\in\mathcal M$ conmuta con $U(\mathcal M)$, entonces $UXU^*=X$, luego $X=E(X)$ es central; el
promedio de Haar es el análogo continuo del promedio sobre el grupo del capítulo 2, ejercicio 3. Si,
en cambio, $U$ es de Haar en $U(D)$, la misma demostración (con $C=1$) da
$\mathbb E[UXU^*]=\frac{\operatorname{Tr}X}{D}1$: se pierden los pesos de los sectores, y $UXU^*$ ni
siquiera está en $\mathcal M$.

> **Significado físico.** Una unitaria aleatoria de Haar del álgebra es una evolución aleatoria que
> respeta la regla de superselección. Promediar sobre ella desordena todo dentro de cada sector y
> conserva exactamente la parte clásica: el peso de cada carga. Este es el "twirl" de la
> información cuántica, y es la esperanza condicional sobre el centro del capítulo 2.

Para los momentos de segundo orden necesitamos el conmutante de $\{V\otimes V\}$ sobre
$\mathbb C^k\otimes\mathbb C^k$. Sea $F$ el **intercambio** (*swap*), $F(x\otimes y)=y\otimes x$;
es autoadjunto y unitario, y
$\operatorname{Tr}1=k^2$, $\operatorname{Tr}F=\sum_{i,j}\langle e_i\otimes e_j,e_j\otimes e_i\rangle=k$, $F^2=1$.

**Lema 5.8.** Sea $k\ge2$. Una matriz $X$ sobre $\mathbb C^k\otimes\mathbb C^k$ conmuta con
$V\otimes V$ para toda $V\in U(k)$ si y solo si $X=a1+bF$ con $a,b\in\mathbb C$; además, $1$ y $F$
son linealmente independientes.

*Demostración.* $F(V\otimes V)(x\otimes y)=Vy\otimes Vx=(V\otimes V)F(x\otimes y)$, luego $1$ y $F$
conmutan con todos los $V\otimes V$. Recíprocamente, supongamos que $X$ conmuta con ellos, y
escribamos $X_{ij,pq}=\langle e_i\otimes e_j,X(e_p\otimes e_q)\rangle$.
*Unitarias diagonales:* $\Delta\otimes\Delta$ multiplica $e_i\otimes e_j$ por
$e^{i(\theta_i+\theta_j)}$, luego $(e^{i(\theta_i+\theta_j)}-e^{i(\theta_p+\theta_q)})X_{ij,pq}=0$.
Si los multiconjuntos $\{i,j\}$ y $\{p,q\}$ son distintos, algún índice $m$ aparece
$n\in\{\pm1,\pm2\}$ veces más en el primero; con $\theta_m=\pi/2$ y todos los demás $\theta=0$, el
corchete vale $e^{i(\theta_p+\theta_q)}(e^{in\pi/2}-1)\neq0$, luego $X_{ij,pq}=0$. Por tanto
$X(e_i\otimes e_i)=a_ie_i\otimes e_i$ y, para $i\neq j$,
$X(e_i\otimes e_j)=b_{ij}\,e_i\otimes e_j+c_{ij}\,e_j\otimes e_i$.
*Permutaciones:* aplicando $X(P_\sigma\otimes P_\sigma)=(P_\sigma\otimes P_\sigma)X$ a $e_i\otimes e_j$
se obtiene $a_{\sigma(i)}=a_i$, $b_{\sigma(i)\sigma(j)}=b_{ij}$, $c_{\sigma(i)\sigma(j)}=c_{ij}$. Las
permutaciones actúan transitivamente sobre los pares ordenados de índices distintos, así que
$a_i=a'$, $b_{ij}=b'$, $c_{ij}=c'$; por tanto, $X=b'1+c'F$ sobre el subespacio generado por los
$e_i\otimes e_j$ con $i\ne j$.
*Una rotación:* sea $R$ la unitaria con $Re_1=(e_1+e_2)/\sqrt2$, $Re_2=(-e_1+e_2)/\sqrt2$, $Re_m=e_m$ para $m\ge3$.
Entonces $(R\otimes R)(e_1\otimes e_1)=\frac12(e_1\otimes e_1+e_1\otimes e_2+e_2\otimes e_1+e_2\otimes e_2)$, y
$X(R\otimes R)(e_1\otimes e_1)=\frac12\big(a'e_1\otimes e_1+(b'+c')(e_1\otimes e_2+e_2\otimes e_1)+a'e_2\otimes e_2\big)$
debe ser igual a $(R\otimes R)X(e_1\otimes e_1)=a'(R\otimes R)(e_1\otimes e_1)$. Comparando los
coeficientes de $e_1\otimes e_2$ se obtiene $a'=b'+c'$, que es el valor de $b'1+c'F$ sobre
$e_i\otimes e_i$. Así, $X=b'1+c'F$. Por último, $1$ y $F$ llevan $e_1\otimes e_2$ a los vectores
independientes $e_1\otimes e_2$ y $e_2\otimes e_1$. $\square$

**Corolario 5.9 (dualidad de Schur–Weyl de orden dos).** Sean $k\ge2$ y $P_\pm=\frac12(1\pm F)$.
Son las proyecciones sobre los subespacios simétrico y antisimétrico
$\mathrm{Sym}=\{x:Fx=x\}$ y $\mathrm{Alt}=\{x:Fx=-x\}$, de dimensiones $k(k+1)/2$ y $k(k-1)/2$.
La representación $V\mapsto V\otimes V$ de $U(k)$ es irreducible sobre $\mathrm{Sym}$ y sobre
$\mathrm{Alt}$, y las dos representaciones no son equivalentes.

*Demostración.* $F=F^*$ y $F^2=1$ hacen de $P_\pm$ proyecciones con $P_++P_-=1$ e imágenes
$\mathrm{Sym}$, $\mathrm{Alt}$; sus dimensiones son $\operatorname{Tr}P_\pm=(k^2\pm k)/2$. El
subespacio $\mathcal A$ generado por $\{V\otimes V\}$ es una $*$-álgebra ($1\otimes1$, productos
$VW\otimes VW$, adjuntos $V^*\otimes V^*$) con conmutante
$\mathcal A'=\operatorname{span}\{1,F\}=\operatorname{span}\{P_+,P_-\}$ por el lema 5.8. Sea
$K\subseteq\mathrm{Sym}$ invariante bajo todos los $V\otimes V$. Por el lema 2.1, su proyección está
en $\mathcal A'$, luego es $xP_++yP_-$ con $x,y\in\{0,1\}$ (es una proyección); se anula sobre
$\mathrm{Alt}$, así que $y=0$, y $K=0$ o $K=\mathrm{Sym}$. El mismo argumento vale para
$\mathrm{Alt}$. Si $T:\mathrm{Sym}\to\mathrm{Alt}$ entrelaza las dos representaciones, entonces
$TP_+$ (visto como operador sobre $\mathbb C^k\otimes\mathbb C^k$) conmuta con $\mathcal A$, luego
$TP_+=xP_++yP_-$; se anula sobre $\mathrm{Alt}$, así que $y=0$, y su imagen está contenida en
$\mathrm{Alt}\cap\mathrm{Sym}=0$, luego $T=0$. $\square$

**Teorema 5.10 (fórmula de Weingarten de orden dos).** Sean $k\ge2$, $U$ distribuida según la
medida de Haar en $U(k)$ y $X$ una matriz sobre $\mathbb C^k\otimes\mathbb C^k$. Entonces
$$
\mathbb E\,[(U\otimes U)X(U\otimes U)^*]=\alpha\,1+\beta F,\qquad
\alpha=\frac{\operatorname{Tr}X-\operatorname{Tr}(FX)/k}{k^2-1},\quad
\beta=\frac{\operatorname{Tr}(FX)-\operatorname{Tr}X/k}{k^2-1}.
$$
Equivalentemente, con $\mathrm{Wg}(e)=\frac1{k^2-1}$, $\mathrm{Wg}(\tau)=-\frac1{k(k^2-1)}$ sobre
el grupo $S_2=\{e,\tau\}$,
$$
\mathbb E\big[U_{i_1j_1}U_{i_2j_2}\overline{U_{i'_1j'_1}U_{i'_2j'_2}}\big]
=\sum_{\sigma,\pi\in S_2}\delta_{i_1i'_{\sigma(1)}}\delta_{i_2i'_{\sigma(2)}}\,
\delta_{j_1j'_{\pi(1)}}\delta_{j_2j'_{\pi(2)}}\,\mathrm{Wg}(\sigma\pi^{-1}).
$$

*Demostración.* Sea $T=\mathbb E[(U\otimes U)X(U\otimes U)^*]$. Como en la proposición 5.7, $T$
conmuta con todos los $V\otimes V$, luego $T=\alpha1+\beta F$ por el lema 5.8. Tomando trazas,
$\operatorname{Tr}T=\operatorname{Tr}X$ y, como $F$ conmuta con $U\otimes U$,
$\operatorname{Tr}(FT)=\mathbb E\operatorname{Tr}((U\otimes U)^*F(U\otimes U)X)=\operatorname{Tr}(FX)$.
Esto da $k^2\alpha+k\beta=\operatorname{Tr}X$ y $k\alpha+k^2\beta=\operatorname{Tr}(FX)$, un sistema
lineal de determinante $k^2(k^2-1)\neq0$, cuya solución es la indicada. Para las entradas, tomemos
$X=|e_{j_1}\otimes e_{j_2}\rangle\langle e_{j'_1}\otimes e_{j'_2}|$: entonces
$\langle e_{i_1}\otimes e_{i_2},T\,e_{i'_1}\otimes e_{i'_2}\rangle$ es el primer miembro,
$\operatorname{Tr}X=\delta_{j_1j'_1}\delta_{j_2j'_2}$ (el término $\pi=e$),
$\operatorname{Tr}(FX)=\delta_{j_1j'_2}\delta_{j_2j'_1}$ (el término $\pi=\tau$), y $1$, $F$ aportan
$\delta_{i_1i'_1}\delta_{i_2i'_2}$ (el término $\sigma=e$) y $\delta_{i_1i'_2}\delta_{i_2i'_1}$
($\sigma=\tau$). Agrupando términos,
$\alpha=\sum_\pi\mathrm{Wg}(e\pi^{-1})[\pi]$ y $\beta=\sum_\pi\mathrm{Wg}(\tau\pi^{-1})[\pi]$,
donde $[\pi]$ denota las deltas en $j$ de $\pi$. $\square$

La matriz $\big(\mathrm{Wg}(\sigma^{-1}\pi)\big)_{\sigma,\pi}$ es la inversa de la matriz de Gram
$\begin{pmatrix}k^2&k\\k&k^2\end{pmatrix}$ de $1$ y $F$ respecto de
$\langle A,B\rangle=\operatorname{Tr}(A^*B)$; esta es la forma que se generaliza.

**Corolario 5.11.** Para $k\ge2$:
$$
\mathbb E|U_{11}|^4=\frac2{k(k+1)},\qquad
\mathbb E|U_{11}|^2|U_{12}|^2=\frac1{k(k+1)},\qquad
\mathbb E\big[U_{11}U_{22}\overline{U_{12}U_{21}}\big]=-\frac1{(k-1)k(k+1)}.
$$

*Demostración.* En el teorema 5.10, (i) todos los índices iguales a $1$: todo producto de deltas
vale $1$, así que el momento es $2(\mathrm{Wg}(e)+\mathrm{Wg}(\tau))=2\frac{k-1}{k(k^2-1)}$.
(ii) $i=(1,1)$, $j=(1,2)$, $i'=(1,1)$, $j'=(1,2)$: solo sobrevive $\pi=e$, lo que da
$\mathrm{Wg}(e)+\mathrm{Wg}(\tau)$. (iii) $i=(1,2)$, $j=(1,2)$, $i'=(1,2)$, $j'=(2,1)$: solo
sobreviven $\sigma=e$, $\pi=\tau$, lo que da $\mathrm{Wg}(\tau)$. $\square$

El último momento es negativo: las entradas de una unitaria de Haar no son independientes, y las
correlaciones son de orden relativo $1/k$. Para el álgebra, las entradas de bloques distintos son
independientes (corolario 5.4), de modo que los momentos mixtos entre sectores se factorizan en
momentos de bloques individuales.

Para momentos de cualquier orden se aplica el mismo razonamiento una vez conocido el conmutante de
$\{V^{\otimes n}\}$; está generado por las permutaciones de los $n$ factores tensoriales (dualidad
de Schur–Weyl; véase R. Goodman, N. R. Wallach, *Symmetry, Representations, and Invariants*,
Springer, 2009, capítulo 9). Enunciamos el resultado sin demostración.

**Teorema 5.12 (fórmula de Weingarten; enunciado sin demostración).** Sean $n\ge1$ y $U$ de Haar
en $U(k)$. Entonces
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\pi\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\pi(m)}}\;\mathrm{Wg}(\sigma\pi^{-1},k),
$$
donde $\mathrm{Wg}(\cdot,k)$ es una función sobre $S_n$ que solo depende del tipo de ciclos. Para
$n\le k$, la matriz $\big(\mathrm{Wg}(\sigma^{-1}\pi,k)\big)_{\sigma,\pi}$ es la inversa de la
matriz de Gram $\big(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)}\big)_{\sigma,\pi}$; para $n>k$ esta
matriz es singular y se toma su pseudoinversa. Los momentos con distinto número de factores $U$ y
$\overline U$ se anulan. (D. Weingarten, J. Math. Phys. 19 (1978), 999–1001; B. Collins, Int. Math.
Res. Not. 2003, 953–982; B. Collins, P. Śniady, Comm. Math. Phys. 264 (2006), 773–795.)

Comprobamos la proposición 5.7 frente a `channels.center_expectation`, el conmutante del lema 5.8
mediante álgebra lineal exacta, y el teorema 5.10 y el corolario 5.11 mediante Monte Carlo con el
bloque $U(3)$ de las muestras anteriores:

```python
X = alg.from_blocks([torch.randn(k, k, dtype=dt) for k in alg.k_factors])      # a fixed element
T = (U @ X @ U.adjoint()).matrix                                    # B samples of U X U*
mean, err = T.mean(0), T.std(0) / math.sqrt(B)
EX = channels.center_expectation(alg)(X).matrix[0]
random = err > 1e-12                                                 # entries that are not constant
z = ((mean - EX).abs()[random] / err[random]).max().item()
print(f"twirl over U(M) vs E(X): largest deviation {z:.2f} standard errors")
assert ((mean - EX).abs() <= 5 * err + 1e-12).all()                 # Proposition 5.7

def commutant_dim(ops, tol=1e-9):
    """Dimension of {X : XA = AX for all A in ops}, as the null space of X -> (AX - XA)_A."""
    d = ops[0].shape[0]
    eye = torch.eye(d, dtype=dt)
    ops = [A.contiguous() for A in ops]                              # torch.kron needs contiguous input
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return int((torch.linalg.svdvals(L) < tol).sum())

k = 3
F = torch.zeros(k * k, k * k, dtype=dt)
for i in range(k):
    for j in range(k):
        F[j * k + i, i * k + j] = 1                                  # F (e_i (x) e_j) = e_j (x) e_i
VV = [torch.kron(W, W) for W in alg.random_unitary(k, batch_size=3).contiguous()]
assert all(close(F @ W, W @ F) for W in VV)
assert commutant_dim(VV) == 2                                        # Lemma 5.8: span{1, F}

u3 = u[2]                                                            # Haar on U(3)
mc(u3[:, 0, 0].abs() ** 4, 2 / (k * (k + 1)), "E|U_11|^4")
mc(u3[:, 0, 0].abs() ** 2 * u3[:, 0, 1].abs() ** 2, 1 / (k * (k + 1)), "E|U_11|^2 |U_12|^2")
mc((u3[:, 0, 0] * u3[:, 1, 1] * (u3[:, 0, 1] * u3[:, 1, 0]).conj()).real,
   -1 / ((k - 1) * k * (k + 1)), "E U_11 U_22 conj(U_12 U_21)")

n_s = 20_000                                                         # Theorem 5.10 for a random X
X2 = torch.randn(k * k, k * k, dtype=dt)
trX, trFX = X2.diagonal().sum(), (F @ X2).diagonal().sum()
alpha, beta = (trX - trFX / k) / (k ** 2 - 1), (trFX - trX / k) / (k ** 2 - 1)
UU = torch.einsum('bij,bkl->bikjl', u3[:n_s], u3[:n_s]).reshape(n_s, k * k, k * k)   # U (x) U
T2 = UU @ X2 @ UU.mH
mean, err = T2.mean(0), T2.std(0) / math.sqrt(n_s)
assert ((mean - alpha * torch.eye(k * k, dtype=dt) - beta * F).abs() <= 5 * err).all()
```

## 5.3 Estadística de niveles

Sea $H$ una matriz autoadjunta aleatoria $n\times n$ con valores propios
$\lambda_1\le\dots\le\lambda_n$. Sus **espaciamientos** son $s_i=\lambda_{i+1}-\lambda_i$. Hay dos
conjuntos estándar: el **GOE**, $H=(G+G^T)/2$ con $G$ real de entradas independientes $N(0,1)$, y
el **GUE**, $H=(G+G^*)/\sqrt2$ con $G$ de Ginibre (sección 1.7). En ambos, las entradas diagonales
son $N(0,1)$, las partes real e imaginaria de las entradas por encima de la diagonal son
$N(0,\tfrac12)$, y la densidad de $H$ es proporcional a $e^{-\operatorname{Tr}H^2/2}$.

**Proposición 5.13 (conjetura de Wigner).** Sean $n=2$, $s=\lambda_2-\lambda_1$ y
$S=s/\mathbb Es$. Para el GOE ($\beta=1$) y el GUE ($\beta=2$), $S$ tiene las densidades
$$
p_1(S)=\frac\pi2\,S\,e^{-\pi S^2/4},\qquad p_2(S)=\frac{32}{\pi^2}\,S^2\,e^{-4S^2/\pi},
$$
con funciones de distribución $P(S\le t)=1-e^{-\pi t^2/4}$ y
$P(S\le t)=\operatorname{erf}(2t/\sqrt\pi)-\frac{4t}\pi e^{-4t^2/\pi}$. En particular,
$P(S\le\varepsilon)$ es del orden de $\varepsilon^{\beta+1}$ cuando $\varepsilon\to0$.

*Demostración.* Para $H=\begin{pmatrix}a&b\\\bar b&c\end{pmatrix}$ los valores propios son
$\frac{a+c}2\pm\frac12\sqrt{(a-c)^2+4|b|^2}$, luego $s=\sqrt{(a-c)^2+4|b|^2}$. *GOE:* $x=a-c$ e
$y=2b$ son variables $N(0,2)$ independientes y, en coordenadas polares,
$P(s\le t)=\int_{x^2+y^2\le t^2}\frac{e^{-(x^2+y^2)/4}}{4\pi}\,dx\,dy=\int_0^t\frac r2e^{-r^2/4}dr=1-e^{-t^2/4}$.
Entonces $\mathbb Es=\int_0^\infty P(s>t)\,dt=\int_0^\infty e^{-t^2/4}dt=\sqrt\pi$, y
$P(S\le t)=P(s\le\sqrt\pi t)=1-e^{-\pi t^2/4}$; basta derivar. *GUE:* $x=a-c$,
$y=2\operatorname{Re}b$, $y'=2\operatorname{Im}b$ son $N(0,2)$ independientes y, en coordenadas
esféricas, la densidad de $s=\sqrt{x^2+y^2+y'^2}$ es
$4\pi r^2(4\pi)^{-3/2}e^{-r^2/4}=\frac{r^2e^{-r^2/4}}{2\sqrt\pi}$. Su media es
$\frac1{2\sqrt\pi}\int_0^\infty r^3e^{-r^2/4}dr=\frac{8}{2\sqrt\pi}=\frac4{\sqrt\pi}$. La sustitución
$r=4S/\sqrt\pi$ da $p_2$, e integrando $p_2$ por partes
($\int x^2e^{-ax^2}dx=-\frac{x}{2a}e^{-ax^2}+\frac1{2a}\int e^{-ax^2}dx$ con $a=4/\pi$) se obtiene la
función de distribución. El comportamiento para $\varepsilon$ pequeño se sigue de
$p_\beta(S)\sim cS^\beta$. $\square$

El exponente cuenta dimensiones: una matriz simétrica real $2\times2$ es degenerada si y solo si se
cumplen dos condiciones reales ($a=c$, $b=0$), y una hermítica compleja si y solo si se cumplen
tres. Para todo $n$, la densidad conjunta de los valores propios es proporcional a
$\prod_{i<j}|\lambda_i-\lambda_j|^\beta e^{-\sum_i\lambda_i^2/2}$, y cuando $n\to\infty$ la
distribución de los espaciamientos en el interior del espectro converge a un límite (Gaudin–Mehta)
que está cerca de la conjetura, pero no coincide con ella (M. L. Mehta, *Random Matrices*, 3.ª ed.,
Elsevier, 2004; se enuncia sin demostración).

> **Significado físico.** Wigner propuso que los niveles de energía de los sistemas cuánticos
> complejos (núcleos pesados) se comportan estadísticamente como los valores propios de un
> hamiltoniano aleatorio con la simetría adecuada: simétrico real si hay simetría de inversión
> temporal ($\beta=1$), hermítico complejo en otro caso ($\beta=2$). Los niveles se *repelen*
> entonces: dos niveles de la misma simetría casi nunca se acercan (von Neumann y Wigner, 1929).
> La repulsión es la firma típica del caos cuántico, y los niveles independientes (sin repulsión)
> la de los sistemas integrables o localizados (conjeturas de Bohigas–Giannoni–Schmit y de
> Berry–Tabor).

Los espaciamientos dependen de la densidad local de niveles, que varía a lo largo del espectro;
compararlos con $p_\beta$ requiere un **desdoblamiento** (*unfolding*), el cambio de variables
$\lambda\mapsto\bar N(\lambda)$ mediante una aproximación suave $\bar N$ de la función de conteo,
que hace que el espaciamiento medio sea $1$. La elección de $\bar N$ es una fuente de sesgo. El
cociente de espaciamientos consecutivos lo evita.

**Definición.** El **cociente de espaciamientos** es
$r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})\in[0,1]$ (V. Oganesyan, D. A. Huse, Phys. Rev. B 75
(2007), 155111). No cambia bajo $\lambda\mapsto a\lambda+b$, $a>0$, así que no necesita
desdoblamiento allí donde la densidad es aproximadamente constante a lo largo de dos
espaciamientos.

**Proposición 5.14 (niveles independientes).** Sean $s,s'$ independientes con la misma
distribución exponencial (como dos espaciamientos consecutivos de un proceso de Poisson). Entonces
$r=\min(s,s')/\max(s,s')$ cumple $P(r\le t)=\frac{2t}{1+t}$ en $[0,1]$, tiene densidad
$\frac2{(1+t)^2}$ y
$$
\mathbb Er=2\ln2-1\approx0.3863 .
$$

*Demostración.* Por escalado, la intensidad es $1$. Para $x=s/s'$, $P(x\le t)=\mathbb E\,P(s\le ts'\mid s')=\mathbb E(1-e^{-ts'})=1-\frac1{1+t}=\frac t{1+t}$.
Para $t<1$, los sucesos $\{x\le t\}$ y $\{1/x\le t\}$ son disjuntos, $x$ y $1/x$ tienen la misma
ley ($s,s'$ son intercambiables), y $\{r\le t\}$ es su unión, luego $P(r\le t)=\frac{2t}{1+t}$.
Entonces
$\mathbb Er=\int_0^1\frac{2t}{(1+t)^2}dt=2\Big[\ln(1+t)+\frac1{1+t}\Big]_0^1=2\ln2-1$. $\square$

Lo mismo vale exactamente para $n$ niveles uniformes independientes en un intervalo, porque sus
espaciamientos tienen la ley de $(e_1,\dots,e_{n+1})/\sum_je_j$ para exponenciales independientes
$e_j$, de modo que cada cociente $s_i/s_{i+1}=e_i/e_{i+1}$ (W. Feller, *An Introduction to
Probability Theory and Its Applications*, vol. II, Wiley, 1971, capítulo I; se enuncia sin
demostración). Para el GOE y el GUE en el límite $n\to\infty$, $\langle r\rangle=0.5307$ y
$0.5996$ (valores numéricos; los análogos $3\times3$ de la conjetura de Wigner dan
$4-2\sqrt3\approx0.5359$ y $2\sqrt3/\pi-\frac12\approx0.6027$; Y. Y. Atas, E. Bogomolny,
O. Giraud, G. Roux, Phys. Rev. Lett. 110 (2013), 084101; se enuncia sin demostración).

La función `SpinChain.level_spacing_ratio` (en `chains`) devuelve el cociente medio de cada fila de
niveles. Comprobamos las proposiciones 5.13 y 5.14 y el valor del GOE:

```python
ratio = SpinChain.level_spacing_ratio
n2 = 200_000
alg_r = TypeIAlgebra([2], [2], complex_valued=False, precision='double')     # real symmetric 2 x 2
alg_c = TypeIAlgebra([2], [2], precision='double')
G = torch.randn(n2, 1, 2, 2, dtype=torch.float64)
Gc = torch.randn(n2, 1, 2, 2, dtype=dt)
goe2 = alg_r.operator((G + G.mT) / 2, is_self_adjoint=True)
gue2 = alg_c.operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
cdf = {1: lambda t: 1 - math.exp(-math.pi * t ** 2 / 4),
       2: lambda t: math.erf(2 * t / math.sqrt(math.pi)) - 4 * t / math.pi * math.exp(-4 * t ** 2 / math.pi)}
for beta, H, mean_s in ((1, goe2, math.sqrt(math.pi)), (2, gue2, 4 / math.sqrt(math.pi))):
    lam = H.eigenvalues()[0]                                         # (n2, 2), ascending
    S = (lam[:, 1] - lam[:, 0]) / mean_s
    mc(S, 1.0, f"beta={beta}: E S")
    for t in (0.25, 1.0, 2.0):
        mc((S <= t).double(), cdf[beta](t), f"beta={beta}: P(S <= {t})")

levels = torch.rand(4000, 50, dtype=torch.float64)                   # independent levels
mc(ratio(levels), 2 * math.log(2) - 1, "<r>, independent uniform levels")

n = 64
alg64 = TypeIAlgebra([n], [n], complex_valued=False, precision='double')
G = torch.randn(800, 1, n, n, dtype=torch.float64)
lam = alg64.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()[0]
mc(ratio(lam[:, n // 4: 3 * n // 4]), 0.5307, "<r>, GOE(64), middle half of the spectrum")
```

## 5.4 Por qué no deben mezclarse los sectores

Sea $H=\bigoplus_cH_c\in\mathcal M$ un elemento autoadjunto aleatorio con **bloques
independientes**, por ejemplo bloques GOE. Como operador sobre $\mathbb C^D$, su espectro es la
unión de los espectros de los bloques. Dentro de un bloque los niveles se repelen; mostramos que
los niveles de bloques distintos no lo hacen.

Supongamos que cada bloque tiene una **densidad de niveles** $\rho_c$ continua y acotada:
$\mathbb E\sum_i\varphi(\lambda_i(H_c))=\int\varphi\,\rho_c$ para toda $\varphi$ de Borel acotada
(de modo que $\int\rho_c=k_c$). Para el GOE y el GUE esto se sigue de la densidad conjunta citada
más arriba.

**Teorema 5.15 (no hay repulsión entre sectores).** Para $c\neq d$ y toda
$\psi:\mathbb R^2\to\mathbb R$ de Borel acotada,
$$
\mathbb E\sum_{i,j}\psi\big(\lambda_i(H_c),\lambda_j(H_d)\big)=\iint\psi(x,y)\,\rho_c(x)\rho_d(y)\,dx\,dy .
$$
En consecuencia, si $N_\varepsilon$ es el número de pares de niveles de bloques distintos a
distancia $<\varepsilon$, entonces $\mathbb EN_\varepsilon/(2\varepsilon)\to\sum_{c<d}\int\rho_c\rho_d$
cuando $\varepsilon\to0$, y el menor espaciamiento $\delta$ del espectro fusionado cumple
$$
P(\delta<\varepsilon)\ \ge\ \frac{\mathbb EN_\varepsilon}{\sum_{c<d}k_ck_d}
=\frac{2\varepsilon\sum_{c<d}\int\rho_c\rho_d}{\sum_{c<d}k_ck_d}+o(\varepsilon).
$$

*Demostración.* Los valores propios ordenados son funciones continuas de la matriz, así que
$\Lambda_c=(\lambda_i(H_c))_i$ y $\Lambda_d$ son vectores aleatorios independientes. Por Fubini
(condicionando a $\Lambda_d$),
$\mathbb E\sum_{i,j}\psi(\Lambda_{c,i},\Lambda_{d,j})=\mathbb E\sum_j g(\Lambda_{d,j})$ con
$g(y)=\mathbb E\sum_i\psi(\Lambda_{c,i},y)=\int\psi(x,y)\rho_c(x)\,dx$, que es acotada y de Borel;
basta aplicar la definición de $\rho_d$ a $g$. Con $\psi(x,y)=1\{|x-y|<\varepsilon\}$,
$$
\frac{\mathbb EN_\varepsilon}{2\varepsilon}=\sum_{c<d}\int\rho_c(x)\Big(\frac1{2\varepsilon}\int_{x-\varepsilon}^{x+\varepsilon}\rho_d(y)\,dy\Big)dx ;
$$
el paréntesis tiende a $\rho_d(x)$ por continuidad y está acotado por $\sup\rho_d$, y $\rho_c$ es
integrable, así que la convergencia dominada da el límite. Por último,
$N_\varepsilon\le\sum_{c<d}k_ck_d$, luego
$\mathbb EN_\varepsilon\le\big(\sum_{c<d}k_ck_d\big)P(N_\varepsilon\ge1)$, y $N_\varepsilon\ge1$
implica $\delta<\varepsilon$. $\square$

Así pues, para dos bloques con espectros solapados, $P(\delta<\varepsilon)$ es al menos del orden
de $\varepsilon$, mientras que dentro de un bloque GOE $2\times2$ vale
$1-e^{-\varepsilon^2/4}\le\varepsilon^2/4$ (proposición 5.13, espaciamiento sin normalizar), y para
un bloque GOE $n\times n$ es del orden de $\varepsilon^2$ (por la densidad conjunta; no se demuestra
aquí). La prueba de repulsión de niveles falla para el espectro fusionado aunque cada bloque sea un
GOE perfecto.

*Muchos sectores (una deducción, no una demostración).* Para espectros estacionarios, la
probabilidad $\mathcal G(s)$ de no encontrar ningún nivel en un intervalo de longitud $s$ determina
la densidad de espaciamientos, $p(s)=\mathcal G''(s)$ para espaciamiento medio unidad (Mehta,
*op. cit.*). Los espectros independientes se superponen multiplicando estas probabilidades de
hueco. Si $C$ sectores aportan cada uno una fracción $1/C$ de los niveles y cada uno tiene la
probabilidad de hueco $\mathcal G_1$ en sus propias unidades desdobladas, entonces
$\mathcal G(s)=\mathcal G_1(s/C)^C$, con $\mathcal G_1(x)=1-x+\int_0^x(x-t)p_1(t)\,dt$. Derivando dos
veces en $0$,
$$
p(0)=1-\frac1C+\frac{p_1(0)}C=1-\frac1C ,
$$
el análogo del teorema 5.15 para espectros desdoblados, y $\mathcal G_1(s/C)^C\to e^{-s}$, la
probabilidad de hueco de un proceso de Poisson: cuando $C\to\infty$, el espectro fusionado se
parece al de niveles independientes (N. Rosenzweig, C. E. Porter, Phys. Rev. 120 (1960), 1698–1714;
M. V. Berry, M. Robnik, J. Phys. A 17 (1984), 2413–2421; las versiones rigurosas del límite de
Poisson son teoremas límite para superposiciones de procesos puntuales dispersos, véase
B. Grigelionis, Theory Probab. Appl. 8 (1963)).

En la biblioteca, los bloques están separados desde el principio: `eigenvalues()` devuelve un
espectro por sector, y mezclarlos requiere un `torch.cat` explícito. Ponemos $C$ bloques GOE
independientes de tamaño $64$ en una misma álgebra y comparamos el cociente medio por sector con el
del espectro fusionado:

```python
k, merged_r = 64, {}
for C in (1, 2, 4, 8, 16):
    alg_C = TypeIAlgebra([k] * C, [k] * C, complex_valued=False, precision='double')
    G = torch.randn(800 // C, C, k, k, dtype=torch.float64)
    lam = alg_C.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()    # one spectrum per sector
    per_sector = torch.stack([ratio(l[:, k // 4: 3 * k // 4]) for l in lam])
    merged = torch.cat(lam, dim=-1).sort(dim=-1).values                       # the spectrum on C^D
    r_m = ratio(merged[:, k * C // 4: 3 * k * C // 4])
    merged_r[C] = r_m.mean().item()
    print(f"C = {C:2d}: <r> per sector {per_sector.mean():.4f}, merged {merged_r[C]:.4f}"
          f" +- {r_m.std() / math.sqrt(len(r_m)):.4f}")
    assert abs(per_sector.mean() - 0.5307) < 0.01                    # every sector is GOE
assert merged_r[1] > merged_r[2] > merged_r[4] > merged_r[16]
assert abs(merged_r[16] - (2 * math.log(2) - 1)) < 0.01              # close to independent levels
```

Cada sector tiene $\langle r\rangle\approx0.53$; el espectro fusionado baja a $0.42$ con dos
sectores y a $0.39$, el valor de los niveles independientes, con $8$ y $16$.

Lo mismo ocurre en un modelo físico. La cadena de Heisenberg con campos aleatorios (capítulo 9)
conserva la magnetización, así que su hamiltoniano está en $\bigoplus_{N=0}^LM_{\binom LN}$; para
desorden débil, los niveles dentro de cada sector tienen estadística GOE (numéricamente; A. Pal,
D. A. Huse, Phys. Rev. B 82 (2010), 174411). Sus bloques *no* son independientes —todos dependen de
los mismos campos aleatorios—, así que el teorema 5.15 no se aplica tal como está enunciado; el
experimento muestra el mismo efecto.

```python
chain = SpinChain(10, complex_valued=False)                         # sectors N = 0, ..., 10
lam = chain.random_field_heisenberg(W=0.5, batch_size=20).eigenvalues()
k5 = chain.sector_dim(5)                                            # 252
r_sector = ratio(lam[5][:, k5 // 4: 3 * k5 // 4]).mean().item()
spectrum = torch.cat(lam, dim=-1).sort(dim=-1).values               # all 1024 levels mixed
r_all = ratio(spectrum[:, 256:768]).mean().item()
print(f"Heisenberg chain, L = 10: <r> in sector N = 5: {r_sector:.4f}, all sectors mixed: {r_all:.4f}")
assert r_sector > 0.49 and r_all < 0.42
```

Dentro del sector $N=5$ el cociente es $0.52$ (un valor de tamaño finito ligeramente inferior a
$0.5307$); el espectro con todos los sectores mezclados da $0.39$.

> **Significado físico.** Los niveles con números cuánticos conservados distintos no interactúan,
> así que se cruzan libremente; solo se repelen los niveles de la misma simetría. Un espectro que
> mezcla sectores de simetría muestra, por tanto, una estadística de tipo Poisson incluso cuando la
> dinámica en cada sector es completamente caótica, y una simetría no resuelta puede confundirse
> con integrabilidad o localización. La regla del oficio —resolver todas las simetrías antes de
> calcular la estadística de niveles— equivale a decir que hay que calcular en el álgebra
> $\bigoplus_cM_{k_c}$, no en $M_D$.

## 5.5 Con qué comparar

Unitarias de Haar: `scipy.stats.unitary_group` (la QR corregida del teorema 5.3) y `ortho_group`,
las unitarias aleatorias de Haar de QuTiP, y cualquier QR de `numpy`/`torch` seguida de la
corrección de fases. Cálculo de Weingarten: los paquetes simbólicos IntU (Mathematica;
Z. Puchała, J. A. Miszczak), RTNI (Mathematica y Python; M. Fukuda, R. König, I. Nechita, J. Phys.
A 52 (2019), 425303) y el paquete de Python Haarpy calculan momentos de cualquier orden de forma
simbólica; las fórmulas demostradas en este capítulo son su caso de orden dos. La estadística de
niveles suele calcularse con unas pocas líneas de `numpy`, y los códigos de diagonalización exacta
como QuSpin proporcionan los espectros sector a sector. Lo que añade la biblioteca es el grupo
unitario de un álgebra como un único objeto: bloques de Haar independientes de todos los tamaños
generados en llamadas por lotes, el twirl identificado con la esperanza condicional
`channels.center_expectation`, y espectros devueltos por sector, de modo que mezclar sectores nunca
es la opción por defecto.

## Ejercicios

1. Sea $U$ de Haar en $U(\mathcal M)$ y $\operatorname{Tr}U=\sum_c\operatorname{Tr}U_c$. Demuestre
   que $\mathbb E|\operatorname{Tr}U|^2=C$ y que $\mathbb E|\operatorname{Tr}U|^4=\sum_cm(k_c)+2C(C-1)$
   con $m(1)=1$, $m(k)=2$ para $k\ge2$; compárelo con $1$ y $2$ para una unitaria de Haar de
   $U(D)$, $D\ge2$. Compruébelo para $\mathcal M=M_1\oplus M_2\oplus M_3$.
2. Usando el teorema 5.12, demuestre que $\mathbb E|U_{11}|^{2n}=\frac{n!}{k(k+1)\cdots(k+n-1)}$
   para $n\le k$. Calcule numéricamente la matriz de Weingarten de orden $3$ para $k=3$ y compruebe
   que $\mathbb E|U_{11}|^6=\frac1{10}$.
3. El **conjunto ortogonal circular** (`measure='coe'`) es $S=W^TW$ con $W$ de Haar en $U(k)$.
   Demuestre que $S$ es unitaria y simétrica, que su ley es invariante bajo $S\mapsto V^TSV$ para
   toda $V\in U(k)$, y que la densidad media de sus fases propias es uniforme. Compruebe
   numéricamente que sus fases propias tienen $\langle r\rangle\approx0.53$, como el GOE.

<details>
<summary>Soluciones</summary>

1. Escribamos $t_c=\operatorname{Tr}U_c$; los $t_c$ son independientes (corolario 5.4).
   Multiplicar $U_c$ por $e^{i\theta}$ no cambia su ley y multiplica $t_c$ por $e^{i\theta}$; así,
   un momento $\mathbb E[t_at_b\bar t_c\bar t_d]$ se anula salvo que cada sector aparezca tantas
   veces en $(a,b)$ como en $(c,d)$, es decir, $\{a,b\}=\{c,d\}$ como multiconjuntos, y del mismo
   modo $\mathbb E[t_a\bar t_b]=0$ para $a\neq b$. Por la proposición 5.7,
   $\mathbb E[(U_c)_{ii}\overline{(U_c)_{jj}}]=\langle e_i,\mathbb E[U_ce_ie_j^*U_c^*]e_j\rangle=\delta_{ij}/k_c$,
   luego $\mathbb E|t_c|^2=1$ y $\mathbb E|\operatorname{Tr}U|^2=C$. En $\mathbb E|\sum_at_a|^4$, los
   términos con $a=b$ dan $\sum_a\mathbb E|t_a|^4$, y los términos con $a\neq b$ dan
   $2\sum_{a\neq b}\mathbb E|t_a|^2\,\mathbb E|t_b|^2=2C(C-1)$ (dos elecciones de $(c,d)$). Para
   $k=1$, $|t|=1$. Para $k\ge2$,
   $\mathbb E|t|^4=\sum_{i,j,i',j'}\mathbb E[U_{ii}U_{jj}\overline{U_{i'i'}U_{j'j'}}]$; en el
   teorema 5.10, los dos términos con $\sigma=\pi$ imponen $\{i',j'\}=\{i,j\}$ de una manera y dan
   cada uno $k^2\mathrm{Wg}(e)$ tras sumar, y los dos con $\sigma\neq\pi$ imponen $i=j=i'=j'$ y dan
   cada uno $k\,\mathrm{Wg}(\tau)$: en total, $2\big(\frac{k^2}{k^2-1}-\frac1{k^2-1}\big)=2$. El caso
   $C=1$, $k=D$ da $1$ y $2$ para $U(D)$. Para $M_1\oplus M_2\oplus M_3$:
   $\mathbb E|\operatorname{Tr}U|^4=1+2+2+12=17$, frente a $2$ para una unitaria de Haar de $U(6)$.

```python
tr = U.trace                                                        # Tr U = sum_c Tr U_c, shape (B,)
mc(tr.abs() ** 2, 3.0, "E|Tr U|^2")
mc(tr.abs() ** 4, 17.0, "E|Tr U|^4")
```

2. Con todos los índices iguales a $1$, toda delta vale $1$, así que
   $\mathbb E|U_{11}|^{2n}=\sum_{\sigma,\pi}\mathrm{Wg}(\sigma\pi^{-1})=\mathbf 1^TW\mathbf 1$, donde
   $W=\Gamma^{-1}=(\mathrm{Wg}(\sigma^{-1}\pi))$ y $\Gamma=(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)})$.
   Aquí hemos usado que $\mathrm{Wg}$ solo depende del tipo de ciclos, y que
   $\sigma\pi^{-1}=\sigma(\pi^{-1}\sigma)\sigma^{-1}$ es conjugada de
   $\pi^{-1}\sigma=(\sigma^{-1}\pi)^{-1}$, que tiene el tipo de ciclos de $\sigma^{-1}\pi$. Cada
   suma por filas de $\Gamma$ es $c_n(k)=\sum_{\rho\in S_n}k^{\#\mathrm{cycles}(\rho)}$, luego
   $\Gamma\mathbf 1=c_n(k)\mathbf 1$, $W\mathbf 1=\mathbf 1/c_n(k)$ y
   $\mathbf 1^TW\mathbf 1=n!/c_n(k)$. Por último, $c_n(x)=x(x+1)\cdots(x+n-1)$ por inducción: una
   permutación de $\{1,\dots,n\}$ se obtiene de una de $\{1,\dots,n-1\}$ o bien añadiendo $n$ como
   punto fijo (un ciclo más, factor $x$), o bien insertando $n$ tras uno de los $n-1$ elementos en
   su ciclo ($n-1$ maneras, mismo número de ciclos); esto es una biyección, luego
   $c_n(x)=(x+n-1)c_{n-1}(x)$.

```python
import itertools

def n_cycles(p):                                                    # number of cycles of a permutation
    seen, count = set(), 0
    for i in range(len(p)):
        count += i not in seen
        while i not in seen:
            seen.add(i)
            i = p[i]
    return count

perms = list(itertools.permutations(range(3)))                      # S_3; s.index(j) = s^{-1}(j)
Gram = torch.tensor([[3.0 ** n_cycles([s.index(p[i]) for i in range(3)]) for p in perms] for s in perms],
                    dtype=torch.float64)                            # k^{#cycles(s^{-1} p)}, k = 3
Wg3 = torch.linalg.inv(Gram)                                        # the Weingarten matrix, n = k = 3
assert abs(Wg3.sum().item() - 6 / (3 * 4 * 5)) < 1e-12
mc(u3[:, 0, 0].abs() ** 6, 0.1, "E|U_11|^6")
```

3. $S^*S=W^*\bar WW^TW=1$, ya que $W^T$ es unitaria, y $S^T=W^TW=S$. Para $V\in U(k)$, $WV$ tiene
   la ley de $W$ (invariancia por la derecha, lema 5.1), y $(WV)^T(WV)=V^TSV$. Con
   $V=e^{i\theta/2}1$, $e^{i\theta}S$ tiene la ley de $S$, así que la ley del conjunto de fases
   propias es invariante por rotaciones, y también lo es la densidad media, que es, por tanto,
   uniforme. En particular, los espaciamientos no necesitan desdoblamiento. La estadística local del
   COE coincide con la del GOE para $k$ grande (F. J. Dyson, J. Math. Phys. 3 (1962), 140–156;
   Mehta, *op. cit.*; se enuncia sin demostración).

```python
alg_s = TypeIAlgebra([32], [32], precision='double')
S = alg_s.random_unitary_operator(400, measure='coe').matrix[:, 0]
assert close(S, S.mT) and close(S @ S.mH, torch.eye(32, dtype=dt))
phases = torch.linalg.eigvals(S).angle().sort(dim=-1).values
s = torch.diff(torch.cat([phases, phases[:, :1] + 2 * math.pi], dim=-1), dim=-1)   # around the circle
r = torch.minimum(s, s.roll(-1, -1)) / torch.maximum(s, s.roll(-1, -1))
mc(r.mean(-1), 0.5307, "<r>, COE(32) eigenphases")
```

</details>

**Siguiente:** el [capítulo 6](06_quantum_channels.md) trata los canales cuánticos: las
representaciones de Kraus y de Stinespring, los duales respecto de una traza y las aplicaciones
entre sectores. La esperanza condicional $E$, que aquí aparece como un promedio de Haar, es el
primer ejemplo.
