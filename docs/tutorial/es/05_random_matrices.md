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
condicional sobre el centro, la fórmula de Weingarten de orden dos, la conjetura de Wigner, el
cociente de espaciamientos y su valor exacto para matrices gaussianas $3\times3$, y por qué la
estadística de niveles calculada sobre un espectro con sectores mezclados se parece a la de niveles
independientes.

En todo el capítulo, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ con $D=\sum_ck_c$, $P_c$ son las
proyecciones centrales minimales y $E$ es la esperanza condicional sobre el centro de la
proposición 2.8. $U(k)$ es el grupo de las matrices unitarias $k\times k$ y
$U(\mathcal M)=\{U\in\mathcal M:U^*U=1\}$ es el grupo unitario de $\mathcal M$. Las variables
aleatorias viven en un espacio de probabilidad con esperanza $\mathbb E$; la esperanza de una matriz
aleatoria se toma entrada a entrada. Una variable *gaussiana compleja estándar* tiene partes real e
imaginaria independientes $\sim N(0,\tfrac12)$, de modo que $\mathbb E|z|^2=1$ (esto es lo que
genera `torch.randn` para un dtype complejo, según su documentación; el primer bloque de código
comprueba los segundos momentos); una *matriz de Ginibre* tiene entradas independientes gaussianas
complejas estándar.

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
continuas son iguales (corolario A.2: las funciones continuas $\max(0,1-n\operatorname{dist}(\cdot,K))$
decrecen hacia la indicatriz de un cerrado $K$, y los cerrados determinan la medida por el lema de
Dynkin, lema A.1). Si $\mu_1,\mu_2$ son invariantes por la izquierda, entonces
$\nu(B)=\mu_2(B^{-1})$ es invariante por la derecha, ya que $(Bg)^{-1}=g^{-1}B^{-1}$; luego
$\mu_1=\nu=\mu_2$, y $\mu_1$ es invariante por la derecha y por inversión. $\square$

**Definición.** La **medida de Haar** de $G$ es su medida de probabilidad de Borel invariante por
la izquierda. Un elemento aleatorio con esta ley está **distribuido según la medida de Haar**.

El lema 5.1 no dice nada sobre la existencia. El curso solo necesita las medidas de Haar de $U(k)$,
$U(1)$, $U(\mathcal M)$, $SU(k)$, $O(k)$ y de los grupos diagonales $T$ del lema C.2, y cada una de
ellas se construye explícitamente: $U(k)$ en el teorema 5.3, $U(\mathcal M)$ en el corolario 5.4,
$SU(k)$ a continuación de este, $O(k)$ en el teorema B.17 y $T$ en el lema C.2 (véase la
sección B.7).

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
`torch.linalg.qr`), entonces $Q=Q'\Lambda$ con $\Lambda=\operatorname{diag}(R'_{ii}/|R'_{ii}|)$
está distribuida según la medida de Haar.

*Demostración.* *Invertibilidad.* $|\det Z|^2$ es un polinomio real en las $2k^2$ coordenadas
reales de $Z$, no idénticamente nulo (vale $1$ en $Z=1$), y $Z$ tiene la densidad
$\pi^{-k^2}e^{-\|Z\|_2^2}$ (proposición A.16(b)); luego $P(\det Z=0)=0$ por el corolario A.7(b).

*Invariancia de la ley gaussiana.* Para $V\in U(k)$, la aplicación $Z\mapsto VZ$ es lineal real y
conserva $\|Z\|_2$ (sección 1.4), así que es una aplicación ortogonal de
$M_k(\mathbb C)\cong\mathbb R^{2k^2}$, que conserva la medida de Lebesgue (lema A.3(b)) y la
densidad: $VZ$ tiene la misma ley que $Z$ (proposición A.16(b)).

*Equivariancia.* Si $Z=QR$, entonces $VZ=(VQ)R$ es la descomposición del lema 5.2 de $VZ$, luego
$Q(VZ)=VQ(Z)$. Como $Q$ es continua sobre las matrices invertibles (lema 5.2), $Q(Z)$ es un
elemento aleatorio de $U(k)$ y, para todo boreliano $B\subseteq U(k)$,
$P(VQ(Z)\in B)=P(Q(VZ)\in B)=P(Q(Z)\in B)$: la ley de $Q(Z)$ es invariante por la izquierda y,
por el lema 5.1, es la medida de Haar (esto prueba también que la medida de Haar de $U(k)$
existe). Por último, $Q'R'=(Q'\Lambda)(\Lambda^*R')$, y $\Lambda^*R'$ es triangular superior con
diagonal $|R'_{ii}|>0$, luego $Q'\Lambda=Q(Z)$ por unicidad. $\square$

La misma demostración con matrices gaussianas reales da la medida de Haar del grupo ortogonal
$O(k)$, con la corrección $\operatorname{diag}(\operatorname{sgn}R'_{ii})$; esto se desarrolla en el
teorema B.17. Sin la corrección $\Lambda$, la ley de $Q'$ depende del convenio de signos de la
rutina QR. Para la rutina a la que llama `torch.linalg.qr` (LAPACK), el código siguiente encuentra
$\mathbb E\operatorname{Re}\operatorname{Tr}Q'$ a muchos errores estándar de $0$, mientras que una
unitaria de Haar tiene $\mathbb EU=0$ (corolario 5.5): así que este $Q'$ no está distribuido según
la medida de Haar. Esta es una afirmación sobre esa rutina, comprobada por el código, no un teorema
(la receta con $\Lambda$ se debe a F. Mezzadri, *How to generate random matrices from the classical
compact groups*, Notices AMS 54 (2007), 592–604).

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

El grupo $U(\mathcal M)$ es mucho más pequeño que $U(D)$: para $C\ge2$, una unitaria de Haar $U$ de
$U(D)$ pertenece a $\mathcal M$ con probabilidad cero. En efecto, por el teorema 5.3 y el lema 5.1,
$U$ tiene la ley del factor $Q$ de una matriz de Ginibre $Z$ de tamaño $D\times D$, y por la
construcción de Gram–Schmidt del lema 5.2 la primera columna de $Q$ es $z/\|z\|$ para la primera
columna $z$ de $Z$, un vector gaussiano complejo estándar de $\mathbb C^D$ con densidad
(proposición A.16(a)). Para cada $i$, $\{z_i=0\}$ es un subespacio propio de $\mathbb C^D$, luego
$P(z_i=0)=0$ (corolario A.7(c)), y todas las entradas de la primera columna son no nulas casi
seguramente. Un elemento de $\mathcal M$ tiene entradas nulas en la primera columna fuera de las
filas del primer bloque, y estas filas existen cuando $C\ge2$.

La medida de Haar de $SU(k)$ se obtiene con el mismo espíritu: si $U$ es de Haar en $U(k)$ y
$\theta=\arg\det U\in(-\pi,\pi]$ (una función de Borel de $U$), entonces $V=e^{-i\theta/k}U\in SU(k)$,
y para $W\in SU(k)$, $\det(WU)=\det U$ da $V(WU)=WV(U)$; como $WU$ tiene la ley de $U$, la ley de
$V$ es invariante por la izquierda en $SU(k)$ y, por tanto, es la medida de Haar de $SU(k)$ por el
lema 5.1.

**Corolario 5.5 (primeros momentos de las entradas).** Si $U$ está distribuida según la medida de
Haar en $U(k)$, entonces $\mathbb EU_{ij}=0$ y $\mathbb E|U_{ij}|^2=1/k$.

*Demostración.* $-1\in U(k)$, luego $-U$ tiene la ley de $U$ y $\mathbb EU=-\mathbb EU$. Para una
matriz de permutación $P$, $PU$ tiene la ley de $U$, así que $\mathbb E|U_{ij}|^2$ no depende de
$i$; como $\sum_i|U_{ij}|^2=1$, vale $1/k$. $\square$

La biblioteca genera unitarias de Haar de $U(\mathcal M)$ con `random_unitary_operator`, mediante
una única descomposición QR por lotes por cada tamaño de bloque, para todos los sectores y todas
las muestras a la vez (en una GPU, mediante Gram–Schmidt, que da la misma $Q$ por la unicidad del
lema 5.2). La medida `'haar'` es el teorema 5.3 en cada bloque; para un álgebra compleja,
`'haar_su'` multiplica por $e^{-i\theta/k}$ como antes; `'coe'` es el conjunto $S=W^TW$ del
ejercicio 3 (las demás medidas de la biblioteca no se usan en este curso). Comprobamos los
corolarios 5.4 y 5.5, y que el muestreador de la biblioteca es exactamente la QR corregida del
teorema 5.3:

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

z = torch.randn(B, dtype=dt)                                         # torch.randn draws standard complex Gaussians:
mc(z.real ** 2, 0.5, "E (Re z)^2")                                   # Re z, Im z ~ N(0, 1/2), uncorrelated
mc(z.imag ** 2, 0.5, "E (Im z)^2")
mc(z.real * z.imag, 0.0, "E Re z Im z")

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

*Observación numérica.* Sin la corrección, el código encuentra
$\mathbb E\operatorname{Re}\operatorname{Tr}Q'\approx-0.84$ para $k=2$, mientras que una unitaria de
Haar tiene $\mathbb EU=0$.

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
siquiera tiene por qué estar en $\mathcal M$ (para $C\ge2$, $X=P_1$ y la rotación $U$ de ángulo
$\pi/4$ en el plano de $e_1$ y $e_{k_1+1}$, $UP_1U^*$ tiene la entrada $\frac12$ en la posición
$(1,k_1+1)$).

> **Significado físico** (interpretación física, no un teorema de este curso). Una unitaria aleatoria de Haar del álgebra es una evolución aleatoria que
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
de Schur–Weyl, teorema B.6(c)). El apéndice B lo lleva a cabo.

**Teorema 5.12 (fórmula de Weingarten de cualquier orden).** Sean $n\ge1$ y $U$ de Haar en $U(k)$.
Entonces
$$
\mathbb E\big[U_{i_1j_1}\cdots U_{i_nj_n}\overline{U_{i'_1j'_1}\cdots U_{i'_nj'_n}}\big]
=\sum_{\sigma,\pi\in S_n}\prod_{m=1}^n\delta_{i_mi'_{\sigma(m)}}\delta_{j_mj'_{\pi(m)}}\;\mathrm{Wg}(\sigma\pi^{-1},k),
$$
donde $\mathrm{Wg}(\cdot,k)$ es una función sobre $S_n$ que solo depende del tipo de ciclos. Para
$n\le k$, la matriz $\big(\mathrm{Wg}(\sigma^{-1}\pi,k)\big)_{\sigma,\pi}$ es la inversa de la
matriz de Gram $\big(k^{\#\mathrm{cycles}(\sigma^{-1}\pi)}\big)_{\sigma,\pi}$; para $n>k$ esta
matriz es singular y se toma su pseudoinversa. Los momentos con distinto número de factores $U$ y
$\overline U$ se anulan.

*Demostración.* Es el teorema B.13 con $d=k$: la parte (a) da $W_{\sigma\pi}=\mathrm{Wg}(\sigma^{-1}\pi)$
con $\mathrm{Wg}$ una función del tipo de ciclos, la parte (b) la invertibilidad de la matriz de
Gram si y solo si $n\le k$, la parte (d) la fórmula de los momentos y la parte (e) la anulación de
los momentos no equilibrados. $\square$

(La fórmula se debe a D. Weingarten, J. Math. Phys. 19 (1978), 999–1001, y a B. Collins, Int. Math.
Res. Not. 2003, 953–982; véase también B. Collins, P. Śniady, Comm. Math. Phys. 264 (2006), 773–795.)

Comprobamos la proposición 5.7 frente a `channels.center_expectation`, el lema 5.8 mediante álgebra
lineal (el conmutante de tres $V\otimes V$ muestreados, calculado como un núcleo, contiene a $F$ y
tiene dimensión $2$, así que, por el lema 5.8, ya es $\operatorname{span}\{1,F\}$), y el
teorema 5.10 y el corolario 5.11 mediante Monte Carlo con el bloque $U(3)$ de las muestras
anteriores:

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
son $N(0,1)$, las entradas por encima de la diagonal son $N(0,\tfrac12)$ en el GOE y tienen partes
real e imaginaria independientes $N(0,\tfrac12)$ en el GUE, las entradas $h_{ij}$, $i\le j$, son
independientes, y la densidad de $H$ es proporcional a $e^{-\operatorname{Tr}H^2/2}$
(proposición A.17; la proporcionalidad se debe a $\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$).

**Proposición 5.13 (conjetura de Wigner).** Sean $n=2$, $s=\lambda_2-\lambda_1$ y
$S=s/\mathbb Es$. Para el GOE ($\beta=1$) y el GUE ($\beta=2$), $S$ tiene las densidades
$$
p_1(S)=\frac\pi2\,S\,e^{-\pi S^2/4},\qquad p_2(S)=\frac{32}{\pi^2}\,S^2\,e^{-4S^2/\pi},
$$
con funciones de distribución $P(S\le t)=1-e^{-\pi t^2/4}$ y
$P(S\le t)=\operatorname{erf}(2t/\sqrt\pi)-\frac{4t}\pi e^{-4t^2/\pi}$. En particular,
$P(S\le\varepsilon)$ es del orden de $\varepsilon^{\beta+1}$ cuando $\varepsilon\to0$.

*Demostración.* Para $H=\begin{pmatrix}a&b\\\bar b&c\end{pmatrix}$ los valores propios son
$\frac{a+c}2\pm\frac12\sqrt{(a-c)^2+4|b|^2}$, luego $s=\sqrt{(a-c)^2+4|b|^2}$. Por la
proposición A.17, $a$, $c$ y las coordenadas reales de $b$ son independientes, $a,c\sim N(0,1)$ y
cada coordenada real de $b$ es $N(0,\frac12)$; así, $a-c\sim N(0,2)$ (proposición A.14(c)) y
$2\operatorname{Re}b$, $2\operatorname{Im}b$ son $N(0,2)$, todas independientes. *GOE:* $x=a-c$ e
$y=2b$ son variables $N(0,2)$ independientes y, en coordenadas polares (corolario A.19 con $n=2$,
$s_2=2\pi$),
$P(s\le t)=\int_{x^2+y^2\le t^2}\frac{e^{-(x^2+y^2)/4}}{4\pi}\,dx\,dy=\int_0^t\frac r2e^{-r^2/4}dr=1-e^{-t^2/4}$.
Entonces $\mathbb Es=\int_0^\infty P(s>t)\,dt=\int_0^\infty e^{-t^2/4}dt=\sqrt\pi$, y
$P(S\le t)=P(s\le\sqrt\pi t)=1-e^{-\pi t^2/4}$; basta derivar. *GUE:* $x=a-c$,
$y=2\operatorname{Re}b$, $y'=2\operatorname{Im}b$ son $N(0,2)$ independientes y, por el
teorema A.20 ($n=3$, $s_3=4\pi$), la densidad de $s=\sqrt{x^2+y^2+y'^2}$ es
$4\pi r^2(4\pi)^{-3/2}e^{-r^2/4}=\frac{r^2e^{-r^2/4}}{2\sqrt\pi}$. Su media es
$\frac1{2\sqrt\pi}\int_0^\infty r^3e^{-r^2/4}dr=\frac{8}{2\sqrt\pi}=\frac4{\sqrt\pi}$. La sustitución
$r=4S/\sqrt\pi$ da $p_2$, e integrando $p_2$ por partes
($\int x^2e^{-ax^2}dx=-\frac{x}{2a}e^{-ax^2}+\frac1{2a}\int e^{-ax^2}dx$ con $a=4/\pi$) se obtiene la
función de distribución. El comportamiento para $\varepsilon$ pequeño se sigue de
$p_\beta(S)\sim cS^\beta$. $\square$

El exponente cuenta dimensiones: una matriz simétrica real $2\times2$ es degenerada si y solo si se
cumplen dos condiciones reales ($a=c$, $b=0$), y una hermítica compleja si y solo si se cumplen
tres. Para todo $n$, los valores propios ordenados $\lambda_1<\dots<\lambda_n$ del GOE ($\beta=1$)
y del GUE ($\beta=2$) tienen la densidad conjunta
$$
Z_{n,\beta}^{-1}\prod_{i<j}(\lambda_j-\lambda_i)^\beta\,e^{-\sum_i\lambda_i^2/2}
$$
(teorema C.7, demostrado en el apéndice C con la fórmula de integración de Weyl), y el
espaciamiento mínimo $\delta$ cumple $c\,\varepsilon^{\beta+1}\le P(\delta<\varepsilon)\le C\varepsilon^{\beta+1}$
para $0<\varepsilon\le1$, con constantes $0<c\le C$ que dependen de $n$ y de $\beta$
(proposición C.9). El límite $n\to\infty$ no se trata en este curso (véase M. L. Mehta, *Random
Matrices*, 3.ª ed., Elsevier, 2004).

> **Significado físico** (interpretación física, no un teorema de este curso). Wigner propuso que
> los niveles de energía de los sistemas cuánticos complejos (núcleos pesados) se comportan
> estadísticamente como los valores propios de un hamiltoniano aleatorio con la simetría adecuada:
> simétrico real si hay simetría de inversión temporal ($\beta=1$), hermítico complejo en otro caso
> ($\beta=2$). Los niveles se *repelen* entonces: dos niveles de la misma simetría casi nunca se
> acercan (von Neumann y Wigner, 1929). La repulsión es la firma típica del caos cuántico, y los
> niveles independientes (sin repulsión) la de los sistemas integrables o localizados; estas son
> las conjeturas de Bohigas–Giannoni–Schmit y de Berry–Tabor de la física, que ni se demuestran ni
> se usan en este curso.

Los espaciamientos dependen de la densidad local de niveles, que varía a lo largo del espectro;
compararlos con $p_\beta$ requiere un **desdoblamiento** (*unfolding*), el cambio de variables
$\lambda\mapsto\bar N(\lambda)$ mediante una aproximación suave $\bar N$ de la función de conteo,
que hace que el espaciamiento medio sea $1$. El resultado depende de la elección de $\bar N$. El
cociente de espaciamientos consecutivos evita esta elección.

**Definición.** El **cociente de espaciamientos** es
$r_i=\min(s_i,s_{i+1})/\max(s_i,s_{i+1})\in[0,1]$ (V. Oganesyan, D. A. Huse, Phys. Rev. B 75
(2007), 155111). No cambia bajo $\lambda\mapsto a\lambda+b$, $a>0$, así que no depende del
espaciamiento medio local; por eso se usa sin desdoblamiento.

Un cociente necesita tres niveles, así que las matrices gaussianas más pequeñas con un cociente de
espaciamientos son las $3\times3$; para ellas, el cociente medio puede calcularse exactamente, como
se hizo con la conjetura de Wigner para $2\times2$.

**Proposición 5.14 (cocientes medios de espaciamientos).** (a) *Niveles independientes.* Sean
$s,s'$ independientes con la misma distribución exponencial (como dos espaciamientos consecutivos de
los niveles $e_0+\dots+e_{k-1}$ construidos a partir de exponenciales independientes $e_j$,
sección A.7). Entonces $r=\min(s,s')/\max(s,s')$ cumple $P(r\le t)=\frac{2t}{1+t}$ en $[0,1]$,
tiene densidad $\frac2{(1+t)^2}$ y
$$
\mathbb Er=2\ln2-1\approx0.3863 .
$$
(b) *GOE y GUE $3\times3$.* Sea $H$ una matriz $3\times3$ del GOE ($\beta=1$) o del GUE
($\beta=2$) con valores propios $\lambda_1<\lambda_2<\lambda_3$, espaciamientos
$s_1=\lambda_2-\lambda_1$, $s_2=\lambda_3-\lambda_2$, y $r=\min(s_1,s_2)/\max(s_1,s_2)$. Entonces
$x=s_2/s_1$ tiene una densidad proporcional a
$$
f_\beta(x)=\frac{(x+x^2)^\beta}{(1+x+x^2)^{1+3\beta/2}}\qquad\text{on }(0,\infty),
$$
$r$ tiene la densidad $f_\beta/I_0$ sobre $(0,1)$, donde $I_0=\int_0^1f_\beta$, y
$$
\mathbb Er=4-2\sqrt3\approx0.5359\quad(\text{GOE}),\qquad
\mathbb Er=\frac{2\sqrt3}\pi-\frac12\approx0.6027\quad(\text{GUE}).
$$

*Demostración.* (a) Por escalado, la intensidad es $1$. Para $x=s/s'$, por la independencia y el
teorema de Tonelli,
$P(x\le t)=\int_0^\infty e^{-y}(1-e^{-ty})\,dy=1-\frac1{1+t}=\frac t{1+t}$.
Para $t<1$, los sucesos $\{x\le t\}$ y $\{1/x\le t\}$ son disjuntos, $x$ y $1/x$ tienen la misma
ley ($s,s'$ son intercambiables), y $\{r\le t\}$ es su unión, luego $P(r\le t)=\frac{2t}{1+t}$.
Entonces
$\mathbb Er=\int_0^1\frac{2t}{(1+t)^2}dt=2\Big[\ln(1+t)+\frac1{1+t}\Big]_0^1=2\ln2-1$.

(b) *La ley de $x$.* Por el teorema C.7, el vector $\lambda=(\lambda_1,\lambda_2,\lambda_3)$ tiene la
densidad $Z^{-1}\Delta(\lambda)^\beta e^{-|\lambda|^2/2}$ sobre $C=\{\lambda_1<\lambda_2<\lambda_3\}$,
con $\Delta(\lambda)=\prod_{i<j}(\lambda_j-\lambda_i)$. La aplicación lineal
$$
(m,s_1,s_2)\mapsto\Big(m-\tfrac{2s_1+s_2}3,\ m+\tfrac{s_1-s_2}3,\ m+\tfrac{s_1+2s_2}3\Big)
$$
tiene determinante $1$ y lleva $\mathbb R\times(0,\infty)^2$ sobre $C$; su inversa es
$m=\frac13\sum_i\lambda_i$, $s_1=\lambda_2-\lambda_1$, $s_2=\lambda_3-\lambda_2$. En estas variables,
$\Delta=s_1s_2(s_1+s_2)$ y $|\lambda|^2=3m^2+\frac23(s_1^2+s_1s_2+s_2^2)$ (desarrollando). Por el
lema A.3(c) y el teorema de Fubini (integrando en $m$), $(s_1,s_2)$ tiene la densidad
$c\,\big(s_1s_2(s_1+s_2)\big)^\beta e^{-(s_1^2+s_1s_2+s_2^2)/3}$ sobre $(0,\infty)^2$ para una
constante $c>0$. Sea $q(x)=1+x+x^2$. Para todo boreliano $A\subseteq(0,\infty)$, sustituyendo
$s_2=xs_1$ en la integral interior y usando el teorema de Tonelli,
$$
P(x\in A)=c\int_0^\infty\!\!\int_A s_1^{3\beta+1}(x+x^2)^\beta e^{-s_1^2q(x)/3}\,dx\,ds_1
=c\,\gamma_\beta\int_A(x+x^2)^\beta\Big(\frac3{q(x)}\Big)^{1+3\beta/2}dx ,
$$
donde la sustitución $s_1=(3/q(x))^{1/2}u$ da el factor, con
$\gamma_\beta=\int_0^\infty u^{3\beta+1}e^{-u^2}du\in(0,\infty)$. Así, $x$ tiene la densidad
$f_\beta/\int_0^\infty f_\beta$.

*Reducción a $[0,1]$.* $f_\beta(1/x)=x^2f_\beta(x)$ (multiplicando numerador y denominador por
potencias de $x$). Con $x=1/y$,
$\int_{1/A}f_\beta(x)\,dx=\int_Af_\beta(1/y)y^{-2}dy=\int_Af_\beta$ para todo boreliano
$A\subseteq(0,1)$, donde $1/A=\{1/y:y\in A\}$; en particular, $\int_0^\infty f_\beta=2I_0$ con
$I_j=\int_0^1x^jf_\beta(x)\,dx$. Como $r=\min(x,1/x)$ y $P(x=1)=0$,
$\{r\in A\}=\{x\in A\}\cup\{x\in1/A\}$ (unión disjunta), luego
$$
P(r\in A)=\frac{2\int_Af_\beta}{2I_0}=\int_A\frac{f_\beta}{I_0},\qquad \mathbb Er=\frac{I_1}{I_0} .
$$
Como $x+x^2=q-1$, $f_\beta=g_\beta(q)$ con $g_1(u)=u^{-3/2}-u^{-5/2}$ y $g_2(u)=(u-1)^2u^{-4}$.
Escribamos $x=\frac12(q'(x)-1)$ con $q'=2x+1$ y sustituyamos $u=q(x)$ (de $q(0)=1$ a $q(1)=3$):
$$
I_1=\frac12\int_1^3g_\beta(u)\,du-\frac12I_0,\qquad \mathbb Er=\frac{\int_1^3g_\beta}{2I_0}-\frac12 .
$$
Aquí $\int_1^3g_1=\big[-2u^{-1/2}+\frac23u^{-3/2}\big]_1^3=\frac43-\frac{16}{9\sqrt3}$ y
$\int_1^3g_2=\big[-u^{-1}+u^{-2}-\frac13u^{-3}\big]_1^3=\frac8{81}$.

*Las integrales $I_0$.* Sea $L_a=\int_0^1q^{-a}dx$. De $(2x+1)^2=4q-3$ se obtiene
$\frac d{dx}\big[(2x+1)q^{-a}\big]=2q^{-a}-a(2x+1)^2q^{-a-1}=(2-4a)q^{-a}+3aq^{-a-1}$ e, integrando
sobre $[0,1]$,
$$
3a\,L_{a+1}=3^{1-a}-1+(4a-2)L_a .
$$
GOE: $f_1=q^{-3/2}-q^{-5/2}$. Con $a=\frac12$, $L_{3/2}=\frac23(\sqrt3-1)$; con $a=\frac32$,
$L_{5/2}=\frac29\big(\frac1{\sqrt3}-1+4L_{3/2}\big)=\frac{18\sqrt3-22}{27}$. Así,
$I_0=L_{3/2}-L_{5/2}=\frac4{27}$ y
$\mathbb Er=\frac{27}8\big(\frac43-\frac{16}{9\sqrt3}\big)-\frac12=\frac92-2\sqrt3-\frac12=4-2\sqrt3$.
GUE: $f_2=q^{-2}-2q^{-3}+q^{-4}$. Como $q=(x+\frac12)^2+\frac34$,
$L_1=\frac2{\sqrt3}\big[\arctan\frac{2x+1}{\sqrt3}\big]_0^1=\frac2{\sqrt3}\big(\frac\pi3-\frac\pi6\big)=\frac\pi{3\sqrt3}$, y la
recurrencia con $a=1,2,3$ da $L_2=\frac23L_1=\frac{2\pi}{9\sqrt3}$, $L_3=L_2-\frac19$, $L_4=\frac19\big(10L_3-\frac89\big)$. Así,
$I_0=L_2-2L_3+L_4=\frac{2\pi}{81\sqrt3}$ y $\mathbb Er=\frac{81\sqrt3}{4\pi}\cdot\frac8{81}-\frac12=\frac{2\sqrt3}\pi-\frac12$. $\square$

Los valores de (b) fueron hallados por Y. Y. Atas, E. Bogomolny, O. Giraud y G. Roux (Phys. Rev.
Lett. 110 (2013), 084101). Para $n$ mayor, el cociente medio no se calcula en este curso; el código
siguiente lo *observa* numéricamente para el GOE con $n=64$.

La parte (a) vale exactamente para $n$ niveles uniformes independientes en un intervalo: sus $n+1$
espaciamientos (incluidos los dos huecos de los extremos) tienen la ley de
$(e_0,\dots,e_n)/\sum_je_j$ para exponenciales independientes $e_j$, de modo que cada cociente
$s_i/s_{i+1}$ de espaciamientos consecutivos tiene la ley de $e_i/e_{i+1}$, y cada cociente de
espaciamientos $r_i$ tiene la ley de $r$ en (a) (teorema A.25 y corolario A.26(b),(d)).

La función `SpinChain.level_spacing_ratio` (en `chains`) devuelve el cociente medio de cada fila de
niveles. Comprobamos las proposiciones 5.13 y 5.14 y registramos el valor del GOE para $n=64$:

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
mc(ratio(levels), 2 * math.log(2) - 1, "<r>, independent uniform levels")      # Proposition 5.14(a)

n3 = 100_000                                                         # Proposition 5.14(b): 3 x 3
G = torch.randn(n3, 1, 3, 3, dtype=torch.float64)
Gc = torch.randn(n3, 1, 3, 3, dtype=dt)
goe3 = TypeIAlgebra([3], [3], complex_valued=False, precision='double').operator((G + G.mT) / 2, is_self_adjoint=True)
gue3 = TypeIAlgebra([3], [3], precision='double').operator((Gc + Gc.mH) / math.sqrt(2), is_self_adjoint=True)
h = 1e-4
xs = (torch.arange(int(1 / h), dtype=torch.float64) + 0.5) * h      # midpoint rule on [0, 1]
for beta, H, exact in ((1, goe3, 4 - 2 * math.sqrt(3)), (2, gue3, 2 * math.sqrt(3) / math.pi - 0.5)):
    f = (xs + xs ** 2) ** beta / (1 + xs + xs ** 2) ** (1 + 1.5 * beta)
    assert abs((xs * f).sum() / f.sum() - exact) < 1e-7              # E r = I_1 / I_0
    r3 = ratio(H.eigenvalues()[0])                                   # one ratio per matrix
    mc(r3, exact, f"beta={beta}: <r>, 3 x 3")
    mc((r3 <= 0.5).double(), (f * (xs <= 0.5)).sum().item() / f.sum().item(), f"beta={beta}: P(r <= 1/2), 3 x 3")

n = 64                                                               # numerical observation, no theorem
alg64 = TypeIAlgebra([n], [n], complex_valued=False, precision='double')
G = torch.randn(800, 1, n, n, dtype=torch.float64)
lam = alg64.operator((G + G.mT) / 2, is_self_adjoint=True).eigenvalues()[0]
r_goe = ratio(lam[:, n // 4: 3 * n // 4])                            # 800 samples of the mean ratio
print(f"<r>, GOE(64), middle half of the spectrum: {r_goe.mean():.4f} +- {r_goe.std() / math.sqrt(800):.4f}")

def same_mean(a, b, name):
    """Two independent Monte Carlo samples of the same law: means agree within 5 combined standard errors."""
    ea, eb = a.std().item() / math.sqrt(a.numel()), b.std().item() / math.sqrt(b.numel())
    print(f"{name}: {a.mean():.4f} vs {b.mean():.4f} +- {math.hypot(ea, eb):.4f}")
    assert abs(a.mean().item() - b.mean().item()) <= 5 * math.hypot(ea, eb)
```

*Observación numérica.* El código encuentra $\langle r\rangle\approx0.53$ para la mitad central del
espectro del GOE con $n=64$, por debajo del valor exacto $3\times3$, $4-2\sqrt3\approx0.536$, y muy
por encima del valor $2\ln2-1\approx0.386$ de los niveles independientes. Este número es un
resultado del cálculo, no un teorema; más abajo sirve de referencia para los bloques GOE del mismo
tamaño $64$, cuyos cocientes de espaciamientos tienen exactamente la misma ley.

## 5.4 Por qué no deben mezclarse los sectores

Sea $H=\bigoplus_cH_c\in\mathcal M$ un elemento autoadjunto aleatorio con **bloques
independientes**, por ejemplo bloques GOE. Como operador sobre $\mathbb C^D$, su espectro es la
unión de los espectros de los bloques. Dentro de un bloque los niveles se repelen; mostramos que
los niveles de bloques distintos no lo hacen.

Supongamos que cada bloque tiene una **densidad de niveles** $\rho_c$ continua y acotada:
$\mathbb E\sum_i\varphi(\lambda_i(H_c))=\int\varphi\,\rho_c$ para toda $\varphi$ de Borel acotada
(de modo que $\int\rho_c=k_c$). Para bloques GOE y GUE esto se cumple por el corolario C.8, que
deduce $\rho_c$ de la densidad conjunta del teorema C.7.

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

*Demostración.* Los valores propios ordenados son funciones continuas de la matriz
(corolario A.11), así que $\Lambda_c=(\lambda_i(H_c))_i$ y $\Lambda_d$ son vectores aleatorios
independientes. Por el teorema de Fubini (su ley conjunta es el producto de sus leyes),
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

Así pues, para dos bloques GOE o GUE, $P(\delta<\varepsilon)$ es al menos del orden de
$\varepsilon$: sus densidades de niveles son positivas en todas partes (en el corolario C.8,
$q(x,\cdot)>0$ fuera del conjunto nulo donde coinciden dos coordenadas de $(x,y)$), luego
$\int\rho_c\rho_d>0$. Dentro de un bloque GOE $2\times2$ vale
$1-e^{-\varepsilon^2/4}\le\varepsilon^2/4$ (proposición 5.13, espaciamiento sin normalizar), y
dentro de un bloque GOE $n\times n$ es a lo sumo una constante por $\varepsilon^2$
(proposición C.9). La prueba de repulsión de niveles falla para el espectro fusionado aunque cada
bloque sea un GOE perfecto.

Con más sectores el efecto es más fuerte. No demostramos ningún teorema límite para muchos
sectores; el experimento siguiente es una observación numérica. En la biblioteca, los bloques están
separados desde el principio: `eigenvalues()` devuelve un espectro por sector, y mezclarlos requiere
un `torch.cat` explícito. Ponemos $C$ bloques GOE independientes de tamaño $64$ en una misma álgebra
y comparamos el cociente medio por sector con el del espectro fusionado. Cada sector es una matriz
GOE de tamaño $64$, así que sus cocientes tienen exactamente la ley de `r_goe` de más arriba, y el
código lo comprueba; el comportamiento del espectro fusionado solo se observa:

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
    same_mean(per_sector.flatten(), r_goe, f"C = {C:2d}: sectors vs GOE(64)")   # every sector is GOE(64)
assert merged_r[1] > merged_r[2] > merged_r[4] > merged_r[16]       # numerical observation
assert abs(merged_r[16] - (2 * math.log(2) - 1)) < 0.01              # numerical observation
```

*Observación numérica.* Cada sector tiene $\langle r\rangle\approx0.53$, el valor GOE(64) de la
sección 5.3 (esta coincidencia es la comprobación de un hecho: los sectores son matrices GOE(64)).
El espectro fusionado baja a $0.42$ con dos sectores y a aproximadamente $0.39$, a menos de $0.01$
del valor $2\ln2-1$ de los niveles independientes (proposición 5.14(a)), con $8$ y $16$ sectores.
Que el espectro fusionado se acerque a niveles independientes cuando $C$ crece es lo que observa el
código; no se demuestra en este curso.

Lo mismo ocurre en un modelo físico. La cadena de Heisenberg con campos aleatorios (capítulo 9)
conserva la magnetización, así que su hamiltoniano está en $\bigoplus_{N=0}^LM_{\binom LN}$
(proposición 9.1). Sus bloques *no* son independientes —todos dependen de los mismos campos
aleatorios—, así que el teorema 5.15 no se aplica tal como está enunciado; el experimento siguiente
muestra el mismo efecto, como observación numérica.

```python
chain = SpinChain(10, complex_valued=False)                         # sectors N = 0, ..., 10
lam = chain.random_field_heisenberg(W=0.5, batch_size=20).eigenvalues()
k5 = chain.sector_dim(5)                                            # 252
r_sector = ratio(lam[5][:, k5 // 4: 3 * k5 // 4]).mean().item()
spectrum = torch.cat(lam, dim=-1).sort(dim=-1).values               # all 1024 levels mixed
r_all = ratio(spectrum[:, 256:768]).mean().item()
print(f"Heisenberg chain, L = 10: <r> in sector N = 5: {r_sector:.4f}, all sectors mixed: {r_all:.4f}")
assert r_sector > 0.49 and r_all < 0.42                             # numerical observation
```

*Observación numérica.* Dentro del sector $N=5$ el código encuentra $\langle r\rangle\approx0.53$,
cerca del valor GOE(64) de más arriba, y el espectro con todos los sectores mezclados da $0.39$,
cerca de $2\ln2-1$.

> **Significado físico** (interpretación física, no un teorema de este curso). Los niveles con
> números cuánticos conservados distintos no interactúan, así que se cruzan libremente; solo se
> repelen los niveles de la misma simetría. Un espectro que mezcla sectores de simetría muestra,
> por tanto, una estadística de tipo Poisson incluso cuando la dinámica en cada sector es
> completamente caótica, y una simetría no resuelta puede confundirse con integrabilidad o
> localización. La regla del oficio —resolver todas las simetrías antes de calcular la estadística
> de niveles— equivale a decir que hay que calcular en el álgebra $\bigoplus_cM_{k_c}$, no en $M_D$.

## 5.5 Con qué comparar

Unitarias de Haar: `scipy.stats.unitary_group` y `ortho_group`, las unitarias aleatorias de Haar
de QuTiP, y cualquier QR de `numpy`/`torch` seguida de la corrección de fases. Cálculo de
Weingarten: los paquetes simbólicos IntU (Mathematica; Z. Puchała, J. A. Miszczak), RTNI
(Mathematica y Python; M. Fukuda, R. König, I. Nechita, J. Phys. A 52 (2019), 425303) y el paquete
de Python Haarpy calculan momentos de cualquier orden de forma simbólica; el teorema 5.10 es su
caso de orden dos, y el teorema 5.12 (demostrado en el apéndice B), la fórmula general que
implementan. La estadística de niveles suele calcularse con unas pocas líneas de `numpy`, y los
códigos de diagonalización exacta como QuSpin proporcionan los espectros sector a sector. Lo que
añade la biblioteca es el grupo unitario de un álgebra como un único objeto: bloques de Haar
independientes de todos los tamaños generados en llamadas por lotes, el twirl identificado con la
esperanza condicional `channels.center_expectation`, y espectros devueltos por sector, de modo que
mezclar sectores nunca es la opción por defecto.

## Ejercicios

1. Sea $U$ de Haar en $U(\mathcal M)$ y $\operatorname{Tr}U=\sum_c\operatorname{Tr}U_c$. Demuestre
   que $\mathbb E|\operatorname{Tr}U|^2=C$ y que $\mathbb E|\operatorname{Tr}U|^4=\sum_cm(k_c)+2C(C-1)$
   con $m(1)=1$, $m(k)=2$ para $k\ge2$; compárelo con $1$ y $2$ para una unitaria de Haar de
   $U(D)$, $D\ge2$. Compruébelo para $\mathcal M=M_1\oplus M_2\oplus M_3$.
2. Usando el teorema 5.12, demuestre que $\mathbb E|U_{11}|^{2n}=\frac{n!}{k(k+1)\cdots(k+n-1)}$
   para $n\le k$ (el corolario B.14 lo extiende a $n>k$). Calcule numéricamente la matriz de
   Weingarten de orden $3$ para $k=3$ y compruebe que $\mathbb E|U_{11}|^6=\frac1{10}$.
3. El **conjunto ortogonal circular** (`measure='coe'`) es $S=W^TW$ con $W$ de Haar en $U(k)$.
   Demuestre que $S$ es unitaria y simétrica, y que su ley es invariante bajo $S\mapsto V^TSV$ para
   toda $V\in U(k)$. Sean $\lambda_1,\dots,\lambda_k$ los valores propios de $S$ (con multiplicidad)
   y $N_B(S)=\#\{j:\lambda_j\in B\}$ para los borelianos $B$ de la circunferencia unidad
   $\mathbb T$. Demuestre que $N_B(S)$ es una variable aleatoria, que $\nu(B)=\frac1k\mathbb EN_B(S)$
   es la longitud de arco normalizada (la densidad media de las fases propias es uniforme) y que
   $\mathbb E\operatorname{Tr}S^m=0$ para todo entero $m\neq0$. Calcule $\langle r\rangle$ para las
   fases propias del COE con $k=32$ y compárelo con el valor del GOE de la sección 5.3 (una
   observación numérica).

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
   $V=e^{i\alpha/2}1$: $e^{i\alpha}S$ tiene la ley de $S$ para todo $\alpha$ real.

   *Medibilidad.* Por el teorema 1.3, $S=Y\operatorname{diag}(\lambda_1,\dots,\lambda_k)Y^*$ con $Y$
   unitaria. Si $v$ es un vector unitario con $Sv=\mu v$, escribir $v=\sum_jc_jYe_j$ da
   $c_j(\lambda_j-\mu)=0$, luego $v$ está en $E_K=\operatorname{span}\{Ye_j:\lambda_j\in K\}$ siempre
   que $\mu\in K$; y $\dim E_K=N_K(S)$. Por tanto, para $K\subseteq\mathbb T$ y $m\ge1$:
   $N_K(S)\ge m$ si y solo si $S$ tiene $m$ vectores propios ortonormales con valores propios en $K$.
   Si $K$ es cerrado, el conjunto $F_{K,m}$ de las unitarias con esta propiedad es cerrado en $U(k)$:
   si $S_l\to S$ y $S_l$ tiene vectores propios ortonormales $v_i^{(l)}$ con valores propios
   $\mu_i^{(l)}\in K$ ($i\le m$), una subsucesión cumple $v_i^{(l)}\to v_i$ y
   $\mu_i^{(l)}\to\mu_i\in K$ para todo $i\le m$ (Bolzano–Weierstrass en el compacto
   $(\text{esfera unidad}\times K)^m$), y los límites son ortonormales con $Sv_i=\mu_iv_i$.
   Sean $0\le\theta_1(S)\le\dots\le\theta_k(S)<2\pi$ las fases propias, de modo que los valores
   propios son los $e^{i\theta_j(S)}$ (corolario 1.4: $|\lambda_j|=1$). Para $0\le t<2\pi$, el arco
   $K_t=\{e^{i\varphi}:0\le\varphi\le t\}$ es cerrado y $\{\theta_j\le t\}=\{N_{K_t}\ge j\}=F_{K_t,j}$,
   así que cada $\theta_j$ es una función de Borel sobre $U(k)$. Por tanto,
   $N_B(S)=\sum_j1_B(e^{i\theta_j(S)})$ es una variable aleatoria, y
   $\nu(B)=\frac1k\sum_jP(e^{i\theta_j(S)}\in B)$ es una medida de probabilidad de Borel sobre
   $\mathbb T$, el promedio de las leyes de los $k$ valores propios.

   *Uniformidad.* Los valores propios de $e^{i\alpha}S$ son los $e^{i\alpha}\lambda_j$, luego
   $N_B(e^{i\alpha}S)=N_{e^{-i\alpha}B}(S)$ y, como $e^{i\alpha}S$ tiene la ley de $S$,
   $\nu(e^{-i\alpha}B)=\frac1k\mathbb EN_B(e^{i\alpha}S)=\nu(B)$. Así, $\nu$ es invariante por
   rotaciones y, por la proposición A.22 (con $N=1$), es la longitud de arco normalizada:
   $\mathbb EN_B(S)=k\,|B|/2\pi$ para todo arco de longitud $|B|$. Del mismo modo,
   $\operatorname{Tr}S^m=\sum_j\lambda_j^m$ (también para $m<0$, pues $S^{-1}=S^*$) es una función
   continua y acotada de $S$, y
   $\mathbb E\operatorname{Tr}S^m=\mathbb E\operatorname{Tr}(e^{i\alpha}S)^m=e^{im\alpha}\,\mathbb E\operatorname{Tr}S^m$
   para todo $\alpha$, así que se anula para $m\neq0$: son $k$ veces los coeficientes de Fourier de
   $\nu$. En particular, la densidad media de las fases propias es constante, y los espaciamientos
   alrededor de la circunferencia no necesitan desdoblamiento.

   *Observación numérica.* El código encuentra $\langle r\rangle\approx0.53$ para el COE con
   $k=32$, a menos de $0.015$ del valor GOE(64) de la sección 5.3. Esta coincidencia es un resultado
   del cálculo; ningún teorema de este curso relaciona la estadística de espaciamientos del COE con
   la del GOE (sobre los conjuntos circulares, véase F. J. Dyson, J. Math. Phys. 3 (1962), 140–156).

```python
alg_s = TypeIAlgebra([32], [32], precision='double')
torch.manual_seed(2)
S = alg_s.random_unitary_operator(400, measure='coe').matrix[:, 0]
torch.manual_seed(2)
W = alg_s.random_unitary(32, batch_size=400)                         # same random numbers
assert close(S, W.mT @ W)                                            # the library's COE is W^T W
assert close(S, S.mT) and close(S @ S.mH, torch.eye(32, dtype=dt))
for m in (1, 2, 3):                                                  # E Tr S^m = 0
    trm = torch.linalg.matrix_power(S, m).diagonal(dim1=-2, dim2=-1).sum(-1)
    mc(trm.real, 0.0, f"E Re Tr S^{m}")
    mc(trm.imag, 0.0, f"E Im Tr S^{m}")
phases = torch.linalg.eigvals(S).angle().sort(dim=-1).values         # in (-pi, pi]
mc(((phases > 0) & (phases <= 1)).sum(-1).double(), 32 / (2 * math.pi), "E #eigenphases in (0, 1]")
s = torch.diff(torch.cat([phases, phases[:, :1] + 2 * math.pi], dim=-1), dim=-1)   # around the circle
r = torch.minimum(s, s.roll(-1, -1)) / torch.maximum(s, s.roll(-1, -1))
r_coe = r.mean(-1)
print(f"<r>, COE(32) eigenphases: {r_coe.mean():.4f}; GOE(64), Section 5.3: {r_goe.mean():.4f}")
assert abs(r_coe.mean() - r_goe.mean()) < 0.015                      # numerical observation
```

</details>

**Siguiente:** el [capítulo 6](06_quantum_channels.md) trata los canales cuánticos: las
representaciones de Kraus y de Stinespring, los duales respecto de una traza y las aplicaciones
entre sectores. La esperanza condicional $E$, que aquí aparece como un promedio de Haar, es el
primer ejemplo.
