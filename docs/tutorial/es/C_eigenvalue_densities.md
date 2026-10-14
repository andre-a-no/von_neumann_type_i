[English](../C_eigenvalue_densities.md) · **Español** · [Français](../fr/C_eigenvalue_densities.md) · [Deutsch](../de/C_eigenvalue_densities.md) · [中文](../zh/C_eigenvalue_densities.md) · [日本語](../ja/C_eigenvalue_densities.md) · [Русский](../ru/C_eigenvalue_densities.md)

# Apéndice C. Densidades de valores propios y la fórmula de Page

El capítulo 5 utiliza la densidad conjunta de los valores propios de los conjuntos gaussianos, y el capítulo 7
utiliza la fórmula de Page para la entropía de entrelazamiento media de un estado puro aleatorio. Este apéndice
demuestra ambas. El camino es: la fórmula de integración de Weyl (el jacobiano de la «diagonalización»), los
conjuntos gaussianos, la función de un punto de una densidad determinantal mediante polinomios ortogonales, el
conjunto de Wishart–Laguerre y la medida inducida sobre los estados reducidos y, por último, un cálculo exacto con
polinomios de Laguerre.

**Aprenderá** por qué los valores propios de una matriz aleatoria invariante se repelen como
$\prod_{i<j}|\lambda_i-\lambda_j|^\beta$, por qué los vectores propios están entonces distribuidos según la medida de
Haar y son independientes de los valores propios, cómo se integra un determinante de polinomios ortonormales, y cómo
$\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}$ se deduce de un puñado de integrales.

Además de los requisitos previos del curso (teoría de la medida, la fórmula del cambio de variables en
$\mathbb R^N$ para aplicaciones $C^1$ inyectivas con derivada invertible, los teoremas de la función inversa y de la
función implícita), utilizamos del capítulo 5 el lema 5.1 (unicidad e invariancia por la derecha de la medida de
Haar), el lema 5.2 y el teorema 5.3 (existencia de la medida de Haar en $U(n)$ mediante la QR gaussiana), del apéndice B
el teorema B.17 (lo mismo para $O(n)$), del apéndice A la proposición A.16(a) (invariancia unitaria de los vectores
gaussianos complejos) y el teorema A.20(b) (la ley de $|z|^2$), y del capítulo 7 la
representación $\rho_A=GG^*/\operatorname{Tr}GG^*$ de la demostración del teorema 7.11.

**Notación.** $\beta\in\{1,2\}$ y, en consecuencia, $\mathbb F=\mathbb R$ o $\mathbb C$. $\mathcal H$ es el
espacio de las matrices $n\times n$ reales simétricas ($\beta=1$) o complejas hermíticas ($\beta=2$), con la
medida de Lebesgue $dH$ en las coordenadas $h_{ii}$, $\operatorname{Re}h_{ij}$ y (para $\beta=2$)
$\operatorname{Im}h_{ij}$, $i<j$; su dimensión es $n+\beta n(n-1)/2$. $K=O(n)$ o $U(n)$, $T$ es el
grupo de los elementos diagonales de $K$ (signos, o fases), $\mathfrak k=\{X\in M_n(\mathbb F):X^*=-X\}$,
$\mathfrak t$ su parte diagonal y $\mathfrak m$ su parte fuera de la diagonal (diagonal nula), con
coordenadas $\operatorname{Re}X_{ij}$, $\operatorname{Im}X_{ij}$, $i<j$. Para $\lambda\in\mathbb R^n$,
$\Lambda=\operatorname{diag}(\lambda)$ y $\Delta(\lambda)=\prod_{i<j}(\lambda_j-\lambda_i)$; $C=\{\lambda:\lambda_1<\dots<\lambda_n\}$,
sobre el cual $\Delta>0$. Para $H\in\mathcal H$, $\lambda(H)$ es el vector de sus valores propios en orden
creciente. $\|X\|_2^2=\operatorname{Tr}X^*X$. Para $V\in K$ la aplicación $H\mapsto VHV^*$ es lineal sobre $\mathcal H$
y conserva $\|\cdot\|_2$; una isometría lineal tiene $|\det|=1$, así que **conserva $dH$**.

## C.1 Espectros simples y marcos propios

**Lema C.1.** (a) El conjunto $\Omega\subseteq\mathcal H$ de las matrices con $n$ valores propios distintos es
abierto, y $\mathcal H\setminus\Omega$ es un conjunto de medida de Lebesgue nula. (b) Sobre $\Omega$ la aplicación $H\mapsto\lambda(H)$ es
$C^\infty$, y las proyecciones espectrales $P_j(H)=\prod_{k\neq j}\frac{H-\lambda_k(H)}{\lambda_j(H)-\lambda_k(H)}$ son continuas.

*Demostración.* (a) Sea $A_{ij}=\lambda_j^{i-1}$ con $\lambda=\lambda(H)$. Entonces $\det A=\Delta(\lambda)$ (Vandermonde) y
$(AA^T)_{ij}=\sum_k\lambda_k^{i+j-2}=\operatorname{Tr}H^{i+j-2}$, así que
$\operatorname{disc}(H)=\det\big(\operatorname{Tr}H^{i+j-2}\big)_{i,j=1}^n=\Delta(\lambda(H))^2$ es un polinomio
en las coordenadas de $H$, no nulo en $\operatorname{diag}(1,\dots,n)$, y $\Omega=\{\operatorname{disc}\neq0\}$ es abierto.
El conjunto de ceros de un polinomio no nulo $p$ sobre $\mathbb R^N$ es nulo, por inducción sobre $N$: para $N=1$ es
finito; en otro caso se escribe $p=\sum_jp_j(x')x_N^j$ con algún $p_{j_0}\neq0$; para $x'$ fuera del conjunto de ceros de
$p_{j_0}$ (nulo por inducción) $x_N\mapsto p(x',x_N)$ tiene un número finito de ceros, y el teorema de Tonelli aplicado
a la indicatriz del conjunto de ceros da medida $0$.
(b) $\chi(H,t)=\det(t-H)=\prod_i(t-\lambda_i)$ es un polinomio real en $(H,t)$. En $H_0\in\Omega$,
$\partial_t\chi(H_0,\lambda_j(H_0))=\prod_{k\neq j}(\lambda_j-\lambda_k)\neq0$, así que el teorema de la función implícita da
funciones $C^\infty$ $r_j$ cerca de $H_0$ con $\chi(H,r_j(H))=0$ y $r_j(H_0)=\lambda_j(H_0)$. Por continuidad
$r_1<\dots<r_n$ cerca de $H_0$: son $n$ raíces distintas de un polinomio de grado $n$, luego todas ellas, en
orden, y $r_j=\lambda_j$. Si $H=\sum_k\lambda_kQ_k$ es la descomposición espectral, la interpolación de Lagrange da
$P_j(H)=Q_j$; la continuidad se sigue de la de $\lambda$. $\square$

Para $H\in\Omega$, un **marco propio** es una $V\in K$ con $H=V\Lambda V^*$, $\Lambda=\operatorname{diag}\lambda(H)$. Si $V,V'$
son marcos propios, $M=V'^*V$ conmuta con $\Lambda$, así que $(\lambda_i-\lambda_j)M_{ij}=0$ y $M\in T$; recíprocamente, $VD$
es un marco propio para toda $D\in T$. Un **marco propio de Borel** $F(H)$ se obtiene tomando como $j$-ésima columna
$P_j(H)e_m/\|P_j(H)e_m\|$ con $m$ el menor índice tal que $P_j(H)e_m\neq0$: es continuo sobre cada uno de los
borelianos en los que estos índices son fijos, sus columnas son vectores propios ortonormales, y es real cuando $H$ lo es.

**Lema C.2 (medidas de Haar).** $U(n)$, $O(n)$ y $T$ tienen medidas de Haar, y $\int_Kf(VW)\,dV=\int_Kf\,dV$
para todo $W\in K$ y toda $f$ de Borel acotada.

*Demostración.* $U(n)$: teorema 5.3. $O(n)$: teorema B.17 (la construcción QR gaussiana real: una matriz gaussiana real es
invertible casi seguramente, y el factor ortogonal de su descomposición QR está distribuido según la medida de Haar). $T$: la medida uniforme sobre $\{\pm1\}^n$, o el producto de las longitudes de arco normalizadas sobre las
circunferencias, es invariante por la multiplicación por un elemento fijo. La invariancia por la derecha es el lema 5.1. $\square$

**Lema C.3 (marcos propios de conjuntos invariantes).** Sea $H$ un elemento aleatorio de $\mathcal H$ con
$P(H\in\Omega)=1$ cuya ley es invariante por $H\mapsto VHV^*$ para toda $V\in K$. Sean $D$ uniforme sobre $T$ e
independiente de $H$, y $\tilde U=F(H)D$. Entonces $\tilde U$ es un marco propio de $H$, está distribuido según la medida de Haar en
$K$, y es independiente de $\lambda(H)$.

*Demostración.* Fijemos $V\in K$. Tanto $F(VHV^*)$ como $VF(H)$ son marcos propios de $VHV^*$, así que
$F(VHV^*)=VF(H)E(H)$ con $E(H)=(VF(H))^*F(VHV^*)\in T$ función de Borel de $H$. Para funciones de Borel acotadas
$\varphi$ sobre $\mathbb R^n$ y $\psi$ sobre $K$, como $VHV^*$ tiene la ley de $H$ y $\lambda(VHV^*)=\lambda(H)$,
$$
\mathbb E\big[\varphi(\lambda(H))\psi(\tilde U)\big]=\mathbb E\big[\varphi(\lambda(H))\psi(VF(H)E(H)D)\big]
=\mathbb E\big[\varphi(\lambda(H))\psi(V\tilde U)\big],
$$
donde el segundo paso integra primero sobre $D$ (Fubini): para $H$ fijo, $E(H)D$ es uniforme sobre $T$. Sea ahora $V$
distribuida según la medida de Haar e independiente de todo lo demás; por Fubini y el lema C.2,
$\mathbb E[\varphi(\lambda(H))\psi(\tilde U)]=\mathbb E\big[\varphi(\lambda(H))\int_K\psi(V\tilde U)\,dV\big]
=\mathbb E[\varphi(\lambda(H))]\int_K\psi\,dV$. Con $\varphi=1$, $\tilde U$ es de Haar; la forma de producto es la
independencia. $\square$

## C.2 La fórmula de integración de Weyl

Para $\varepsilon>0$ sean $\mathfrak m_\varepsilon=\{X\in\mathfrak m:\|X\|_2<\varepsilon\}$ y
$$
\Psi(X,\lambda)=e^X\Lambda e^{-X},\qquad (X,\lambda)\in\mathfrak m_\varepsilon\times C .
$$
Para $X\in\mathfrak k$, $(e^X)^*=e^{-X}=(e^X)^{-1}$, así que $e^X\in K$.

**Lema C.4 (una carta).** Existen $\varepsilon>0$ y una función continua $j:\overline{\mathfrak m_\varepsilon}\to(0,\infty)$
con $j(0)=1$ tales que: (a) $|\det D\Psi(X,\lambda)|=\Delta(\lambda)^\beta j(X)$; (b) $\Psi$ es inyectiva; (c) para
$H\in\Omega$, $H\in\Psi(\mathfrak m_\varepsilon\times C)$ si y solo si $F(H)\in S=\{e^XD:X\in\mathfrak m_\varepsilon,D\in T\}$, y
$S$ es un boreliano con $ST=S$.

*Demostración.* (a) $\Psi$ es $C^\infty$ (series de potencias). Para $Y\in\mathfrak m$ pongamos $A_X(Y)=e^{-X}\frac{d}{dt}\big|_0e^{X+tY}$.
Derivando $W(t)^*W(t)=1$ para $W(t)=e^{X+tY}\in K$ se obtiene $A_X(Y)\in\mathfrak k$, y con
$\frac{d}{dt}W^{-1}=-W^{-1}W'W^{-1}$, para $g=e^X$,
$$
D\Psi(X,\lambda)(Y,\mu)=g\big([A_X(Y),\Lambda]+\operatorname{diag}(\mu)\big)g^{-1}.
$$
Las matrices diagonales conmutan con $\Lambda$, así que $[A_X(Y),\Lambda]=[B_XY,\Lambda]$ con $B_X=\pi\circ A_X:\mathfrak m\to\mathfrak m$ y
$\pi$ la proyección de $\mathfrak k$ sobre $\mathfrak m$ paralela a $\mathfrak t$. Así $D\Psi=\operatorname{Ad}_g\circ L_\lambda\circ(B_X\times1)$
con $L_\lambda(Z,\mu)=[Z,\Lambda]+\operatorname{diag}(\mu)$. Como $([Z,\Lambda])_{ij}=Z_{ij}(\lambda_j-\lambda_i)$, en las coordenadas
$L_\lambda$ multiplica $\operatorname{Re}Z_{ij}$ (y $\operatorname{Im}Z_{ij}$) por $\lambda_j-\lambda_i$ y envía $\mu$ a la diagonal:
$|\det L_\lambda|=\Delta(\lambda)^\beta$. $\operatorname{Ad}_g$ conserva $dH$. Pongamos $j(X)=|\det B_X|$: es continua, y
$B_0=1$ porque $\frac d{dt}e^{tY}|_0=Y$; se elige $\varepsilon$ con $j>0$ sobre $\overline{\mathfrak m_\varepsilon}$.
(b) Si $\Psi(X,\lambda)=\Psi(X',\lambda')=H$, entonces $\lambda=\lambda'=\lambda(H)$, y $e^X,e^{X'}$ son marcos propios de $H\in\Omega$,
así que $e^X=e^{X'}D$ con $D\in T$. Sea $G(X,t)=e^Xe^t$ sobre $\mathfrak m\times\mathfrak t$ (para $\beta=1$, $\mathfrak t=\{0\}$); $DG(0,0)(Y,s)=Y+s$ es inyectiva,
y como $Y$ está fuera de la diagonal y $s$ es diagonal, $\|Y+s\|_2^2=\|Y\|_2^2+\|s\|_2^2$, así que $\|DG(0,0)v\|=\|v\|$ para la norma producto; pongamos $m=1$. Se elige $\delta>0$ con $\|DG(z)-DG(0)\|\le m/2$ para $\|z\|<\delta$; para $u,v$ en esta
bola, $G(u)-G(v)=\int_0^1DG(v+s(u-v))(u-v)\,ds$ da $\|G(u)-G(v)\|\ge\frac m2\|u-v\|$, así que $G$ es inyectiva en ella.
($\delta$ depende solo de $n,\beta$; $\varepsilon$ se elige después.)
Ahora bien, $\|D-1\|_2=\|e^{X'}(D-1)\|_2=\|e^X-e^{X'}\|_2\le\|e^X-1\|_2+\|e^{X'}-1\|_2\le2(e^\varepsilon-1)$, por la serie de potencias y
$\|AB\|_2\le\|A\|_2\|B\|_2$. Para $\beta=1$, $D=\operatorname{diag}(\pm1)$ y una entrada $-1$ daría
$\|D-1\|_2\ge2$; para $\varepsilon<\log2$, $D=1=e^t$ con $t=0$. Para $\beta=2$, $D=e^t$ con $t=i\operatorname{diag}(\theta)$,
$\theta_j\in(-\pi,\pi]$, y $|\theta_j|\le\frac\pi2|e^{i\theta_j}-1|$ (pues $|e^{i\theta}-1|=2|\sin\frac\theta2|$ y $\sin x\ge2x/\pi$ sobre $[0,\pi/2]$), así que
$\|t\|_2\le\frac\pi2\|D-1\|_2\le\pi(e^\varepsilon-1)$. Para $\varepsilon$ pequeño ($\varepsilon^2+\pi^2(e^\varepsilon-1)^2<\delta^2$), $(X,0)$ y $(X',t)$ están en la $\delta$-bola y
$G(X,0)=e^X=e^{X'}D=G(X',t)$, así que $X=X'$.
(c) Si $H=\Psi(X,\lambda')$, entonces $e^X$ es un marco propio, así que $F(H)\in e^XT\subseteq S$. Si $F(H)=e^XD$, entonces
$H=F(H)\Lambda F(H)^*=\Psi(X,\lambda(H))$. $S=G(\mathfrak m_\varepsilon\times\mathfrak t)$ (o $e^{\mathfrak m_\varepsilon}T$, con $T$ finito) es una
imagen continua de un conjunto $\sigma$-compacto, luego es $\sigma$-compacto y de Borel; $ST=S$ es evidente. $\square$

**Teorema C.5 (fórmula de integración de Weyl).** Existe una constante $c_{n,\beta}>0$ tal que para toda función de Borel
$f:\mathcal H\to[0,\infty]$,
$$
\int_{\mathcal H}f(H)\,dH=c_{n,\beta}\int_C\Delta(\lambda)^\beta\Big(\int_Kf(V\Lambda V^*)\,dV\Big)d\lambda .
$$
Equivalentemente, $\Phi(V,\lambda)=V\Lambda V^*$ transporta la medida $c_{n,\beta}\,dV\otimes\Delta(\lambda)^\beta d\lambda$ sobre $K\times C$
a la medida de Lebesgue sobre $\mathcal H$. El factor $\Delta(\lambda)^\beta$ es el jacobiano de la carta $\Psi$ en $X=0$
(lema C.4(a)); la constante $c_{n,\beta}$ refleja la normalización de la medida de Haar y nunca se necesita
explícitamente más abajo, pues toda densidad obtenida de la fórmula se normaliza por su masa total.

*Demostración.* *Conjunto de referencia.* Sea $H_0$ con densidad $p_0(H)=Z_0^{-1}e^{-\operatorname{Tr}H^2/2}$; como
$\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$, $Z_0$ es un producto finito de integrales gaussianas. $p_0$ y $dH$
son invariantes por $H\mapsto VHV^*$, luego también lo es la ley de $H_0$; $P(H_0\in\Omega)=1$ por el lema C.1. Por el lema C.3,
$\tilde U$ es de Haar, independiente de $\lambda=\lambda(H_0)$, y $H_0=\tilde U\Lambda\tilde U^*$.
*La ley de $\lambda$.* Sea $A\subseteq C$ un boreliano. Por el lema C.4(c) y $ST=S$,
$\{H_0\in\Psi(\mathfrak m_\varepsilon\times A)\}=\{\lambda\in A,\ \tilde U\in S\}$ salvo un suceso nulo, cuya probabilidad es
$P(\lambda\in A)\,h$ con $h=\int_K1_S\,dV$. Por otra parte, $\Psi$ es inyectiva y $C^1$ con jacobiano no nulo sobre el
abierto $\mathfrak m_\varepsilon\times C$, y $p_0(\Psi(X,\lambda))=p_0(\Lambda)$; la fórmula del cambio de variables y el lema C.4(a) dan
$$
P\big(H_0\in\Psi(\mathfrak m_\varepsilon\times A)\big)=\int_{\mathfrak m_\varepsilon}\!\int_Ap_0(\Lambda)\Delta(\lambda)^\beta j(X)\,d\lambda\,dX
=\kappa\int_Ap_0(\Lambda)\Delta(\lambda)^\beta\,d\lambda,\qquad \kappa=\int_{\mathfrak m_\varepsilon}j>0 .
$$
Con $A=C$, $h=\kappa\int_Cp_0\Delta^\beta>0$. Por tanto $\lambda$ tiene la densidad $c\,p_0(\Lambda)\Delta(\lambda)^\beta$ sobre $C$, $c=\kappa/h$.
*$f$ general.* Como $p_0>0$, $\int f\,dH=\mathbb E\big[(f/p_0)(\tilde U\Lambda\tilde U^*)\big]$. La ley conjunta de $(\tilde U,\lambda)$ es
$dV\otimes c\,p_0(\Lambda)\Delta^\beta d\lambda$, y $p_0(V\Lambda V^*)=p_0(\Lambda)$; el teorema de Tonelli da la fórmula con $c_{n,\beta}=c$.
La constante queda determinada por la propia fórmula (tómese $f=p_0$). $\square$

**Corolario C.6 (densidades invariantes).** Sea $H$ con una densidad $p$ sobre $\mathcal H$ tal que $p(VHV^*)=p(H)$ para toda
$V\in K$. Entonces $\lambda(H)$ tiene la densidad $c_{n,\beta}\,p(\Lambda)\Delta(\lambda)^\beta$ sobre $C$; el marco propio $\tilde U$ del lema C.3
está distribuido según la medida de Haar y es independiente de $\lambda(H)$; y para toda función de Borel simétrica $g\ge0$ sobre $\mathbb R^n$,
$$
\mathbb E\,g(\lambda(H))=\frac{c_{n,\beta}}{n!}\int_{\mathbb R^n}g(x)\,p(\operatorname{diag}x)\,|\Delta(x)|^\beta\,dx .
$$

*Demostración.* $P(H\in\Omega)=1$ por el lema C.1. Para un boreliano $A\subseteq C$ se aplica el teorema C.5 a $f=1_A(\lambda(\cdot))\,p$, para la cual
$f(V\Lambda V^*)=1_A(\lambda)p(\Lambda)$. La ley de $H$ es invariante, así que se aplica el lema C.3. Salvo el conjunto nulo en el que dos
coordenadas coinciden, $\mathbb R^n$ es la unión disjunta de las $n!$ imágenes de $C$ por las permutaciones de coordenadas, y
$g$, $|\Delta|$ y $x\mapsto p(\operatorname{diag}x)$ son invariantes por permutaciones (una matriz de permutación está en $K$). $\square$

La comprobación siguiente deriva $\Psi$ numéricamente para $n=2,3$ y ambos $\beta$: en $X=0$ el jacobiano es $\Delta(\lambda)^\beta$, y
en un $X\neq0$ fijo su cociente por $\Delta(\lambda)^\beta$ no depende de $\lambda$ (lema C.4(a)).

```python
import math
import torch
from fractions import Fraction
from torch_vn_algebra import TypeIAlgebra, states

torch.manual_seed(0)
dt, f64 = torch.complex128, torch.float64

def mc(samples, exact, name):
    """Monte Carlo mean of real samples with its standard error; must agree within 5 standard errors."""
    mean, err = samples.mean().item(), samples.std().item() / math.sqrt(samples.numel())
    print(f"{name}: {mean:.4f} +- {err:.4f}   (theory {exact:.4f})")
    assert abs(mean - exact) <= 5 * err + 1e-12
    return mean

def vandermonde(lam):
    """Delta(lam) = prod_{i<j} (lam_j - lam_i) along the last axis."""
    n = lam.shape[-1]
    return math.prod(lam[..., j] - lam[..., i] for i in range(n) for j in range(i + 1, n))

def coords(H, beta):
    """Coordinates of dH: h_ii, Re h_ij (i < j) and, for beta = 2, Im h_ij."""
    iu = torch.triu_indices(H.shape[-1], H.shape[-1], 1)
    off = H[..., iu[0], iu[1]]
    return torch.cat([H.diagonal(dim1=-2, dim2=-1).real, off.real] + ([off.imag] if beta == 2 else []), dim=-1)

def from_m(c, n, beta):
    """The element X of m (X* = -X, zero diagonal) with coordinates Re X_ij (, Im X_ij), i < j."""
    iu = torch.triu_indices(n, n, 1)
    p = iu.shape[1]
    X = torch.zeros(n, n, dtype=dt)
    X[iu[0], iu[1]] = c[:p] + 1j * c[p:] if beta == 2 else c[:p].to(dt)
    return X - X.mH

def Psi(p, n, beta):                                       # (X, lambda) -> e^X Lambda e^{-X}
    d = beta * n * (n - 1) // 2
    E = torch.linalg.matrix_exp(from_m(p[:d], n, beta))
    return coords(E @ torch.diag(p[d:]).to(dt) @ E.mH, beta)

def jacobian(p, n, beta, h=1e-6):
    J = torch.stack([(Psi(p + h * e, n, beta) - Psi(p - h * e, n, beta)) / (2 * h)
                     for e in torch.eye(len(p), dtype=f64)], dim=1)
    return torch.linalg.det(J).abs().item()

for n in (2, 3):
    for beta in (1, 2):
        d = beta * n * (n - 1) // 2
        X = 0.4 * torch.randn(d, dtype=f64)
        ratios = []
        for trial in range(3):
            lam = torch.randn(n, dtype=f64).sort().values
            D = vandermonde(lam).item() ** beta
            assert abs(jacobian(torch.cat([torch.zeros(d, dtype=f64), lam]), n, beta) / D - 1) < 1e-6
            ratios.append(jacobian(torch.cat([X, lam]), n, beta) / D)
        print(f"n = {n}, beta = {beta}: j(X) = {ratios[0]:.6f} for a random X")
        assert max(ratios) - min(ratios) < 1e-6              # j(X) does not depend on lambda
```

## C.3 Los conjuntos gaussianos

**Teorema C.7 (GOE y GUE).** Sea $H=(G+G^T)/2$ con $G$ real de entradas independientes $N(0,1)$ ($\beta=1$, GOE), o
$H=(G+G^*)/\sqrt2$ con $G$ de Ginibre ($\beta=2$, GUE), como en la sección 5.3. Entonces $H$ tiene la densidad
$(2\pi)^{-n/2}\pi^{-\beta n(n-1)/4}e^{-\operatorname{Tr}H^2/2}$ sobre $\mathcal H$, sus valores propios tienen la densidad
$$
Z_{n,\beta}^{-1}\,\Delta(\lambda)^\beta\,e^{-\sum_i\lambda_i^2/2}\quad\text{on }C,\qquad Z_{n,\beta}=\int_C\Delta^\beta e^{-|\lambda|^2/2}\,d\lambda\in(0,\infty),
$$
y el marco propio $\tilde U$ está distribuido según la medida de Haar en $K$ y es independiente de los valores propios.

*Demostración.* GOE: $h_{ii}=G_{ii}\sim N(0,1)$ y $h_{ij}=(G_{ij}+G_{ji})/2\sim N(0,\frac12)$ para $i<j$. GUE: $h_{ii}=\sqrt2\operatorname{Re}G_{ii}\sim N(0,1)$,
$\operatorname{Re}h_{ij}=(\operatorname{Re}G_{ij}+\operatorname{Re}G_{ji})/\sqrt2$ y $\operatorname{Im}h_{ij}=(\operatorname{Im}G_{ij}-\operatorname{Im}G_{ji})/\sqrt2$ son $N(0,\frac12)$.
Coordenadas distintas son funciones de conjuntos disjuntos de variables independientes, luego son independientes, y la
densidad es $\prod_i(2\pi)^{-1/2}e^{-h_{ii}^2/2}\prod_{i<j}\pi^{-\beta/2}e^{-|h_{ij}|^2}$, que es la enunciada porque
$\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$. Es invariante, y se aplica el corolario C.6; la constante queda fijada por la
normalización. $Z_{n,\beta}<\infty$ porque $|\lambda_j-\lambda_i|\le(1+|\lambda_i|)(1+|\lambda_j|)$, así que
$|\Delta(\lambda)|\le\prod_i(1+|\lambda_i|)^{n-1}$. $\square$

**Corolario C.8 (densidad de niveles).** Para el GOE o el GUE de tamaño $n$, sean $q(x)=(n!Z_{n,\beta})^{-1}|\Delta(x)|^\beta e^{-|x|^2/2}$ sobre
$\mathbb R^n$ y $\rho(x)=n\int_{\mathbb R^{n-1}}q(x,y)\,dy$. Entonces $\mathbb E\sum_i\varphi(\lambda_i(H))=\int\varphi\rho$ para toda función de Borel acotada
$\varphi$, $\int\rho=n$, y $\rho$ es continua con $0\le\rho(x)\le C_n(1+|x|)^{\beta(n-1)}e^{-x^2/2}$. En particular, la
hipótesis del teorema 5.15 se cumple para bloques GOE o GUE independientes.

*Demostración.* Por el corolario C.6 con $g(x)=\sum_i\varphi(x_i)$ (para $\varphi\ge0$; en general se aplica a $\varphi^\pm=\max(\pm\varphi,0)$) y la simetría de $q$,
$\mathbb E\sum_i\varphi(\lambda_i)=n\int\varphi(x_1)q(x)\,dx=\int\varphi\rho$ (Fubini). Por la cota de $|\Delta|$,
$q(x,y)\le(n!Z_{n,\beta})^{-1}(1+|x|)^{\beta(n-1)}e^{-x^2/2}\prod_{i}(1+|y_i|)^{\beta(n-1)}e^{-y_i^2/2}$, lo que da la cota tras
integrar sobre $y$ y, como $\sup_x(1+|x|)^{\beta(n-1)}e^{-x^2/2}<\infty$, una mayorante integrable independiente de $x$: la continuidad
se sigue por convergencia dominada. $\square$

**Proposición C.9 (espaciamiento mínimo).** Para el GOE o el GUE de tamaño $n\ge2$ sea $\delta=\min_i(\lambda_{i+1}-\lambda_i)$. Existen
constantes $0<c\le C$ (que dependen de $n,\beta$) con $c\,\varepsilon^{\beta+1}\le P(\delta<\varepsilon)\le C\varepsilon^{\beta+1}$ para $0<\varepsilon\le1$.

*Demostración.* *Cota superior.* $P(\delta<\varepsilon)\le\sum_iP(\lambda_{i+1}-\lambda_i<\varepsilon)$. Sobre $\{0<u=x_{i+1}-x_i<\varepsilon\}$ se acota el factor del
par $(i,i+1)$ por $u^\beta$, los demás como arriba, $1+|x_{i+1}|\le2(1+|x_i|)$ y $e^{-x_{i+1}^2/2}\le1$; sustituyendo
$x_{i+1}=x_i+u$,
$$
P(\lambda_{i+1}-\lambda_i<\varepsilon)\le\frac{2^{\beta(n-1)}}{Z_{n,\beta}}\int_0^\varepsilon u^\beta du\int(1+|x|)^{2\beta(n-1)}e^{-\frac{x^2}2}dx
\Big(\int(1+|x|)^{\beta(n-1)}e^{-\frac{x^2}2}dx\Big)^{n-2}.
$$
*Cota inferior.* Sea $R_\varepsilon=\{x:0<x_1<1,\ 0<x_2-x_1<\varepsilon,\ 3k<x_k<3k+1\ (3\le k\le n)\}\subseteq C$. Sobre $R_\varepsilon$,
$\delta\le x_2-x_1$, todo par distinto de $(1,2)$ está a distancia $\ge1$, y $|x|^2\le M$ con $M$ independiente de $\varepsilon\le1$; así que
$P(\delta<\varepsilon)\ge P(\lambda\in R_\varepsilon)\ge Z_{n,\beta}^{-1}e^{-M/2}\int_0^\varepsilon u^\beta du$. $\square$

Para el GOE este es el orden $\varepsilon^2$ citado al final de la sección 5.4. Comprobamos el teorema C.7 para $n=3$: la densidad
de niveles del GOE se calcula a partir de la densidad conjunta integrando sobre los otros dos valores propios en una malla, y se compara
con el número medio de valores propios en cada intervalo; comprobamos también que los vectores propios no están correlacionados con los
valores propios ($|\tilde U_{11}|^2$ no depende de las fases $D$, así que puede leerse en cualquier marco propio devuelto por
`eigh`). Para $U\in O(3)$ de Haar, $\mathbb E|U_{11}|^2=\frac13$: la multiplicación por la izquierda por una matriz de permutación permuta las
entradas de la primera columna y conserva la medida de Haar, así que las tres $\mathbb E|U_{i1}|^2$ son iguales, y suman $1$.

```python
B, n = 200_000, 3
G = torch.randn(B, n, n, dtype=f64)
goe_lam, goe_vec = torch.linalg.eigh((G + G.mT) / 2)                    # Section 5.3 normalisation
edges = torch.arange(-4.0, 4.01, 1.0, dtype=f64)

h = 0.04                                                                # grid for the other two eigenvalues
xs = torch.arange(-6, 6 + h / 2, h, dtype=f64)
Y1, Y2 = torch.meshgrid(xs, xs, indexing='ij')
joint = lambda x1, x2, x3: (vandermonde(torch.stack([x1, x2, x3], -1)).abs()
                            * torch.exp(-(x1 ** 2 + x2 ** 2 + x3 ** 2) / 2))
marg = torch.stack([joint(x + 0 * Y1, Y1, Y2).sum() * h * h for x in xs])   # int |Delta| e^{-|x|^2/2} dy dz
rho_goe = 3 * marg / (marg.sum() * h)                                   # Corollary C.8, normalised by n! Z
for lo, hi in zip(edges[:-1], edges[1:]):
    sel = (xs >= lo - 1e-9) & (xs <= hi + 1e-9)
    count = ((goe_lam >= lo) & (goe_lam < hi)).sum(-1).double()
    mc(count, torch.trapezoid(rho_goe[sel], xs[sel]).item(), f"GOE(3): E #eigenvalues in [{lo:.0f}, {hi:.0f})")

u11 = goe_vec[:, 0, 0] ** 2                                             # |U_11|^2, Haar: mean 1/3
mc(u11, 1 / 3, "GOE(3): E |U_11|^2")
mc((u11 - 1 / 3) * goe_lam[:, 0], 0.0, "GOE(3): E (|U_11|^2 - 1/3) lambda_1")   # independence (Theorem C.7)
```

## C.4 Polinomios ortogonales y la función de un punto

**Lema C.10 (función de un punto).** Sean $I\subseteq\mathbb R$ un intervalo y $w\ge0$ una función de Borel sobre $I$, positiva sobre un
conjunto de medida positiva, con $\int_I|x|^mw<\infty$ para todo $m$. Sean $p_0,\dots,p_{n-1}$ polinomios reales, $\deg p_k=k$
con coeficiente principal $\kappa_k\neq0$, ortonormales en $L^2(w)$ (existen: se aplica Gram–Schmidt a $1,x,\dots,x^{n-1}$, que son
independientes en $L^2(w)$ porque un polinomio no nulo tiene un número finito de ceros; se puede elegir cualquier signo para cada $p_k$). Sean $q=Z_n^{-1}\Delta(x)^2\prod_iw(x_i)$ sobre
$I^n$ con $Z_n=\int_{I^n}\Delta^2\prod_iw(x_i)\,dx$. Entonces $Z_n=n!\prod_k\kappa_k^{-2}$ y
$$
n\int_{I^{n-1}}q(x,y)\,dy=K_n(x)\,w(x),\qquad K_n(x)=\sum_{k=0}^{n-1}p_k(x)^2=v(x)^TM^{-1}v(x),
$$
donde $v(x)=(1,x,\dots,x^{n-1})$ y $M_{ij}=\int_Ix^{i+j}w$, $0\le i,j<n$.

*Demostración.* $\Delta(x)=\det(x_j^{i-1})_{i,j}$. Sumar a cada fila múltiplos de las filas situadas encima sustituye $x_j^{i-1}$ por
$p_{i-1}(x_j)/\kappa_{i-1}$, así que $\Delta(x)=(\prod_k\kappa_k)^{-1}\det(p_{i-1}(x_j))$ y
$$
\Delta(x)^2=\Big(\prod_k\kappa_k\Big)^{-2}\sum_{\sigma,\tau\in S_n}\operatorname{sgn}\sigma\operatorname{sgn}\tau\prod_{j=1}^np_{\sigma(j)-1}(x_j)\,p_{\tau(j)-1}(x_j).
$$
Integrando sobre $x_2,\dots,x_n$ respecto de $\prod_{j\ge2}w(x_j)$, el factor $j$ da $\delta_{\sigma(j)\tau(j)}$, así que solo sobrevive
$\sigma=\tau$, y el resultado es $(\prod_k\kappa_k)^{-2}\sum_\sigma p_{\sigma(1)-1}(x_1)^2=(\prod_k\kappa_k)^{-2}(n-1)!\,K_n(x_1)$.
Integrando respecto de $w(x_1)$ se obtiene $Z_n=(\prod_k\kappa_k)^{-2}n!$, y dividiendo se obtiene la fórmula. Por último, si
$p_k(x)=\sum_iC_{ki}x^i$, la ortonormalidad se escribe $CMC^T=1$ con $C$ invertible, así que $M^{-1}=C^TC$ y
$v^TM^{-1}v=|Cv|^2=\sum_kp_k^2$. $\square$

Por el corolario C.6, los valores propios no ordenados del GUE tienen la densidad $q$ del lema C.10 con $w(x)=e^{-x^2/2}$ sobre
$\mathbb R$, así que la densidad de niveles del GUE es $K_n(x)e^{-x^2/2}$. Lo comprobamos para $n=3$ con los momentos
$\int x^{2m}e^{-x^2/2}dx=\sqrt{2\pi}(2m-1)(2m-3)\cdots1$:

```python
def one_point(x, n, moment, w):
    """K_n(x) w(x) = v(x)^T M^{-1} v(x) w(x), M the matrix of moments (Lemma C.10)."""
    M = torch.tensor([[moment(i + j) for j in range(n)] for i in range(n)], dtype=f64)
    v = torch.stack([x ** i for i in range(n)], dim=-1)
    return (v @ torch.linalg.inv(M) * v).sum(-1) * w(x)

gauss_moment = lambda m: math.sqrt(2 * math.pi) * math.prod(range(m - 1, 0, -2)) if m % 2 == 0 else 0.0
Gc = torch.randn(B, n, n, dtype=dt)
gue_lam = torch.linalg.eigvalsh((Gc + Gc.mH) / math.sqrt(2))
for lo, hi in zip(edges[:-1], edges[1:]):
    x = torch.linspace(lo, hi, 401, dtype=f64)
    exact = torch.trapezoid(one_point(x, n, gauss_moment, lambda t: torch.exp(-t ** 2 / 2)), x).item()
    mc(((gue_lam >= lo) & (gue_lam < hi)).sum(-1).double(), exact, f"GUE(3): E #eigenvalues in [{lo:.0f}, {hi:.0f})")
```

## C.5 Matrices de Wishart y la medida inducida

Sean $a\le b$, $G$ una matriz de Ginibre $a\times b$ (entradas independientes con partes real e imaginaria independientes $N(0,\frac12)$,
densidad $\pi^{-ab}e^{-\operatorname{Tr}GG^*}$) y $W=GG^*$. Sea $T_+$ el conjunto de las matrices $a\times a$ triangulares inferiores con
diagonal positiva, un abierto de $\mathbb R^{a^2}$ en las coordenadas $L_{ii}$, $\operatorname{Re}L_{ij}$, $\operatorname{Im}L_{ij}$ ($i>j$), y
$P_a$ el conjunto de las matrices definidas positivas, abierto en las matrices hermíticas $a\times a$.

**Lema C.11 (descomposición de Bartlett).** Casi seguramente las filas $r_1,\dots,r_a$ de $G$ son linealmente independientes, y
Gram–Schmidt da $G=LV$ con $L\in T_+$ y $V$ de filas ortonormales, así que $W=LL^*$. Las entradas $L_{ij}$, $i\ge j$, son
independientes, $L_{ij}$ ($i>j$) es gaussiana compleja estándar y $L_{ii}^2$ tiene la densidad $x^{b-i}e^{-x}/(b-i)!$ sobre $(0,\infty)$.
Por tanto $L$ tiene la densidad $c_L\prod_iL_{ii}^{2(b-i)+1}e^{-\operatorname{Tr}LL^*}$ sobre $T_+$.

*Demostración.* Sea $\mathcal F_{i-1}=\sigma(r_1,\dots,r_{i-1})$ y supongamos (inducción) que $r_1,\dots,r_{i-1}$ son c.s. independientes, con
vectores de Gram–Schmidt $v_1,\dots,v_{i-1}$. Sea $u_1,\dots,u_b$ la base ortonormal obtenida por Gram–Schmidt a partir de
$v_1,\dots,v_{i-1},e_1,\dots,e_b$ (descartando los vectores dependientes); $u_k=v_k$ para $k<i$, y los $u_k$ son funciones de Borel de
$r_1,\dots,r_{i-1}$. Pongamos $c_k=\langle u_k,r_i\rangle$. Como $r_i$ es independiente de $\mathcal F_{i-1}$ y su densidad $\pi^{-b}e^{-\|r\|^2}$ es
invariante por las unitarias (proposición A.16(a)), condicionalmente a $\mathcal F_{i-1}$ el vector $c$ es gaussiano complejo estándar (Fubini). Ahora bien,
$L_{ij}=c_j$ para $j<i$ y $r_i-\sum_{j<i}c_jv_j=\sum_{k\ge i}c_ku_k$, así que $L_{ii}^2=\sum_{k=i}^b|c_k|^2>0$ c.s. La ley condicional de
la fila $i$ de $L$ dado $\mathcal F_{i-1}$ es, por tanto, una ley fija, así que la fila $i$ es independiente de las filas anteriores, con las coordenadas
$c_1,\dots,c_{i-1}$ y $L_{ii}$ independientes; por inducción todas las filas de $L$ son independientes. Por el teorema A.20(b) con $N=m=b-i+1$,
$L_{ii}^2=\|(c_i,\dots,c_b)\|^2$ tiene la densidad $s^{m-1}e^{-s}/(m-1)!$. La densidad de $L_{ii}=\sqrt s$ es entonces
$2t^{2(b-i)+1}e^{-t^2}/(b-i)!$, y $\operatorname{Tr}LL^*=\sum_{i\ge j}|L_{ij}|^2$. $\square$

**Lema C.12 (aplicación de Cholesky).** $\chi(L)=LL^*$ es una biyección $C^\infty$ de $T_+$ sobre $P_a$ con
$|\det D\chi(L)|=2^a\prod_{i=1}^aL_{ii}^{2(a-i)+1}$.

*Demostración.* *Sobreyectiva:* para $W\in P_a$ se aplica el lema 5.2 a $W^{1/2}=QR$; entonces $W=R^*R$ y $L=R^*\in T_+$. *Inyectiva:* si
$LL^*=MM^*$, entonces $M^{-1}L=M^*(L^*)^{-1}$ es triangular inferior y superior, luego diagonal con diagonal positiva, y
unitaria, luego $1$. *Jacobiano:* para $i\ge j$, $W_{ij}=L_{ij}L_{jj}+\sum_{k<j}L_{ik}\overline{L_{jk}}$. Ordenemos las coordenadas de $L$ y del
triángulo inferior de $W$ columna por columna y, dentro de la columna $j$, la diagonal primero. Entonces $W_{ij}$ depende solo de
$L_{ij}$ y de coordenadas anteriores, así que $D\chi$ es triangular, con entradas diagonales $2L_{jj}$ ($\partial W_{jj}/\partial L_{jj}$) y,
para cada $i>j$, dos veces $L_{jj}$ (partes real e imaginaria de $W_{ij}$). El producto es $\prod_j2L_{jj}\cdot L_{jj}^{2(a-j)}$. $\square$

**Teorema C.13 (conjuntos de Wishart y de Laguerre).** $W=GG^*$ tiene la densidad $c_W\det(W)^{b-a}e^{-\operatorname{Tr}W}$ sobre $P_a$
(nula fuera), y sus valores propios tienen la densidad
$$
Z_{a,b}^{-1}\,\Delta(\lambda)^2\prod_{i=1}^a\lambda_i^{\,b-a}e^{-\lambda_i}\quad\text{on }C_+=\{0<\lambda_1<\dots<\lambda_a\}.
$$

*Demostración.* Para $g\ge0$ de Borel, el lema C.11 y el cambio de variables $W=\chi(L)$ (lema C.12) dan
$\mathbb E\,g(W)=\int_{P_a}g(W)\,c_L\prod_iL_{ii}^{2(b-i)+1}e^{-\operatorname{Tr}W}\big/\big(2^a\prod_iL_{ii}^{2(a-i)+1}\big)dW$ con
$L=\chi^{-1}(W)$, y $\prod_iL_{ii}^{2(b-a)}=\det(W)^{b-a}$. Esta densidad es invariante por $W\mapsto VWV^*$, $V\in U(a)$, y
$\Lambda\in P_a$ si y solo si $\lambda_1>0$; se aplica el corolario C.6 con $\beta=2$. $\square$

**Teorema C.14 (medida inducida).** Sean $\psi$ un estado puro aleatorio de Haar de $\mathbb C^a\otimes\mathbb C^b$ (sección 7.5),
$\rho_A=\operatorname{Tr}_2|\psi\rangle\langle\psi|=W/s$ con $s=\operatorname{Tr}W$ (demostración del teorema 7.11), y
$\mu_1<\dots<\mu_a$ los valores propios de $\rho_A$. Entonces $s$ y $\mu$ son independientes, $s$ tiene la densidad
$s^{ab-1}e^{-s}/(ab-1)!$, y $(\mu_1,\dots,\mu_{a-1})$ tiene la densidad
$$
Z'^{-1}_{a,b}\,\Delta(\mu)^2\prod_{i=1}^a\mu_i^{\,b-a}\quad\text{on }D_a=\{0<\mu_1<\dots<\mu_{a-1}<\mu_a\},\qquad\mu_a=1-\sum_{i<a}\mu_i ,
$$
que es la densidad «$\propto\delta(1-\sum_i\mu_i)\,\Delta(\mu)^2\prod_i\mu_i^{b-a}$» de la literatura física.

*Demostración.* $\Theta(s,\mu_1,\dots,\mu_{a-1})=s\,(\mu_1,\dots,\mu_{a-1},1-\sum_{i<a}\mu_i)$ es una biyección $C^\infty$ de $(0,\infty)\times D_a$ sobre
$C_+$, con inversa $\lambda\mapsto(\sum\lambda_i,\lambda_1/\sum\lambda_i,\dots)$. Sumar todas las filas de $D\Theta$ a la última la convierte en
$(1,0,\dots,0)$, así que $|\det D\Theta|=s^{a-1}$. Como $\Delta$ es homogéneo de grado $a(a-1)/2$, la densidad de $(s,\mu)$ es
$$
Z_{a,b}^{-1}\,s^{a(a-1)+a(b-a)+a-1}e^{-s}\cdot\Delta(\mu)^2\prod_i\mu_i^{b-a}=Z_{a,b}^{-1}\,s^{ab-1}e^{-s}\cdot\Delta(\mu)^2\prod_i\mu_i^{b-a},
$$
un producto, lo que demuestra la independencia e identifica ambos factores tras la normalización. $\square$

Comprobamos el lema C.11 sobre el factor de Cholesky calculado por `torch.linalg.cholesky` (que es la única $L\in T_+$ por el
lema C.12), y el teorema C.14 para $a=2$, donde $\mu_1\in(0,\frac12)$ tiene la densidad $\propto(1-2t)^2(t(1-t))^{b-2}$:

```python
a, b, B = 3, 5, 100_000
G = torch.randn(B, a, b, dtype=dt)                        # E|G_ij|^2 = 1
W = G @ G.mH
L = torch.linalg.cholesky(W)                              # W = L L*, L in T_+
for i in range(a):                                        # 0-based i: L_ii^2 ~ Gamma(b - i) (Lemma C.11, index i + 1)
    mc(L[:, i, i].real ** 2, b - i, f"E L_{i + 1}{i + 1}^2")
    mc(L[:, i, i].real ** 4, (b - i) * (b - i + 1), f"E L_{i + 1}{i + 1}^4")
mc(L[:, 2, 0].abs() ** 2 * L[:, 1, 1].real ** 2, 1.0 * (b - 1), "E |L_31|^2 L_22^2")   # independence

for b2 in (2, 4):
    rho = states.random_density_matrix(TypeIAlgebra([2], [2], precision='double'), 100_000, rank=b2)
    mu1 = rho.eigenvalues()[0][:, 0]                      # smallest eigenvalue of rho_A
    t = torch.linspace(0, 0.5, 20_001, dtype=f64)
    dens = (1 - 2 * t) ** 2 * (t * (1 - t)) ** (b2 - 2)
    dens = dens / torch.trapezoid(dens, t)
    for lo, hi in ((0.0, 0.1), (0.1, 0.25), (0.25, 0.4)):
        sel = (t >= lo) & (t <= hi)
        mc(((mu1 >= lo) & (mu1 < hi)).double(), torch.trapezoid(dens[sel], t[sel]).item(),
           f"a = 2, b = {b2}: P({lo} <= mu_1 < {hi})")
```

## C.6 Polinomios de Laguerre y la fórmula de Page

**Lema C.15 (polinomios de Laguerre).** Para enteros $\alpha,k\ge0$ sea
$L_k^{(\alpha)}(x)=\sum_{j=0}^k(-1)^j\binom{k+\alpha}{k-j}\frac{x^j}{j!}$, y $L_{-1}^{(\alpha)}=0$. Entonces:

(a) $L_k^{(\alpha)}(x)\,x^\alpha e^{-x}=\frac1{k!}\frac{d^k}{dx^k}\big(x^{k+\alpha}e^{-x}\big)$; $L_k^{(\alpha)}$ tiene grado $k$ y coeficiente principal $(-1)^k/k!$.

(b) $\frac{d}{dx}L_k^{(\alpha)}=-L_{k-1}^{(\alpha+1)}$, luego $\frac{d^i}{dx^i}L_k^{(\alpha)}=(-1)^iL_{k-i}^{(\alpha+i)}$ para $i\le k$; y $L_k^{(\alpha)}=L_k^{(\alpha+1)}-L_{k-1}^{(\alpha+1)}$.

(c) $\int_0^\infty x^mL_k^{(\alpha)}x^\alpha e^{-x}dx=0$ para $m<k$, y $\int_0^\infty L_k^{(\alpha)}L_l^{(\alpha)}x^\alpha e^{-x}dx=\delta_{kl}\frac{(k+\alpha)!}{k!}$.

En consecuencia, por el lema C.10 y el teorema C.13, los valores propios de $W$ tienen la densidad de niveles
($\mathbb E\sum_i\varphi(\lambda_i(W))=\int_0^\infty\varphi\,\rho_{a,b}$ para toda $\varphi\ge0$ de Borel, y para toda $\varphi$ con $\int|\varphi|\rho_{a,b}<\infty$)
$$
\rho_{a,b}(x)=\sum_{k=0}^{a-1}\frac{k!}{(k+\alpha)!}L_k^{(\alpha)}(x)^2\,x^\alpha e^{-x},\qquad\alpha=b-a .
$$

*Demostración.* (a) Por Leibniz, $\frac{d^k}{dx^k}(x^{k+\alpha}e^{-x})=\sum_j\binom kj(-1)^j\frac{(k+\alpha)!}{(\alpha+j)!}x^{\alpha+j}e^{-x}$ (con $j$ derivadas
sobre $e^{-x}$), y $\frac1{k!}\binom kj\frac{(k+\alpha)!}{(\alpha+j)!}=\binom{k+\alpha}{k-j}\frac1{j!}$. (b) Derivando término a término,
$-\sum_{i=0}^{k-1}(-1)^i\binom{k+\alpha}{k-1-i}\frac{x^i}{i!}=-L_{k-1}^{(\alpha+1)}$; la segunda identidad es la regla de Pascal
$\binom{k+\alpha+1}{k-j}-\binom{k+\alpha}{k-1-j}=\binom{k+\alpha}{k-j}$. (c) Se integra (a) por partes $k$ veces contra $x^m$: los
términos de borde se anulan, porque $\frac{d^r}{dx^r}(x^{k+\alpha}e^{-x})$ tiene el factor $x^{k+\alpha-r}$, $r<k$, y decae en $\infty$. El
resultado es $\frac{(-1)^k}{k!}\int\frac{d^kx^m}{dx^k}x^{k+\alpha}e^{-x}dx$, que vale $0$ para $m<k$ y $(-1)^k(k+\alpha)!$ para $m=k$. Para $l<k$,
$L_l^{(\alpha)}$ es combinación de $x^m$, $m<k$; para $l=k$ solo cuenta el coeficiente principal $(-1)^k/k!$. Así,
$p_k=\big(k!/(k+\alpha)!\big)^{1/2}L_k^{(\alpha)}$ son ortonormales para $w(x)=x^\alpha e^{-x}$ sobre $(0,\infty)$, y los valores propios no ordenados
de $W$ tienen la densidad $q$ del lema C.10 (con $n=a$) por el corolario C.6, de modo que, como en la demostración del corolario C.8,
$\mathbb E\sum_i\varphi(\lambda_i)=a\int\varphi(x_1)q(x)\,dx=\int\varphi\,K_aw=\int\varphi\,\rho_{a,b}$ para toda $\varphi\ge0$ de Borel, y por linealidad para
$\varphi$ con $\int|\varphi|\rho_{a,b}<\infty$. $\square$

**Lema C.16 (integrales logarítmicas).** Sean $I_m=\int_0^\infty x^m\log x\,e^{-x}dx$ y $H_m=\sum_{i=1}^m\frac1i$. Entonces
$I_m=m!\,(I_0+H_m)$. Para enteros $\alpha\ge0$ y $0\le k\le m$,
$$
A^{(\alpha)}_{k,m}:=\int_0^\infty\log x\,L_k^{(\alpha)}(x)L_m^{(\alpha)}(x)\,x^\alpha e^{-x}dx=
\begin{cases}I_{k+\alpha}/k!, & m=k,\\ -\dfrac{(k+\alpha)!}{k!\,(m-k)}, & m>k.\end{cases}
$$
($I_0$ se cancela en el teorema C.17, así que su valor no es necesario.)

*Demostración.* Para $m\ge1$, la integración por partes con $(-e^{-x})'=e^{-x}$ da $I_m=mI_{m-1}+(m-1)!$, es decir,
$I_m/m!=I_{m-1}/(m-1)!+1/m$. Para $A_{k,m}$, se escribe $L_m^{(\alpha)}x^\alpha e^{-x}=\frac1{m!}g^{(m)}$ con $g=x^{m+\alpha}e^{-x}$
(lema C.15(a)) y $f=\log x\,L_k^{(\alpha)}$, y se integra por partes $m$ veces. Los términos de borde
$f^{(i)}g^{(m-1-i)}$ se anulan en $\infty$ y en $0$, donde $f^{(i)}=O(|\log x|+x^{-i})$ y $g^{(m-1-i)}=O(x^{\alpha+1+i})$. Así,
$A_{k,m}=\frac{(-1)^m}{m!}\int_0^\infty f^{(m)}g\,dx$ (el integrando es $O(x^\alpha(1+|\log x|))$ en $0$) con
$f^{(m)}=\sum_{i\le k}\binom mi(\log x)^{(m-i)}\frac{d^i}{dx^i}L_k^{(\alpha)}$. Para $m-i\ge1$, $(\log x)^{(m-i)}=(-1)^{m-i-1}(m-i-1)!\,x^{i-m}$,
y por el lema C.15(b,c)
$$
\int_0^\infty x^{i-m}\frac{d^iL_k^{(\alpha)}}{dx^i}\,x^{m+\alpha}e^{-x}dx=(-1)^i\int_0^\infty L_{k-i}^{(\alpha+i)}\,x^{\alpha+i}e^{-x}dx=(-1)^i(\alpha+i)!\,\delta_{ik}.
$$
Si $m>k$, solo contribuye $i=k$: $A_{k,m}=\frac{(-1)^m}{m!}\binom mk(-1)^{m-k-1}(m-k-1)!(-1)^k(\alpha+k)!=-\frac{(k+\alpha)!}{k!(m-k)}$.
Si $m=k$, los términos $i<k$ se anulan y el término $i=k$ es $\log x\cdot(-1)^k$, lo que da $A_{k,k}=\frac{(-1)^k}{k!}(-1)^kI_{k+\alpha}$. $\square$

**Teorema C.17 (fórmula de Page).** Para un estado puro aleatorio de Haar de $\mathbb C^a\otimes\mathbb C^b$ con $a\le b$,
$$
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
$$

*Demostración.* *Reducción.* Con la notación del teorema C.14, $\lambda_i=s\mu_i$ y $\sum_i\mu_i=1$, así que
$\sum_i\lambda_i\log\lambda_i=s\log s-s\,S(\rho_A)$. Todos los términos son integrables ($|x\log x|\le1+x^2$, y $\mathbb Es^2<\infty$,
$0\le S(\rho_A)\le a/e$ pues $0\le-x\log x\le1/e$ sobre $[0,1]$), $\mathbb E\sum_i\lambda_i\log\lambda_i=\int x\log x\,\rho_{a,b}(x)\,dx$ por el lema C.15, y $s$ es independiente de
$S(\rho_A)$, así que $\mathbb E\sum_i\lambda_i\log\lambda_i=\mathbb E[s\log s]-\mathbb E[s]\,\mathbb E\,S(\rho_A)$. De la densidad de $s$,
$\mathbb Es=ab$ y $\mathbb E[s\log s]=I_{ab}/(ab-1)!=ab\,(I_0+H_{ab})$ (lema C.16). Por tanto
$$
\mathbb E\,S(\rho_A)=I_0+H_{ab}-\frac1{ab}\int_0^\infty x\log x\,\rho_{a,b}(x)\,dx .
$$
*La integral.* Por el lema C.15(b), $L_k^{(\alpha)}=L_k^{(\alpha+1)}-L_{k-1}^{(\alpha+1)}$, así que con $x\cdot x^\alpha=x^{\alpha+1}$ y el lema C.16
para $\alpha+1$,
$$
\frac{k!}{(k+\alpha)!}\int_0^\infty x\log x\,L_k^{(\alpha)}(x)^2x^\alpha e^{-x}dx=\frac{k!}{(k+\alpha)!}\Big(A^{(\alpha+1)}_{k,k}-2A^{(\alpha+1)}_{k-1,k}+A^{(\alpha+1)}_{k-1,k-1}\Big)
$$
$$
=(k+\alpha+1)(I_0+H_{k+\alpha+1})+2k+k\,(I_0+H_{k+\alpha}),
$$
también para $k=0$. Como $\sum_{k=0}^{a-1}(2k+\alpha+1)=a(a-1)+a(b-a+1)=ab$,
$$
\int_0^\infty x\log x\,\rho_{a,b}\,dx=ab\,I_0+\Sigma,\qquad \Sigma=\sum_{k=0}^{a-1}\big[(k+\alpha+1)H_{k+\alpha+1}+kH_{k+\alpha}+2k\big].
$$
*La suma.* Pongamos $m=\alpha+1=b-a+1$ y $j=k+m\in\{m,\dots,b\}$. Como $(j-m)H_{j-1}=(j-m)H_j-1+\frac mj$,
$\Sigma=\sum_{j=m}^b(2j-m)H_j-a+m(H_b-H_{m-1})+a(a-1)$. Con $F(j)=j(j+1-m)$ se tiene $F(j)-F(j-1)=2j-m$, $F(m-1)=0$,
$F(b)=ab$, y la suma por partes da
$$
\sum_{j=m}^b(2j-m)H_j=abH_b-\sum_{j=m}^{b-1}\frac{F(j)}{j+1}=abH_b-\sum_{j=m}^{b-1}\Big(j-m+\frac m{j+1}\Big)=abH_b-\frac{(a-1)(a-2)}2-m(H_b-H_m).
$$
Como $m(H_m-H_{m-1})=1$, $\Sigma=abH_b-\frac{(a-1)(a-2)}2+1-a+a(a-1)=abH_b+\frac{a(a-1)}2$. Por tanto
$\mathbb E\,S(\rho_A)=H_{ab}-H_b-\frac{a-1}{2b}$. $\square$

**Corolario C.18.** Para $a\le b$,
$\log a-\frac{a-1}{ab+1}-\frac{a-1}{2b}\le\mathbb E\,S(\rho_A)\le\log a-\frac{a-1}{2b}$; en particular
$\big|\mathbb E\,S(\rho_A)-\log a+\frac a{2b}\big|\le\frac1{2b}$. Para dos cúbits, $\mathbb E\,S(\rho_A)=\frac13$.

*Demostración.* Como $1/x$ es decreciente, $\int_k^{k+1}\frac{dx}x\le\frac1k\le\int_{k-1}^k\frac{dx}x$; sumando sobre $b<k\le ab$,
$\log\frac{ab+1}{b+1}\le H_{ab}-H_b\le\log a$, y $\log\frac{ab+1}{b+1}=\log a-\log\big(1+\frac{a-1}{ab+1}\big)\ge\log a-\frac{a-1}{ab+1}$.
Por último, $\frac{a-1}{ab+1}\le\frac1b$, y $a=b=2$ da $\frac13+\frac14-\frac14$. $\square$

La fórmula fue conjeturada por D. N. Page (Phys. Rev. Lett. 71, 1291, 1993) y demostrada por S. K. Foong y S. Kanno
(Phys. Rev. Lett. 72, 1148, 1994), J. Sánchez-Ruiz (Phys. Rev. E 52, 5653, 1995) y S. Sen (Phys. Rev. Lett. 77, 1, 1996);
el cálculo anterior sigue la reducción de Sen al conjunto de Laguerre, pero suma la densidad de un punto término a término
en lugar de usar la fórmula de Christoffel–Darboux. La medida inducida es la de K. Życzkowski y H.-J. Sommers
(J. Phys. A 34, 7111, 2001); la descomposición del lema C.11 se debe a M. S. Bartlett (1933), y el teorema C.5 a H. Weyl.

Comprobamos la densidad de niveles de Laguerre frente a un histograma, los valores $A^{(\alpha)}_{k,m}$ del lema C.16 y el resultado
intermedio $\int x\log x\,\rho_{a,b}=ab\,I_0+\Sigma$ por cuadratura (con el propio $I_0$ calculado por cuadratura), la identidad para
$\Sigma$ en aritmética exacta, la fórmula de Page frente a una integración directa de la
densidad inducida para $a=2$ (valores racionales exactos) y frente a Monte Carlo, y el corolario C.18:

```python
def laguerre(k, al, x):
    return sum((-1) ** j * math.comb(k + al, k - j) * x ** j / math.factorial(j) for j in range(k + 1))

def rho_laguerre(x, a, b):                                # Lemma C.15
    al = b - a
    return sum(math.factorial(k) / math.factorial(k + al) * laguerre(k, al, x) ** 2 for k in range(a)) * x ** al * torch.exp(-x)

lam = torch.linalg.eigvalsh(W)                            # a = 3, b = 5 from above
x = torch.linspace(0, 80, 160_001, dtype=f64)
r = rho_laguerre(x, a, b)
assert abs(torch.trapezoid(r, x).item() - a) < 1e-8 and abs(torch.trapezoid(x * r, x).item() - a * b) < 1e-8
for lo, hi in ((0, 1), (1, 2), (2, 4), (4, 8), (8, 12), (12, 20)):
    sel = (x >= lo) & (x <= hi)
    mc(((lam >= lo) & (lam < hi)).sum(-1).double(), torch.trapezoid(r[sel], x[sel]).item(),
       f"Laguerre(3, 5): E #eigenvalues in [{lo}, {hi})")

def Hn(m):
    return sum(Fraction(1, i) for i in range(1, m + 1))

def page(a, b):
    return Hn(a * b) - Hn(b) - Fraction(a - 1, 2 * b)

vv = torch.linspace(0, 12, 1_200_001, dtype=f64)[1:]      # x = v^2 smooths log x at 0
def log_integral(g):                                      # int_0^oo log(x) g(x) dx
    xv = vv ** 2
    return torch.trapezoid(xv.log() * g(xv) * 2 * vv, vv).item()
I0 = log_integral(lambda t: torch.exp(-t))
print(f"I_0 = {I0:.10f} (quadrature)")
for al in range(3):                                       # Lemma C.16
    for k in range(4):
        for m in range(k, 5):
            A = log_integral(lambda t: laguerre(k, al, t) * laguerre(m, al, t) * t ** al * torch.exp(-t))
            Im = math.factorial(k + al) * (I0 + float(Hn(k + al)))
            assert abs(A - (Im / math.factorial(k) if m == k else -math.factorial(k + al) / (math.factorial(k) * (m - k)))) < 1e-6
for a_, b_ in ((2, 2), (2, 3), (3, 3), (3, 5), (4, 6)):
    al = b_ - a_
    Sigma = sum((k + al + 1) * Hn(k + al + 1) + k * Hn(k + al) + 2 * k for k in range(a_))
    assert Sigma == a_ * b_ * Hn(b_) + Fraction(a_ * (a_ - 1), 2)            # proof of Theorem C.17, exactly
    integral = log_integral(lambda t: t * rho_laguerre(t, a_, b_))
    assert abs(integral - (a_ * b_ * I0 + float(Sigma))) < 1e-6

print("a = 2:", ", ".join(f"b = {b_}: {page(2, b_)}" for b_ in range(2, 7)))
assert page(2, 2) == Fraction(1, 3)
for b_ in range(2, 7):                                    # a = 2: integrate the density of Theorem C.14
    v = torch.linspace(0, math.sqrt(0.5), 200_001, dtype=f64)[1:]
    t = v ** 2                                            # mu_1 = v^2 smooths the logarithm at 0
    dens = (1 - 2 * t) ** 2 * (t * (1 - t)) ** (b_ - 2) * 2 * v
    S = -(t * t.log() + (1 - t) * (1 - t).log())
    assert abs((torch.trapezoid(dens * S, v) / torch.trapezoid(dens, v)).item() - float(page(2, b_))) < 1e-9

for a_, b_ in ((2, 2), (2, 3), (3, 3), (3, 5)):
    rho = states.random_density_matrix(TypeIAlgebra([a_], [a_], precision='double'), 40_000, rank=b_)
    mc(rho.von_neumann_entropy(), float(page(a_, b_)), f"E S(rho_A), a = {a_}, b = {b_}")   # Theorem C.17

for a_ in range(1, 12):
    for b_ in range(a_, 40):
        p = float(page(a_, b_))
        assert math.log(a_) - (a_ - 1) / (a_ * b_ + 1) - (a_ - 1) / (2 * b_) - 1e-12 <= p <= math.log(a_) - (a_ - 1) / (2 * b_) + 1e-12
        assert abs(p - math.log(a_) + a_ / (2 * b_)) <= 1 / (2 * b_) + 1e-12  # Corollary C.18
```

**Volver:** el [capítulo 5](05_random_matrices.md) utiliza el teorema C.7, el corolario C.8 y la proposición C.9;
el [capítulo 7](07_tensor_products.md) utiliza los teoremas C.14 y C.17 y el corolario C.18.
