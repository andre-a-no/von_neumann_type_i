[English](../13_jordan_algebras.md) · **Español** · [Français](../fr/13_jordan_algebras.md) · [Deutsch](../de/13_jordan_algebras.md) · [中文](../zh/13_jordan_algebras.md) · [日本語](../ja/13_jordan_algebras.md) · [Русский](../ru/13_jordan_algebras.md)

# 13. Álgebras de Jordan: por qué la mecánica cuántica usa *-álgebras

El producto $xy$ de dos observables no es un observable salvo que $x$ e $y$ conmuten, pero el producto simetrizado
$x\circ y=\frac12(xy+yx)$ siempre lo es. En 1933 P. Jordan propuso tomar los observables con este producto como la
estructura básica de la mecánica cuántica, y Jordan, von Neumann y Wigner clasificaron las álgebras de dimensión finita
que aparecen. Este capítulo desarrolla la teoría en dimensión finita, con demostraciones completas cuando son cortas, y
se pregunta qué se le escapa a la estructura de Jordan. La respuesta es el signo del conmutador: el álgebra de Jordan
de un sistema determina el conmutador solo salvo signo (sección 13.5), pero los sistemas compuestos necesitan el signo,
y este decide si un estado está entrelazado.

**Aprenderá** el producto de Jordan y la representación cuadrática $P(y)x=yxy$, que, a diferencia de $\circ$,
conserva la positividad; los factores de espín, su teoría espectral y el cono de Lorentz; subálgebras de Jordan de
matrices hermíticas que no son partes autoadjuntas de $*$-álgebras; la clasificación de Jordan–von Neumann–Wigner; por
qué la traspuesta es un automorfismo de Jordan y no un $*$-automorfismo, y qué hace sobre un sistema compuesto; y una
desigualdad de trazas del capítulo 11 escrita en lenguaje de Jordan y demostrada en un factor de espín.

En todo el capítulo, $\operatorname{Herm}_n(\mathbb K)$ es el espacio vectorial real de las matrices hermíticas
$n\times n$ sobre $\mathbb K=\mathbb R$ o $\mathbb C$, $M_{sa}$ es la parte autoadjunta de
$\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$, $X,Y,Z$ son las matrices de Pauli y $\sigma=(X,Y,Z)$, de modo que
$u\cdot\sigma=u_1X+u_2Y+u_3Z$ para $u\in\mathbb R^3$. $[a,b]=ab-ba$, y $|u|$ es la norma euclídea de $u\in\mathbb R^m$.

## 13.1 Los observables y el producto de Jordan

**Definición.** Un **álgebra de Jordan** es un espacio vectorial real $J$ con un producto bilineal y conmutativo
$x\circ y$ que satisface la **identidad de Jordan** $(x^2\circ y)\circ x=x^2\circ(y\circ x)$, donde $x^2=x\circ x$. Es
**especial** si es isomorfa a un subespacio de un álgebra asociativa cerrado bajo $\frac12(xy+yx)$. Sobre $J$,
$L(y)x=y\circ x$, y la **representación cuadrática** es $P(y)=2L(y)^2-L(y^2)$.

**Lema 13.1.** Sea $\mathcal A$ un álgebra asociativa real o compleja y $x\circ y=\frac12(xy+yx)$. Entonces
$\circ$ satisface la identidad de Jordan, y las potencias de Jordan $x^{\circ1}=x$, $x^{\circ(k+1)}=x\circ x^{\circ k}$
coinciden con las potencias asociativas $x^k$. En particular, $M_{sa}$, $\operatorname{Herm}_n(\mathbb C)$ y
$\operatorname{Herm}_n(\mathbb R)$ son álgebras de Jordan especiales.

*Demostración.* Desarrollando, $4(x^2\circ y)\circ x=x^2yx+yx^3+x^3y+xyx^2$ y
$4x^2\circ(y\circ x)=x^2yx+x^3y+yx^3+xyx^2$, los mismos cuatro términos. Por inducción,
$x\circ x^k=\frac12(x^{k+1}+x^{k+1})=x^{k+1}$. Si $x,y$ son hermíticas, $(xy+yx)^*=yx+xy$, de modo que los tres
espacios son cerrados bajo $\circ$; para $x,y$ simétricas reales el producto es real. $\square$

**Proposición 13.2.** Sean $x,y,z\in M_{sa}$.

1. $P(y)x=yxy$, y se cumple la **fórmula fundamental** $P(P(y)x)=P(y)P(x)P(y)$.
2. $\operatorname{Tr}P(y)x=\operatorname{Tr}y^2x$.
3. Si $x\ge0$, entonces $P(y)x\ge0$. Si $y$ es invertible, $P(y)$ es una biyección del cono positivo sobre sí mismo,
   con inversa $P(y^{-1})$.
4. El producto de Jordan de elementos positivos no tiene por qué ser positivo: para $p=\begin{pmatrix}1&0\\0&0\end{pmatrix}$ y
   $q=\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}$, $\det(p\circ q)=-\frac1{16}$.

*Demostración.* (1) $2y\circ(y\circ x)=\frac12(y^2x+2yxy+xy^2)$ e $y^2\circ x=\frac12(y^2x+xy^2)$; réstense. Entonces
$P(P(y)x)z=yxy\,z\,yxy=P(y)P(x)P(y)z$ por asociatividad. (2) Ciclicidad de la traza. (3)
$yxy=(x^{1/2}y)^*(x^{1/2}y)\ge0$ (corolario 1.7), y $P(y)P(y^{-1})=P(y^{-1})P(y)=\mathrm{id}$. (4)
$p\circ q=\frac14\begin{pmatrix}2&1\\1&0\end{pmatrix}$. $\square$

Escribiendo $P(y)$ como $2L(y)^2-L(y^2)$, cada identidad de (1) se expresa solo con $\circ$, así que vale en toda álgebra de Jordan especial; la fórmula
fundamental vale en toda álgebra de Jordan (K. McCrimmon, *A Taste of Jordan Algebras*, Springer 2004, que no
reproducimos). El módulo `torch_vn_algebra.jordan` implementa estas aplicaciones sobre operadores de un `TypeIAlgebra`,
de forma perezosa y en lotes sobre todos los sectores:

```python
import math
import torch
from torch_vn_algebra import TypeIAlgebra, SpinFactor
from torch_vn_algebra import jordan as jd

torch.manual_seed(0)
dt = torch.complex128
jp, P = jd.jordan_product, jd.quadratic_representation

def close(A, B, tol=1e-10):
    A = A.matrix if hasattr(A, 'matrix') else A
    B = B.matrix if hasattr(B, 'matrix') else B
    return (A - B).abs().max().item() < tol

alg = TypeIAlgebra([2, 3], [2, 3], precision='double')           # M_2 (+) M_3
def sa(B=4):
    return alg.operator_from_eigenvalues(lambda d: torch.randn(B, d), batch_size=B, force_self_adjoint=True)
x, y, z = sa(), sa(), sa()
x2 = jp(x, x)
assert close(jp(jp(x2, y), x), jp(x2, jp(y, x)))                  # Jordan identity (Lemma 13.1)
assert close(jd.jordan_power(x, 5), x @ x @ x @ x @ x)            # power associativity
assert close(P(y, x), 2 * jp(y, jp(y, x)) - jp(jp(y, y), x))      # P(y) = 2L(y)^2 - L(y^2)
assert close(P(P(y, x), z), P(y, P(x, P(y, z))))                  # fundamental formula
assert close(P(y, x).trace, (jp(y, y) @ x).trace)                 # Tr P(y)x = Tr y^2 x
pos = alg.operator_from_eigenvalues(lambda d: torch.rand(4, d), batch_size=4, force_positive=True)
assert all((e > -1e-12).all() for e in P(y, pos).eigenvalues())  # Proposition 13.2(3)

M2 = TypeIAlgebra([2], [2], precision='double')
p = M2.operator(torch.tensor([[[[1.0, 0], [0, 0]]]]), is_self_adjoint=True, is_positive=True)
q = M2.operator(torch.tensor([[[[0.5, 0.5], [0.5, 0.5]]]]), is_self_adjoint=True, is_positive=True)
print(jp(p, q).eigenvalues()[0])                                  # (1 -+ sqrt 2)/4: not positive
assert abs(torch.linalg.det(jp(p, q).matrix[0, 0]).real.item() + 1 / 16) < 1e-14
```

> **Significado físico.** $x\circ y$ es el producto simetrizado de observables: la única combinación
> $\alpha xy+\beta yx$ que es simétrica en $x,y$ y se reduce a $x^2$ para $y=x$, y es de nuevo un observable. La
> representación cuadrática es la **regla de Lüders**: medir una pregunta de sí o no $p$ y quedarse con el resultado "sí"
> transforma una densidad $\rho$ en $P(p)\rho/\operatorname{Tr}P(p)\rho=p\rho p/\operatorname{Tr}(p\rho)$. La
> proposición 13.2(4) dice que $\circ$ por sí solo no tiene tal significado operacional para elementos positivos.

## 13.2 Factores de espín y el cono de Lorentz

**Definición.** Para $m\ge1$, el **factor de espín** $V_m=\mathbb R\oplus\mathbb R^m$ tiene elementos $x=(a,u)$, el
producto $(a,u)\circ(b,v)=(ab+u\cdot v,\;av+bu)$ y la unidad $1=(1,0)$. Su **traza** es $\operatorname{tr}x=2a$.

El producto es bilineal y conmutativo; la identidad de Jordan se sigue de una representación matricial.

**Lema 13.3 (simetrías que anticonmutan).** Sean $p=\lfloor m/2\rfloor$ y $D=2^p$ ($D=2$ si $m=1$). Sobre
$(\mathbb C^2)^{\otimes p}$ pongamos $s_{2j-1}=Z^{\otimes(j-1)}\otimes X\otimes1$, $s_{2j}=Z^{\otimes(j-1)}\otimes Y\otimes1$
para $j\le p$, y $s_m=Z^{\otimes p}$ si $m$ es impar ($s_1=Z$ si $m=1$). Entonces $s_i=s_i^*$ y
$s_is_j+s_js_i=2\delta_{ij}1$. La aplicación $\phi(a,u)=a1+\sum_iu_is_i$ es inyectiva y satisface
$\phi(x\circ y)=\phi(x)\circ\phi(y)$; por tanto, $V_m$ es un álgebra de Jordan especial. Para $m\ge2$, las $s_i$
generan $M_D(\mathbb C)$ como $*$-álgebra.

*Demostración.* Las matrices de Pauli tienen cuadrado $1$ y anticonmutan dos a dos; lo mismo hace $Z$ con $X$ y con
$Y$. Dos cadenas distintas $s_i,s_j$ difieren exactamente en un sitio de forma anticonmutante: para $s_{2j-1},s_{2j}$
en el sitio $j$ ($X$ frente a $Y$); para cadenas con índices $j<j'$, en el sitio $j$ ($X$ o $Y$ frente a $Z$),
mientras que en los demás sitios llevan factores iguales o un factor $1$; lo mismo vale para $s_m=Z^{\otimes p}$. Por
tanto, $s_is_j=-s_js_i$ para $i\ne j$, y $s_i^2=1$. Así, $(u\cdot s)\circ(v\cdot s)=\sum_{ij}u_iv_j\,s_i\circ s_j=(u\cdot v)1$,
lo que da $\phi(x)\circ\phi(y)=\phi(x\circ y)$. Cada $s_i$ tiene traza nula (es un producto tensorial con un factor
de traza nula) y $\operatorname{Tr}(s_is_j)=\operatorname{Tr}(s_i\circ s_j)=D\delta_{ij}$, de modo que $1,s_1,\dots,s_m$
son ortogonales para $\operatorname{Tr}(ab)$ y $\phi$ es inyectiva. La identidad de Jordan vale en
$\operatorname{Herm}_D$ (lema 13.1) y $\phi$ la transporta de vuelta. Por último,
$s_{2j-1}s_{2j}=1\otimes\cdots\otimes XY\otimes\cdots\otimes1=iZ_j$ ($Z$ en el sitio $j$), de modo que la
$*$-álgebra generada contiene $Z_j$, luego $X_j=Z_1\cdots Z_{j-1}s_{2j-1}$ e $Y_j$, y por tanto todos los productos
tensoriales de matrices de Pauli, una base de $M_D(\mathbb C)$. $\square$

Salvo signos, $s_{2j-1}$ y $s_{2j}$ son los operadores de Majorana $c_j+c_j^*$ e $i(c_j-c_j^*)$ de los fermiones de
Jordan–Wigner $c_j$ de la sección 9.4.

**Teorema 13.4 (teoría espectral de $V_m$).** Sean $x=(a,u)\in V_m$ y $\lambda_\pm=a\pm|u|$.

1. Los idempotentes de $V_m$ son $0$, $1$ y $c_w=\frac12(1,w)$ con $|w|=1$; y $c_w\circ c_{w'}=0$ si y solo si
   $w'=-w$.
2. Si $u\ne0$, $x=\lambda_+e_++\lambda_-e_-$ con $e_\pm=c_{\pm u/|u|}$; si $u=0$, $x=a1$. Si $x=\mu_1c_w+\mu_2c_{-w}$,
   entonces $\{\mu_1,\mu_2\}=\{\lambda_+,\lambda_-\}$. Llamamos a $\lambda_\pm$ los **valores propios** de $x$.
3. Para $f:\{\lambda_+,\lambda_-\}\to\mathbb R$ pongamos $f(x)=f(\lambda_+)e_++f(\lambda_-)e_-$ ($f(a)1$ si $u=0$).
   Entonces $x^{\circ k}=\lambda_+^ke_++\lambda_-^ke_-$, $\phi(f(x))=f(\phi(x))$, los valores propios de $\phi(x)$ son
   $\lambda_\pm$, cada uno con multiplicidad $D/2$, y $\operatorname{Tr}\phi(x)=\frac D2\operatorname{tr}x$. Además,
   $\operatorname{tr}x=\lambda_++\lambda_-$ y $\det x:=\lambda_+\lambda_-=a^2-|u|^2$.
4. El cono de los cuadrados $Q=\{x^2:x\in V_m\}$ es el **cono de Lorentz** $\{(a,u):a\ge|u|\}=\{x:\lambda_-\ge0\}=\{x:\phi(x)\ge0\}$.
   Es **autodual**: $y\in Q$ si y solo si $\operatorname{tr}(x\circ y)\ge0$ para todo $x\in Q$.
5. $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$; por tanto, $P(y)Q\subseteq Q$ para todo $y\in V_m$.

*Demostración.* (1) $(a,u)^2=(a^2+|u|^2,2au)$ es igual a $(a,u)$ si y solo si $u=0$, $a\in\{0,1\}$, o $a=\frac12$,
$|u|=\frac12$. $c_w\circ c_{w'}=\frac14(1+w\cdot w',\,w+w')$ se anula si y solo si $w'=-w$. (2)
$\lambda_+e_++\lambda_-e_-=(a,0)+|u|(0,u/|u|)=x$, y $\mu_1c_w+\mu_2c_{-w}=(\frac{\mu_1+\mu_2}2,\frac{\mu_1-\mu_2}2w)$, luego
$\mu_1+\mu_2=2a$ y $|\mu_1-\mu_2|=2|u|$. (3) Por (1), $(\alpha e_++\beta e_-)\circ(\gamma e_++\delta e_-)=\alpha\gamma
e_++\beta\delta e_-$, lo que da las potencias por inducción. Con $w=u/|u|$, $(w\cdot s)^2=1$, así que
$\phi(e_\pm)=\frac12(1\pm w\cdot s)$ son proyecciones ortogonales de suma $1$ y producto $0$, cada una de traza $D/2$;
por tanto, $\phi(x)=\lambda_+\phi(e_+)+\lambda_-\phi(e_-)$ es una descomposición espectral, y
$f(\phi(x))=f(\lambda_+)\phi(e_+)+f(\lambda_-)\phi(e_-)=\phi(f(x))$ (teorema 1.5). (4) $x^2$ tiene valores propios
$\lambda_\pm^2\ge0$; recíprocamente, si $\lambda_-\ge0$, entonces $x=(\sqrt x)^2$. $\lambda_-\ge0$ si y solo si
$a\ge|u|$, si y solo si $\phi(x)\ge0$ por (3). Para la autodualidad, $\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$. Si
$a\ge|u|$ y $b\ge|v|$, entonces $ab\ge|u||v|\ge-u\cdot v$. Si $b<|v|$, entonces $x=(|v|,-v)\in Q$ da
$2(|v|b-|v|^2)<0$ (para $v=0$, $b<0$, tómese $x=1$). (5) $P(y)x=2y\circ(y\circ x)-y^2\circ x$ y $\phi$ conserva
$\circ$, así que $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ por la proposición 13.2(1); aplíquense (4) y la
proposición 13.2(3). $\square$

**Proposición 13.5.** (a) Para $m=3$, $\phi(a,u)=a1+u\cdot\sigma$ es un isomorfismo de Jordan
$V_3\cong\operatorname{Herm}_2(\mathbb C)$. (b) Los automorfismos de Jordan de $V_m$ son las aplicaciones
$(a,u)\mapsto(a,Ou)$, $O\in O(m)$. (c) Para $m\ge2$, $\phi(V_m)$ tiene dimensión $m+1$, y la parte autoadjunta de la
$*$-álgebra que genera tiene dimensión $D^2=4^{\lfloor m/2\rfloor}$; solo son iguales para $m=3$.

*Demostración.* (a) Para $m=3$, $(s_1,s_2,s_3)=(X,Y,Z)$, y $\phi$ es inyectiva entre espacios de dimensión $4$. (b) Las
aplicaciones $(a,u)\mapsto(a,Ou)$ conservan $u\cdot v$ y, por tanto, $\circ$. Recíprocamente, sea $\alpha$ un
automorfismo de Jordan. Entonces $\alpha(1)$ es una unidad de $V_m$, así que $\alpha(1)=1$. Por la fórmula de $x^2$,
$x^2=1$ si y solo si $x=\pm1$ o $x=(0,w)$ con $|w|=1$. Para $|w|=1$, $\alpha(0,w)^2=1$ y
$\alpha(0,w)\ne\pm1=\alpha(\pm1)$, luego $\alpha(0,w)=(0,w')$ con $|w'|=1$. Escribiendo $\alpha(0,u)=(\ell(u),Ou)$ con
$\ell,O$ lineales, esto dice que $\ell=0$ y $|Ou|=|u|$, es decir, $O\in O(m)$. (c) Lema 13.3;
$4^{\lfloor m/2\rfloor}>m+1$ para $m=2$ y $m\ge4$. $\square$

```python
V = SpinFactor(5)
xs, ys = V.random(1000), V.random(1000)                          # (1000, 6): batches of (a, u)
x2 = V.square(xs)
assert (V.product(V.product(x2, ys), xs) - V.product(x2, V.product(ys, xs))).abs().max() < 1e-12
lam, e = V.spectral_decomposition(xs)                            # (lambda-, lambda+), (e-, e+)
assert (lam[:, 1, None] * e[:, 1] + lam[:, 0, None] * e[:, 0] - xs).abs().max() < 1e-12
Xm, Ym, Dm = V.to_matrix(xs), V.to_matrix(ys), V.representation_dim  # phi into Herm_4(C)
assert (V.to_matrix(V.product(xs, ys)) - 0.5 * (Xm @ Ym + Ym @ Xm)).abs().max() < 1e-12
assert (torch.linalg.eigvalsh(Xm) - lam.repeat_interleave(Dm // 2, dim=-1)).abs().max() < 1e-12
assert (V.to_matrix(V.exp(xs)) - torch.linalg.matrix_exp(Xm)).abs().max() < 1e-9       # phi(f(x)) = f(phi(x))
cone = V.random(1000, cone=True)                                 # squares
assert V.in_cone(cone).all() and (torch.linalg.eigvalsh(V.to_matrix(cone)) > -1e-12).all()
assert (V.square(V.sqrt(cone)) - cone).abs().max() < 1e-10
O, _ = torch.linalg.qr(torch.randn(5, 5, dtype=torch.float64))  # Proposition 13.5(b)
rot = lambda w: torch.cat([w[..., :1], w[..., 1:] @ O.T], dim=-1)
assert (rot(V.product(xs, ys)) - V.product(rot(xs), rot(ys))).abs().max() < 1e-12

for m in range(1, 7):                                            # Proposition 13.5(c)
    r = jd.compare_with_star_algebra(jd.anticommuting_symmetries(m))
    print(f"m = {m}: dim V_m = {r['jordan_dim']}, dim of the self-adjoint part of the *-algebra = {r['star_dim']}")
    assert r['jordan_dim'] == m + 1 and r['equal'] == (m in (1, 3))
```

Para $m=1$, $s_1=Z$ genera las matrices diagonales, de dimensión $2=m+1$, y $V_1\cong\mathbb R\oplus\mathbb R$ no es
simple. Las dimensiones impresas son $3,4$ para $m=2$; $5,16$ para $m=4$; $7,64$ para $m=6$.

> **Significado físico.** Un estado de $V_m$ (una forma lineal positiva y normalizada) es, por autodualidad,
> $x\mapsto\operatorname{tr}(\rho\circ x)$ con $\rho=\frac12(1,r)\in Q$, $|r|\le1$: el espacio de estados es una bola
> de dimensión $m$. Para $m=3$ es la bola de Bloch de un qubit (proposición 13.5(a)). El cono de Lorentz es el cono de
> luz futuro de la relatividad especial en $1+m$ dimensiones, y la proposición 13.5(b) dice que los automorfismos de
> Jordan son las aplicaciones ortogonales del espacio que fijan el eje temporal.

Los conos $Q$ son los conjuntos de restricciones de la **programación cónica de segundo orden** (SOCP): optimización
lineal sobre productos de conos de Lorentz. Los métodos de punto interior para SOCP usan exactamente el teorema 13.4:
valores propios, $\det$ y $P(y)$ para el escalado (F. Alizadeh, D. Goldfarb, *Second-order cone programming*, Math.
Program. 95 (2003) 3–51).

## 13.3 Subálgebras de Jordan que no son *-subálgebras

Una subálgebra de Jordan $J\subseteq\operatorname{Herm}_n(\mathbb C)$ genera una $*$-álgebra compleja $\mathcal A$ con
$J\subseteq\mathcal A_{sa}$. La proposición 13.5(c) dio ejemplos con $J\ne\mathcal A_{sa}$; he aquí dos más.

**Proposición 13.6.** (a) $\operatorname{Herm}_n(\mathbb R)$ es una subálgebra de Jordan de
$\operatorname{Herm}_n(\mathbb C)$ de dimensión $n(n+1)/2$; la $*$-álgebra que genera es $M_n(\mathbb C)$, con parte
autoadjunta de dimensión $n^2$. (b) Para $x,y\in\operatorname{Herm}_2(\mathbb C)$ existe una unitaria $U$ tal que
$UxU^*$ y $UyU^*$ son simétricas reales. El álgebra de Jordan generada por $1,x,y$ es $\operatorname{span}\{1,x,y\}$; si
$x,y,1$ son linealmente independientes, es $U^*\operatorname{Herm}_2(\mathbb R)U\cong V_2$, de dimensión $3$, mientras
que la $*$-álgebra generada por $x,y$ es $M_2(\mathbb C)$, de dimensión $4$.

*Demostración.* (a) Lema 13.1. Para $i\ne j$, $E_{ij}=E_{ii}(E_{ij}+E_{ji})$ es un producto de matrices simétricas
reales, así que el álgebra generada contiene todas las unidades matriciales. (b) Escribamos $x=a1+u\cdot\sigma$,
$y=b1+v\cdot\sigma$ y elijamos un vector unitario $n\perp u,v$. Las matrices $n\cdot\sigma$ e $Y$ son hermíticas, de
traza nula y de cuadrado $1$, de modo que ambas tienen los valores propios $1,-1$, y $U(n\cdot\sigma)U^*=Y$ para alguna
unitaria $U$ (teorema 1.3). La conjugación por $U$ conserva $\operatorname{Tr}$ y $\circ$, así que lleva
$u\cdot\sigma$ a $(Ru)\cdot\sigma$ con $R$ lineal y
$Ru\cdot Rv=\frac12\operatorname{Tr}\big((Ru\cdot\sigma)(Rv\cdot\sigma)\big)=u\cdot v$: $R\in O(3)$, y $Rn=e_2$. Por
tanto, $Ru,Rv\perp e_2$, y $UxU^*=a1+(Ru)_1X+(Ru)_3Z$ es real, ya que $X,Z$ lo son; análogamente $UyU^*$. Además,
$x\circ y=(ab+u\cdot v)1+av\cdot\sigma+bu\cdot\sigma\in\operatorname{span}\{1,x,y\}$, y lo mismo $x^2,y^2$, de modo que
el subespacio generado es cerrado. Si tiene dimensión $3$, su conjugado por $U$ es todo
$\operatorname{Herm}_2(\mathbb R)=\operatorname{span}\{1,X,Z\}$, que es $\cong V_2$ por el lema 13.3 ($X,Z$
anticonmutan), y genera $M_2(\mathbb C)$ por (a). $\square$

```python
for n in range(2, 6):                                            # Proposition 13.6(a), a random real pair
    G = torch.randn(2, n, n, dtype=torch.float64)
    r = jd.compare_with_star_algebra(G + G.mT)
    assert (r['jordan_dim'], r['star_dim']) == (n * (n + 1) // 2, n * n)

H = torch.randn(3, 2, 2, dtype=dt)
H = H + H.mH                                                     # three random Hermitian 2 x 2 matrices
r2, r3 = jd.compare_with_star_algebra(H[:2]), jd.compare_with_star_algebra(H)
print("two generators:", r2['jordan_dim'], r2['star_dim'], "  three:", r3['jordan_dim'], r3['star_dim'])
assert (r2['jordan_dim'], r2['star_dim'], r3['jordan_dim']) == (3, 4, 4)

V3 = SpinFactor(3)
u, v = V3.from_matrix(H[0])[1:], V3.from_matrix(H[1])[1:]        # Bloch vectors
nrm = torch.linalg.cross(u, v)
nrm = nrm / nrm.norm()
_, W1 = torch.linalg.eigh(V3.to_matrix(torch.cat([torch.zeros(1, dtype=torch.float64), nrm])))
_, W2 = torch.linalg.eigh(jd.anticommuting_symmetries(3)[1])     # Y
Um = W2 @ W1.mH                                                  # U (n . sigma) U* = Y
assert all((Um @ A @ Um.mH).imag.abs().max() < 1e-12 for A in H[:2])
```

Tres matrices hermíticas $2\times2$ genéricas generan todo $\operatorname{Herm}_2(\mathbb C)$: dimensión $4$ en ambas
columnas. Por tanto, en un bloque $M_2$ de $\mathcal M$, dos observables nunca revelan la estructura compleja; un
tercero puede hacerlo.

## 13.4 La clasificación

Un álgebra de Jordan es **formalmente real** si $\sum_ix_i^2=0$ implica que todos los $x_i=0$; en dimensión finita
esto equivale a ser **euclídea**, es decir, a tener un producto escalar con
$\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$ [FK94]. $M_{sa}$ (con $\operatorname{Tr}(xy)$) y $V_m$ (con
$\operatorname{tr}(x\circ y)$) son euclídeas.

**Teorema 13.7 (Jordan–von Neumann–Wigner; enunciado sin demostración).** Toda álgebra de Jordan formalmente real de
dimensión finita es suma directa de álgebras simples, y toda álgebra simple es isomorfa exactamente a una de las
siguientes: $\mathbb R$; $\operatorname{Herm}_r(\mathbb K)$ para $\mathbb K=\mathbb R,\mathbb C,\mathbb H$ y $r\ge3$;
los factores de espín $V_m$, $m\ge2$; el **álgebra de Albert** $\operatorname{Herm}_3(\mathbb O)$ de las matrices
hermíticas $3\times3$ sobre los octoniones (dimensión $27$). En rango $2$,
$\operatorname{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}$ para $\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O$.
(P. Jordan, J. von Neumann, E. Wigner, *On an algebraic generalization of the quantum mechanical formalism*, Ann. of
Math. 35 (1934) 29–64; una demostración moderna está en [FK94, capítulo V].)

Todas ellas salvo el álgebra de Albert son especiales: $\operatorname{Herm}_r(\mathbb H)$ se sumerge en
$\operatorname{Herm}_{2r}(\mathbb C)$ escribiendo los cuaterniones como matrices complejas $2\times2$, y $V_m$ por el
lema 13.3. El álgebra de Albert es **excepcional**: no es isomorfa a ninguna subálgebra de Jordan de un álgebra
asociativa (A. A. Albert, *On a certain algebra of quantum mechanics*, Ann. of Math. 35 (1934) 65–73). Una razón que
puede enunciarse brevemente: C. M. Glennie encontró identidades polinómicas que valen en toda álgebra de Jordan especial
pero fallan en $\operatorname{Herm}_3(\mathbb O)$ (*Some identities valid in special Jordan algebras but not valid in
all Jordan algebras*, Pacific J. Math. 16 (1966) 47–59; véase McCrimmon, *A Taste of Jordan Algebras*). Los octoniones
no son asociativos, y $\operatorname{Herm}_r(\mathbb O)$ satisface la identidad de Jordan solo para $r\le3$, razón por
la cual la lista se detiene ahí. Así, $\operatorname{Herm}_3(\mathbb O)$ no puede representarse fielmente mediante
matrices reales o complejas (ni dentro de ninguna álgebra asociativa) con el producto $\frac12(xy+yx)$; el software
tiene que trabajar con entradas octoniónicas.

**Teorema 13.8 (Koecher–Vinberg; enunciado sin demostración).** Los interiores de los conos de cuadrados de las
álgebras de Jordan euclídeas son exactamente los **conos simétricos**: conos convexos abiertos en un espacio euclídeo
que son autoduales y homogéneos (sus automorfismos lineales actúan transitivamente). [FK94, capítulo III.]

El teorema 13.4(4) y el ejercicio 1 demuestran ambas propiedades para el cono de Lorentz. Así, los conos de la
programación semidefinida (productos de $\operatorname{Herm}_n(\mathbb R)^+$ o $\operatorname{Herm}_n(\mathbb C)^+$), de
SOCP y de la programación lineal ($\mathbb R_+^n$, el cono de $\mathbb R^n=\mathbb R\oplus\dots\oplus\mathbb R$) son
todos de esta forma, y los métodos de punto interior pueden escribirse de una vez para todas las álgebras de Jordan
euclídeas.

[FK94]: J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, Oxford University Press, 1994.

## 13.5 Aplicaciones de Jordan frente a *-aplicaciones

Un **$*$-homomorfismo de Jordan** entre álgebras de von Neumann es una aplicación $\mathbb C$-lineal $\Phi$ con
$\Phi(x^*)=\Phi(x)^*$ y $\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$; sobre las partes autoadjuntas es un homomorfismo de
Jordan, y todo homomorfismo de Jordan entre partes autoadjuntas se extiende a uno de ellos mediante
$\Phi(x+iy)=\Phi(x)+i\Phi(y)$.

**Proposición 13.9 (la traspuesta).** Sea $T(A)=\bigoplus_cA_c^T$ sobre $\mathcal M$. (a) $T$ es un
$*$-automorfismo de Jordan unital con $T(AB)=T(B)T(A)$; conserva los espectros, de modo que $T$ y $T^{-1}=T$ son
positivas. (b) Sobre $M_k$, $k\ge2$, la matriz de Choi $\sum_{ij}E_{ij}\otimes T(E_{ij})$ es el intercambio
$F(e_a\otimes e_b)=e_b\otimes e_a$, con valores propios $+1$ (multiplicidad $k(k+1)/2$) y $-1$ (multiplicidad
$k(k-1)/2$). Así, $T$ no es completamente positiva.

*Demostración.* (a) $(A^T)^*=(A^*)^T$, $(AB)^T=B^TA^T$, luego $(A\circ B)^T=B^T\circ A^T=A^T\circ B^T$;
$\det(A^T-\lambda)=\det(A-\lambda)$. (b) $\sum_{ij}E_{ij}e_a\otimes E_{ji}e_b=e_b\otimes e_a$. $F^2=1$, y $F$ vale $+1$
sobre los tensores simétricos y $-1$ sobre los antisimétricos, que generan todo el espacio. Por el teorema 6.2 (y la
proposición 6.1), $T$ no es CP. $\square$

Así, un automorfismo de Jordan no tiene por qué ser un $*$-automorfismo. Salvo este fenómeno, la positividad determina
la estructura de Jordan:

**Teorema 13.10 (Kadison; enunciado sin demostración).** Sea $\Phi$ una biyección lineal entre C\*-álgebras unitales
con $\Phi(1)=1$ tal que $\Phi$ y $\Phi^{-1}$ son positivas. Entonces $\Phi$ es un $*$-isomorfismo de Jordan. (R. V.
Kadison, *A generalized Schwarz inequality and algebraic invariants for operator algebras*, Ann. of Math. 56 (1952)
494–503, corolario 5. En *Isometries of operator algebras*, Ann. of Math. 54 (1951) 325–338, había demostrado que una
isometría lineal $\Phi$ de una C\*-álgebra unital sobre otra es $\Phi(1)$, una unitaria, por un $*$-isomorfismo de
Jordan.)

El recíproco es inmediato: un $*$-isomorfismo de Jordan lleva los cuadrados $h^2=h\circ h$ de elementos autoadjuntos
a cuadrados y, por tanto, el cono positivo sobre el cono positivo. Para los factores $M_n$ se conocen exactamente los
automorfismos de Jordan: todo $*$-automorfismo de Jordan de $M_n(\mathbb C)$ es $A\mapsto UAU^*$ o $A\mapsto UA^TU^*$
con $U$ unitaria (N. Jacobson, C. E. Rickart, *Jordan homomorphisms of rings*, Trans. Amer. Math. Soc. 69 (1950)
479–502; I. N. Herstein, *Jordan homomorphisms*, Trans. Amer. Math. Soc. 81 (1956) 331–341). Demostramos el caso
$n=2$.

**Proposición 13.11.** Todo automorfismo de Jordan de $\operatorname{Herm}_2(\mathbb C)$ es $A\mapsto UAU^*$ o
$A\mapsto UA^TU^*$ con $U\in U(2)$.

*Demostración.* Por la proposición 13.5(a,b), un automorfismo de Jordan actúa como
$a1+u\cdot\sigma\mapsto a1+(Ou)\cdot\sigma$ con $O\in O(3)$. Escribamos $R_U\in O(3)$ para la matriz de
$A\mapsto UAU^*$ en este sentido; $R_{UV}=R_UR_V$. (i) $\det R_U=1$:
$U=W\operatorname{diag}(e^{i\theta_1},e^{i\theta_2})W^*$ está unida a $1$ por el camino
$U_t=W\operatorname{diag}(e^{it\theta_1},e^{it\theta_2})W^*$, y $\det R_{U_t}\in\{\pm1\}$ depende continuamente de $t$.
(ii) Todo $R\in SO(3)$ es algún $R_U$: como en la proposición 13.6, existe $U_1$ con
$U_1ZU_1^*=(Re_3)\cdot\sigma$, es decir, $R_{U_1}e_3=Re_3$. Entonces $R_{U_1}^{-1}R\in SO(3)$ fija $e_3$, así que es
una rotación de cierto ángulo $\beta$ alrededor de $e_3$. Un cálculo directo da, para
$U_2=\operatorname{diag}(e^{-i\beta/2},e^{i\beta/2})$, $U_2XU_2^*=\cos\beta\,X+\sin\beta\,Y$ y
$U_2YU_2^*=-\sin\beta\,X+\cos\beta\,Y$, de modo que $R=R_{U_1U_2}$. (iii) La traspuesta fija $X,Z$ y lleva
$Y\mapsto-Y$: $O_T=\operatorname{diag}(1,-1,1)$, $\det O_T=-1$. Por tanto, $O=R_U$ si $\det O=1$, y $OO_T=R_U$, es
decir, el automorfismo es $A\mapsto UA^TU^*$, si $\det O=-1$. $\square$

**Sistemas compuestos.** La estructura de Jordan de un factor no puede distinguir $T$ de la identidad. Para dos
sistemas esta ambigüedad se vuelve visible:

**Proposición 13.12.** Sean $x,x'\in\operatorname{Herm}_n(\mathbb C)$, $y,y'\in\operatorname{Herm}_m(\mathbb C)$. En
$\operatorname{Herm}_{nm}(\mathbb C)$, que actúa sobre $\mathbb C^n\otimes\mathbb C^m$,
$$
(x\otimes y)\circ(x'\otimes y')=(x\circ x')\otimes(y\circ y')+\tfrac14[x,x']\otimes[y,y'].
$$
Sean $\theta=\mathrm{id}\otimes T$ (la **traspuesta parcial**) y $A\circ_\theta B=\theta(\theta(A)\circ\theta(B))$.
Entonces $\circ_\theta$ es un producto de Jordan sobre $\operatorname{Herm}_{nm}(\mathbb C)$, isomorfo a $\circ$ vía
$\theta$, que coincide con $\circ$ sobre todos los pares de observables locales $x\otimes1$, $1\otimes y$, y
$(x\otimes y)\circ_\theta(x'\otimes y')=(x\circ x')\otimes(y\circ y')-\frac14[x,x']\otimes[y,y']$. Su cono de cuadrados
es $\theta(\operatorname{Herm}_{nm}^+)$. Para $n=m=2$ y $\Omega=e_1\otimes e_1+e_2\otimes e_2$, el estado
$\omega=\frac12|\Omega\rangle\langle\Omega|\ge0$ no está en este cono: $\theta(\omega)=\frac12F$ tiene el valor propio
$-\frac12$.

*Demostración.* Con $j=x\circ x'$, $c=\frac12[x,x']$, $k=y\circ y'$, $d=\frac12[y,y']$: $xx'=j+c$, $x'x=j-c$, y
análogamente para $y$, de modo que
$\frac12(xx'\otimes yy'+x'x\otimes y'y)=\frac12\big((j+c)\otimes(k+d)+(j-c)\otimes(k-d)\big)=j\otimes k+c\otimes d$.
$\theta$ es una involución $\mathbb R$-lineal de $\operatorname{Herm}_{nm}$: lleva $x\otimes y$ a $x\otimes y^T$, y los
$x\otimes y$ generan $\operatorname{Herm}_{nm}$ sobre $\mathbb R$ (los productos $h_\alpha\otimes g_\beta$ de bases
reales de $\operatorname{Herm}_n,\operatorname{Herm}_m$ forman una base compleja de $M_{nm}$, y una matriz hermítica
tiene coeficientes reales en ella). Así, $\circ_\theta$ es un producto de Jordan y $\theta$ un isomorfismo. Como
$(y^T\circ y'^T)^T=y\circ y'$ y $[y^T,y'^T]^T=-[y,y']$, aplicar $\theta$ a la fórmula para
$(x\otimes y^T)\circ(x'\otimes y'^T)$ da la fórmula para $\circ_\theta$. Sobre observables locales el término del
conmutador se anula ($[1,1]=0$, $[x,1]=0$). Los cuadrados son $A\circ_\theta A=\theta(\theta(A)^2)$, y $\omega$ está
en $\theta(\operatorname{Herm}_{nm}^+)$ si y solo si $\theta(\omega)\ge0$, ya que $\theta^2=\mathrm{id}$. Por último,
$\theta(\omega)=\frac12\sum_{ij}E_{ij}\otimes E_{ji}=\frac12F$ (proposición 13.9(b); es el cálculo de la demostración
de la proposición 6.1). $\square$

Los dos productos son álgebras de Jordan igual de buenas, coinciden en todo lo que puede medirse sobre un subsistema
cada vez, y discrepan en qué operadores conjuntos son positivos, es decir, en cuáles son los estados del sistema
compuesto. Lo que los distingue es el signo del término del conmutador, y este signo lo fija la $*$-álgebra
$M_n\otimes M_m$ y no las álgebras de Jordan de los factores. ($\circ_\theta$ es el producto de Jordan del producto
asociativo $A\star B=\theta(\theta(A)\theta(B))$, para el cual $(x\otimes y)\star(x'\otimes y')=xx'\otimes y'y$: el
álgebra $M_n\otimes M_m^{\rm op}$.) El cono $\theta(\operatorname{Herm}_{nm}^+)$ es lo que da el **criterio PPT**
(A. Peres, Phys. Rev. Lett. 77 (1996) 1413): un estado separable $\rho=\sum_kp_ka_k\otimes b_k$, $p_k\ge0$,
$a_k,b_k\ge0$, está en ambos conos, ya que $\theta(\rho)=\sum_kp_ka_k\otimes b_k^T$ con $b_k^T\ge0$ (proposición
13.9(a)) y $a\otimes b=(a^{1/2}\otimes b^{1/2})^2\ge0$. Así, $\omega$, que no está en
$\theta(\operatorname{Herm}_4^+)$, está entrelazado. Para resultados sobre qué compuestos de álgebras de Jordan existen,
véanse H. Hanche-Olsen, *On the structure and tensor products of JC-algebras*, Canad. J. Math. 35 (1983) 1059–1074, y
H. Barnum, M. Graydon, A. Wilce, *Composites and categories of Euclidean Jordan algebras*, Quantum 4 (2020) 359.

**Dinámica.** Cabría esperar que la evolución temporal obligue a usar el producto asociativo. No es así:

**Proposición 13.13.** Sean $H\in\operatorname{Herm}_n(\mathbb C)$ y $\delta(z)=i[H,z]$. (a) $\delta$ lleva
$\operatorname{Herm}_n$ en sí mismo y es una derivación de $\circ$; $e^{t\delta}(z)=e^{itH}ze^{-itH}$ es un
automorfismo de Jordan. (b) $[L(x),L(y)]z=\frac14[[x,y],z]$ para $x,y,z\in\operatorname{Herm}_n$. (c)
$\delta=4\sum_k[L(x_k),L(y_k)]$ para ciertos $x_k,y_k\in\operatorname{Herm}_n$: la dinámica hamiltoniana está
generada por el producto de Jordan por sí solo.

*Demostración.* (a) $(i[H,z])^*=-i[z,H]=i[H,z]$; $[H,\cdot]$ es una derivación del producto asociativo y, por tanto,
de $\circ$. (b) $4x\circ(y\circ z)=xyz+xzy+yzx+zyx$; réstese lo mismo con $x,y$ intercambiados. (c) Por (b) basta
escribir $iH-\frac{\operatorname{Tr}H}n i1$ (que tiene los mismos conmutadores que $iH$) como suma de $[x_k,y_k]$.
Para $i\ne j$ y $a=E_{ij}+E_{ji}$, $b=i(E_{ij}-E_{ji})$: $[a,E_{ii}]=E_{ji}-E_{ij}$,
$[b,E_{ii}]=-i(E_{ij}+E_{ji})$ y $[a,b]=2i(E_{jj}-E_{ii})$; estas generan sobre $\mathbb R$ las matrices
antihermíticas de traza nula. $\square$

Lo que el álgebra de Jordan no contiene es la **correspondencia** $H\mapsto i[H,\cdot]$ entre un observable y la
dinámica que genera. En $\operatorname{Herm}_n(\mathbb R)$, por ejemplo, las derivaciones
$[L(x),L(y)]=\frac14[[x,y],\cdot]$ provienen de matrices $[x,y]$ reales antisimétricas, que no son $i$ por una matriz
simétrica real. E. M. Alfsen y F. W. Shultz demostraron que una JBW-álgebra (una versión de dimensión infinita de un
álgebra de Jordan formalmente real; en dimensión finita ambas nociones coinciden) es la parte autoadjunta de un álgebra
de von Neumann si y solo si admite tal "correspondencia dinámica", y caracterizaron geométricamente los espacios de
estados de las C\*-álgebras (*Geometry of State Spaces of Operator Algebras*, Birkhäuser, 2003; enunciado sin
demostración). Una correspondencia no es única: sobre $\operatorname{Herm}_n(\mathbb C)$, tanto $H\mapsto i[H,\cdot]$
como $H\mapsto-i[H,\cdot]$ sirven; la traspuesta las intercambia, ya que $T(i[H,T(z)])=-i[H^T,z]$, y pertenecen a $M_n$
y a su álgebra opuesta: los dos signos de la proposición 13.12.

```python
x, y = sa(), sa()
T = jd.transpose_map
assert close(T(jp(x, y)), jp(T(x), T(y))) and close(T(x @ y), T(y) @ T(x))   # Proposition 13.9(a)
assert not close(T(x @ y), T(x) @ T(y))

k = 3
Mk = TypeIAlgebra([k], [k], precision='double')
units = torch.eye(k * k, dtype=dt).reshape(k * k, 1, k, k)                     # E_ij, batch index i*k + j
TE = T(Mk.operator(units)).matrix[:, 0]
choi = TE.reshape(k, k, k, k).permute(0, 2, 1, 3).reshape(k * k, k * k)       # block (i, j) = T(E_ij)
swap = torch.eye(k * k, dtype=dt).reshape(k, k, k, k).permute(0, 1, 3, 2).reshape(k * k, k * k)
w = torch.linalg.eigvalsh(choi)
assert torch.equal(choi, swap) and (w < 0).sum() == k * (k - 1) // 2          # Proposition 13.9(b)

sig = jd.anticommuting_symmetries(3)                                           # X, Y, Z
O_T = torch.stack([V3.from_matrix(s.T)[1:] for s in sig], dim=1)               # transpose on Bloch vectors
assert torch.allclose(O_T, torch.diag(torch.tensor([1.0, -1.0, 1.0], dtype=torch.float64)))

def ptranspose(R, n, m):                                                       # theta = id (x) T
    return R.reshape(n, m, n, m).permute(0, 3, 2, 1).reshape(n * m, n * m)
def herm(n):
    A = torch.randn(n, n, dtype=dt)
    return A + A.mH
jm = lambda A, B: 0.5 * (A @ B + B @ A)
cm = lambda A, B: A @ B - B @ A
a, a2, b, b2 = herm(2), herm(2), herm(3), herm(3)
local = torch.kron(jm(a, a2), jm(b, b2))
assert close(jm(torch.kron(a, b), torch.kron(a2, b2)), local + 0.25 * torch.kron(cm(a, a2), cm(b, b2)))
theta = lambda R: ptranspose(R, 2, 3)
assert close(theta(jm(theta(torch.kron(a, b)), theta(torch.kron(a2, b2)))),
             local - 0.25 * torch.kron(cm(a, a2), cm(b, b2)))                  # Proposition 13.12
Om = torch.tensor([1.0, 0, 0, 1.0], dtype=dt) / math.sqrt(2)
omega = torch.outer(Om, Om.conj())                                             # Bell state
print("spectrum of the partial transpose of the Bell state:", torch.linalg.eigvalsh(ptranspose(omega, 2, 2)))
assert abs(torch.linalg.eigvalsh(ptranspose(omega, 2, 2))[0] + 0.5) < 1e-12

Hm, z1, z2 = herm(3), herm(3), herm(3)
delta = lambda A: 1j * cm(Hm, A)
assert close(delta(jm(z1, z2)), jm(delta(z1), z2) + jm(z1, delta(z2)))         # Proposition 13.13(a)
assert close(jm(z1, jm(z2, Hm)) - jm(z2, jm(z1, Hm)), 0.25 * cm(cm(z1, z2), Hm))   # (b)
```

> **Significado físico.** Para una partícula sin espín, la traspuesta (la conjugación compleja en una base real) es la
> **inversión temporal**. Es una simetría de los observables de un único sistema, pero es antiunitaria y, aplicada
> solo a un subsistema, no es una operación física: la proposición 13.12 muestra que convierte el estado de Bell en un
> operador con un valor propio negativo. La mecánica cuántica usa $*$-álgebras complejas porque los sistemas compuestos
> necesitan el conmutador con su signo, no solo el producto simetrizado.

## 13.6 De vuelta a las desigualdades de trazas

La proposición 11.9 (Exp. 3) afirma $\operatorname{Tr}|YXY|\le\operatorname{Tr}(Y|X|Y)$ para $X=X^*$ e $Y\ge0$ (aquí
$X,Y$ son operadores de $\mathcal M$, no matrices de Pauli). Su demostración nunca usa $Y\ge0$: $YX_\pm Y=T^*T$ con
$T=X_\pm^{1/2}Y$ solo necesita $Y=Y^*$, y el lema 11.8 se aplica a $YX_+Y$ e $YX_-Y$. Así que vale para todo $Y=Y^*$
y, en lenguaje de Jordan, se lee
$$
\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|,
$$
lo que tiene sentido en toda álgebra de Jordan euclídea. La demostramos en $V_m$ con un argumento que usa solo tres
propiedades.

**Proposición 13.14 (Exp. 3 en un factor de espín).** Para $x,y\in V_m$,
$\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|$, con igualdad si $x\in Q$ o $-x\in Q$.

*Demostración.* Sean $x_\pm=f_\pm(x)$ con $f_+(t)=\max(t,0)$, $f_-(t)=\max(-t,0)$ (teorema 13.4(3)): $x_\pm\in Q$,
$x=x_+-x_-$ y $|x|=x_++x_-$. Pongamos $A=P(y)x_+$ y $B=P(y)x_-$; están en $Q$ por el teorema 13.4(5), y
$P(y)x=A-B$. Sean $P(y)x=\mu_1c_1+\mu_2c_2$ una descomposición espectral (teorema 13.4(2); $c_1=c_w$, $c_2=c_{-w}$, de
modo que $c_1+c_2=1$ y $c_1\circ c_2=0$), $c$ la suma de los $c_i$ con $\mu_i>0$, y
$s=2c-1=\sum_i\operatorname{sgn}'(\mu_i)c_i$ con $\operatorname{sgn}'(t)=1$ para $t>0$ y $-1$ en otro caso. Entonces
$s\circ P(y)x=\sum_i|\mu_i|c_i=|P(y)x|$ (por el teorema 13.4(2,3), esto no depende de la elección de $w$ cuando
$\mu_1=\mu_2$), y $c,1-c\in Q$ (los idempotentes son cuadrados). Por la mitad fácil de la autodualidad (teorema
13.4(4): $\operatorname{tr}(a\circ b)\ge0$ para $a,b\in Q$), aplicada a $\operatorname{tr}(A\circ(1-c))$ y a
$\operatorname{tr}(B\circ c)$,
$$
\operatorname{tr}|P(y)x|=\operatorname{tr}(A\circ c)-\operatorname{tr}(A\circ(1-c))-\operatorname{tr}(B\circ c)+\operatorname{tr}(B\circ(1-c))
\le\operatorname{tr}A+\operatorname{tr}B=\operatorname{tr}P(y)|x|.
$$
Si $x\in Q$, entonces $x_-=0$, $B=0$ y $|P(y)x|=A$, de modo que ambos lados valen $\operatorname{tr}A$; análogamente
para $-x\in Q$. $\square$

La demostración usa solo: (i) $P(y)Q\subseteq Q$ para todo $y$; (ii) el teorema espectral: todo $z$ es
$\sum_i\mu_ic_i$ con idempotentes ortogonales $c_i$ de suma $1$, y $f(z)=\sum_if(\mu_i)c_i$ (de modo que $x_\pm$ y los
idempotentes son cuadrados); (iii) $\operatorname{tr}$ es lineal,
$\operatorname{tr}(\sum_i\mu_ic_i)=\sum_i\mu_i\operatorname{tr}c_i$ con $\operatorname{tr}c_i\ge0$, y
$\operatorname{tr}(a\circ b)\ge0$ para $a,b\in Q$. En un álgebra de Jordan euclídea de rango $r$, (ii) vale con $r$
idempotentes primitivos $c_i$ (un marco de Jordan), la traza es $\operatorname{tr}z=\sum_i\mu_i$, la forma
$\operatorname{tr}(x\circ y)$ es asociativa, $\operatorname{tr}((x\circ y)\circ z)=\operatorname{tr}(x\circ(y\circ z))$, y
(i) se cumple [FK94, capítulos II–III]; no demostramos estos hechos. La última parte de (iii) se sigue entonces de
(i): desarrollando $P(a^{1/2})=2L(a^{1/2})^2-L(a)$ y usando la asociatividad,
$\operatorname{tr}(a\circ b)=\operatorname{tr}P(a^{1/2})b$, la traza de un elemento de $Q$, luego $\ge0$. Así, la
desigualdad vale en toda álgebra de Jordan euclídea, incluidas $\operatorname{Herm}_r(\mathbb H)$ y el álgebra de
Albert, que no son isomorfas como álgebras de Jordan a partes autoadjuntas de $*$-álgebras complejas. En $M_{sa}$, el
mismo argumento, con $\operatorname{Tr}$ y $P(y)x=yxy$, vuelve a demostrar la proposición 11.9. En $V_m$, la
desigualdad se sigue también de la proposición 11.9 a través de $\phi$, ya que $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ y
$\phi(|x|)=|\phi(x)|$; el código comprueba ambas cosas:

```python
for m in (2, 4, 7):
    Vm = SpinFactor(m)
    xs, ys = Vm.random(20_000), Vm.random(20_000)
    lhs = Vm.trace(Vm.apply_function(Vm.quadratic_representation(ys, xs), torch.abs))
    rhs = Vm.trace(Vm.quadratic_representation(ys, Vm.apply_function(xs, torch.abs)))
    definite = Vm.determinant(xs) >= 0                                 # x or -x in the cone
    print(f"m = {m}: min gap {(rhs - lhs).min().item():.1e}, max gap {(rhs - lhs).max().item():.1f}, "
          f"definite x: {definite.double().mean().item():.3f}")
    assert (rhs - lhs >= -1e-12 * rhs).all()                            # Proposition 13.14
    assert ((rhs - lhs)[definite].abs() <= 1e-12 * rhs[definite]).all()
    Xs, Ys = Vm.to_matrix(xs[:200]), Vm.to_matrix(ys[:200])
    tn = torch.linalg.svdvals(Ys @ Xs @ Ys).sum(-1)                     # Tr |phi(y) phi(x) phi(y)|
    assert (tn - Vm.representation_dim / 2 * lhs[:200]).abs().max() < 1e-10 * tn.max()
```

La diferencia $\operatorname{tr}P(y)|x|-\operatorname{tr}|P(y)x|$ se anula para $x$ definido y es positiva para un $x$
indefinido típico (los máximos impresos son del orden de $10^1$–$10^2$ para estas muestras gaussianas).

## 13.7 Con qué comparar

Hasta donde sabemos, no existe ninguna biblioteca numérica de uso extendido para álgebras de Jordan; las operaciones
de este capítulo suelen escribirse a mano en términos de matrices. El álgebra computacional cubre el lado algebraico:
`JordanAlgebra` de SageMath construye álgebras de Jordan especiales, álgebras de Jordan de formas bilineales simétricas
(que incluyen los factores de espín) y, a partir de un álgebra de octoniones, el álgebra de Jordan excepcional de
dimensión 27; Magma también tiene constructores de álgebras de Jordan, incluida la excepcional. Son herramientas
algebraicas exactas, más que numéricas y en lotes. Los conos de Lorentz son estándar en optimización: los solvers de
SOCP (ECOS, SCS, Clarabel, MOSEK, accesibles mediante la restricción `SOC` de cvxpy) usan la descomposición espectral de
$V_m$ (para proyectar sobre el cono) o su representación cuadrática (para el escalado de Nesterov–Todd de los métodos de
punto interior), y los solvers semidefinidos hacen lo mismo para $\operatorname{Herm}_n$. La biblioteca añade productos
de Jordan, la representación cuadrática y la traspuesta como operadores en lotes y diferenciables sobre
$\bigoplus_cM_{k_c}$; factores de espín con su cálculo funcional y su representación matricial; y un cálculo
aleatorizado del álgebra de Jordan generada por unas matrices dadas, comparada con la $*$-álgebra que generan.

## Ejercicios

1. (**Homogeneidad del cono de Lorentz.**) Sea $y\in V_m$ con $\det y\ne0$. Demuestre que $y^{-1}=(a,-u)/\det y$
   satisface $y\circ y^{-1}=1$, que $P(y)P(y^{-1})=\mathrm{id}$, que $P(y)$ lleva el interior de $Q$ sobre sí mismo, y
   que $P(x^{1/2})1=x$ para $x$ en el interior. Deduzca que el interior de $Q$ es homogéneo.
2. (**Estados isótropos.**) Para $p\in[0,1]$ sea $\rho_p=p\,\omega+(1-p)\frac14$ sobre $\mathbb C^2\otimes\mathbb C^2$
   (proposición 13.12). Calcule el espectro de $\theta(\rho_p)$ y demuestre que $\rho_p$ está entrelazado para
   $p>\frac13$.
3. (**Productos de puntos de la frontera.**) Sean $|w|=|w'|=1$ en $\mathbb R^m$. Demuestre que $(1,w)\circ(1,w')\in Q$
   si y solo si $w'=w$ o $w'=-w$: el producto de Jordan de dos rayos de la frontera del cono de Lorentz casi nunca está
   en el cono.
4. (**Mecánica cuántica real.**) Demuestre que, para $x,y\in\operatorname{Herm}_n(\mathbb R)$, la derivación
   $[L(x),L(y)]$ de $\operatorname{Herm}_n(\mathbb R)$ es $z\mapsto[K,z]$ con $K$ real antisimétrica, y que para
   $n\ge2$ tal aplicación nunca es $z\mapsto i[H,z]$ con $H\in\operatorname{Herm}_n(\mathbb R)$, salvo que sea $0$.

<details>
<summary>Soluciones</summary>

1. $(a,u)\circ(a,-u)=(a^2-|u|^2,0)$, lo que da $y\circ y^{-1}=1$. $\phi(y)$ y $\phi(y^{-1})$ conmutan (ambas están en
   $\operatorname{span}\{1,u\cdot s\}$), así que $\phi(y)\phi(y^{-1})=\phi(y)\circ\phi(y^{-1})=1$ y, por el teorema
   13.4(5), $\phi(P(y)P(y^{-1})x)=\phi(x)$; $\phi$ es inyectiva. $P(y)$ lleva $Q$ en $Q$, y lo mismo hace su inversa
   $P(y^{-1})$, luego $P(y)Q=Q$; al ser un homeomorfismo lineal, lleva el interior sobre el interior. Por último,
   $P(y)1=2y\circ y-y^2=y^2$, así que $P(x^{1/2})1=x$, y $P(x^{1/2})$ es invertible cuando $x$ es interior
   ($\det x^{1/2}=\sqrt{\lambda_+\lambda_-}>0$). Dos puntos interiores cualesquiera están unidos por
   $P(x'^{1/2})P(x^{1/2})^{-1}$.

```python
V4 = SpinFactor(4)
xi = V4.random(100, cone=True)                                   # interior points (almost surely)
r = V4.sqrt(xi)
assert (V4.quadratic_representation(r, V4.unit((100,))) - xi).abs().max() < 1e-10
yv = V4.random(100)
back = V4.quadratic_representation(V4.inverse(yv), V4.quadratic_representation(yv, xi))
assert (back - xi).abs().max() < 1e-8 * xi.abs().max()
assert V4.in_cone(V4.quadratic_representation(yv, xi)).all()
```

2. $\theta(\rho_p)=\frac p2F+\frac{1-p}4$ por la proposición 13.12. $F$ tiene valores propios $1$ (tres veces) y $-1$
   (una vez), de modo que el espectro es $\frac{1+p}4$ (tres veces) y $\frac{1-3p}4$, que es negativo si y solo si
   $p>\frac13$; entonces $\rho_p$ no es separable por el criterio PPT. (Para sistemas $2\times2$ también vale el
   recíproco: M., P. y R. Horodecki, Phys. Lett. A 223 (1996) 1; no lo demostramos.)

```python
for pp in (0.2, 1 / 3, 0.5):
    rho_p = pp * omega + (1 - pp) * torch.eye(4, dtype=dt) / 4
    lam_min = torch.linalg.eigvalsh(ptranspose(rho_p, 2, 2))[0].item()
    assert abs(lam_min - (1 - 3 * pp) / 4) < 1e-12
```

3. $(1,w)\circ(1,w')=(1+t,\,w+w')$ con $t=w\cdot w'\in[-1,1]$ y $|w+w'|=\sqrt{2+2t}$. Está en $Q$ si y solo si
   $1+t\ge\sqrt{2(1+t)}$, es decir, si y solo si $t=-1$ (el producto es $0$) o $1+t\ge2$, es decir, $t=1$.

```python
V3c = SpinFactor(3)
wv = torch.nn.functional.normalize(torch.randn(2, 1000, 3, dtype=torch.float64), dim=-1)
prod = V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[1]))
assert not V3c.in_cone(prod).any()
assert V3c.in_cone(V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[0]))).all()
```

4. Por la proposición 13.13(b), $[L(x),L(y)]z=[K,z]$ con $K=\frac14[x,y]$, y $K^T=\frac14[y^T,x^T]=-K$ es real. Para
   $H$ y $z$ reales, $i[H,z]$ es puramente imaginaria, mientras que $[K,z]$ es real, así que solo coinciden si ambas se
   anulan. Así, en la mecánica cuántica real los generadores de la dinámica, $K$, no son observables.

</details>

**Cómo seguir.** Este es el último capítulo del curso. La referencia estándar para la teoría en dimensión finita es
[FK94]; K. McCrimmon, *A Taste of Jordan Algebras* (Springer 2004), cubre el álgebra y la historia, y H. Hanche-Olsen y
E. Størmer, *Jordan Operator Algebras* (Pitman 1984), y Alfsen–Shultz (2003), la teoría en dimensión infinita y su
relación con las C\*-álgebras. Para el papel de los sistemas compuestos en las reconstrucciones de la teoría cuántica a
partir de axiomas operacionales, véanse los artículos de Barnum y Wilce citados más arriba. Para desigualdades de trazas
más allá del capítulo 11, E. Carlen, *Trace inequalities and quantum entropy* (Contemp. Math. 529, 2010). La
representación cuadrática y el cono de Lorentz reaparecen en la optimización cónica: Alizadeh–Goldfarb (2003) y los
métodos de punto interior jordanianos que allí se revisan.
