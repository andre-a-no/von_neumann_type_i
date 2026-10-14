[English](../13_jordan_algebras.md) · **Español** · [Français](../fr/13_jordan_algebras.md) · [Deutsch](../de/13_jordan_algebras.md) · [中文](../zh/13_jordan_algebras.md) · [日本語](../ja/13_jordan_algebras.md) · [Русский](../ru/13_jordan_algebras.md)

# 13. Álgebras de Jordan: por qué la mecánica cuántica usa *-álgebras

El producto $xy$ de dos observables no es un observable salvo que $x$ e $y$ conmuten, pero el producto simetrizado
$x\circ y=\frac12(xy+yx)$ siempre lo es. En 1933 P. Jordan propuso tomar los observables con este producto como la
estructura básica de la mecánica cuántica, y Jordan, von Neumann y Wigner clasificaron las álgebras de dimensión finita
que aparecen. Este capítulo desarrolla la teoría en dimensión finita, con demostraciones completas aquí o, para las
tres largas (la clasificación, el teorema de Koecher–Vinberg, el teorema de Kadison y los automorfismos de Jordan de
$\bigoplus_cM_{k_c}$), en los apéndices E y F, y se pregunta qué se le escapa a la estructura de Jordan. La respuesta es
el signo del conmutador: el álgebra de Jordan de un sistema determina el conmutador solo salvo signo (sección 13.5),
pero los sistemas compuestos necesitan el signo: decide qué operadores conjuntos son positivos, y el test PPT
construido a partir de él detecta el entrelazamiento (proposición 13.12).

**Aprenderá** el producto de Jordan y la representación cuadrática $P(y)x=yxy$, que, a diferencia de $\circ$,
conserva la positividad; los factores de espín, su teoría espectral y el cono de Lorentz; subálgebras de Jordan de
matrices hermíticas que no son partes autoadjuntas de $*$-álgebras; la clasificación de Jordan–von Neumann–Wigner; por
qué la traspuesta es un automorfismo de Jordan y no un $*$-automorfismo, y qué hace sobre un sistema compuesto; y una
desigualdad de trazas del capítulo 11 escrita en lenguaje de Jordan y demostrada en un factor de espín.

En todo el capítulo, $\operatorname{Herm}_n(\mathbb K)$ es el espacio vectorial real de las matrices hermíticas
$n\times n$ sobre $\mathbb K=\mathbb R$ o $\mathbb C$ (la sección 13.4 usa también los cuaterniones $\mathbb H$ y los
octoniones $\mathbb O$, construidos en la sección E.6, y $\operatorname{Herm}_n$ sobre ellos, sección E.7), $M_{sa}$ es
la parte autoadjunta de
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

Escribiendo $P(y)$ como $2L(y)^2-L(y^2)$, cada identidad de (1) se expresa solo con $\circ$. La demostración de (1)
usa solo la asociatividad del producto matricial, así que en toda álgebra de Jordan especial $J\subseteq\mathcal A$
($\mathcal A$ asociativa) da $P(y)x=yxy$ y la fórmula fundamental $P(P(y)x)=P(y)P(x)P(y)$ sobre $J$. (Esto cubre todas
las álgebras en las que el curso la usa; la fórmula fundamental no se usa para álgebras no especiales.)

Entre los productos $\alpha xy+\beta yx$ ($\alpha,\beta\in\mathbb R$) sobre $\operatorname{Herm}_n(\mathbb C)$, $n\ge2$,
$\frac12(xy+yx)$ es el único que es simétrico en $x,y$ e igual a $x^2$ para $y=x$: la simetría significa
$(\alpha-\beta)(xy-yx)=0$ para todos los $x,y$, y $XZ\ne ZX$ (en el bloque $2\times2$ superior izquierdo) obliga a
$\alpha=\beta$; entonces $(\alpha+\beta)1=1$ para $x=y=1$ da $\alpha=\beta=\frac12$. El módulo `torch_vn_algebra.jordan`
implementa estas aplicaciones sobre operadores de un `TypeIAlgebra`, de forma perezosa y en lotes sobre todos los
sectores:

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

> **Significado físico** (interpretación física, no un teorema de este curso). $x\circ y$ es el producto simetrizado
> de observables (la única combinación de este tipo, como se ha visto arriba), y es de nuevo un observable. La
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
Jordan–Wigner $c_j$ de la sección 9.4, con el sitio $j-1$ del capítulo 9 como $j$-ésimo factor tensorial: allí
$s^a=\frac12\sigma_a$, $n=s^z+\frac12$ y $s^\pm=s^x\pm is^y$, de modo que $(-1)^n=1-2n=-Z$ y $s^-=\frac12(X-iY)$, y
por tanto $c_j=(-1)^{j-1}Z^{\otimes(j-1)}\otimes s^-\otimes1$. Como $s^-+s^+=X$ e $i(s^--s^+)=i(-iY)=Y$,
$$
s_{2j-1}=(-1)^{j-1}(c_j+c_j^*),\qquad s_{2j}=(-1)^{j-1}\,i(c_j-c_j^*).
$$

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

pq = 3                                                           # Majorana operators of Section 9.4, p = 3 sites
Xp, Yp, Zp = jd.anticommuting_symmetries(3)
sminus, I2 = 0.5 * (Xp - 1j * Yp), torch.eye(2, dtype=dt)
parity = I2 - 2 * (0.5 * Zp + 0.5 * I2)                          # (-1)^n with n = s^z + 1/2
def kron_all(fs):
    out = torch.ones(1, 1, dtype=dt)
    for f_ in fs:
        out = torch.kron(out, f_)
    return out
s6 = jd.anticommuting_symmetries(2 * pq)
for j in range(1, pq + 1):
    c = kron_all([parity] * (j - 1) + [sminus] + [I2] * (pq - j))
    sign = (-1) ** (j - 1)
    assert close(s6[2 * j - 2], sign * (c + c.mH)) and close(s6[2 * j - 1], sign * 1j * (c - c.mH))
```

Para $m=1$, $s_1=Z$ genera las matrices diagonales, de dimensión $2=m+1$, y $V_1\cong\mathbb R\oplus\mathbb R$ no es
simple. Las dimensiones impresas son $3,4$ para $m=2$; $5,16$ para $m=4$; $7,64$ para $m=6$.

**Estados de $V_m$.** Un **estado** de $V_m$ es una forma lineal $\omega$ con $\omega(Q)\ge0$ y $\omega(1)=1$. Los
estados son exactamente las aplicaciones $x\mapsto\operatorname{tr}(\rho\circ x)$ con $\rho=\frac12(1,r)$, $|r|\le1$;
así, el espacio de estados es una bola de dimensión $m$. *Demostración.*
$\langle x,y\rangle=\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$ es un producto escalar sobre $V_m$, de modo que toda
$\omega$ lineal es $\operatorname{tr}(\rho\circ\cdot)$ para exactamente un $\rho=(b,v)$ (representación de Riesz en
dimensión finita: en una base ortonormal, $\rho=\sum_k\omega(e_k)e_k$). Por la autodualidad del teorema 13.4(4),
$\omega(Q)\ge0$ si y solo si $\rho\in Q$, es decir, $b\ge|v|$; y $\omega(1)=\operatorname{tr}\rho=2b=1$. Así,
$\rho=\frac12(1,r)$ con $r=2v$, $|r|\le1$; recíprocamente, todo $\rho$ de este tipo da un estado. $\square$

> **Significado físico** (interpretación física, no un teorema de este curso). Para $m=3$ el espacio de estados de
> $V_m$ es la bola de Bloch de un qubit (proposición 13.5(a)). El cono de Lorentz es el cono de
> luz futuro de la relatividad especial en $1+m$ dimensiones, y la proposición 13.5(b) dice que los automorfismos de
> Jordan son las aplicaciones ortogonales del espacio que fijan el eje temporal.

Los conos $Q$ son los conjuntos de restricciones de la **programación cónica de segundo orden** (SOCP): optimización
lineal sobre productos de conos de Lorentz. Las magnitudes del teorema 13.4 (valores propios, $\det$ y $P(y)$) son
aquellas en las que suelen formularse los métodos de punto interior para SOCP; estos algoritmos no se tratan en este
curso (lectura adicional: F. Alizadeh, D. Goldfarb, *Second-order cone programming*, Math. Program. 95 (2003) 3–51).

## 13.3 Subálgebras de Jordan que no son *-subálgebras

Una subálgebra de Jordan $J\subseteq\operatorname{Herm}_n(\mathbb C)$ genera una $*$-álgebra compleja $\mathcal A$ con
$J\subseteq\mathcal A_{sa}$. La proposición 13.5(c) dio ejemplos con $J\ne\mathcal A_{sa}$; he aquí dos más.

**Proposición 13.6.** (a) $\operatorname{Herm}_n(\mathbb R)$ es una subálgebra de Jordan de
$\operatorname{Herm}_n(\mathbb C)$ de dimensión $n(n+1)/2$; la $*$-álgebra que genera es $M_n(\mathbb C)$, con parte
autoadjunta de dimensión $n^2$. (b) Para $x,y\in\operatorname{Herm}_2(\mathbb C)$ existe una unitaria $U$ tal que
$UxU^*$ y $UyU^*$ son simétricas reales. El álgebra de Jordan generada por $1,x,y$ es $\operatorname{span}\{1,x,y\}$; si
$x,y,1$ son linealmente independientes, es $U^*\operatorname{Herm}_2(\mathbb R)U\cong V_2$, de dimensión $3$, mientras
que la $*$-álgebra generada por $1,x,y$ es $M_2(\mathbb C)$, de dimensión $4$.

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
$\operatorname{Herm}_2(\mathbb R)=\operatorname{span}\{1,X,Z\}$, que es $\cong V_2$ vía $(a,u)\mapsto a1+u_1X+u_2Z$
(la demostración del lema 13.3 usa solo que las $s_i$ son simetrías hermíticas que anticonmutan, y $X,Z$ lo son), y
genera $M_2(\mathbb C)$ por (a). $\square$

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

Un álgebra de Jordan de dimensión finita es **formalmente real** si $\sum_ix_i^2=0$ implica que todos los $x_i=0$, y
**euclídea** si tiene una unidad $1$ y un producto escalar con $\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$.
Ambas nociones coinciden (teorema E.5(4); la dirección «euclídea $\Rightarrow$ formalmente real» ocupa una línea:
$0=\sum_i\langle x_i^2,1\rangle=\sum_i\langle x_i,x_i\rangle$). $M_{sa}$ con $\langle x,y\rangle=\operatorname{Tr}(xy)$
es euclídea: $\operatorname{Tr}(x^2)=\sum_{ij}|x_{ij}|^2$ y, por ciclicidad,
$\operatorname{Tr}((x\circ y)z)=\frac12\operatorname{Tr}(xyz+yxz)=\operatorname{Tr}(y(x\circ z))$. También lo es $V_m$
con $\langle x,y\rangle=\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$: para $z=(c,w)$,
$\operatorname{tr}((x\circ y)\circ z)=2(abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w)$ es simétrico en $x,y,z$.

**Teorema 13.7 (Jordan–von Neumann–Wigner).** Toda álgebra de Jordan formalmente real de dimensión finita es suma
directa de ideales simples, y toda álgebra simple es isomorfa exactamente a una de las siguientes: $\mathbb R$;
$\operatorname{Herm}_r(\mathbb K)$ para $\mathbb K=\mathbb R,\mathbb C,\mathbb H$ y $r\ge3$; los factores de espín
$V_m$, $m\ge2$; el **álgebra de Albert** $\operatorname{Herm}_3(\mathbb O)$ de las matrices hermíticas $3\times3$ sobre
los octoniones (dimensión $27$). Además, $\operatorname{Herm}_1(\mathbb K)\cong\mathbb R$ y
$\operatorname{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}$ para $\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O$.

*Demostración.* El apéndice E lo demuestra por completo: la descomposición es el teorema E.9, la lista y «exactamente
una» son el teorema E.27(1),(2), y los isomorfismos en rango uno y dos son el teorema E.24(3). $\square$

(P. Jordan, J. von Neumann, E. Wigner, *On an algebraic generalization of the quantum mechanical formalism*, Ann. of
Math. 35 (1934) 29–64; la demostración del apéndice E sigue [FK94, capítulos II–V].)

Todas ellas salvo el álgebra de Albert son especiales: $\mathbb R=\operatorname{Herm}_1(\mathbb R)$ y
$\operatorname{Herm}_r(\mathbb K)$, $\mathbb K=\mathbb R,\mathbb C$, por el lema 13.1, $\operatorname{Herm}_r(\mathbb H)$
porque escribir los cuaterniones como matrices complejas $2\times2$ la sumerge en $\operatorname{Herm}_{2r}(\mathbb C)$
(teorema E.24(1)), y $V_m$ por el lema 13.3. El álgebra de Albert es un álgebra de Jordan (teorema E.24(2)) y es
**excepcional**: no es isomorfa a ninguna subálgebra de Jordan de un álgebra asociativa real o compleja (teorema E.26,
un argumento del tipo usado por A. A. Albert, *On a certain algebra of quantum mechanics*, Ann. of Math. 35 (1934)
65–73: una inmersión de Jordan en un álgebra asociativa haría asociativo el producto de los octoniones, que no lo es,
lema E.15(2)). Para $r\ge4$, $\operatorname{Herm}_r(\mathbb O)$ no satisface la identidad de Jordan (corolario E.25),
y un álgebra simple con un marco de Jordan de $r\ge4$ elementos es $\operatorname{Herm}_r(D)$ con un álgebra de
composición asociativa $D\ne\mathbb O$ (teorema E.22), razón por la cual los octoniones solo aparecen para $r\le3$. Así,
$\operatorname{Herm}_3(\mathbb O)$ no puede representarse fielmente mediante matrices reales o complejas (ni dentro de
ninguna álgebra asociativa) con el producto $\frac12(xy+yx)$; un cálculo en ella tiene que usar el producto de los
octoniones (o sus constantes de estructura) en lugar de un producto matricial.

Sea $\Omega$ un cono convexo abierto no vacío en un espacio euclídeo $V$ (un espacio real con producto escalar de
dimensión finita $\ne\{0\}$). Su **cono dual** es
$\Omega^*=\{y:\langle x,y\rangle>0\text{ para todo }x\in\bar\Omega\setminus\{0\}\}$; $\Omega$ es **autodual** si
$\Omega^*=\Omega$, **homogéneo** si las aplicaciones lineales $g$ con $g\Omega=\Omega$ actúan transitivamente sobre
$\Omega$, y **simétrico** si es ambas cosas.

**Teorema 13.8 (Koecher–Vinberg).** Un cono convexo abierto no vacío $\Omega$ en un espacio euclídeo $V$ es simétrico si
y solo si existe un producto $\circ$ sobre $V$ tal que $(V,\circ,\langle\cdot,\cdot\rangle)$ es un álgebra de Jordan
euclídea y $\Omega$ es el interior de su cono de cuadrados.

*Demostración.* Teorema F.12 («si»: teoremas F.5(2) y F.6(2)). $\square$

(M. Koecher, 1957; E. B. Vinberg, 1960; la demostración del apéndice F sigue [FK94, capítulos I–III].) El teorema
13.4(4) demuestra la forma cerrada de la autodualidad para el cono de Lorentz ($y\in Q$ si y solo si
$\operatorname{tr}(x\circ y)\ge0$ para todo $x\in Q$), y el ejercicio 1 su homogeneidad. El cono de cuadrados de una
suma directa de álgebras de Jordan euclídeas (con la suma ortogonal de los productos escalares) es el producto de sus
conos de cuadrados, ya que $(x_1,\dots,x_k)^2=(x_1^2,\dots,x_k^2)$. Así, por la dirección «si», los interiores de los
conos de la programación semidefinida (productos de $\operatorname{Herm}_n(\mathbb R)^+$ o
$\operatorname{Herm}_n(\mathbb C)^+$, que son los conos de cuadrados por el teorema 1.6; para $x\ge0$ real la raíz
cuadrada es real, ya que su conjugada compleja es también una raíz cuadrada positiva de $x$), de SOCP (productos de
conos de Lorentz) y de la programación lineal ($\mathbb R_+^n$, el cono de cuadrados de
$\mathbb R^n=\mathbb R\oplus\dots\oplus\mathbb R$) son conos simétricos.

[FK94]: J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, Oxford University Press, 1994.

## 13.5 Aplicaciones de Jordan frente a *-aplicaciones

Un **$*$-homomorfismo de Jordan** entre álgebras de von Neumann es una aplicación $\mathbb C$-lineal $\Phi$ con
$\Phi(x^*)=\Phi(x)^*$ y $\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$; sobre las partes autoadjuntas es un homomorfismo de
Jordan, y todo homomorfismo de Jordan entre partes autoadjuntas se extiende a uno de ellos mediante
$\Phi(x+iy)=\Phi(x)+i\Phi(y)$ ($x,y$ autoadjuntos). (Esta $\Phi$ es $\mathbb C$-lineal y satisface
$\Phi(x-iy)=\Phi(x+iy)^*$; y como $\circ$ es $\mathbb C$-bilineal y lleva pares de elementos autoadjuntos a elementos
autoadjuntos, $(x+iy)\circ(x'+iy')=(x\circ x'-y\circ y')+i(x\circ y'+y\circ x')$ se lleva a
$\Phi(x+iy)\circ\Phi(x'+iy')$.)

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

**Teorema 13.10 (Kadison, dimensión finita).** Sean $\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$ y
$\mathcal N=\bigoplus_dM_{k'_d}(\mathbb C)$, y sea $\Phi:\mathcal M\to\mathcal N$ una biyección lineal con $\Phi(1)=1$
tal que $\Phi$ y $\Phi^{-1}$ son positivas. Entonces $\Phi$ es un $*$-isomorfismo de Jordan. Recíprocamente, todo
$*$-isomorfismo de Jordan $\mathcal M\to\mathcal N$ es unital, y tanto él como su inverso son positivos.

*Demostración.* Teorema F.15. (La idea: $\Phi$ lleva el conjunto convexo $\{0\le a\le1\}$ sobre el conjunto
correspondiente de $\mathcal N$ y, por tanto, los puntos extremos, que son las proyecciones, sobre proyecciones;
conserva $p+q\le1$, lo que para proyecciones significa $pq=0$; así lleva descomposiciones espectrales a
descomposiciones espectrales y $\Phi(x^2)=\Phi(x)^2$.) $\square$

R. V. Kadison demostró el teorema para todas las C\*-álgebras unitales (*A generalized Schwarz inequality and algebraic
invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503); en este curso solo se demuestra y se usa el caso
de dimensión finita. Para los factores $M_n$ se conocen exactamente los automorfismos de Jordan: todo
$*$-automorfismo de Jordan de $M_n(\mathbb C)$ es $A\mapsto UAU^*$ o $A\mapsto UA^TU^*$ con $U$ unitaria (teorema
F.17, para todo $n$; el resultado se debe a N. Jacobson, C. E. Rickart, Trans. Amer. Math. Soc. 69 (1950) 479–502, y
a I. N. Herstein, Trans. Amer. Math. Soc. 81 (1956) 331–341). Las primeras aplicaciones son multiplicativas y las
segundas antimultiplicativas, de modo que un $*$-automorfismo de Jordan de $M_n$ o bien conserva todos los
conmutadores o bien los invierte todos: $\Phi([x,y])=\pm[\Phi(x),\Phi(y)]$. Sobre $\bigoplus_cM_{k_c}$, un
$*$-isomorfismo de Jordan permuta bloques del mismo tamaño y actúa sobre cada bloque de una de las dos maneras
(teorema F.18); es unital y positivo en ambos sentidos (teorema 13.10) y, por el corolario F.19, es completamente
positivo si y solo si no traspone ningún bloque de tamaño $\ge2$, si y solo si es un $*$-isomorfismo. He aquí una
demostración directa del caso $n=2$ en el lenguaje de los factores de espín.

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
(A. Peres, Phys. Rev. Lett. 77 (1996) 1413). Se dice que un estado $\rho$ de $M_n\otimes M_m$ es **separable** si
$\rho=\sum_kp_ka_k\otimes b_k$ (una suma finita) con $p_k\ge0$, $a_k,b_k\ge0$, y **entrelazado** en caso contrario. Un
estado separable está en ambos conos, ya que $\theta(\rho)=\sum_kp_ka_k\otimes b_k^T$ con $b_k^T\ge0$ (proposición
13.9(a)) y $a\otimes b=(a^{1/2}\otimes b^{1/2})^2\ge0$. Así, $\omega$, que no está en
$\theta(\operatorname{Herm}_4^+)$, está entrelazado. (Lecturas adicionales sobre compuestos de álgebras de Jordan, que
no se usan aquí: H. Hanche-Olsen, *On the structure and tensor products of JC-algebras*, Canad. J. Math. 35 (1983)
1059–1074; H. Barnum, M. Graydon, A. Wilce, *Composites and categories of Euclidean Jordan algebras*, Quantum 4 (2020)
359.)

**Dinámica.** Cabría esperar que la evolución temporal obligue a usar el producto asociativo. No es así:

**Proposición 13.13.** Sean $H\in\operatorname{Herm}_n(\mathbb C)$ y $\delta(z)=i[H,z]$. (a) $\delta$ lleva
$\operatorname{Herm}_n$ en sí mismo y es una derivación de $\circ$; $e^{t\delta}(z)=e^{itH}ze^{-itH}$ es un
automorfismo de Jordan. (b) $[L(x),L(y)]z=\frac14[[x,y],z]$ para $x,y,z\in\operatorname{Herm}_n$. (c)
$\delta=4\sum_k[L(x_k),L(y_k)]$ para ciertos $x_k,y_k\in\operatorname{Herm}_n$: la dinámica hamiltoniana está
generada por el producto de Jordan por sí solo.

*Demostración.* (a) $(i[H,z])^*=-i[z,H]=i[H,z]$; $[H,\cdot]$ es una derivación del producto asociativo y, por tanto,
de $\circ$. $z(t)=e^{itH}ze^{-itH}$ satisface $z'=i[H,z]$ y $z(0)=z$ (proposición 8.2(3) con $\mathcal M=M_n$),
así que $z(t)=e^{t\delta}(z)$ por la unicidad del lema 8.1; la conjugación por la unitaria $e^{itH}$ conserva productos
y adjuntos y, por tanto, $\circ$ y $\operatorname{Herm}_n$.
(b) $4x\circ(y\circ z)=xyz+xzy+yzx+zyx$; réstese lo mismo con $x,y$ intercambiados. (c) Por (b) basta
escribir $iH-\frac{\operatorname{Tr}H}n i1$ (que tiene los mismos conmutadores que $iH$) como suma de $[x_k,y_k]$.
Para $i\ne j$ y $a=E_{ij}+E_{ji}$, $b=i(E_{ij}-E_{ji})$: $[a,E_{ii}]=E_{ji}-E_{ij}$,
$[b,E_{ii}]=-i(E_{ij}+E_{ji})$ y $[a,b]=2i(E_{jj}-E_{ii})$; estas generan sobre $\mathbb R$ las matrices
antihermíticas de traza nula. $\square$

Lo que el álgebra de Jordan no contiene es la **correspondencia** $H\mapsto i[H,\cdot]$ entre un observable y la
dinámica que genera. En $\operatorname{Herm}_n(\mathbb R)$, por ejemplo, las derivaciones
$[L(x),L(y)]=\frac14[[x,y],\cdot]$ provienen de matrices $[x,y]$ reales antisimétricas, que no son $i$ por una matriz
simétrica real (ejercicio 4). La correspondencia no está determinada por $\circ$: sobre
$\operatorname{Herm}_n(\mathbb C)$, tanto $H\mapsto i[H,\cdot]$ como $H\mapsto-i[H,\cdot]$ envían todo observable a una
derivación de $\circ$ (proposición 13.13(a)), y la traspuesta, un automorfismo de Jordan, las intercambia, ya que
$T(i[H,T(z)])=-i[H^T,z]$. La segunda es $i$ por el conmutador del álgebra opuesta ($a\cdot_{\rm op}b=ba$):
$-i[H,z]=i(H\cdot_{\rm op}z-z\cdot_{\rm op}H)$ — los dos signos de la proposición 13.12. (Lecturas adicionales sobre
tales «correspondencias dinámicas», que no se usan aquí: E. M. Alfsen, F. W. Shultz, *Geometry of State Spaces of
Operator Algebras*, Birkhäuser, 2003.)

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

> **Significado físico** (interpretación física, no un teorema de este curso). Para una partícula sin espín, la
> traspuesta (la conjugación compleja en una base real) es la **inversión temporal**. Es una simetría de los observables de un único sistema, pero es antiunitaria y, aplicada
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
lo que tiene sentido en toda álgebra de Jordan euclídea una vez fijada una traza. La demostramos en $V_m$ con un argumento que usa solo tres
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
idempotentes son cuadrados); (iii) $\operatorname{tr}$ es lineal, $\operatorname{tr}(a\circ1)=\operatorname{tr}a$, y
$\operatorname{tr}(a\circ b)\ge0$ para $a,b\in Q$. Las tres valen en toda álgebra de Jordan euclídea $V$ con unidad $e$
para la forma lineal $\tau(x)=\langle x,e\rangle$, para la cual $\tau(a\circ b)=\langle b,a\circ e\rangle=\langle a,b\rangle$:
(ii) es el teorema F.3, (i) es el teorema F.6(4), y $\langle a,b\rangle\ge0$ para $a,b\in Q$ es el teorema F.5(1). El
corolario F.7 lleva a cabo el argumento: $\tau(|P(y)x|)\le\tau(P(y)|x|)$ en toda álgebra de Jordan euclídea, con
igualdad si $x\in Q$ o $-x\in Q$. En $M_{sa}$ con $\langle x,y\rangle=\operatorname{Tr}(xy)$, $\tau=\operatorname{Tr}$ y
$P(y)x=yxy$; los idempotentes de $(M_{sa},\circ)$ son las proyecciones ($c\circ c=c^2$), y para proyecciones
$c\circ c'=0$ si y solo si $cc'=0$ (lema F.14(2)), así que, por la unicidad del teorema F.3, su descomposición de $x$ es
la descomposición espectral en proyecciones propias y $|x|$ es el valor absoluto matricial. Así, esto vuelve a demostrar
la proposición 11.9 para todo $Y=Y^*$. En $V_m$ con $\langle x,y\rangle=\operatorname{tr}(x\circ y)$,
$\tau=\operatorname{tr}$, la descomposición del teorema F.3 es la del teorema 13.4(2) (con $x=a1$ si $u=0$), y esto es
la proposición 13.14. En un álgebra de Jordan euclídea simple $J$ de rango $r$, la traza
$\operatorname{tr}x=\frac r{\dim J}\operatorname{Tr}L(x)$ del corolario E.28 cumple
$\operatorname{tr}(\sum_i\mu_ic_i)=\sum_i\mu_i$ para todo marco de Jordan $(c_i)$ (teorema E.7), y
$\langle x,y\rangle=\operatorname{tr}(x\circ y)$ es un producto escalar asociativo; para él $\tau=\operatorname{tr}$, así
que la desigualdad vale para esta traza. Esto se aplica a $\operatorname{Herm}_r(\mathbb H)$ y al álgebra de Albert, que
son álgebras de Jordan euclídeas simples (teorema E.24). Para $r\ge3$ no son isomorfas como álgebras de Jordan a
$M_{sa}$ para ningún $\mathcal M=\bigoplus_cM_{k_c}$: $M_{sa}=\bigoplus_c\operatorname{Herm}_{k_c}(\mathbb C)$ es una
suma directa de ideales, así que solo es simple con un único bloque, y $\operatorname{Herm}_k(\mathbb C)$
($\cong\mathbb R$ para $k=1$, $\cong V_3$ para $k=2$) no es isomorfa a $\operatorname{Herm}_r(\mathbb H)$, $r\ge3$, ni a
$\operatorname{Herm}_3(\mathbb O)$, por la unicidad del teorema 13.7. En $V_m$, la desigualdad se sigue también de la
proposición 11.9 (extendida más arriba a $Y=Y^*$) a través de $\phi$, ya que $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ y
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
    wx, Ux = torch.linalg.eigh(Xs)
    absX = Ux @ torch.diag_embed(wx.abs().to(dt)) @ Ux.mH                  # |phi(x)|
    tr_rhs = (Ys @ absX @ Ys).diagonal(dim1=-2, dim2=-1).sum(-1).real    # Tr phi(y) |phi(x)| phi(y)
    assert (torch.linalg.eigvalsh(Ys)[:, 0] < 0).any()                  # some phi(y) are indefinite
    assert (tr_rhs - tn >= -1e-10 * tr_rhs).all()                       # Proposition 11.9 for Y = Y*
```

La diferencia $\operatorname{tr}P(y)|x|-\operatorname{tr}|P(y)x|$ se anula para $x$ definido (proposición 13.14).
Observación numérica (no un teorema): para estas muestras gaussianas, las diferencias máximas impresas son del orden de
$10^1$–$10^2$, así que la desigualdad está lejos de ser una igualdad para $x$ indefinido.

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
   ($\det x^{1/2}=\sqrt{\lambda_+\lambda_-}>0$; el interior de $Q$ es $\{a>|u|\}=\{\lambda_->0\}$, ya que $a-|u|$ es
   continua y $(a-\varepsilon,u)\notin Q$ cuando $a=|u|$). Dos puntos interiores cualesquiera están unidos por
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
   $p>\frac13$; entonces $\rho_p$ no es separable por el criterio PPT.

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

**Cómo seguir.** Este es el último capítulo del curso; los apéndices E y F completan sus demostraciones. Como lectura
adicional, la referencia estándar para la teoría en dimensión finita es [FK94]; K. McCrimmon, *A Taste of Jordan Algebras* (Springer 2004), cubre el álgebra y la historia, y H. Hanche-Olsen y
E. Størmer, *Jordan Operator Algebras* (Pitman 1984), y Alfsen–Shultz (2003), la teoría en dimensión infinita y su
relación con las C\*-álgebras. Para el papel de los sistemas compuestos en las reconstrucciones de la teoría cuántica a
partir de axiomas operacionales, véanse los artículos de Barnum y Wilce citados más arriba. Para desigualdades de trazas
más allá del capítulo 11, E. Carlen, *Trace inequalities and quantum entropy* (Contemp. Math. 529, 2010). La
representación cuadrática y el cono de Lorentz reaparecen en la optimización cónica: Alizadeh–Goldfarb (2003) y los
métodos de punto interior jordanianos que allí se revisan.
