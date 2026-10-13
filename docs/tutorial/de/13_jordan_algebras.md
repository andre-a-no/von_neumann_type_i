[English](../13_jordan_algebras.md) · [Español](../es/13_jordan_algebras.md) · [Français](../fr/13_jordan_algebras.md) · **Deutsch** · [中文](../zh/13_jordan_algebras.md) · [日本語](../ja/13_jordan_algebras.md) · [Русский](../ru/13_jordan_algebras.md)

# 13. Jordan-Algebren: warum die Quantenmechanik *-Algebren verwendet

Das Produkt $xy$ zweier Observablen ist keine Observable, es sei denn, $x$ und $y$ kommutieren, das symmetrisierte
Produkt $x\circ y=\frac12(xy+yx)$ dagegen immer. 1933 schlug P. Jordan vor, die Observablen mit diesem Produkt als
Grundstruktur der Quantenmechanik zu nehmen, und Jordan, von Neumann und Wigner klassifizierten die
endlichdimensionalen Algebren, die dabei auftreten. Dieses Kapitel entwickelt die endlichdimensionale Theorie, mit
vollständigen Beweisen, wo diese kurz sind, und fragt, was der Jordan-Struktur entgeht. Die Antwort ist das Vorzeichen
des Kommutators: Die Jordan-Algebra eines Systems bestimmt den Kommutator nur bis auf das Vorzeichen (Abschnitt 13.5),
zusammengesetzte Systeme brauchen aber das Vorzeichen, und es entscheidet darüber, ob ein Zustand verschränkt ist.

**Sie lernen** das Jordan-Produkt und die quadratische Darstellung $P(y)x=yxy$ kennen, die im Gegensatz zu $\circ$
die Positivität erhält; Spinfaktoren, ihre Spektraltheorie und den Lorentz-Kegel; Jordan-Unteralgebren hermitescher
Matrizen, die keine selbstadjungierten Teile von $*$-Algebren sind; die Klassifikation von Jordan–von Neumann–Wigner;
warum die Transposition ein Jordan-Automorphismus und kein $*$-Automorphismus ist und was sie auf einem
zusammengesetzten System bewirkt; und eine Spurungleichung aus Kapitel 11, in Jordan-Sprache geschrieben und in einem
Spinfaktor bewiesen.

Im Folgenden ist $\operatorname{Herm}_n(\mathbb K)$ der reelle Vektorraum der hermiteschen $n\times n$-Matrizen über
$\mathbb K=\mathbb R$ oder $\mathbb C$, $M_{sa}$ der selbstadjungierte Teil von $\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$,
$X,Y,Z$ sind die Pauli-Matrizen und $\sigma=(X,Y,Z)$, also $u\cdot\sigma=u_1X+u_2Y+u_3Z$ für $u\in\mathbb R^3$.
$[a,b]=ab-ba$, und $|u|$ ist die euklidische Norm von $u\in\mathbb R^m$.

## 13.1 Observablen und das Jordan-Produkt

**Definition.** Eine **Jordan-Algebra** ist ein reeller Vektorraum $J$ mit einem bilinearen, kommutativen Produkt
$x\circ y$, das die **Jordan-Identität** $(x^2\circ y)\circ x=x^2\circ(y\circ x)$ erfüllt, wobei $x^2=x\circ x$. Sie
heißt **speziell**, wenn sie zu einem Unterraum einer assoziativen Algebra isomorph ist, der unter $\frac12(xy+yx)$
abgeschlossen ist. Auf $J$ sei $L(y)x=y\circ x$; die **quadratische Darstellung** ist $P(y)=2L(y)^2-L(y^2)$.

**Lemma 13.1.** Sei $\mathcal A$ eine assoziative reelle oder komplexe Algebra und $x\circ y=\frac12(xy+yx)$. Dann
erfüllt $\circ$ die Jordan-Identität, und die Jordan-Potenzen $x^{\circ1}=x$, $x^{\circ(k+1)}=x\circ x^{\circ k}$
stimmen mit den assoziativen Potenzen $x^k$ überein. Insbesondere sind $M_{sa}$, $\operatorname{Herm}_n(\mathbb C)$
und $\operatorname{Herm}_n(\mathbb R)$ spezielle Jordan-Algebren.

*Beweis.* Ausmultiplizieren ergibt $4(x^2\circ y)\circ x=x^2yx+yx^3+x^3y+xyx^2$ und
$4x^2\circ(y\circ x)=x^2yx+x^3y+yx^3+xyx^2$, dieselben vier Terme. Durch Induktion ist
$x\circ x^k=\frac12(x^{k+1}+x^{k+1})=x^{k+1}$. Sind $x,y$ hermitesch, so ist $(xy+yx)^*=yx+xy$, also sind die drei Räume
unter $\circ$ abgeschlossen; für reelle symmetrische $x,y$ ist das Produkt reell. $\square$

**Proposition 13.2.** Seien $x,y,z\in M_{sa}$.

1. $P(y)x=yxy$, und es gilt die **Fundamentalformel** $P(P(y)x)=P(y)P(x)P(y)$.
2. $\operatorname{Tr}P(y)x=\operatorname{Tr}y^2x$.
3. Ist $x\ge0$, so ist $P(y)x\ge0$. Ist $y$ invertierbar, so ist $P(y)$ eine Bijektion des positiven Kegels auf sich,
   mit Inverser $P(y^{-1})$.
4. Das Jordan-Produkt positiver Elemente muss nicht positiv sein: Für $p=\begin{pmatrix}1&0\\0&0\end{pmatrix}$ und
   $q=\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}$ ist $\det(p\circ q)=-\frac1{16}$.

*Beweis.* (1) $2y\circ(y\circ x)=\frac12(y^2x+2yxy+xy^2)$ und $y^2\circ x=\frac12(y^2x+xy^2)$; man subtrahiere. Dann
ist $P(P(y)x)z=yxy\,z\,yxy=P(y)P(x)P(y)z$ wegen der Assoziativität. (2) Zyklizität der Spur. (3)
$yxy=(x^{1/2}y)^*(x^{1/2}y)\ge0$ (Korollar 1.7), und $P(y)P(y^{-1})=P(y^{-1})P(y)=\mathrm{id}$. (4)
$p\circ q=\frac14\begin{pmatrix}2&1\\1&0\end{pmatrix}$. $\square$

Schreibt man $P(y)$ als $2L(y)^2-L(y^2)$, so ist jede Identität in (1) allein mit $\circ$ ausgedrückt und gilt daher in jeder speziellen Jordan-Algebra; die
Fundamentalformel gilt in jeder Jordan-Algebra (K. McCrimmon, *A Taste of Jordan Algebras*, Springer 2004; wir geben
den Beweis nicht wieder). Das Modul `torch_vn_algebra.jordan` implementiert diese Abbildungen auf Operatoren einer
`TypeIAlgebra`, verzögert ausgewertet und gebündelt über alle Sektoren:

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

> **Physikalische Bedeutung.** $x\circ y$ ist das symmetrisierte Produkt von Observablen: die einzige Kombination
> $\alpha xy+\beta yx$, die in $x,y$ symmetrisch ist und sich für $y=x$ auf $x^2$ reduziert, und sie ist wieder eine
> Observable. Die quadratische Darstellung ist die **Lüders-Regel**: Misst man eine Ja/Nein-Frage $p$ und behält das
> Ergebnis „ja“, so geht eine Dichte $\rho$ in $P(p)\rho/\operatorname{Tr}P(p)\rho=p\rho p/\operatorname{Tr}(p\rho)$
> über. Proposition 13.2(4) besagt, dass $\circ$ selbst für positive Elemente keine solche operationelle Bedeutung hat.

## 13.2 Spinfaktoren und der Lorentz-Kegel

**Definition.** Für $m\ge1$ hat der **Spinfaktor** $V_m=\mathbb R\oplus\mathbb R^m$ die Elemente $x=(a,u)$, das
Produkt $(a,u)\circ(b,v)=(ab+u\cdot v,\;av+bu)$ und die Eins $1=(1,0)$. Seine **Spur** ist $\operatorname{tr}x=2a$.

Das Produkt ist bilinear und kommutativ; die Jordan-Identität folgt aus einer Matrixdarstellung.

**Lemma 13.3 (antikommutierende Symmetrien).** Seien $p=\lfloor m/2\rfloor$ und $D=2^p$ ($D=2$ für $m=1$). Auf
$(\mathbb C^2)^{\otimes p}$ setze man $s_{2j-1}=Z^{\otimes(j-1)}\otimes X\otimes1$, $s_{2j}=Z^{\otimes(j-1)}\otimes Y\otimes1$
für $j\le p$ und $s_m=Z^{\otimes p}$, falls $m$ ungerade ist ($s_1=Z$ für $m=1$). Dann gilt $s_i=s_i^*$ und
$s_is_j+s_js_i=2\delta_{ij}1$. Die Abbildung $\phi(a,u)=a1+\sum_iu_is_i$ ist injektiv und erfüllt
$\phi(x\circ y)=\phi(x)\circ\phi(y)$; folglich ist $V_m$ eine spezielle Jordan-Algebra. Für $m\ge2$ erzeugen die
$s_i$ die $*$-Algebra $M_D(\mathbb C)$.

*Beweis.* Pauli-Matrizen quadrieren zu $1$ und antikommutieren paarweise; dasselbe gilt für $Z$ mit $X$ und $Y$. Zwei
verschiedene Ketten $s_i,s_j$ unterscheiden sich an genau einem Platz in antikommutierender Weise: für
$s_{2j-1},s_{2j}$ am Platz $j$ ($X$ gegen $Y$); für Ketten mit Indizes $j<j'$ am Platz $j$ ($X$ oder $Y$ gegen $Z$),
während sie an den übrigen Plätzen gleiche Faktoren oder einen Faktor $1$ tragen; dasselbe gilt für
$s_m=Z^{\otimes p}$. Also ist $s_is_j=-s_js_i$ für $i\ne j$ und $s_i^2=1$. Somit ist
$(u\cdot s)\circ(v\cdot s)=\sum_{ij}u_iv_j\,s_i\circ s_j=(u\cdot v)1$, was $\phi(x)\circ\phi(y)=\phi(x\circ y)$ ergibt.
Jedes $s_i$ ist spurfrei (ein Tensorprodukt mit einem spurfreien Faktor), und
$\operatorname{Tr}(s_is_j)=\operatorname{Tr}(s_i\circ s_j)=D\delta_{ij}$, also sind $1,s_1,\dots,s_m$ orthogonal
bezüglich $\operatorname{Tr}(ab)$, und $\phi$ ist injektiv. Die Jordan-Identität gilt in $\operatorname{Herm}_D$
(Lemma 13.1) und wird durch $\phi$ zurücktransportiert. Schließlich ist
$s_{2j-1}s_{2j}=1\otimes\cdots\otimes XY\otimes\cdots\otimes1=iZ_j$ ($Z$ am Platz $j$), also enthält die erzeugte
$*$-Algebra $Z_j$, dann $X_j=Z_1\cdots Z_{j-1}s_{2j-1}$ und $Y_j$, also alle Tensorprodukte von Pauli-Matrizen, eine
Basis von $M_D(\mathbb C)$. $\square$

Bis auf Vorzeichen sind $s_{2j-1}$ und $s_{2j}$ die Majorana-Operatoren $c_j+c_j^*$ und $i(c_j-c_j^*)$ der
Jordan-Wigner-Fermionen $c_j$ aus Abschnitt 9.4.

**Satz 13.4 (Spektraltheorie von $V_m$).** Seien $x=(a,u)\in V_m$ und $\lambda_\pm=a\pm|u|$.

1. Die Idempotente von $V_m$ sind $0$, $1$ und $c_w=\frac12(1,w)$ mit $|w|=1$; und $c_w\circ c_{w'}=0$ genau dann,
   wenn $w'=-w$.
2. Ist $u\ne0$, so ist $x=\lambda_+e_++\lambda_-e_-$ mit $e_\pm=c_{\pm u/|u|}$; ist $u=0$, so ist $x=a1$. Ist
   $x=\mu_1c_w+\mu_2c_{-w}$, so gilt $\{\mu_1,\mu_2\}=\{\lambda_+,\lambda_-\}$. Wir nennen $\lambda_\pm$ die
   **Eigenwerte** von $x$.
3. Für $f:\{\lambda_+,\lambda_-\}\to\mathbb R$ setze man $f(x)=f(\lambda_+)e_++f(\lambda_-)e_-$ ($f(a)1$ für $u=0$).
   Dann ist $x^{\circ k}=\lambda_+^ke_++\lambda_-^ke_-$, $\phi(f(x))=f(\phi(x))$, die Eigenwerte von $\phi(x)$ sind
   $\lambda_\pm$, jeweils mit Vielfachheit $D/2$, und $\operatorname{Tr}\phi(x)=\frac D2\operatorname{tr}x$. Außerdem
   ist $\operatorname{tr}x=\lambda_++\lambda_-$ und $\det x:=\lambda_+\lambda_-=a^2-|u|^2$.
4. Der Kegel der Quadrate $Q=\{x^2:x\in V_m\}$ ist der **Lorentz-Kegel** $\{(a,u):a\ge|u|\}=\{x:\lambda_-\ge0\}=\{x:\phi(x)\ge0\}$.
   Er ist **selbstdual**: $y\in Q$ genau dann, wenn $\operatorname{tr}(x\circ y)\ge0$ für alle $x\in Q$.
5. $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$; folglich ist $P(y)Q\subseteq Q$ für jedes $y\in V_m$.

*Beweis.* (1) $(a,u)^2=(a^2+|u|^2,2au)$ ist genau dann gleich $(a,u)$, wenn $u=0$, $a\in\{0,1\}$, oder $a=\frac12$,
$|u|=\frac12$. $c_w\circ c_{w'}=\frac14(1+w\cdot w',\,w+w')$ verschwindet genau dann, wenn $w'=-w$. (2)
$\lambda_+e_++\lambda_-e_-=(a,0)+|u|(0,u/|u|)=x$, und $\mu_1c_w+\mu_2c_{-w}=(\frac{\mu_1+\mu_2}2,\frac{\mu_1-\mu_2}2w)$,
also $\mu_1+\mu_2=2a$ und $|\mu_1-\mu_2|=2|u|$. (3) Nach (1) ist $(\alpha e_++\beta e_-)\circ(\gamma e_++\delta e_-)=\alpha\gamma
e_++\beta\delta e_-$, was durch Induktion die Potenzen liefert. Mit $w=u/|u|$ ist $(w\cdot s)^2=1$, also sind
$\phi(e_\pm)=\frac12(1\pm w\cdot s)$ Orthogonalprojektionen mit Summe $1$ und Produkt $0$, jede mit Spur $D/2$; somit
ist $\phi(x)=\lambda_+\phi(e_+)+\lambda_-\phi(e_-)$ eine Spektralzerlegung, und
$f(\phi(x))=f(\lambda_+)\phi(e_+)+f(\lambda_-)\phi(e_-)=\phi(f(x))$ (Satz 1.5). (4) $x^2$ hat die Eigenwerte
$\lambda_\pm^2\ge0$; ist umgekehrt $\lambda_-\ge0$, so ist $x=(\sqrt x)^2$. Es ist $\lambda_-\ge0$ genau dann, wenn
$a\ge|u|$, und nach (3) genau dann, wenn $\phi(x)\ge0$. Zur Selbstdualität: $\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$.
Ist $a\ge|u|$ und $b\ge|v|$, so ist $ab\ge|u||v|\ge-u\cdot v$. Ist $b<|v|$, so liefert $x=(|v|,-v)\in Q$ den Wert
$2(|v|b-|v|^2)<0$ (für $v=0$, $b<0$ nehme man $x=1$). (5) $P(y)x=2y\circ(y\circ x)-y^2\circ x$, und $\phi$ erhält
$\circ$, also ist $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ nach Proposition 13.2(1); man wende (4) und Proposition 13.2(3)
an. $\square$

**Proposition 13.5.** (a) Für $m=3$ ist $\phi(a,u)=a1+u\cdot\sigma$ ein Jordan-Isomorphismus
$V_3\cong\operatorname{Herm}_2(\mathbb C)$. (b) Die Jordan-Automorphismen von $V_m$ sind die Abbildungen
$(a,u)\mapsto(a,Ou)$, $O\in O(m)$. (c) Für $m\ge2$ hat $\phi(V_m)$ die Dimension $m+1$, und der selbstadjungierte Teil
der von ihm erzeugten $*$-Algebra hat die Dimension $D^2=4^{\lfloor m/2\rfloor}$; gleich sind sie nur für $m=3$.

*Beweis.* (a) Für $m=3$ ist $(s_1,s_2,s_3)=(X,Y,Z)$, und $\phi$ ist injektiv zwischen Räumen der Dimension $4$. (b)
Die Abbildungen $(a,u)\mapsto(a,Ou)$ erhalten $u\cdot v$, also $\circ$. Sei umgekehrt $\alpha$ ein
Jordan-Automorphismus. Dann ist $\alpha(1)$ eine Eins von $V_m$, also $\alpha(1)=1$. Nach der Formel für $x^2$ gilt
$x^2=1$ genau dann, wenn $x=\pm1$ oder $x=(0,w)$ mit $|w|=1$. Für $|w|=1$ ist $\alpha(0,w)^2=1$ und
$\alpha(0,w)\ne\pm1=\alpha(\pm1)$, also $\alpha(0,w)=(0,w')$ mit $|w'|=1$. Schreibt man $\alpha(0,u)=(\ell(u),Ou)$ mit
linearen $\ell,O$, so besagt dies $\ell=0$ und $|Ou|=|u|$, d. h. $O\in O(m)$. (c) Lemma 13.3;
$4^{\lfloor m/2\rfloor}>m+1$ für $m=2$ und $m\ge4$. $\square$

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

Für $m=1$ erzeugt $s_1=Z$ die Diagonalmatrizen, der Dimension $2=m+1$, und $V_1\cong\mathbb R\oplus\mathbb R$ ist
nicht einfach. Die ausgegebenen Dimensionen sind $3,4$ für $m=2$; $5,16$ für $m=4$; $7,64$ für $m=6$.

> **Physikalische Bedeutung.** Ein Zustand von $V_m$ (ein positives normiertes lineares Funktional) ist wegen der
> Selbstdualität von der Form $x\mapsto\operatorname{tr}(\rho\circ x)$ mit $\rho=\frac12(1,r)\in Q$, $|r|\le1$: Der
> Zustandsraum ist eine Kugel der Dimension $m$. Für $m=3$ ist es die Bloch-Kugel eines Qubits (Proposition
> 13.5(a)). Der Lorentz-Kegel ist der Zukunftslichtkegel der speziellen Relativitätstheorie in $1+m$ Dimensionen, und
> Proposition 13.5(b) besagt, dass die Jordan-Automorphismen die orthogonalen Abbildungen des Raumes sind, die die
> Zeitachse festhalten.

Die Kegel $Q$ sind die Nebenbedingungsmengen der **Kegelprogrammierung zweiter Ordnung** (SOCP): lineare
Optimierung über Produkten von Lorentz-Kegeln. Innere-Punkte-Verfahren für SOCP verwenden genau Satz 13.4:
Eigenwerte, $\det$ und $P(y)$ für die Skalierung (F. Alizadeh, D. Goldfarb, *Second-order cone programming*, Math.
Program. 95 (2003) 3–51).

## 13.3 Jordan-Unteralgebren, die keine *-Unteralgebren sind

Eine Jordan-Unteralgebra $J\subseteq\operatorname{Herm}_n(\mathbb C)$ erzeugt eine komplexe $*$-Algebra $\mathcal A$
mit $J\subseteq\mathcal A_{sa}$. Proposition 13.5(c) lieferte Beispiele mit $J\ne\mathcal A_{sa}$; hier sind zwei
weitere.

**Proposition 13.6.** (a) $\operatorname{Herm}_n(\mathbb R)$ ist eine Jordan-Unteralgebra von
$\operatorname{Herm}_n(\mathbb C)$ der Dimension $n(n+1)/2$; die von ihr erzeugte $*$-Algebra ist $M_n(\mathbb C)$,
deren selbstadjungierter Teil die Dimension $n^2$ hat. (b) Zu $x,y\in\operatorname{Herm}_2(\mathbb C)$ gibt es eine
unitäre Matrix $U$, sodass $UxU^*$ und $UyU^*$ reell symmetrisch sind. Die von $1,x,y$ erzeugte Jordan-Algebra ist
$\operatorname{span}\{1,x,y\}$; sind $x,y,1$ linear unabhängig, so ist sie $U^*\operatorname{Herm}_2(\mathbb R)U\cong V_2$,
der Dimension $3$, während die von $x,y$ erzeugte $*$-Algebra $M_2(\mathbb C)$ ist, der Dimension $4$.

*Beweis.* (a) Lemma 13.1. Für $i\ne j$ ist $E_{ij}=E_{ii}(E_{ij}+E_{ji})$ ein Produkt reeller symmetrischer Matrizen,
also enthält die erzeugte Algebra alle Matrixeinheiten. (b) Man schreibe $x=a1+u\cdot\sigma$, $y=b1+v\cdot\sigma$ und
wähle einen Einheitsvektor $n\perp u,v$. Die Matrizen $n\cdot\sigma$ und $Y$ sind hermitesch, spurfrei und quadrieren
zu $1$, also haben beide die Eigenwerte $1,-1$, und $U(n\cdot\sigma)U^*=Y$ für ein unitäres $U$ (Satz 1.3). Konjugation
mit $U$ erhält $\operatorname{Tr}$ und $\circ$, bildet also $u\cdot\sigma$ auf $(Ru)\cdot\sigma$ ab, mit linearem $R$
und $Ru\cdot Rv=\frac12\operatorname{Tr}\big((Ru\cdot\sigma)(Rv\cdot\sigma)\big)=u\cdot v$: $R\in O(3)$, und $Rn=e_2$.
Daher ist $Ru,Rv\perp e_2$, und $UxU^*=a1+(Ru)_1X+(Ru)_3Z$ ist reell, da $X,Z$ es sind; ebenso $UyU^*$. Ferner ist
$x\circ y=(ab+u\cdot v)1+av\cdot\sigma+bu\cdot\sigma\in\operatorname{span}\{1,x,y\}$, und ebenso $x^2,y^2$, also ist der
Spann abgeschlossen. Hat er die Dimension $3$, so ist sein Konjugiertes unter $U$ ganz
$\operatorname{Herm}_2(\mathbb R)=\operatorname{span}\{1,X,Z\}$, das nach Lemma 13.3 $\cong V_2$ ist ($X,Z$
antikommutieren) und nach (a) $M_2(\mathbb C)$ erzeugt. $\square$

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

Drei generische hermitesche $2\times2$-Matrizen erzeugen ganz $\operatorname{Herm}_2(\mathbb C)$: Dimension $4$ in
beiden Spalten. In einem Block $M_2$ von $\mathcal M$ verraten zwei Observablen die komplexe Struktur also nie; eine
dritte kann es.

## 13.4 Die Klassifikation

Eine Jordan-Algebra heißt **formal reell**, wenn aus $\sum_ix_i^2=0$ folgt, dass alle $x_i=0$ sind; in endlicher
Dimension ist das gleichbedeutend damit, **euklidisch** zu sein, d. h. ein Skalarprodukt mit
$\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$ zu besitzen [FK94]. $M_{sa}$ (mit $\operatorname{Tr}(xy)$) und
$V_m$ (mit $\operatorname{tr}(x\circ y)$) sind euklidisch.

**Satz 13.7 (Jordan–von Neumann–Wigner; ohne Beweis angegeben).** Jede endlichdimensionale formal reelle
Jordan-Algebra ist eine direkte Summe einfacher Jordan-Algebren, und jede einfache ist zu genau einer der folgenden
isomorph: $\mathbb R$; $\operatorname{Herm}_r(\mathbb K)$ für $\mathbb K=\mathbb R,\mathbb C,\mathbb H$ und $r\ge3$;
den Spinfaktoren $V_m$, $m\ge2$; der **Albert-Algebra** $\operatorname{Herm}_3(\mathbb O)$ der hermiteschen
$3\times3$-Matrizen über den Oktonionen (Dimension $27$). Im Rang $2$ gilt
$\operatorname{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}$ für $\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O$.
(P. Jordan, J. von Neumann, E. Wigner, *On an algebraic generalization of the quantum mechanical formalism*, Ann. of
Math. 35 (1934) 29–64; ein moderner Beweis steht in [FK94, Kapitel V].)

Alle außer der Albert-Algebra sind speziell: $\operatorname{Herm}_r(\mathbb H)$ bettet sich in
$\operatorname{Herm}_{2r}(\mathbb C)$ ein, indem man Quaternionen als komplexe $2\times2$-Matrizen schreibt, und $V_m$
nach Lemma 13.3. Die Albert-Algebra ist **exzeptionell**: Sie ist zu keiner Jordan-Unteralgebra einer assoziativen
Algebra isomorph (A. A. Albert, *On a certain algebra of quantum mechanics*, Ann. of Math. 35 (1934) 65–73). Ein Grund,
der sich kurz angeben lässt: C. M. Glennie fand polynomiale Identitäten, die in jeder speziellen Jordan-Algebra
gelten, in $\operatorname{Herm}_3(\mathbb O)$ aber nicht (*Some identities valid in special Jordan algebras but not
valid in all Jordan algebras*, Pacific J. Math. 16 (1966) 47–59; siehe McCrimmon, *A Taste of Jordan Algebras*). Die
Oktonionen sind nicht assoziativ, und $\operatorname{Herm}_r(\mathbb O)$ erfüllt die Jordan-Identität nur für
$r\le3$, weshalb die Liste dort endet. Daher lässt sich $\operatorname{Herm}_3(\mathbb O)$ nicht treu durch reelle
oder komplexe Matrizen (oder innerhalb irgendeiner assoziativen Algebra) mit dem Produkt $\frac12(xy+yx)$ darstellen;
Software muss mit oktonionischen Einträgen arbeiten.

**Satz 13.8 (Koecher–Vinberg; ohne Beweis angegeben).** Die Inneren der Kegel der Quadrate euklidischer
Jordan-Algebren sind genau die **symmetrischen Kegel**: offene konvexe Kegel in einem euklidischen Raum, die selbstdual
und homogen sind (ihre linearen Automorphismen wirken transitiv). [FK94, Kapitel III.]

Satz 13.4(4) und Aufgabe 1 beweisen beide Eigenschaften für den Lorentz-Kegel. Somit sind die Kegel der
semidefiniten Programmierung (Produkte von $\operatorname{Herm}_n(\mathbb R)^+$ oder $\operatorname{Herm}_n(\mathbb C)^+$),
der SOCP und der linearen Programmierung ($\mathbb R_+^n$, der Kegel von $\mathbb R^n=\mathbb R\oplus\dots\oplus\mathbb R$)
alle von dieser Form, und Innere-Punkte-Verfahren lassen sich ein für alle Mal für alle euklidischen Jordan-Algebren
formulieren.

[FK94]: J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, Oxford University Press, 1994.

## 13.5 Jordan-Abbildungen und *-Abbildungen

Ein **Jordan-$*$-Homomorphismus** von von-Neumann-Algebren ist eine komplex-lineare Abbildung $\Phi$ mit
$\Phi(x^*)=\Phi(x)^*$ und $\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$; auf den selbstadjungierten Teilen ist er ein
Jordan-Homomorphismus, und jeder Jordan-Homomorphismus selbstadjungierter Teile setzt sich durch
$\Phi(x+iy)=\Phi(x)+i\Phi(y)$ zu einem solchen fort.

**Proposition 13.9 (die Transposition).** Sei $T(A)=\bigoplus_cA_c^T$ auf $\mathcal M$. (a) $T$ ist ein unitaler
Jordan-$*$-Automorphismus mit $T(AB)=T(B)T(A)$; er erhält Spektren, also sind $T$ und $T^{-1}=T$ positiv. (b) Auf
$M_k$, $k\ge2$, ist die Choi-Matrix $\sum_{ij}E_{ij}\otimes T(E_{ij})$ der Swap $F(e_a\otimes e_b)=e_b\otimes e_a$,
mit den Eigenwerten $+1$ (Vielfachheit $k(k+1)/2$) und $-1$ (Vielfachheit $k(k-1)/2$). Also ist $T$ nicht vollständig
positiv.

*Beweis.* (a) $(A^T)^*=(A^*)^T$, $(AB)^T=B^TA^T$, also $(A\circ B)^T=B^T\circ A^T=A^T\circ B^T$;
$\det(A^T-\lambda)=\det(A-\lambda)$. (b) $\sum_{ij}E_{ij}e_a\otimes E_{ji}e_b=e_b\otimes e_a$. $F^2=1$, und $F$ ist $+1$
auf den symmetrischen und $-1$ auf den antisymmetrischen Tensoren, die zusammen den Raum aufspannen. Nach Satz 6.2
(und Proposition 6.1) ist $T$ nicht CP. $\square$

Ein Jordan-Automorphismus muss also kein $*$-Automorphismus sein. Bis auf dieses Phänomen bestimmt die Positivität die
Jordan-Struktur:

**Satz 13.10 (Kadison; ohne Beweis angegeben).** Sei $\Phi$ eine lineare Bijektion zwischen unitalen C\*-Algebren mit
$\Phi(1)=1$, sodass $\Phi$ und $\Phi^{-1}$ positiv sind. Dann ist $\Phi$ ein Jordan-$*$-Isomorphismus. (R. V. Kadison,
*A generalized Schwarz inequality and algebraic invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503,
Korollar 5. In *Isometries of operator algebras*, Ann. of Math. 54 (1951) 325–338, hatte er gezeigt, dass eine lineare
Isometrie $\Phi$ einer unitalen C\*-Algebra auf eine andere gleich $\Phi(1)$, einem Unitären, mal einem
Jordan-$*$-Isomorphismus ist.)

Die Umkehrung ist unmittelbar: Ein Jordan-$*$-Isomorphismus bildet Quadrate $h^2=h\circ h$ selbstadjungierter Elemente
auf Quadrate ab, also den positiven Kegel auf den positiven Kegel. Für die Faktoren $M_n$ sind die
Jordan-Automorphismen genau bekannt: Jeder Jordan-$*$-Automorphismus von $M_n(\mathbb C)$ ist $A\mapsto UAU^*$ oder
$A\mapsto UA^TU^*$ mit unitärem $U$ (N. Jacobson, C. E. Rickart, *Jordan homomorphisms of rings*, Trans. Amer. Math.
Soc. 69 (1950) 479–502; I. N. Herstein, *Jordan homomorphisms*, Trans. Amer. Math. Soc. 81 (1956) 331–341). Wir
beweisen den Fall $n=2$.

**Proposition 13.11.** Jeder Jordan-Automorphismus von $\operatorname{Herm}_2(\mathbb C)$ ist $A\mapsto UAU^*$ oder
$A\mapsto UA^TU^*$ mit $U\in U(2)$.

*Beweis.* Nach Proposition 13.5(a,b) wirkt ein Jordan-Automorphismus als $a1+u\cdot\sigma\mapsto a1+(Ou)\cdot\sigma$
mit $O\in O(3)$. Man schreibe $R_U\in O(3)$ für die Matrix von $A\mapsto UAU^*$ in diesem Sinne; $R_{UV}=R_UR_V$. (i)
$\det R_U=1$: $U=W\operatorname{diag}(e^{i\theta_1},e^{i\theta_2})W^*$ ist durch den Weg
$U_t=W\operatorname{diag}(e^{it\theta_1},e^{it\theta_2})W^*$ mit $1$ verbunden, und $\det R_{U_t}\in\{\pm1\}$ hängt
stetig von $t$ ab. (ii) Jedes $R\in SO(3)$ ist ein $R_U$: Wie in Proposition 13.6 gibt es $U_1$ mit
$U_1ZU_1^*=(Re_3)\cdot\sigma$, d. h. $R_{U_1}e_3=Re_3$. Dann hält $R_{U_1}^{-1}R\in SO(3)$ den Vektor $e_3$ fest, ist
also eine Drehung um einen Winkel $\beta$ um $e_3$. Eine direkte Rechnung ergibt für
$U_2=\operatorname{diag}(e^{-i\beta/2},e^{i\beta/2})$, dass $U_2XU_2^*=\cos\beta\,X+\sin\beta\,Y$ und
$U_2YU_2^*=-\sin\beta\,X+\cos\beta\,Y$, also $R=R_{U_1U_2}$. (iii) Die Transposition hält $X,Z$ fest und bildet
$Y\mapsto-Y$ ab: $O_T=\operatorname{diag}(1,-1,1)$, $\det O_T=-1$. Daher ist $O=R_U$, falls $\det O=1$, und $OO_T=R_U$,
d. h. der Automorphismus ist $A\mapsto UA^TU^*$, falls $\det O=-1$. $\square$

**Zusammengesetzte Systeme.** Die Jordan-Struktur eines Faktors kann $T$ nicht von der Identität unterscheiden. Bei
zwei Systemen wird diese Mehrdeutigkeit sichtbar:

**Proposition 13.12.** Seien $x,x'\in\operatorname{Herm}_n(\mathbb C)$, $y,y'\in\operatorname{Herm}_m(\mathbb C)$. In
$\operatorname{Herm}_{nm}(\mathbb C)$, wirkend auf $\mathbb C^n\otimes\mathbb C^m$, gilt
$$
(x\otimes y)\circ(x'\otimes y')=(x\circ x')\otimes(y\circ y')+\tfrac14[x,x']\otimes[y,y'].
$$
Seien $\theta=\mathrm{id}\otimes T$ (die **partielle Transposition**) und
$A\circ_\theta B=\theta(\theta(A)\circ\theta(B))$. Dann ist $\circ_\theta$ ein Jordan-Produkt auf
$\operatorname{Herm}_{nm}(\mathbb C)$, über $\theta$ isomorph zu $\circ$, das mit $\circ$ auf allen Paaren lokaler
Observablen $x\otimes1$, $1\otimes y$ übereinstimmt, und
$(x\otimes y)\circ_\theta(x'\otimes y')=(x\circ x')\otimes(y\circ y')-\frac14[x,x']\otimes[y,y']$. Sein Kegel der
Quadrate ist $\theta(\operatorname{Herm}_{nm}^+)$. Für $n=m=2$ und $\Omega=e_1\otimes e_1+e_2\otimes e_2$ liegt der
Zustand $\omega=\frac12|\Omega\rangle\langle\Omega|\ge0$ nicht in diesem Kegel: $\theta(\omega)=\frac12F$ hat den
Eigenwert $-\frac12$.

*Beweis.* Mit $j=x\circ x'$, $c=\frac12[x,x']$, $k=y\circ y'$, $d=\frac12[y,y']$ ist $xx'=j+c$, $x'x=j-c$, und ebenso
für $y$, also
$\frac12(xx'\otimes yy'+x'x\otimes y'y)=\frac12\big((j+c)\otimes(k+d)+(j-c)\otimes(k-d)\big)=j\otimes k+c\otimes d$.
$\theta$ ist eine reell-lineare Involution von $\operatorname{Herm}_{nm}$: Sie bildet $x\otimes y$ auf $x\otimes y^T$
ab, und die $x\otimes y$ spannen $\operatorname{Herm}_{nm}$ über $\mathbb R$ auf (die Produkte
$h_\alpha\otimes g_\beta$ reeller Basen von $\operatorname{Herm}_n,\operatorname{Herm}_m$ bilden eine komplexe Basis
von $M_{nm}$, und eine hermitesche Matrix hat darin reelle Koeffizienten). Also ist $\circ_\theta$ ein Jordan-Produkt
und $\theta$ ein Isomorphismus. Da $(y^T\circ y'^T)^T=y\circ y'$ und $[y^T,y'^T]^T=-[y,y']$, ergibt die Anwendung von
$\theta$ auf die Formel für $(x\otimes y^T)\circ(x'\otimes y'^T)$ die Formel für $\circ_\theta$. Auf lokalen
Observablen verschwindet der Kommutatorterm ($[1,1]=0$, $[x,1]=0$). Die Quadrate sind
$A\circ_\theta A=\theta(\theta(A)^2)$, und $\omega$ liegt genau dann in $\theta(\operatorname{Herm}_{nm}^+)$, wenn
$\theta(\omega)\ge0$, da $\theta^2=\mathrm{id}$. Schließlich ist
$\theta(\omega)=\frac12\sum_{ij}E_{ij}\otimes E_{ji}=\frac12F$ (Proposition 13.9(b); das ist die Rechnung aus dem
Beweis von Proposition 6.1). $\square$

Die beiden Produkte sind gleich gute Jordan-Algebren, sie stimmen in allem überein, was sich an jeweils einem
Teilsystem messen lässt, und sie unterscheiden sich darin, welche gemeinsamen Operatoren positiv sind, also darin,
was die Zustände des zusammengesetzten Systems sind. Was sie unterscheidet, ist das Vorzeichen des Kommutatorterms,
und dieses Vorzeichen wird durch die $*$-Algebra $M_n\otimes M_m$ festgelegt und nicht durch die Jordan-Algebren der
Faktoren. ($\circ_\theta$ ist das Jordan-Produkt des assoziativen Produkts $A\star B=\theta(\theta(A)\theta(B))$, für
das $(x\otimes y)\star(x'\otimes y')=xx'\otimes y'y$ gilt: die Algebra $M_n\otimes M_m^{\rm op}$.) Der Kegel
$\theta(\operatorname{Herm}_{nm}^+)$ ist es, der das **PPT-Kriterium** liefert (A. Peres, Phys. Rev. Lett. 77 (1996)
1413): Ein separabler Zustand $\rho=\sum_kp_ka_k\otimes b_k$, $p_k\ge0$, $a_k,b_k\ge0$, liegt in beiden Kegeln, da
$\theta(\rho)=\sum_kp_ka_k\otimes b_k^T$ mit $b_k^T\ge0$ (Proposition 13.9(a)) und
$a\otimes b=(a^{1/2}\otimes b^{1/2})^2\ge0$. Also ist $\omega$, das nicht in $\theta(\operatorname{Herm}_4^+)$ liegt,
verschränkt. Zu Ergebnissen darüber, welche Zusammensetzungen von Jordan-Algebren existieren, siehe H. Hanche-Olsen,
*On the structure and tensor products of JC-algebras*, Canad. J. Math. 35 (1983) 1059–1074, und H. Barnum,
M. Graydon, A. Wilce, *Composites and categories of Euclidean Jordan algebras*, Quantum 4 (2020) 359.

**Dynamik.** Man könnte hoffen, dass die Zeitentwicklung das assoziative Produkt erzwingt. Das tut sie nicht:

**Proposition 13.13.** Seien $H\in\operatorname{Herm}_n(\mathbb C)$ und $\delta(z)=i[H,z]$. (a) $\delta$ bildet
$\operatorname{Herm}_n$ in sich ab und ist eine Derivation von $\circ$; $e^{t\delta}(z)=e^{itH}ze^{-itH}$ ist ein
Jordan-Automorphismus. (b) $[L(x),L(y)]z=\frac14[[x,y],z]$ für $x,y,z\in\operatorname{Herm}_n$. (c)
$\delta=4\sum_k[L(x_k),L(y_k)]$ für gewisse $x_k,y_k\in\operatorname{Herm}_n$: Die Hamiltonsche Dynamik wird allein
vom Jordan-Produkt erzeugt.

*Beweis.* (a) $(i[H,z])^*=-i[z,H]=i[H,z]$; $[H,\cdot]$ ist eine Derivation des assoziativen Produkts, also auch von
$\circ$. (b) $4x\circ(y\circ z)=xyz+xzy+yzx+zyx$; man subtrahiere denselben Ausdruck mit vertauschten $x,y$. (c) Nach
(b) genügt es, $iH-\frac{\operatorname{Tr}H}n i1$ (das dieselben Kommutatoren hat wie $iH$) als Summe von $[x_k,y_k]$
zu schreiben. Für $i\ne j$ und $a=E_{ij}+E_{ji}$, $b=i(E_{ij}-E_{ji})$ gilt $[a,E_{ii}]=E_{ji}-E_{ij}$,
$[b,E_{ii}]=-i(E_{ij}+E_{ji})$ und $[a,b]=2i(E_{jj}-E_{ii})$; diese spannen über $\mathbb R$ die spurfreien
antihermiteschen Matrizen auf. $\square$

Was die Jordan-Algebra nicht enthält, ist die **Zuordnung** $H\mapsto i[H,\cdot]$ zwischen einer Observablen und der
Dynamik, die sie erzeugt. In $\operatorname{Herm}_n(\mathbb R)$ etwa stammen die Derivationen
$[L(x),L(y)]=\frac14[[x,y],\cdot]$ von reellen antisymmetrischen $[x,y]$, die nicht $i$ mal eine reelle symmetrische
Matrix sind. E. M. Alfsen und F. W. Shultz haben gezeigt, dass eine JBW-Algebra (eine unendlichdimensionale Version
einer formal reellen Jordan-Algebra; in endlicher Dimension stimmen die Begriffe überein) genau dann der
selbstadjungierte Teil einer von-Neumann-Algebra ist, wenn sie eine solche „dynamische Korrespondenz“ zulässt, und
haben die Zustandsräume von C\*-Algebren geometrisch charakterisiert (*Geometry of State Spaces of Operator Algebras*,
Birkhäuser, 2003; ohne Beweis angegeben). Eine solche Korrespondenz ist nicht eindeutig: Auf
$\operatorname{Herm}_n(\mathbb C)$ kommen sowohl $H\mapsto i[H,\cdot]$ als auch $H\mapsto-i[H,\cdot]$ in Frage; die
Transposition vertauscht sie, da $T(i[H,T(z)])=-i[H^T,z]$, und sie gehören zu $M_n$ bzw. zu dessen entgegengesetzter
Algebra – die beiden Vorzeichen aus Proposition 13.12.

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

> **Physikalische Bedeutung.** Für ein spinloses Teilchen ist die Transposition (komplexe Konjugation in einer reellen
> Basis) die **Zeitumkehr**. Sie ist eine Symmetrie der Observablen eines einzelnen Systems, aber sie ist antiunitär,
> und auf nur ein Teilsystem angewandt ist sie keine physikalische Operation: Proposition 13.12 zeigt, dass sie den
> Bell-Zustand in einen Operator mit einem negativen Eigenwert verwandelt. Die Quantenmechanik verwendet komplexe
> $*$-Algebren, weil zusammengesetzte Systeme den Kommutator mit seinem Vorzeichen brauchen, nicht nur das
> symmetrisierte Produkt.

## 13.6 Zurück zu den Spurungleichungen

Proposition 11.9 (Exp. 3) besagt $\operatorname{Tr}|YXY|\le\operatorname{Tr}(Y|X|Y)$ für $X=X^*$ und $Y\ge0$ (hier
sind $X,Y$ Operatoren aus $\mathcal M$, keine Pauli-Matrizen). Ihr Beweis verwendet $Y\ge0$ nirgends: $YX_\pm Y=T^*T$
mit $T=X_\pm^{1/2}Y$ braucht nur $Y=Y^*$, und Lemma 11.8 wird auf $YX_+Y$ und $YX_-Y$ angewandt. Sie gilt also für
alle $Y=Y^*$ und lautet in Jordan-Sprache
$$
\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|,
$$
was in jeder euklidischen Jordan-Algebra sinnvoll ist. Wir beweisen sie in $V_m$ mit einem Argument, das nur drei
Eigenschaften verwendet.

**Proposition 13.14 (Exp. 3 in einem Spinfaktor).** Für $x,y\in V_m$ gilt
$\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|$, mit Gleichheit, falls $x\in Q$ oder $-x\in Q$.

*Beweis.* Seien $x_\pm=f_\pm(x)$ mit $f_+(t)=\max(t,0)$, $f_-(t)=\max(-t,0)$ (Satz 13.4(3)): $x_\pm\in Q$,
$x=x_+-x_-$ und $|x|=x_++x_-$. Man setze $A=P(y)x_+$ und $B=P(y)x_-$; sie liegen nach Satz 13.4(5) in $Q$, und
$P(y)x=A-B$. Sei $P(y)x=\mu_1c_1+\mu_2c_2$ eine Spektralzerlegung (Satz 13.4(2); $c_1=c_w$, $c_2=c_{-w}$, also
$c_1+c_2=1$ und $c_1\circ c_2=0$), $c$ die Summe der $c_i$ mit $\mu_i>0$ und
$s=2c-1=\sum_i\operatorname{sgn}'(\mu_i)c_i$ mit $\operatorname{sgn}'(t)=1$ für $t>0$ und $-1$ sonst. Dann ist
$s\circ P(y)x=\sum_i|\mu_i|c_i=|P(y)x|$ (nach Satz 13.4(2,3) hängt das für $\mu_1=\mu_2$ nicht von der Wahl von $w$
ab), und $c,1-c\in Q$ (Idempotente sind Quadrate). Nach der leichten Hälfte der Selbstdualität (Satz 13.4(4):
$\operatorname{tr}(a\circ b)\ge0$ für $a,b\in Q$), angewandt auf $\operatorname{tr}(A\circ(1-c))$ und
$\operatorname{tr}(B\circ c)$, gilt
$$
\operatorname{tr}|P(y)x|=\operatorname{tr}(A\circ c)-\operatorname{tr}(A\circ(1-c))-\operatorname{tr}(B\circ c)+\operatorname{tr}(B\circ(1-c))
\le\operatorname{tr}A+\operatorname{tr}B=\operatorname{tr}P(y)|x|.
$$
Ist $x\in Q$, so ist $x_-=0$, $B=0$ und $|P(y)x|=A$, also sind beide Seiten gleich $\operatorname{tr}A$; ebenso für
$-x\in Q$. $\square$

Der Beweis verwendet nur: (i) $P(y)Q\subseteq Q$ für jedes $y$; (ii) den Spektralsatz: Jedes $z$ ist
$\sum_i\mu_ic_i$ mit orthogonalen Idempotenten $c_i$ der Summe $1$, und $f(z)=\sum_if(\mu_i)c_i$ (sodass $x_\pm$ und
die Idempotente Quadrate sind); (iii) $\operatorname{tr}$ ist linear,
$\operatorname{tr}(\sum_i\mu_ic_i)=\sum_i\mu_i\operatorname{tr}c_i$ mit $\operatorname{tr}c_i\ge0$, und
$\operatorname{tr}(a\circ b)\ge0$ für $a,b\in Q$. In einer euklidischen Jordan-Algebra vom Rang $r$ gilt (ii) mit $r$
primitiven Idempotenten $c_i$ (einem Jordan-Rahmen), die Spur ist $\operatorname{tr}z=\sum_i\mu_i$, die Form
$\operatorname{tr}(x\circ y)$ ist assoziativ, $\operatorname{tr}((x\circ y)\circ z)=\operatorname{tr}(x\circ(y\circ z))$,
und (i) gilt [FK94, Kapitel II–III]; diese Tatsachen beweisen wir nicht. Der letzte Teil von (iii) folgt dann aus (i):
Entwickelt man $P(a^{1/2})=2L(a^{1/2})^2-L(a)$ und verwendet die Assoziativität, so ist
$\operatorname{tr}(a\circ b)=\operatorname{tr}P(a^{1/2})b$, die Spur eines Elements von $Q$, also $\ge0$. Die
Ungleichung gilt also in jeder euklidischen Jordan-Algebra – einschließlich $\operatorname{Herm}_r(\mathbb H)$ und der
Albert-Algebra, die nicht Jordan-isomorph zu selbstadjungierten Teilen komplexer $*$-Algebren sind. In $M_{sa}$
beweist dasselbe Argument, mit $\operatorname{Tr}$ und $P(y)x=yxy$, Proposition 11.9 erneut. In $V_m$ folgt die
Ungleichung auch aus Proposition 11.9 über $\phi$, da $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ und $\phi(|x|)=|\phi(x)|$;
der Code prüft beides:

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

Die Lücke $\operatorname{tr}P(y)|x|-\operatorname{tr}|P(y)x|$ verschwindet für definite $x$ und ist für typische
indefinite $x$ positiv (die ausgegebenen Maxima liegen für diese Gaußschen Stichproben in der Größenordnung
$10^1$–$10^2$).

## 13.7 Womit man vergleichen sollte

Unseres Wissens gibt es keine verbreitete numerische Bibliothek für Jordan-Algebren; die Operationen dieses Kapitels
werden gewöhnlich von Hand mit Matrizen geschrieben. Computeralgebra deckt die algebraische Seite ab: SageMaths
`JordanAlgebra` konstruiert spezielle Jordan-Algebren, Jordan-Algebren symmetrischer Bilinearformen (zu denen die
Spinfaktoren gehören) und, ausgehend von einer Oktonionenalgebra, die 27-dimensionale exzeptionelle Jordan-Algebra;
auch Magma hat Konstruktoren für Jordan-Algebren, einschließlich der exzeptionellen. Das sind exakte algebraische
Werkzeuge, keine gebündelten numerischen. Die Lorentz-Kegel sind Standard in der Optimierung: SOCP-Löser (ECOS, SCS,
Clarabel, MOSEK, zugänglich über die `SOC`-Nebenbedingung von cvxpy) verwenden die Spektralzerlegung von $V_m$ (zur
Projektion auf den Kegel) oder seine quadratische Darstellung (für die Nesterov-Todd-Skalierung von
Innere-Punkte-Verfahren), und semidefinite Löser tun dasselbe für $\operatorname{Herm}_n$. Die Bibliothek fügt
Jordan-Produkte, die quadratische Darstellung und die Transposition als gebündelte, differenzierbare Operatoren auf
$\bigoplus_cM_{k_c}$ hinzu; Spinfaktoren mit ihrem Funktionalkalkül und ihrer Matrixdarstellung; und eine
randomisierte Berechnung der von gegebenen Matrizen erzeugten Jordan-Algebra, verglichen mit der von ihnen erzeugten
$*$-Algebra.

## Aufgaben

1. (**Homogenität des Lorentz-Kegels.**) Sei $y\in V_m$ mit $\det y\ne0$. Man zeige, dass $y^{-1}=(a,-u)/\det y$ die
   Gleichung $y\circ y^{-1}=1$ erfüllt, dass $P(y)P(y^{-1})=\mathrm{id}$, dass $P(y)$ das Innere von $Q$ auf sich
   abbildet und dass $P(x^{1/2})1=x$ für $x$ im Inneren. Man folgere, dass das Innere von $Q$ homogen ist.
2. (**Isotrope Zustände.**) Für $p\in[0,1]$ sei $\rho_p=p\,\omega+(1-p)\frac14$ auf $\mathbb C^2\otimes\mathbb C^2$
   (Proposition 13.12). Man berechne das Spektrum von $\theta(\rho_p)$ und zeige, dass $\rho_p$ für $p>\frac13$
   verschränkt ist.
3. (**Produkte von Randpunkten.**) Seien $|w|=|w'|=1$ in $\mathbb R^m$. Man zeige, dass $(1,w)\circ(1,w')\in Q$ genau
   dann gilt, wenn $w'=w$ oder $w'=-w$: Das Jordan-Produkt zweier Randstrahlen des Lorentz-Kegels liegt fast nie im
   Kegel.
4. (**Reelle Quantenmechanik.**) Man zeige, dass für $x,y\in\operatorname{Herm}_n(\mathbb R)$ die Derivation
   $[L(x),L(y)]$ von $\operatorname{Herm}_n(\mathbb R)$ gleich $z\mapsto[K,z]$ mit reellem antisymmetrischem $K$ ist
   und dass für $n\ge2$ eine solche Abbildung nie von der Form $z\mapsto i[H,z]$ mit
   $H\in\operatorname{Herm}_n(\mathbb R)$ ist, es sei denn, sie ist $0$.

<details>
<summary>Lösungen</summary>

1. $(a,u)\circ(a,-u)=(a^2-|u|^2,0)$, was $y\circ y^{-1}=1$ ergibt. $\phi(y)$ und $\phi(y^{-1})$ kommutieren (beide
   liegen in $\operatorname{span}\{1,u\cdot s\}$), also ist $\phi(y)\phi(y^{-1})=\phi(y)\circ\phi(y^{-1})=1$, und nach
   Satz 13.4(5) ist $\phi(P(y)P(y^{-1})x)=\phi(x)$; $\phi$ ist injektiv. $P(y)$ bildet $Q$ in $Q$ ab, ebenso seine
   Inverse $P(y^{-1})$, also ist $P(y)Q=Q$; als linearer Homöomorphismus bildet es das Innere auf das Innere ab.
   Schließlich ist $P(y)1=2y\circ y-y^2=y^2$, also $P(x^{1/2})1=x$, und $P(x^{1/2})$ ist invertierbar, wenn $x$ im
   Inneren liegt ($\det x^{1/2}=\sqrt{\lambda_+\lambda_-}>0$). Je zwei innere Punkte werden durch
   $P(x'^{1/2})P(x^{1/2})^{-1}$ ineinander überführt.

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

2. $\theta(\rho_p)=\frac p2F+\frac{1-p}4$ nach Proposition 13.12. $F$ hat die Eigenwerte $1$ (dreifach) und $-1$
   (einfach), also besteht das Spektrum aus $\frac{1+p}4$ (dreifach) und $\frac{1-3p}4$, das genau dann negativ ist,
   wenn $p>\frac13$; dann ist $\rho_p$ nach dem PPT-Kriterium nicht separabel. (Für $2\times2$-Systeme gilt auch die
   Umkehrung: M., P. und R. Horodecki, Phys. Lett. A 223 (1996) 1; wir beweisen sie nicht.)

```python
for pp in (0.2, 1 / 3, 0.5):
    rho_p = pp * omega + (1 - pp) * torch.eye(4, dtype=dt) / 4
    lam_min = torch.linalg.eigvalsh(ptranspose(rho_p, 2, 2))[0].item()
    assert abs(lam_min - (1 - 3 * pp) / 4) < 1e-12
```

3. $(1,w)\circ(1,w')=(1+t,\,w+w')$ mit $t=w\cdot w'\in[-1,1]$ und $|w+w'|=\sqrt{2+2t}$. Es liegt genau dann in $Q$,
   wenn $1+t\ge\sqrt{2(1+t)}$, d. h. genau dann, wenn $t=-1$ (das Produkt ist $0$) oder $1+t\ge2$, d. h. $t=1$.

```python
V3c = SpinFactor(3)
wv = torch.nn.functional.normalize(torch.randn(2, 1000, 3, dtype=torch.float64), dim=-1)
prod = V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[1]))
assert not V3c.in_cone(prod).any()
assert V3c.in_cone(V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[0]))).all()
```

4. Nach Proposition 13.13(b) ist $[L(x),L(y)]z=[K,z]$ mit $K=\frac14[x,y]$, und $K^T=\frac14[y^T,x^T]=-K$ ist reell.
   Für reelle $H$ und $z$ ist $i[H,z]$ rein imaginär, während $[K,z]$ reell ist; sie stimmen also nur überein, wenn
   beide verschwinden. In der reellen Quantenmechanik sind die Erzeugenden der Dynamik, $K$, also keine Observablen.

</details>

**Wie es weitergeht.** Dies ist das letzte Kapitel des Kurses. Die Standardreferenz für die endlichdimensionale
Theorie ist [FK94]; K. McCrimmon, *A Taste of Jordan Algebras* (Springer 2004), behandelt die Algebra und die
Geschichte, H. Hanche-Olsen und E. Størmer, *Jordan Operator Algebras* (Pitman 1984), sowie Alfsen–Shultz (2003) die
unendlichdimensionale Theorie und ihre Beziehung zu C\*-Algebren. Zur Rolle zusammengesetzter Systeme bei
Rekonstruktionen der Quantentheorie aus operationellen Axiomen siehe die oben zitierten Arbeiten von Barnum und Wilce.
Zu Spurungleichungen über Kapitel 11 hinaus: E. Carlen, *Trace inequalities and quantum entropy* (Contemp. Math. 529,
2010). Die quadratische Darstellung und der Lorentz-Kegel tauchen in der konischen Optimierung wieder auf:
Alizadeh–Goldfarb (2003) und die dort besprochenen Jordan-algebraischen Innere-Punkte-Verfahren.
