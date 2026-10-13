[English](../E_jordan_classification.md) · [Español](../es/E_jordan_classification.md) · [Français](../fr/E_jordan_classification.md) · **Deutsch** · [中文](../zh/E_jordan_classification.md) · [日本語](../ja/E_jordan_classification.md) · [Русский](../ru/E_jordan_classification.md)

# Anhang E. Klassifikation der euklidischen Jordan-Algebren

Kapitel 13 hat die Jordan–von-Neumann–Wigner-Klassifikation der endlichdimensionalen formal reellen Jordan-Algebren
verwendet. Dieser Anhang beweist sie vollständig, allein aus der Jordan-Identität und endlichdimensionaler linearer
Algebra. Der Weg ist der klassische (J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, Kapitel II–V; K. McCrimmon,
*A Taste of Jordan Algebras*): Potenzassoziativität und ein Spektrallemma liefern Idempotente und eine assoziative,
positive Spurform; die Peirce-Zerlegung bezüglich eines Jordan-Rahmens zerlegt die Algebra in Stücke $J_{ij}$, deren
Produkte strengen Regeln gehorchen; in einer einfachen Algebra vom Rang $r\ge3$ machen diese Regeln den Peirce-Raum
$J_{13}$ zu einer Kompositionsalgebra, die nach dem Satz von Hurwitz $\mathbb R,\mathbb C,\mathbb H$ oder $\mathbb O$
ist, und die ganze Algebra ist dann deren $\operatorname{Herm}_r$. Abschnitt E.6 beweist außerdem den Satz von
Frobenius über reelle assoziative Divisionsalgebren, und Abschnitt E.8 beweist, dass die Albert-Algebra
$\operatorname{Herm}_3(\mathbb O)$ eine Jordan-Algebra ist, die sich in keine assoziative Algebra einbetten lässt.

**Sie lernen**, warum Potenzassoziativität, Positivität der Spurform und die Existenz eines Einselements allein aus der
formalen Realität folgen; die Peirce-Regeln und wie sie die Koordinaten einer einfachen Algebra zwingen, sich wie Zahlen
zu multiplizieren; den Beweis des Satzes von Hurwitz durch Cayley–Dickson-Verdopplung; und warum die Oktonionen nur im
Rang $3$ eine neue einfache Algebra liefern.

**Konventionen.** Alle Algebren sind endlichdimensionale reelle Vektorräume mit einem bilinearen Produkt. Eine
**Jordan-Algebra** ist wie in Abschnitt 13.1 definiert: $x\circ y=y\circ x$ und $(x^2\circ y)\circ x=x^2\circ(y\circ x)$,
$x^2=x\circ x$; ein Einselement wird nicht vorausgesetzt. $L(x)y=x\circ y$, $[A,B]=AB-BA$, und $\operatorname{Tr}$ ist
die Spur einer linearen Abbildung. Eine **Unteralgebra** ist ein unter $\circ$ abgeschlossener Unterraum, ein **Ideal**
ein Unterraum $I$ mit $J\circ I\subseteq I$, und $J$ ist **einfach**, wenn $J\ne0$ und $0$ und $J$ die einzigen Ideale
sind. Ein **Idempotent** ist ein $c$ mit $c^2=c$; Idempotente $c,c'$ heißen **orthogonal**, wenn $c\circ c'=0$. $J$ ist
**formal reell**, wenn aus $\sum_ix_i^2=0$ folgt, dass alle $x_i=0$ sind, und **euklidisch**, wenn sie ein Einselement
$1$ und ein Skalarprodukt mit $\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$ besitzt. Polynomiale Identitäten
werden wie folgt **linearisiert**: Hat $p(t)=\sum_kp_kt^k$ vektorwertige Koeffizienten und ist $p(t)=0$ für alle
reellen $t$, so sind alle $p_k=0$ (man werte an $\deg p+1$ verschiedenen Stellen aus und invertiere die
Vandermonde-Matrix).

## E.1 Linearisierte Jordan-Identitäten und Potenzassoziativität

**Lemma E.1.** In jeder Jordan-Algebra gilt für alle $x,y,z,w$:

1. $[L(x),L(x^2)]=0$ und $[L(x^2),L(z)]+2[L(x\circ z),L(x)]=0$;
2. $L(x^2\circ y)=2L(x\circ y)L(x)+L(x^2)L(y)-2L(x)L(y)L(x)$;
3. $L((x\circ w)\circ y)=L(w\circ y)L(x)+L(x\circ y)L(w)+L(x\circ w)L(y)-L(w)L(y)L(x)-L(x)L(y)L(w)$;
4. $L((x\circ w)\circ y)-L(x\circ(w\circ y))=[L(w),[L(x),L(y)]]$; insbesondere
   $\operatorname{Tr}L((x\circ w)\circ y)=\operatorname{Tr}L(x\circ(w\circ y))$.

*Beweis.* (1) Die Jordan-Identität besagt $L(x)L(x^2)y=L(x^2)L(x)y$. Man ersetze $x$ durch $x+tz$; der Koeffizient von
$t$ ist $[L(z),L(x^2)]+2[L(x),L(x\circ z)]=0$. (2) Man ersetze $x$ durch $x+tz$ in der Jordan-Identität selbst; der
Koeffizient von $t$ ist
$$
2((x\circ z)\circ y)\circ x+(x^2\circ y)\circ z=2(x\circ z)\circ(y\circ x)+x^2\circ(y\circ z).
$$
Als lineare Abbildungen von $z$ gelesen, sind die vier Terme $2L(x)L(y)L(x)z$, $L(x^2\circ y)z$, $2L(x\circ y)L(x)z$
und $L(x^2)L(y)z$. (3) Man ersetze $x$ durch $x+sw$ in (2) und nehme den Koeffizienten von $s$, geteilt durch $2$.
(4) Man schreibe (3) einmal unverändert und einmal mit vertauschten $x$ und $y$. Die ersten drei Terme der rechten
Seiten stimmen überein, und $(y\circ w)\circ x=x\circ(w\circ y)$, also ist die Differenz
$-L(w)L(y)L(x)-L(x)L(y)L(w)+L(w)L(x)L(y)+L(y)L(x)L(w)=[L(w),[L(x),L(y)]]$, ein Kommutator, dessen Spur $0$ ist.
$\square$

**Lemma E.2 (Potenzassoziativität).** Man setze $x^1=x$ und $x^{k+1}=x\circ x^k$. Dann gilt $x^p\circ x^q=x^{p+q}$ für
alle $p,q\ge1$, und die Operatoren $L(x^k)$ vertauschen paarweise. Folglich gilt für reelle Polynome $f,g$ ohne
konstanten Term und $f(x):=\sum_kf_kx^k$, dass $(fg)(x)=f(x)\circ g(x)$, und $\mathbb R[x]=\{f(x)\}$ ist eine
assoziative, kommutative Unteralgebra.

*Beweis.* $L(x)$ und $L(x^2)$ vertauschen nach Lemma E.1(1); sei $\mathcal A$ die von ihnen erzeugte kommutative
Operatoralgebra. Wir zeigen durch Induktion nach $n\ge2$: $x^p\circ x^q=x^{p+q}$, sobald $p+q\le n$, und
$L(x^k)\in\mathcal A$ für $k\le n$. Für $n=2$ ist das die Definition. Es gelte die Aussage für $n$, und sei $p+q=n+1$.
Ist $q=1$, so gilt $x^p\circ x=x^{n+1}$ wegen der Kommutativität. Ist $q\ge2$, so gilt
$x^p\circ x^q=L(x^p)L(x)x^{q-1}=L(x)L(x^p)x^{q-1}=L(x)x^{n}=x^{n+1}$. Weiter liefert Lemma E.1(2) mit $y=x^{n-1}$ und
$x^2\circ x^{n-1}=x^{n+1}$, dass $L(x^{n+1})=2L(x^n)L(x)+L(x^2)L(x^{n-1})-2L(x)^2L(x^{n-1})\in\mathcal A$. Der Rest
folgt aus der Bilinearität. $\square$

## E.2 Idempotente und die Peirce-Zerlegung

**Lemma E.3 (Spektrallemma).** Sei $J$ formal reell und $x\in J$. Es gibt $k\ge0$, verschiedene reelle Zahlen
$\lambda_1,\dots,\lambda_k$ ungleich null und von null verschiedene, paarweise orthogonale Idempotente
$c_1,\dots,c_k\in\mathbb R[x]$ mit
$$
x=\sum_{i=1}^k\lambda_ic_i,\qquad x^n=\sum_{i=1}^k\lambda_i^nc_i\quad(n\ge1),\qquad\mathbb R[x]=\operatorname{span}\{c_1,\dots,c_k\}.
$$
Insbesondere enthält jede von null verschiedene formal reelle Jordan-Algebra ein von null verschiedenes Idempotent.

*Beweis.* Formale Realität mit einem einzigen Term besagt, dass aus $z^2=0$ folgt $z=0$. Sei
$K=\{q\in\mathbb R[t]:(tq)(x)=0\}$. Nach Lemma E.2 ist $K$ ein Ideal von $\mathbb R[t]$: Ist $q\in K$ und $s=s_0+s_1$
mit $s_1(0)=0$, so ist $(tqs)(x)=s_0(tq)(x)+(tq)(x)\circ s_1(x)=0$. Da $x,x^2,\dots,x^{\dim J+1}$ linear abhängig sind,
ist $K\ne0$, also $K=\mu\,\mathbb R[t]$ mit einem normierten $\mu$. (i) $\mu(0)\neq0$: Ist $\mu=th$, so erfüllt
$w=(th)(x)$ die Gleichung $w^2=(t\cdot\mu h)(x)=0$, da $\mu h\in K$; also $w=0$, $h\in K$ und $\mu\mid h$, was aus
Gradgründen unmöglich ist. (ii) $\mu$ hat keinen mehrfachen Faktor: Ist $\mu=g^2h$ mit $\deg g\ge1$, so erfüllt
$w=(tgh)(x)$ die Gleichung $w^2=(t\cdot\mu\cdot th)(x)=0$, also $gh\in K$, aus Gradgründen unmöglich. (iii) $\mu$ hat
keine nichtreelle Nullstelle: Ist $\mu=((t-\alpha)^2+\beta^2)h$ mit $\beta\ne0$, so erfüllen $a=(t(t-\alpha)h)(x)$ und
$b=(\beta th)(x)$ die Gleichung $a^2+b^2=(t\cdot\mu\cdot th)(x)=0$, also $b=0$ und $h\in K$, unmöglich. Daher ist
$\mu=\prod_{i=1}^k(t-\lambda_i)$ mit verschiedenen reellen $\lambda_i$ ungleich null. Man setze
$$
p_i(t)=\frac t{\lambda_i}\prod_{j\ne i}\frac{t-\lambda_j}{\lambda_i-\lambda_j},\qquad c_i=p_i(x).
$$
Das Polynom $p_ip_l-\delta_{il}p_i$ ist $t$ mal ein Polynom, das in jedem $\lambda_j$ verschwindet, also $t$ mal ein
Element von $K$; somit $c_i\circ c_l=\delta_{il}c_i$. Ebenso ist $t-\sum_i\lambda_ip_i=t\big(1-\sum_i\lambda_ip_i/t\big)$,
und die Klammer verschwindet in jedem $\lambda_j$, also $x=\sum_i\lambda_ic_i$. Wäre $c_i=0$, so wäre $p_i/t\in K$, vom
Grad $k-1<\deg\mu$; also $c_i\ne0$. Die Formel für $x^n$ folgt aus der Orthogonalität und zeigt
$\mathbb R[x]\subseteq\operatorname{span}\{c_i\}$.
$\square$

**Proposition E.4 (Peirce-Zerlegung bezüglich eines Idempotents).** Sei $c$ ein Idempotent einer Jordan-Algebra $J$.

1. $L(c)(2L(c)-1)(L(c)-1)=0$. Folglich ist $J=J_1(c)\oplus J_{1/2}(c)\oplus J_0(c)$ mit
   $J_\alpha(c)=\{x:c\circ x=\alpha x\}$, und $c\in J_1(c)$.
2. $[L(c),L(z)]=0$ für $z\in J_0(c)\cup J_1(c)$, und $[L(c),L(y^2)]=0$ für $y\in J_{1/2}(c)$.
3. $J_1\circ J_1\subseteq J_1$, $J_0\circ J_0\subseteq J_0$, $J_1\circ J_0=0$, $(J_1+J_0)\circ J_{1/2}\subseteq J_{1/2}$
   und $J_{1/2}\circ J_{1/2}\subseteq J_1+J_0$ (alles bezüglich $c$).
4. Ist $c'$ ein zu $c$ orthogonales Idempotent, so ist $c+c'$ ein Idempotent, $[L(c),L(c')]=0$ und
   $J_1(c')\subseteq J_0(c)$.

*Beweis.* (1) Lemma E.1(2) mit $x=y=c$ liefert $L(c)=3L(c)^2-2L(c)^3$. Eine lineare Abbildung, die von einem Polynom
mit einfachen Nullstellen annulliert wird, ist diagonalisierbar mit Eigenwerten unter diesen Nullstellen. (2) Lemma
E.1(1) mit $x=c$ liefert $[L(c),L(z)]=2[L(c),L(c\circ z)]$, d. h. $(1-2\alpha)[L(c),L(z)]=0$ für $z\in J_\alpha(c)$.
Lemma E.1(1) mit $x=y$ und $z=c$ liefert
$[L(y^2),L(c)]+2[L(y\circ c),L(y)]=[L(y^2),L(c)]+[L(y),L(y)]=[L(y^2),L(c)]=0$. (3) Nach (2) bildet $L(z)$ jeden
Eigenraum von $L(c)$ in sich ab, wenn $z\in J_0\cup J_1$. Das liefert die erste, zweite und vierte Regel sowie
$J_1\circ J_0\subseteq J_1\cap J_0=0$. Für $y\in J_{1/2}$ vertauscht $L(y^2)$ mit $L(c)$, also ist
$c\circ y^2=y^2\circ c$ ein Eigenvektor von $L(c)$ zum Eigenwert $1$, so wie $c$; schreibt man $y^2=a_1+a_{1/2}+a_0$,
so liegt $c\circ y^2=a_1+\frac12a_{1/2}$ in $J_1$, also $a_{1/2}=0$. Polarisieren ($y\to y+y'$) liefert die letzte
Regel. (4) $(c+c')^2=c+2c\circ c'+c'=c+c'$, und $c'\in J_0(c)$, also vertauschen $L(c)$ und $L(c')$ nach (2) und sind
simultan diagonalisierbar; insbesondere wird $J_1(c')$ von gemeinsamen Eigenvektoren aufgespannt. Auf einem gemeinsamen
Eigenvektor mit den Eigenwerten $\alpha$ für $L(c)$ und $1$ für $L(c')$ hat $L(c+c')$ den Eigenwert
$\alpha+1\in\{0,\frac12,1\}$, also $\alpha=0$. $\square$

## E.3 Die Spurform; euklidisch = formal reell

**Satz E.5.** Sei $J$ formal reell und $\tau(x,y)=\operatorname{Tr}L(x\circ y)$.

1. $\tau$ ist symmetrisch, bilinear und **assoziativ**: $\tau(x\circ y,z)=\tau(x,y\circ z)$.
2. $\tau$ ist positiv definit; für jedes von null verschiedene Idempotent $c$ ist $\tau(c,c)=\operatorname{Tr}L(c)\ge1$.
3. $J$ besitzt ein Einselement.
4. $J$ ist euklidisch, mit dem Skalarprodukt $\tau$. Umgekehrt ist jede euklidische Jordan-Algebra formal reell.

*Beweis.* (1) Lemma E.1(4). (2) Nach Proposition E.4(1) ist
$\operatorname{Tr}L(c)=\dim J_1(c)+\frac12\dim J_{1/2}(c)\ge1$, da $0\ne c\in J_1(c)$. Für $x=\sum_i\lambda_ic_i$ wie in
Lemma E.3 ist $x^2=\sum_i\lambda_i^2c_i$ und $\tau(x,x)=\operatorname{Tr}L(x^2)=\sum_i\lambda_i^2\operatorname{Tr}L(c_i)$,
was $>0$ ist, außer für $k=0$, d. h. $x=0$. (3) Sei $J\neq0$ und $e$ ein von null verschiedenes Idempotent (Lemma E.3)
mit maximalem $\dim J_1(e)$. Ist $J_0(e)\ne0$, so ist es eine formal reelle Unteralgebra (Proposition E.4(3)) und
enthält ein von null verschiedenes Idempotent $f$; $f$ ist orthogonal zu $e$, also ist $e+f$ ein Idempotent,
$J_1(e)\subseteq J_1(e+f)$ (da $f\circ J_1(e)\subseteq J_0\circ J_1=0$) und $f\in J_1(e+f)\setminus J_1(e)$, im
Widerspruch zur Maximalität. Also ist $J_0(e)=0$. Für $y\in J_{1/2}(e)$ ist dann $y^2\in J_1(e)$, und die
Assoziativität liefert
$\frac12\tau(y,y)=\tau(e\circ y,y)=\tau(e,y^2)=\operatorname{Tr}L(e\circ y^2)=\operatorname{Tr}L(y^2)=\tau(y,y)$, also
$y=0$. Daher ist $J=J_1(e)$ und $e$ ein Einselement. (4) Nach (1)–(3) ist $\tau(x\circ y,z)=\tau(y,x\circ z)$.
Umgekehrt: Ist $J$ euklidisch und $\sum_ix_i^2=0$, so ist $0=\sum_i\langle x_i\circ x_i,1\rangle=\sum_i\langle x_i,x_i\rangle$.
$\square$

„Formal reell“ und „euklidisch“ beschreiben also dieselbe Klasse endlichdimensionaler Algebren. Von nun an ist $J$
formal reell, $1$ sein Einselement und $\tau$ seine Spurform. Der Code prüft Lemma E.1(4), Lemma E.3 (die Idempotente
$p_i(x)$, nur mit Jordan-Potenzen berechnet), Proposition E.4 und Satz E.5 in $\operatorname{Herm}_3(\mathbb C)$, gegeben
durch Strukturkonstanten; das zu ihrer Konstruktion verwendete Cayley–Dickson-Produkt `cd` ist das aus Abschnitt E.6.

```python
import math
import torch

torch.manual_seed(0)
dt = torch.float64

def conj(x):                                      # Cayley-Dickson conjugation (a, b)* = (a*, -b) on R^(2^k)
    if x.shape[-1] == 1:
        return x
    h = x.shape[-1] // 2
    return torch.cat([conj(x[..., :h]), -x[..., h:]], -1)

def cd(x, y):                                     # (a, b)(c, d) = (ac - d*b, da + bc*): R, C, H, O for 1, 2, 4, 8
    if x.shape[-1] == 1:
        return x * y
    h = x.shape[-1] // 2
    a, b, c, d = x[..., :h], x[..., h:], y[..., :h], y[..., h:]
    return torch.cat([cd(a, c) - cd(conj(d), b), cd(d, a) + cd(b, conj(c))], -1)

def herm(n, d):
    """Herm_n(D), D = R, C, H, O for d = 1, 2, 4, 8; coordinates: diagonal, then the entries a_ij, i < j."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    def to_mat(v):
        X = torch.zeros(v.shape[:-1] + (n, n, d), dtype=dt)
        X[..., range(n), range(n), 0] = v[..., :n]
        for p, (i, j) in enumerate(pairs):
            X[..., i, j, :] = v[..., n + d * p:n + d * (p + 1)]
            X[..., j, i, :] = conj(X[..., i, j, :])
        return X
    def from_mat(X):
        return torch.cat([X[..., range(n), range(n), 0]] + [X[..., i, j, :] for i, j in pairs], -1)
    def matmul(X, Y):
        return torch.stack([torch.stack([sum(cd(X[..., i, k, :], Y[..., k, j, :]) for k in range(n))
                                         for j in range(n)], -2) for i in range(n)], -3)
    def prod(x, y):
        X, Y = to_mat(x), to_mat(y)
        return from_mat(0.5 * (matmul(X, Y) + matmul(Y, X)))
    return n + d * len(pairs), prod

def spin(m):                                      # V_m: (a, u) o (b, v) = (ab + u.v, av + bu)
    def prod(x, y):
        return torch.cat([x[..., :1] * y[..., :1] + (x[..., 1:] * y[..., 1:]).sum(-1, keepdim=True),
                          x[..., :1] * y[..., 1:] + y[..., :1] * x[..., 1:]], -1)
    return m + 1, prod

def structure(dim, prod):                         # C[a, b] = e_a o e_b
    E = torch.eye(dim, dtype=dt)
    return prod(E[:, None, :].expand(dim, dim, dim), E[None, :, :].expand(dim, dim, dim))

mul = lambda C, x, y: torch.einsum('...a,...b,abc->...c', x, y, C)
Lop = lambda C, x: torch.einsum('...a,abc->...cb', x, C)          # matrix of L(x)
comm = lambda A, B: A @ B - B @ A

C3 = structure(*herm(3, 2))                       # Herm_3(C), dimension 9
x, y, w = torch.randn(3, 9, dtype=dt)
lhs = Lop(C3, mul(C3, mul(C3, x, w), y)) - Lop(C3, mul(C3, x, mul(C3, w, y)))
assert torch.allclose(lhs, comm(Lop(C3, w), comm(Lop(C3, x), Lop(C3, y))))     # Lemma E.1(4)

def jpow(C, x, k):                                # Jordan powers x^1 = x, x^(k+1) = x o x^k
    p = x
    for _ in range(k - 1):
        p = mul(C, x, p)
    return p

def coords(Z):                                    # complex Hermitian 3 x 3 -> coordinates of Herm_3(C)
    return torch.cat([Z.diagonal().real] + [torch.stack([Z[i, j].real, Z[i, j].imag]) for i, j in ((0, 1), (0, 2), (1, 2))])

G = torch.randn(3, 3, dtype=torch.complex128)
Z = G + G.mH
lam = torch.linalg.eigvalsh(Z)                    # distinct and non-zero (almost surely)
xz = coords(Z)
powers = torch.stack([jpow(C3, xz, k) for k in range(1, 4)])  # x, x^2, x^3
cs = []
for i in range(3):                                # Lagrange idempotents p_i(x) of Lemma E.3
    others = [lam[j] for j in range(3) if j != i]
    coef = torch.tensor([others[0] * others[1], -(others[0] + others[1]), 1.0], dtype=dt)
    den = lam[i] * (lam[i] - others[0]) * (lam[i] - others[1])
    cs.append((coef[:, None] * powers).sum(0) / den)
for i in range(3):
    for j in range(3):
        assert torch.allclose(mul(C3, cs[i], cs[j]), cs[i] if i == j else torch.zeros(9, dtype=dt), atol=1e-10)
assert torch.allclose(sum(l * c for l, c in zip(lam, cs)), xz)

c = cs[0] + cs[1]                                 # the idempotent c_1 + c_2: Proposition E.4
Lc = Lop(C3, c)
assert torch.allclose(Lc @ (2 * Lc - torch.eye(9, dtype=dt)) @ (Lc - torch.eye(9, dtype=dt)), torch.zeros(9, 9, dtype=dt), atol=1e-10)
ev = torch.linalg.eigvals(Lc).real.sort().values
print("eigenvalues of L(c):", ev.numpy().round(6))
tau = lambda C: torch.einsum('abc,c->ab', C, torch.einsum('cbb->c', Lop(C, torch.eye(C.shape[0], dtype=dt))))  # Tr L(e_a o e_b)
T3 = tau(C3)
assert torch.linalg.eigvalsh(T3).min() > 0                                     # Theorem E.5(2)
assert torch.allclose(mul(C3, x, y) @ T3 @ w, x @ T3 @ mul(C3, y, w))         # Theorem E.5(1)
```

## E.4 Jordan-Rahmen und der Spektralsatz

Ein von null verschiedenes Idempotent heißt **primitiv**, wenn es nicht Summe zweier von null verschiedener orthogonaler
Idempotente ist. Ein **Jordan-Rahmen** ist eine Familie $(c_1,\dots,c_r)$ paarweise orthogonaler primitiver Idempotente
mit $\sum_ic_i=1$.

**Lemma E.6.** (1) Ist $c$ primitiv, so ist $J_1(c)=\mathbb Rc$. (2) Paarweise orthogonale, von null verschiedene
Idempotente sind linear unabhängig; jedes von null verschiedene Idempotent ist eine Summe paarweise orthogonaler
primitiver Idempotente; und $J$ besitzt einen Jordan-Rahmen.

*Beweis.* (1) $J_1(c)$ ist eine formal reelle Unteralgebra mit Einselement $c$ (Proposition E.4(3)). Für $x\in J_1(c)$
schreibe man $x=\sum_i\lambda_ic_i$ nach Lemma E.3, mit $c_i\in\mathbb R[x]\subseteq J_1(c)$. Dann ist $c-c_i$ ein zu
$c_i$ orthogonales Idempotent ($(c-c_i)^2=c-2c_i+c_i$ und $c_i\circ(c-c_i)=0$), also $c-c_i=0$ wegen der Primitivität.
Verschiedene $c_i$ sind orthogonal, also gibt es höchstens eines, und $x\in\mathbb Rc$. (2) Ist $\sum_ia_if_i=0$, so
liefert das Produkt mit $f_j$, dass $a_jf_j=0$. Ein Idempotent $e$ ist also Summe einer größtmöglichen Anzahl paarweise
orthogonaler, von null verschiedener Idempotente. Wäre eines davon, $f$, gleich $g+g'$ mit orthogonalen, von null
verschiedenen Idempotenten $g,g'$, so wäre $g=f\circ g\in J_1(f)\subseteq J_0(f')$ für jedes andere Mitglied $f'$
(Proposition E.4(4)), ebenso $g'$, und $f$ könnte durch $g,g'$ ersetzt werden. Man nehme $e=1$. $\square$

**Satz E.7 (Spektralsatz).** Zu jedem $x\in J$ gibt es einen Jordan-Rahmen $(c_1,\dots,c_r)$ und reelle Zahlen
$\lambda_1,\dots,\lambda_r$ (nicht notwendig verschieden) mit $x=\sum_i\lambda_ic_i$.

*Beweis.* Man schreibe $x=\sum_{i=1}^k\lambda_ie_i$ nach Lemma E.3 und füge $e_0=1-\sum_ie_i$, $\lambda_0=0$ hinzu;
$e_0$ ist ein zu allen $e_i$ orthogonales Idempotent. Man zerlege jedes von null verschiedene $e_i$ in paarweise
orthogonale primitive Idempotente (Lemma E.6(2)). Ein Stück $f$ von $e_i$ liegt in $J_1(e_i)\subseteq J_0(e_l)$ für
$l\ne i$, und ein Stück von $e_l$ liegt in $J_1(e_l)$; wegen $J_0(e_l)\circ J_1(e_l)=0$ sind alle Stücke paarweise
orthogonal. Ihre Summe ist $\sum_ie_i=1$. $\square$

Dass alle Jordan-Rahmen von $J$ dieselbe Anzahl von Elementen haben, wird in Korollar E.28 bewiesen.

**Satz E.8 (Peirce-Zerlegung bezüglich eines Rahmens).** Sei $(c_1,\dots,c_r)$ ein Jordan-Rahmen, $t_i=\tau(c_i,c_i)$
und, für $i\ne j$, $J_{ij}=J_{ji}=J_{1/2}(c_i)\cap J_{1/2}(c_j)$.

1. Die $L(c_k)$ vertauschen, $J=\bigoplus_i\mathbb Rc_i\oplus\bigoplus_{i<j}J_{ij}$, und
   $c_k\circ x=\frac12(\delta_{ki}+\delta_{kj})x$ für $x\in J_{ij}$.
2. Für verschiedene $i,j,k$ gilt $J_{ij}\circ J_{jk}\subseteq J_{ik}$; und $J_{ij}\circ J_{kl}=0$, falls $\{i,j\}\cap\{k,l\}=\emptyset$.
3. Es gibt eine positiv definite quadratische Form $N$ auf $J_{ij}$ mit $x^2=N(x)(c_i+c_j)$; mit ihrer Polarform
   $\langle\cdot,\cdot\rangle$ gilt $x\circ y=\langle x,y\rangle(c_i+c_j)$ und $\tau(x,y)=(t_i+t_j)\langle x,y\rangle$
   für $x,y\in J_{ij}$. Ist $J_{ij}\ne0$, so ist $t_i=t_j$.

*Beweis.* (1) Die $L(c_k)$ vertauschen nach Proposition E.4(4), sind diagonalisierbar und summieren sich zu
$L(1)=\mathrm{id}$. Ein gemeinsamer Eigenvektor hat Eigenwerte $\alpha_k\in\{0,\frac12,1\}$ mit $\sum_k\alpha_k=1$:
Entweder ist ein $\alpha_i=1$ und die übrigen verschwinden, und dann liegt der Vektor in $J_1(c_i)=\mathbb Rc_i$
(Lemma E.6(1)); oder es ist $\alpha_i=\alpha_j=\frac12$ für ein Paar, und die übrigen verschwinden. Also ist $J_{ij}$
genau der gemeinsame Eigenraum der zweiten Art, und $J$ ist die direkte Summe der gemeinsamen Eigenräume. (2) Sei
$x\in J_{ij}$ und $y\in J_{jk}$. Bezüglich $c_i$ ist $x\in J_{1/2}$ und $y\in J_0$, also $x\circ y\in J_{1/2}(c_i)$;
bezüglich $c_k$ ebenso $x\circ y\in J_{1/2}(c_k)$. Daher ist $x\circ y\in J_{ik}$. Ist stattdessen $y\in J_{kl}$, so
liegt $x\circ y$ in $J_{1/2}(c_p)$ für $p=i,j,k,l$; einen gemeinsamen Eigenvektor mit vier Eigenwerten $\frac12$ gibt
es nicht, also $x\circ y=0$. (3) Für $x,y\in J_{ij}$ ist $x\circ y\in J_0(c_p)+J_1(c_p)$ für $p=i,j$ und
$x\circ y\in J_0(c_k)$ für $k\notin\{i,j\}$ (Proposition E.4(3)); nach (1) ist $x\circ y\in\mathbb Rc_i+\mathbb Rc_j$.
Man schreibe $x^2=\alpha c_i+\beta c_j$. Die Jordan-Identität mit $y=c_i$ liefert $(x^2\circ c_i)\circ x=\frac12\alpha x$
und $x^2\circ(c_i\circ x)=\frac12x^2\circ x=\frac14(\alpha+\beta)x$, also $\alpha=\beta$, falls $x\ne0$; man setze
$N(x)=\alpha$. Da $\tau(c_i,1)=\operatorname{Tr}L(c_i)=t_i$, ist $\tau(x,x)=\tau(x^2,1)=N(x)(t_i+t_j)$; also ist $N$
eine quadratische Form, positiv definit nach Satz E.5(2), und Polarisieren liefert die Formeln für $x\circ y$ und
$\tau$. Schließlich ist $N(x)t_i=\tau(x^2,c_i)=\tau(x,x\circ c_i)=\frac12\tau(x,x)$, und dasselbe gilt mit $j$; also
$t_i=t_j$, falls es ein $x\ne0$ gibt. $\square$

## E.5 Einfache Algebren und Rang zwei

**Satz E.9.** $J$ ist eine direkte Summe $J^{(1)}\oplus\dots\oplus J^{(s)}$ einfacher Ideale mit
$J^{(p)}\circ J^{(q)}=0$ für $p\ne q$. Jedes $J^{(p)}$ ist eine formal reelle Jordan-Algebra.

*Beweis.* Ist $I$ ein Ideal, so ist auch sein orthogonales Komplement $I^\perp$ bezüglich $\tau$ eines:
$\tau(J\circ I^\perp,I)=\tau(I^\perp,J\circ I)=0$. Daher ist $J=I\oplus I^\perp$ und
$I\circ I^\perp\subseteq I\cap I^\perp=0$, also ist $J$ die direkte Summe der Algebren $I$ und $I^\perp$, und ein Ideal
eines Summanden ist ein Ideal von $J$. Sei $I$ ein von null verschiedenes Ideal kleinster Dimension; dann ist $I$
einfach, und man schließt durch Induktion nach der Dimension, angewandt auf die formal reelle Algebra $I^\perp$.
$\square$

**Lemma E.10.** Sei $(c_i)$ ein Jordan-Rahmen, $i,j,k,l$ verschieden, $u\in J_{ij}$, $w\in J_{jk}$ und $z\in J_{kl}$.
Dann gilt (1) $u\circ(u\circ w)=\frac14N(u)w$ und (2) $(u\circ w)\circ z=u\circ(w\circ z)$.

*Beweis.* (1) Man wende Lemma E.1(1) mit $x=u$ auf $c_j$ an. Mit $L(u^2)=N(u)L(c_i+c_j)$, $c_j\circ w=\frac12w$,
$c_j\circ u=\frac12u$ und $c_j\circ(u\circ w)=0$ (da $u\circ w\in J_{ik}$) ist
$[L(u^2),L(w)]c_j=\frac14N(u)w-\frac12N(u)w$ und $2[L(u\circ w),L(u)]c_j=(u\circ w)\circ u$; ihre Summe verschwindet.
(2) Man wende Lemma E.1(3), mit $(u,z,w)$ anstelle von $(x,w,y)$, auf $c_j$ an. Hier ist $u\circ z=0$ (Satz E.8(2)),
$c_j\circ u=\frac12u$ und $c_j\circ z=0$, also $0=\frac12(z\circ w)\circ u-\frac12z\circ(w\circ u)$. $\square$

**Proposition E.11.** Sei $(c_1,\dots,c_r)$ ein Jordan-Rahmen von $J$, und man verbinde $i\neq j$ durch eine Kante,
wenn $J_{ij}\ne0$.

1. $J$ ist genau dann einfach, wenn dieser Graph zusammenhängend ist.
2. Ist $J$ einfach und $r\ge2$, so haben alle $J_{ij}$, $i\ne j$, dieselbe Dimension $d\ge1$, alle $t_i$ sind gleich
   $t=1+\frac{(r-1)d}2$, und $\dim J=r+d\,\frac{r(r-1)}2=rt$.

*Beweis.* (1) Ist der Graph nicht zusammenhängend, so sei $S$ eine Zusammenhangskomponente und $I$ der Spann der $c_i$
und $J_{ij}$ mit $i,j\in S$, $I'$ dasselbe für das Komplement. Dann ist $J=I\oplus I'$ (alle $J_{ij}$ mit $i\in S$,
$j\notin S$ verschwinden), $I\circ I\subseteq I$ und $I\circ I'=0$ nach Satz E.8; also ist $I$ ein Ideal, $\ne0,J$.
Umgekehrt sei der Graph zusammenhängend und $I\ne0$ ein Ideal. $I$ ist invariant unter den vertauschenden,
diagonalisierbaren $L(c_k)$, also die Summe seiner Durchschnitte mit den gemeinsamen Eigenräumen. Ist
$0\ne x\in I\cap J_{ij}$, so ist $c_i\circ x^2=N(x)c_i\in I$; also liegt ein $c_i$ in $I$. Ist $c_i\in I$ und
$0\ne x\in J_{ij}$, so ist $x=2c_i\circ x\in I$ und $c_j=c_j\circ x^2/N(x)\in I$. Wegen des Zusammenhangs liegen alle
$c_j\in I$, also $1\in I$ und $I=J$. (2) Sei $u\in J_{ij}$ mit $N(u)=1$ und $k\notin\{i,j\}$. Nach Satz E.8(2) bildet
$2L(u)$ den Raum $J_{jk}$ nach $J_{ik}$ und $J_{ik}$ nach $J_{jk}$ ab, und nach Lemma E.10(1) (auch mit vertauschten
$i,j$) ist sein Quadrat auf beiden die Identität. Also ist $\dim J_{ik}=\dim J_{jk}$, sobald $J_{ij}\ne0$. Somit
erzwingen Kanten $ij$, $jk$ eine Kante $ik$; ein zusammenhängender Graph ist dann vollständig, und indem man jeweils
einen Index ändert, sieht man, dass alle $\dim J_{ij}$ gleich sind. Schließlich hat $L(c_i)$ den Eigenwert $1$ auf
$\mathbb Rc_i$, $\frac12$ auf den $r-1$ Räumen $J_{ij}$ und $0$ sonst, also $t_i=\operatorname{Tr}L(c_i)=1+\frac{(r-1)d}2$.
$\square$

Im Rest des Abschnitts ist $J$ einfach, $(c_1,\dots,c_r)$ ein Jordan-Rahmen mit $r\ge2$, und für $x\in J_{ij}$,
$y\in J_{jk}$ ($i,j,k$ verschieden) schreiben wir
$$
m(x,y)=2\,x\circ y\in J_{ik}.
$$
Nach Satz E.8(3) und Proposition E.11(2) ist $\tau(x,y)=2t\langle x,y\rangle$ auf jedem $J_{ij}$.

**Lemma E.12.** Für verschiedene $i,j,k$ und $x\in J_{ij}$, $y\in J_{jk}$, $z\in J_{ik}$ gilt:

1. $N(m(x,y))=N(x)N(y)$ und $m(x,m(x,y))=N(x)y$;
2. $\langle m(x,y),z\rangle=\langle y,m(x,z)\rangle=\langle x,m(y,z)\rangle$;
3. $m(m(x,y),v)=m(x,m(y,v))$ für $l\notin\{i,j,k\}$ und $v\in J_{kl}$.

*Beweis.* (2) Die drei Zahlen sind $\tau(x\circ y,z)/t$, $\tau(y,x\circ z)/t$ und $\tau(x,y\circ z)/t$, gleich wegen
Assoziativität und Kommutativität. (1) Die zweite Formel ist Lemma E.10(1). Dann ist
$2tN(m(x,y))=4\tau(x\circ y,x\circ y)=4\tau(y,x\circ(x\circ y))=N(x)\tau(y,y)=2tN(x)N(y)$. (3) Lemma E.10(2). $\square$

**Satz E.13 (Rang eins und zwei).** Sei $J$ einfach mit einem Jordan-Rahmen aus $r$ Elementen. Ist $r=1$, so ist
$J=\mathbb Rc_1\cong\mathbb R$. Ist $r=2$, so ist $J\cong V_{d+1}$ mit $d=\dim J_{12}\ge1$ (Spinfaktoren wie in
Abschnitt 13.2).

*Beweis.* $r=1$: Satz E.8(1). $r=2$: $J=\mathbb R1\oplus W$ mit $W=\mathbb Rf\oplus J_{12}$, $f=c_1-c_2$. Nach Satz E.8
ist $f\circ f=1$, $f\circ x=\frac12x-\frac12x=0$ und $x\circ y=\langle x,y\rangle1$ für $x,y\in J_{12}$. Also ist
$w\circ w'=B(w,w')1$ für $w,w'\in W$, wobei $B(\alpha f+x,\beta f+y)=\alpha\beta+\langle x,y\rangle$ ein Skalarprodukt
ist. Ist $U:W\to\mathbb R^{d+1}$ eine Isometrie, so bildet $a1+w\mapsto(a,Uw)$ das Produkt
$(a1+w)\circ(b1+w')=(ab+B(w,w'))1+aw'+bw$ auf das Produkt von $V_{d+1}$ ab. $\square$

## E.6 Kompositionsalgebren: die Sätze von Hurwitz und Frobenius

Eine **Kompositionsalgebra** ist eine Algebra $A$ mit Einselement $1$ und einer positiv definiten quadratischen Form
$N$, sodass $N(xy)=N(x)N(y)$. Wir schreiben $\langle\cdot,\cdot\rangle$ für die Polarform ($N(x)=\langle x,x\rangle$),
$\operatorname{Re}x=\langle x,1\rangle$ und $\bar x=2\operatorname{Re}(x)1-x$. Da $N(1)=1$ (Lemma E.14(1), dessen Beweis
die Konjugation nicht benutzt), ist dies die Spiegelung an der Geraden $\mathbb R1$, also ist $x\mapsto\bar x$ eine
isometrische Involution und $\operatorname{Re}\bar x=\operatorname{Re}x$. Ein Isomorphismus von Kompositionsalgebren ist
ein isometrischer Algebrenisomorphismus. Der **Assoziator** ist $(x,y,z)=(xy)z-x(yz)$.

**Lemma E.14.** In einer Kompositionsalgebra gilt für alle $x,y,z,w$:

1. $N(1)=1$, $\langle xy,xz\rangle=N(x)\langle y,z\rangle$, $\langle yx,zx\rangle=\langle y,z\rangle N(x)$ und
   $\langle xy,wz\rangle+\langle wy,xz\rangle=2\langle x,w\rangle\langle y,z\rangle$;
2. $\langle xy,z\rangle=\langle y,\bar xz\rangle=\langle x,z\bar y\rangle$;
3. $\bar x(xy)=N(x)y=(yx)\bar x$; folglich $\bar xx=x\bar x=N(x)1$, $x^2-2\operatorname{Re}(x)x+N(x)1=0$,
   $xy+yx=2\operatorname{Re}(x)y+2\operatorname{Re}(y)x-2\langle x,y\rangle1$,
   $\bar x(wy)+\bar w(xy)=2\langle x,w\rangle y$ und $(yx)\bar w+(yw)\bar x=2\langle x,w\rangle y$;
4. $\overline{xy}=\bar y\bar x$, $\operatorname{Re}(xy)=\operatorname{Re}(yx)$ und
   $\operatorname{Re}((xy)z)=\operatorname{Re}(x(yz))$;
5. $x(xy)=x^2y$ und $(yx)x=yx^2$; der Assoziator ist alternierend (er wechselt das Vorzeichen, wenn man zwei Argumente
   vertauscht).

*Beweis.* (1) $N(1)=N(1)^2>0$. Man polarisiere $N(xy)=N(x)N(y)$ in $y$ bzw. in $x$, und dann
$\langle xy,xz\rangle=N(x)\langle y,z\rangle$ in $x$. (2) Die letzte Identität von (1) mit $w=1$ lautet
$\langle xy,z\rangle+\langle y,xz\rangle=2\operatorname{Re}(x)\langle y,z\rangle$, d. h.
$\langle xy,z\rangle=\langle y,\bar xz\rangle$; die Polarisierung von $\langle yx,zx\rangle=\langle y,z\rangle N(x)$ in
$x$, mit einer der beiden Variablen gleich $1$, liefert $\langle yx,z\rangle=\langle y,z\bar x\rangle$. (3) Nach (2)
und (1) ist $\langle\bar x(xy),z\rangle=\langle xy,xz\rangle=N(x)\langle y,z\rangle$ für alle $z$, und ebenso von
rechts. Man setze $y=1$ und substituiere $\bar x=2\operatorname{Re}(x)1-x$; dann polarisiere man. (4) Nach (2) und der
Isometrie der Konjugation ist
$\langle\bar y\bar x,z\rangle=\langle\bar x,yz\rangle=\langle\bar x\bar z,y\rangle=\langle\bar z,xy\rangle=\langle z,\overline{xy}\rangle$.
Weiter ist $\operatorname{Re}(xy)=\langle y,\bar x\rangle=\langle\bar y,x\rangle=\operatorname{Re}(yx)$ und
$\operatorname{Re}((xy)z)=\langle xy,\bar z\rangle=\langle y,\bar x\bar z\rangle=\langle yz,\bar x\rangle=\operatorname{Re}(x(yz))$.
(5) Nach (3) ist $x(xy)=(2\operatorname{Re}(x)-\bar x)(xy)=(2\operatorname{Re}(x)x-N(x))y=x^2y$, und ebenso von rechts.
Also $(x,x,y)=(y,x,x)=0$; durch Polarisieren wechselt der Assoziator das Vorzeichen, wenn man die ersten beiden oder die
letzten beiden Argumente vertauscht, und diese Transpositionen erzeugen alle Permutationen. $\square$

**Die Cayley–Dickson-Konstruktion.** Für eine Algebra $B$ mit Einselement, Involution $b\mapsto\bar b$ und quadratischer
Form $N$ sei $\operatorname{CD}(B)=B\times B$ mit
$$
(a,b)(c,d)=(ac-\bar db,\;da+b\bar c),\qquad N(a,b)=N(a)+N(b),\qquad\overline{(a,b)}=(\bar a,-b).
$$
Ausgehend von $\mathbb R$ ($\bar x=x$, $N(x)=x^2$) definieren wir $\mathbb C=\operatorname{CD}(\mathbb R)$,
$\mathbb H=\operatorname{CD}(\mathbb C)$ und $\mathbb O=\operatorname{CD}(\mathbb H)$, von den Dimensionen $1,2,4,8$,
wobei jede Verdopplung die soeben definierte Involution und Form benutzt (nach Lemma E.15(1) unten sind dies die
Konjugation und die Norm der Kompositionsalgebra). In $\mathbb C$ ist $(a,b)(c,d)=(ac-bd,ad+bc)$: Das sind die
komplexen Zahlen $a+bi$. $B$ ist die Unteralgebra $B\times0$.

**Lemma E.15.** (1) Ist $B$ eine assoziative Kompositionsalgebra, so ist $\operatorname{CD}(B)$ eine Kompositionsalgebra
mit Einselement $(1,0)$ und Konjugation $\overline{(a,b)}=(\bar a,-b)$. (2) $\mathbb C$ ist kommutativ und assoziativ.
$M(a,b)=\begin{pmatrix}a&-b\\\bar b&\bar a\end{pmatrix}$ ist ein injektiver Algebrenhomomorphismus
$\mathbb H\to M_2(\mathbb C)$ mit $M(\bar q)=M(q)^*$; also ist $\mathbb H$ assoziativ; sie ist nicht kommutativ.
$\mathbb O$ ist nicht assoziativ. (3) $\mathbb R,\mathbb C,\mathbb H,\mathbb O$ sind Kompositionsalgebren, und
$\mathbb R,\mathbb C,\mathbb H$ sind Divisionsalgebren: Jedes $x\ne0$ hat das beidseitige Inverse $\bar x/N(x)$, also
folgt aus $xy=0$, dass $x=0$ oder $y=0$.

*Beweis.* (1) $N((a,b)(c,d))=N(a)N(c)+N(d)N(b)+N(d)N(a)+N(b)N(c)-2\langle ac,\bar db\rangle+2\langle da,b\bar c\rangle$,
und nach Lemma E.14(2) ist $\langle ac,\bar db\rangle=\langle d(ac),b\rangle$ und
$\langle da,b\bar c\rangle=\langle(da)c,b\rangle$, was übereinstimmt, wenn $B$ assoziativ ist. $(1,0)$ ist ein
Einselement, und $2\langle(a,b),(1,0)\rangle(1,0)-(a,b)=(\bar a,-b)$. (2) Ausmultiplizieren unter Benutzung der
Kommutativität von $\mathbb C$ liefert $M(a,b)M(c,d)=M(ac-\bar db,da+b\bar c)$; $M$ ist linear und injektiv, und
$M(\bar a,-b)=M(a,b)^*$. In $\mathbb H$ ist $(i,0)(0,1)=(0,i)$, aber $(0,1)(i,0)=(0,-i)$. In
$\mathbb O=\mathbb H\times\mathbb H$ seien $p=((i,0),0)$, $q=((0,1),0)$, $s=(0,1)$; dann ist $(pq)s=(0,(0,i))$ und
$p(qs)=(0,(0,-i))$. (3) Nach (1), (2) und Induktion. Ist $x\ne0$, so ist $x(\bar x/N(x))=(\bar x/N(x))x=1$ nach Lemma
E.14(3); in einer assoziativen Algebra folgt dann aus $xy=0$, dass $y=(\bar x/N(x))(xy)=0$. $\square$

**Lemma E.16 (Verdopplung).** Sei $A$ eine Kompositionsalgebra, $B\subseteq A$ eine Unteralgebra mit $1\in B$ und
$i\in A$ mit $i\perp B$ und $N(i)=1$. Für $a,b,c,d\in B$ gilt $Bi\perp B$, $N(a+bi)=N(a)+N(b)$ und
$$
(a+bi)(c+di)=(ac-\bar db)+(da+b\bar c)\,i .
$$
Also ist $B\oplus Bi$ eine Unteralgebra der Dimension $2\dim B$, und $a+bi\mapsto(h(a),h(b))$ ist ein Isomorphismus
$B\oplus Bi\cong\operatorname{CD}(D)$, sobald $h:B\to D$ ein Isomorphismus von Kompositionsalgebren ist. Außerdem ist
$B$ assoziativ; ist $B\oplus Bi$ assoziativ, so ist $B$ kommutativ; ist $B\oplus Bi$ kommutativ, so ist $B=\mathbb R1$.

*Beweis.* $B$ enthält mit $b$ auch $\bar b=2\operatorname{Re}(b)1-b$. Nach Lemma E.14 ist
$\langle bi,a\rangle=\langle i,\bar ba\rangle=0$, $N(bi)=N(b)$, $\bar i=-i$ und $ci+ic=2\operatorname{Re}(c)i$, d. h.
$ci=i\bar c$, für $c\in B$. Die Identitäten von Lemma E.14(3) liefern der Reihe nach: mit $(x,w,y)=(\bar a,i,b)$ folgt
$a(ib)=i(\bar ab)$; da $ib=\bar bi$ und $i(\bar ab)=\overline{\bar ab}\,i=(\bar ba)i$, liefert das Ersetzen von $b$
durch $\bar b$, dass $a(bi)=(ba)i$. Mit $(y,x,w)=(a,i,\bar c)$ folgt $(ai)c=(a\bar c)i$. Mit $x=i$ folgt $i(iy)=-y$. Mit
$(x,w,y)=(ai,i,\bar d)$, wobei $\overline{ai}=-ai$ (da $ai\perp1$) und $\langle ai,i\rangle=\operatorname{Re}a$ (Lemma
E.14(1)), ergibt sich
$$
-(ai)(i\bar d)-i((ai)\bar d)=2\operatorname{Re}(a)\bar d .
$$
Hier ist $i\bar d=di$ und $i((ai)\bar d)=i((ad)i)=i(i\,\overline{ad})=-\bar d\bar a$, also
$(ai)(di)=\bar d\bar a-2\operatorname{Re}(a)\bar d=-\bar da$. Ausmultiplizieren von $(a+bi)(c+di)$ mit diesen vier
Regeln liefert die Formel, und $N$ und die Produkte entsprechen unter $h$ denen von $\operatorname{CD}(D)$. Da $N$ auf
$B\oplus Bi\subseteq A$ multiplikativ ist, liefert die Rechnung im Beweis von Lemma E.15(1), dass
$\langle d(ac)-(da)c,b\rangle=0$ für alle $b\in B$, also ist $B$ assoziativ. Ist $B\oplus Bi$ assoziativ, so ist
$(ab)i=a(bi)=(ba)i$, also $ab=ba$. Ist sie kommutativ, so ist $ai=ia=\bar ai$, also $a=\bar a\in\mathbb R1$. $\square$

**Satz E.17 (Hurwitz).** Jede Kompositionsalgebra ist zu genau einer der Algebren $\mathbb R,\mathbb C,\mathbb H,\mathbb O$
isomorph.

*Beweis.* Man setze $B_0=\mathbb R1\cong\mathbb R$. Solange $B_k\ne A$, wähle man $i_{k+1}\perp B_k$ mit $N(i_{k+1})=1$
und setze $B_{k+1}=B_k\oplus B_ki_{k+1}$; nach Lemma E.16 hat $B_k\cong\operatorname{CD}^k(\mathbb R)$ die Dimension
$2^k$. Gäbe es $B_4$, so wäre nach Lemma E.16 $B_3$ assoziativ, also $B_2$ kommutativ, also $B_1=\mathbb R1$, was die
Dimension $1\ne2$ hat. Also ist $A=B_k$ für ein $k\le3$. Die vier Algebren haben verschiedene Dimensionen. $\square$

**Lemma E.18.** Sei $A$ eine Kompositionsalgebra. (1) Je zwei Elemente von $A$ liegen in einer assoziativen
Unteralgebra. (2) (**Moufang-Identität**) $(ux)(yu)=u((xy)u)=(u(xy))u$ für alle $u,x,y$.

*Beweis.* (1) Seien $p,q\in A$; wir dürfen $p\notin\mathbb R1$ annehmen (sonst vertausche man sie, und sind beide
reell, nehme man $\mathbb R1$). $B=\operatorname{span}\{1,p\}$ ist eine Unteralgebra ($p^2=2\operatorname{Re}(p)p-N(p)1$),
eine Kompositionsalgebra der Dimension $2$, also $B\cong\mathbb C$ (Satz E.17). Ist $q\in B$, so sind wir fertig;
andernfalls sei $i$ die normierte, zu $B$ orthogonale Komponente von $q$. Nach Lemma E.16 ist
$B\oplus Bi\cong\operatorname{CD}(\mathbb C)=\mathbb H$ assoziativ und enthält $p,q$. (2) Der Assoziator ist
alternierend (Lemma E.14(5)), also $(u,xy,u)=0$, und das ist die zweite Gleichung. Sei $X=(u,x,y)$.
$\operatorname{Re}X=0$ nach Lemma E.14(4), und $\langle u,X\rangle=\langle u\bar y,ux\rangle-N(u)\langle1,xy\rangle
=N(u)(\langle\bar y,x\rangle-\langle\bar x,y\rangle)=0$; also liefert die polarisierte quadratische Identität von Lemma
E.14(3), dass $uX+Xu=2\operatorname{Re}(u)X$, d. h. $\bar uX=Xu$. Ausschreiben der Definitionen zeigt, dass in jeder
Algebra gilt
$$
(wx,y,z)-(w,xy,z)+(w,x,yz)=w(x,y,z)+(w,x,y)z .
$$
Mit $(w,x,y,z)=(u,u,x,y)$ und $(u^2,x,y)=2\operatorname{Re}(u)X$: $(u,ux,y)=2\operatorname{Re}(u)X-uX=\bar uX=Xu$. Mit
$(w,x,y,z)=(u,x,y,u)$, $(ux,y,u)=(u,ux,y)=Xu$ und $(x,y,u)=X$ (zyklische Permutationen): $(u,x,yu)=uX+Xu-Xu=uX$. Daher
ist $(ux)(yu)-u((xy)u)=(u,x,yu)-u(x,y,u)=uX-uX=0$. $\square$

**Satz E.19 (Frobenius).** Jede endlichdimensionale assoziative reelle Algebra $A$ mit Einselement $1\ne0$ und ohne
Nullteiler ist isomorph zu $\mathbb R$, $\mathbb C$ oder $\mathbb H$.

*Beweis.* (i) Für $x\in A$ sind die Potenzen $1,x,x^2,\dots$ linear abhängig; sei $\mu$ das normierte reelle Polynom
kleinsten Grades mit $\mu(x)=0$. Ist $\mu=fg$ mit $\deg f,\deg g\ge1$, so erzwingt $f(x)g(x)=0$, dass $f(x)=0$ oder
$g(x)=0$; also ist $\mu$ irreduzibel, folglich (Fundamentalsatz der Algebra) $\mu=t-\alpha$ oder
$\mu=(t-\alpha)^2+\beta^2$, $\beta>0$. Sei $V=\{v\in A:v^2\in(-\infty,0]\cdot1\}$. Dann ist $A=\mathbb R1+V$ (im
zweiten Fall ist $x-\alpha\in V$), und $\mathbb R1\cap V=0$. (ii) $V$ ist ein Unterraum: Es ist abgeschlossen unter
Skalaren, und für $u,v\in V$ zeigen wir $u+v\in V$. Sind $u,v$ linear abhängig, so ist $u+v$ ein Vielfaches eines von
ihnen; seien sie also unabhängig. Dann sind $1,u,v$ unabhängig: Ist $1=\alpha u+\gamma v$ mit etwa $\alpha\ne0$, so
erzwingt $u^2=(1-2\gamma v+\gamma^2v^2)/\alpha^2\in\mathbb R1$, dass $\gamma=0$ und $u\in\mathbb R1\cap V=0$. Nach (i)
ist $(u\pm v)^2=p_\pm(u\pm v)+q_\pm1$ mit reellen $p_\pm,q_\pm$; Addieren liefert
$2u^2+2v^2=(p_++p_-)u+(p_+-p_-)v+(q_++q_-)1$, also $p_\pm=0$. Ist $q_+>0$, so ist
$(u+v-\sqrt{q_+})(u+v+\sqrt{q_+})=0$ und $u+v\in\mathbb R1$, was unmöglich ist; also $u+v\in V$. (iii) Für $u,v\in V$
definiert $uv+vu=(u+v)^2-u^2-v^2=-2\langle u,v\rangle1$ eine symmetrische Bilinearform mit
$\langle u,u\rangle=-u^2\ge0$, und $\langle u,u\rangle=0$ nur für $u^2=0$, d. h. $u=0$: ein Skalarprodukt. Orthogonale
$u,v\in V$ antikommutieren. (iv) Ist $V=0$, so ist $A=\mathbb R$. Ist $\dim V=1$, so ist $V=\mathbb Ri$ mit $i^2=-1$
und $A\cong\mathbb C$. Ist $\dim V\ge2$, so nehme man orthonormale $i,j\in V$ und $k=ij$. Dann ist
$k^2=ijij=-i^2j^2=-1$, also $k\in V$, $ik+ki=-j+j=0$ und $jk+kj=i-i=0$: $i,j,k$ sind orthonormal, mit $ij=-ji=k$,
$jk=-kj=i$, $ki=-ik=j$. Gäbe es einen Einheitsvektor $l\in V$, der zu $i,j,k$ orthogonal ist, so wäre
$lk=(li)j=-i(lj)=(ij)l=kl$ und $lk=-kl$, also $kl=0$, unmöglich. Daher ist $A=\operatorname{span}\{1,i,j,k\}$. In
$\mathbb H$ haben die Elemente $I=(i,0)$, $J=(0,1)$, $K=(0,i)$ dieselbe Multiplikationstafel (Lemma E.15(2) liefert
$IJ=K=-JI$, und $I^2=J^2=-1$, $K=IJ$ bestimmen den Rest wegen der Assoziativität), also $A\cong\mathbb H$. $\square$

Der Code konstruiert $\mathbb O$ durch Cayley–Dickson-Verdopplung und prüft an zufälligen Oktonionen die Identitäten
von Lemma E.14, die Moufang-Identität von Lemma E.18 und die Nichtassoziativität aus Lemma E.15, sowie den
Homomorphismus $M$:

```python
o = torch.randn(4, 1000, 8, dtype=dt)             # random octonions
a, b, cc, u = o
N = lambda z: (z * z).sum(-1)
ip = lambda p, q: (p * q).sum(-1)
assert torch.allclose(N(cd(a, b)), N(a) * N(b))                               # composition
assert torch.allclose(ip(cd(a, b), cc), ip(b, cd(conj(a), cc)))               # Lemma E.14(2)
assert torch.allclose(ip(cd(a, b), cc), ip(a, cd(cc, conj(b))))
assert torch.allclose(conj(cd(a, b)), cd(conj(b), conj(a)))                   # Lemma E.14(4)
assert torch.allclose(cd(a, cd(a, b)), cd(cd(a, a), b))                       # alternative
assert torch.allclose(cd(cd(u, a), cd(b, u)), cd(u, cd(cd(a, b), u)))         # Moufang
assoc = cd(cd(a, b), cc) - cd(a, cd(b, cc))
assert assoc[..., 0].abs().max() < 1e-12 and assoc.abs().max() > 1          # Re-associative, not associative
e = torch.eye(8, dtype=dt)
print("(e1 e2) e4 =", cd(cd(e[1], e[2]), e[4]).tolist(), "\ne1 (e2 e4) =", cd(e[1], cd(e[2], e[4])).tolist())

def M(q):                                         # H -> M_2(C): (a, b) -> [[a, -b], [b*, a*]]
    za, zb = torch.complex(q[..., 0], q[..., 1]), torch.complex(q[..., 2], q[..., 3])
    return torch.stack([torch.stack([za, -zb], -1), torch.stack([zb.conj(), za.conj()], -1)], -2)
p4, q4 = torch.randn(2, 1000, 4, dtype=dt)
assert torch.allclose(M(cd(p4, q4)), M(p4) @ M(q4)) and torch.allclose(M(conj(p4)), M(p4).mH)
```

## E.7 Koordinatisierung

Für $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$, $d=\dim D$ und $n\ge1$ sei
$\operatorname{Herm}_n(D)=\{X\in M_n(D):X_{ji}=\overline{X_{ij}}\}$ mit $X\circ Y=\frac12(XY+YX)$, wobei
$(XY)_{ik}=\sum_jX_{ij}Y_{jk}$ (für diese Definition wird keine Assoziativität von $D$ benötigt). Diagonaleinträge
erfüllen $\bar x=x$, sind also reell. Seien $E_{ij}$ die Matrixeinheiten und, für $i\ne j$ und $a\in D$,
$a[ij]=aE_{ij}+\bar aE_{ji}$; somit $a[ji]=\bar a[ij]$.

**Lemma E.20.** $\operatorname{Herm}_n(D)$ ist unter $\circ$ abgeschlossen, hat das Einselement $\sum_iE_{ii}$ und die
Dimension $n+d\frac{n(n-1)}2$, und für verschiedene $i,j,k,l$ und $a,b\in D$ gilt
$$
E_{ii}\circ E_{jj}=\delta_{ij}E_{ii},\quad E_{kk}\circ a[ij]=\tfrac12(\delta_{ki}+\delta_{kj})a[ij],\quad
a[ij]\circ b[ij]=\langle a,b\rangle(E_{ii}+E_{jj}),\quad a[ij]\circ b[jk]=\tfrac12(ab)[ik],\quad a[ij]\circ b[kl]=0 .
$$
Außerdem ist $(X|Y)=\operatorname{Re}\operatorname{tr}(XY)$ ein Skalarprodukt mit $(X\circ Y|Z)=(Y|X\circ Z)$.

*Beweis.* Matrizenmultiplikation: $a[ij]\,b[ij]=a\bar bE_{ii}+\bar abE_{jj}$, und
$a\bar b+b\bar a=\bar ab+\bar ba=2\langle a,b\rangle$ (man polarisiere $x\bar x=\bar xx=N(x)$); $a[ij]\,b[jk]=abE_{ik}$
und $b[jk]\,a[ij]=\bar b\bar aE_{ki}=\overline{ab}E_{ki}$; die übrigen Produkte sind unmittelbar. Alle Produkte von
Basiselementen liegen in $\operatorname{Herm}_n(D)$, also auch alle Produkte. Weiter ist
$(X|X)=\sum_{i,j}\operatorname{Re}(X_{ij}\overline{X_{ij}})=\sum_{i,j}N(X_{ij})>0$ für $X\ne0$. Sei
$T(X,Y,Z)=\operatorname{Re}\operatorname{tr}((XY)Z)=\sum_{i,j,k}\operatorname{Re}((X_{ij}Y_{jk})Z_{ki})$. Nach Lemma
E.14(4) ist $\operatorname{Re}\operatorname{tr}(AB)=\operatorname{Re}\operatorname{tr}(BA)$, $T(X,Y,Z)=T(Z,X,Y)$, und
durch Übergang zu Konjugierten ($\operatorname{Re}\bar p=\operatorname{Re}p$, $\overline{(pq)r}=\bar r(\bar q\bar p)$)
unter Benutzung der Hermitizität von $X,Y,Z$ ist $T(X,Y,Z)=\sum\operatorname{Re}(Z_{ik}(Y_{kj}X_{ji}))=T(Z,Y,X)$.
Also ist $T$ symmetrisch, und $(X\circ Y|Z)=\frac12(T(X,Y,Z)+T(Y,X,Z))=T(X,Y,Z)$ ist symmetrisch in $X,Y,Z$. $\square$

**Lemma E.21 (Orientierungen).** Seien $P_{12},P_{13},P_{23}$ Skalarprodukträume, $P_{ji}=P_{ij}$, und für jedes
geordnete Tripel $(i,j,k)$ verschiedener Elemente von $\{1,2,3\}$ sei $m:P_{ij}\times P_{jk}\to P_{ik}$ bilinear,
sodass $m(x,y)=m(y,x)$ (rechts steht die Abbildung für $(k,j,i)$) und
$\langle m(x,y),z\rangle=\langle y,m(x,z)\rangle$ für $x\in P_{ij}$, $y\in P_{jk}$, $z\in P_{ik}$. Sei $D$ eine
Kompositionsalgebra und seien $\iota_{ij}:D\to P_{ij}$ ($i<j$) lineare isometrische Bijektionen mit
$m(\iota_{12}a,\iota_{23}b)=\iota_{13}(ab)$; man setze $\iota_{ji}(a)=\iota_{ij}(\bar a)$. Dann gilt
$m(\iota_{ij}a,\iota_{jk}b)=\iota_{ik}(ab)$ für alle geordneten Tripel und alle $a,b\in D$.

Die Voraussetzungen gelten (i) für drei Indizes eines Jordan-Rahmens eines einfachen $J$, mit den Peirce-Räumen und
$m=2\circ$ (Kommutativität und Lemma E.12(2)); (ii) in $\operatorname{Herm}_3(D)$, mit $P_{ij}=\{a[ij]\}$,
$\langle a[ij],b[ij]\rangle=\langle a,b\rangle$ (wohldefiniert, da die Konjugation isometrisch ist) und $m=2\circ$, denn
$m(a[ij],b[jk])=(ab)[ik]$ und, für $c[ik]$, $\langle ab,c\rangle=\langle b,\bar ac\rangle$ (Lemma E.14(2)).

*Beweis.* Gilt die Behauptung für $(i,j,k)$, so gilt sie für $(k,j,i)$:
$m(\iota_{kj}a,\iota_{ji}b)=m(\iota_{ij}\bar b,\iota_{jk}\bar a)=\iota_{ik}(\bar b\bar a)=\iota_{ik}(\overline{ab})=\iota_{ki}(ab)$.
Es genügt also, $(2,1,3)$ und $(1,3,2)$ zu behandeln, wenn $(1,2,3)$ gegeben ist. Für alle $c\in D$ gilt unter Benutzung
der Adjungiertheit, von $(1,2,3)$ und Lemma E.14(2)
$$
\langle m(\iota_{21}a,\iota_{13}b),\iota_{23}c\rangle=\langle\iota_{13}b,m(\iota_{12}\bar a,\iota_{23}c)\rangle=\langle b,\bar ac\rangle=\langle ab,c\rangle=\langle\iota_{23}(ab),\iota_{23}c\rangle,
$$
$$
\langle m(\iota_{13}a,\iota_{32}b),\iota_{12}c\rangle=\langle\iota_{13}a,m(\iota_{12}c,\iota_{23}\bar b)\rangle=\langle a,c\bar b\rangle=\langle ab,c\rangle=\langle\iota_{12}(ab),\iota_{12}c\rangle,
$$
und $\iota_{23},\iota_{12}$ sind surjektiv. (In der zweiten Zeile wird die Adjungiertheit für $y=\iota_{32}b$ in der
Rolle von $x$ benutzt, nach $m(x,y)=m(y,x)$.) $\square$

**Satz E.22 (Koordinatisierung).** Sei $J$ einfach mit einem Jordan-Rahmen $(c_1,\dots,c_r)$, $r\ge3$. Dann gibt es eine
Kompositionsalgebra $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$ der Dimension $d$ und einen Algebrenisomorphismus
$\Phi:\operatorname{Herm}_r(D)\to J$ mit $\Phi(E_{ii})=c_i$. Ist $r\ge4$, so ist $D$ assoziativ, d. h. $D\ne\mathbb O$.

*Beweis.* *Schritt 1 (ein Dreieck).* Man wähle $e_{12}\in J_{12}$, $e_{23}\in J_{23}$ mit $N=1$ (Proposition E.11). Auf
$D'=J_{13}$ definiere man $p\cdot q=m(m(p,e_{23}),m(e_{12},q))$, unter Benutzung der Produkte
$J_{13}\times J_{32}\to J_{12}$, $J_{21}\times J_{13}\to J_{23}$ und $J_{12}\times J_{23}\to J_{13}$. Nach Lemma
E.12(1) ist $1'=m(e_{12},e_{23})$ ein Einselement: $p\cdot1'=m(m(p,e_{23}),e_{23})=p$ und
$1'\cdot q=m(e_{12},m(e_{12},q))=q$; und $N(p\cdot q)=N(p)N(q)$. Also ist $(D',N)$ eine Kompositionsalgebra, und nach
Satz E.17 gibt es einen Isomorphismus $h:D\to D'$ mit $D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$. Man setze
$\iota_{13}=h$, $\iota_{12}=m(h(\cdot),e_{23})$ und $\iota_{23}=m(e_{12},h(\cdot))$; das sind isometrische Bijektionen
(Lemma E.12(1)), und $m(\iota_{12}a,\iota_{23}b)=h(a)\cdot h(b)=\iota_{13}(ab)$. Nach Lemma E.21 gilt für $k=3$:

> $(\mathrm T_k)$ Es gibt isometrische Bijektionen $\iota_{ij}:D\to J_{ij}$, $i<j\le k$, sodass mit $\iota_{ji}(a)=\iota_{ij}(\bar a)$
> $m(\iota_{ij}a,\iota_{jl}b)=\iota_{il}(ab)$ für alle verschiedenen $i,j,l\le k$ und alle $a,b\in D$ gilt.

*Schritt 2 (Hinzunahme eines Index).* Es gelte $(\mathrm T_{k-1})$ mit $4\le k\le r$. Man wähle $v\in J_{1k}$ mit
$N(v)=1$ und setze $\iota_{ik}(a)=m(\iota_{i1}a,v)$ für $2\le i<k$ und $\iota_{1k}(a)=m(\iota_{12}a,\iota_{2k}1)$,
isometrische Bijektionen nach Lemma E.12(1). Mit Lemma E.12(3) (jedes Mal mit vier verschiedenen Indizes) und
$(\mathrm T_{k-1})$ gilt für verschiedene $i,j\in\{2,\dots,k-1\}$
$$
m(\iota_{ij}a,\iota_{jk}b)=m(\iota_{ij}a,m(\iota_{j1}b,v))=m(m(\iota_{ij}a,\iota_{j1}b),v)=m(\iota_{i1}(ab),v)=\iota_{ik}(ab).
$$
Insbesondere ist $\iota_{2k}(b)=m(\iota_{23}b,\iota_{3k}1)$, und daher
$$
m(\iota_{12}a,\iota_{2k}b)=m(m(\iota_{12}a,\iota_{23}b),\iota_{3k}1)=m(\iota_{13}(ab),\iota_{3k}1)
=m(m(\iota_{12}(ab),\iota_{23}1),\iota_{3k}1)=m(\iota_{12}(ab),\iota_{2k}1)=\iota_{1k}(ab).
$$
Für $3\le j<k$ ist $\iota_{1j}a=m(\iota_{12}a,\iota_{2j}1)$, also
$m(\iota_{1j}a,\iota_{jk}b)=m(\iota_{12}a,m(\iota_{2j}1,\iota_{jk}b))=m(\iota_{12}a,\iota_{2k}b)=\iota_{1k}(ab)$. Somit
gilt die Relation für alle Tripel $i<j<k$, und Lemma E.21, angewandt auf jedes Dreieck $\{i,j,k\}$, liefert
$(\mathrm T_k)$.

*Schritt 3 (Assoziativität).* Ist $r\ge4$, so gilt nach $(\mathrm T_4)$ und Lemma E.12(3)
$\iota_{14}((ab)c)=m(m(\iota_{12}a,\iota_{23}b),\iota_{34}c)=m(\iota_{12}a,m(\iota_{23}b,\iota_{34}c))=\iota_{14}(a(bc))$.

*Schritt 4 (der Isomorphismus).* Mithilfe von $(\mathrm T_r)$ definiere man die lineare Abbildung
$\Phi(\sum_i\alpha_iE_{ii}+\sum_{i<j}a_{ij}[ij])=\sum_i\alpha_ic_i+\sum_{i<j}\iota_{ij}(a_{ij})$; sie ist bijektiv nach
Satz E.8(1), und $\Phi(a[ij])=\iota_{ij}(a)$ für alle $i\ne j$. Nach Satz E.8 und $(\mathrm T_r)$ gilt
$c_i\circ c_j=\delta_{ij}c_i$, $c_k\circ\iota_{ij}(a)=\frac12(\delta_{ki}+\delta_{kj})\iota_{ij}(a)$,
$\iota_{ij}(a)\circ\iota_{ij}(b)=\langle a,b\rangle(c_i+c_j)$, $\iota_{ij}(a)\circ\iota_{jk}(b)=\frac12\iota_{ik}(ab)$ und
$\iota_{ij}(a)\circ\iota_{kl}(b)=0$: die Tafel von Lemma E.20. Wegen der Bilinearität ist $\Phi(X\circ Y)=\Phi(X)\circ\Phi(Y)$.
$\square$

## E.8 Die Albert-Algebra

**Lemma E.23 (Automorphismen von $\operatorname{Herm}_3(D)$).** Sei $D$ eine Kompositionsalgebra und seien
$f_{12},f_{23},f_{13}$ isometrische lineare Bijektionen von $D$ mit $f_{13}(ab)=f_{12}(a)f_{23}(b)$ für alle $a,b$.
Dann ist
$$
\Phi\Big(\sum_i\alpha_iE_{ii}+\sum_{i<j}a_{ij}[ij]\Big)=\sum_i\alpha_iE_{ii}+\sum_{i<j}f_{ij}(a_{ij})[ij]
$$
ein Automorphismus von $(\operatorname{Herm}_3(D),\circ)$. Für jedes $u\in D$ mit $N(u)=1$ ist
$(f_{12},f_{23},f_{13})=(L_u,R_u,B_u)$ mit $L_ua=ua$, $R_ua=au$, $B_ua=uau$ ein solches Tripel.

*Beweis.* Man wende Lemma E.21 im Beispiel (ii) auf $\iota_{ij}(a)=f_{ij}(a)[ij]$ an: Dann ist
$\Phi(a[ij])=\iota_{ij}(a)$ für alle $i\ne j$, und $2\iota_{ij}(a)\circ\iota_{jk}(b)=\iota_{ik}(ab)$ für alle
verschiedenen $i,j,k$. Da $\Phi$ die $E_{ii}$ festlässt und die $f_{ij}$ isometrisch sind, erhält $\Phi$ jedes Produkt
der Tafel von Lemma E.20, also $\circ$. Was das Tripel betrifft, so sind $L_u,R_u,B_u$ isometrisch ($N(ua)=N(a)$) und
bijektiv, und $u(ab)u=(ua)(bu)$ ist die Moufang-Identität (Lemma E.18(2)). $\square$

**Satz E.24.** (1) Für $D\in\{\mathbb R,\mathbb C,\mathbb H\}$ und $n\ge1$ ist $\operatorname{Herm}_n(D)$ eine
euklidische Jordan-Algebra, und $\operatorname{Herm}_n(\mathbb H)$ ist isomorph zu einer Jordan-Unteralgebra von
$\operatorname{Herm}_{2n}(\mathbb C)$. (2) Die **Albert-Algebra** $\operatorname{Herm}_3(\mathbb O)$ (Dimension $27$)
und $\operatorname{Herm}_2(\mathbb O)$ sind euklidische Jordan-Algebren. (3) In jeder dieser Algebren ist
$(E_{11},\dots,E_{nn})$ ein Jordan-Rahmen mit Peirce-Räumen $J_{ij}=\{a[ij]:a\in D\}$ der Dimension $d=\dim D$; die
Algebra ist einfach; $\operatorname{Herm}_1(D)\cong\mathbb R$ und $\operatorname{Herm}_2(D)\cong V_{d+1}$ für
$D=\mathbb R,\mathbb C,\mathbb H,\mathbb O$.

*Beweis.* (1) $M_n(D)$ ist eine assoziative reelle Algebra (Lemma E.15(2)), also ist nach Lemma E.20 und Lemma 13.1
$\operatorname{Herm}_n(D)$ eine Jordan-Algebra, euklidisch nach Lemma E.20. Ersetzt man jeden Eintrag
$X_{ij}\in\mathbb H$ durch den $2\times2$-Block $M(X_{ij})$, so erhält man einen injektiven Algebrenhomomorphismus
$M_n(\mathbb H)\to M_{2n}(\mathbb C)$ (Lemma E.15(2)), der hermitesche Matrizen auf hermitesche Matrizen abbildet, da
$M(\bar q)=M(q)^*$; er erhält also $\circ$.

(2) Sei $J(x,y)=(x^2\circ y)\circ x-x^2\circ(y\circ x)$ auf $\operatorname{Herm}_3(\mathbb O)$; ein Automorphismus
$\Phi$ erfüllt $J(\Phi x,\Phi y)=\Phi J(x,y)$. Sei $x=\sum_i\alpha_iE_{ii}+p[12]+q[23]+s[13]$. Wir setzen drei
Automorphismen der Form aus Lemma E.23 so zusammen, dass die Einträge an den Stellen $(1,3)$ und $(2,3)$ reell werden.
(a) Ist $s\notin\mathbb R$, so schreibe man $s=|s|(\cos\theta+v\sin\theta)$ mit $v\perp1$, $N(v)=1$, sodass $v^2=-1$
(Lemma E.14(3)); für $u=\cos\frac\theta2-v\sin\frac\theta2$ ergibt die Rechnung in der kommutativen Unteralgebra
$\operatorname{span}\{1,v\}\cong\mathbb C$, dass $B_us=u^2s=|s|$. (b) Nun ist $s$ reell. Man wähle $u\perp1$ mit
$N(u)=1$ und $u\perp q$ (möglich, da $\dim\mathbb O=8$); dann ist $B_us=su^2=-s$ reell und
$\operatorname{Re}(qu)=\langle u,\bar q\rangle=0$, da $\bar q\in\operatorname{span}\{1,q\}$. (c) Nun ist $q\perp1$; ist
$q\ne0$, so nehme man $u=q/|q|$: $qu=q^2/|q|=-|q|$ und $B_us=-s$ sind reell. Es gibt also einen Automorphismus $\Phi$,
sodass $\Phi x$ reelle Einträge hat, außer möglicherweise $p'$ an der Stelle $(1,2)$. Da $J(\Phi x,\cdot)$ linear ist,
genügt es, $J(\Phi x,y)=0$ für $y=E_{ii}$ und $y=c[ij]$, $c\in\mathbb O$, zu zeigen. Nach Lemma E.18(1) liegen $p'$ und
$c$ in einer assoziativen Unteralgebra $Q$; dann ist $\Phi x,y\in\operatorname{Herm}_3(Q)$, eine Unteralgebra von
$\operatorname{Herm}_3(\mathbb O)$ innerhalb der assoziativen Algebra $M_3(Q)$, und $J(\Phi x,y)=0$ nach Lemma 13.1.
Daher ist $J(x,\Phi^{-1}y)=\Phi^{-1}J(\Phi x,y)=0$ für alle $y$: Die Jordan-Identität gilt, und Lemma E.20 macht die
Algebra euklidisch. Die Matrizen mit verschwindender dritter Zeile und Spalte bilden eine Unteralgebra (Lemma E.20)
mit Einselement $E_{11}+E_{22}$, und diese ist $\operatorname{Herm}_2(\mathbb O)$.

(3) Nach Lemma E.20 sind die $E_{ii}$ paarweise orthogonale Idempotente mit Summe $1$; $L(E_{kk})$ wirkt auf $E_{jj}$
durch $\delta_{kj}$ und auf $a[ij]$ durch $\frac12(\delta_{ki}+\delta_{kj})$, also $J_1(E_{ii})=\mathbb RE_{ii}$ und $J_{1/2}(E_{ii})\cap
J_{1/2}(E_{jj})=\{a[ij]\}$. $E_{ii}$ ist primitiv: Ist $E_{ii}=g+g'$ mit orthogonalen, von null verschiedenen
Idempotenten, so ist $g=E_{ii}\circ g\in\mathbb RE_{ii}$, also $g=E_{ii}$ und $g'=0$. Der Graph aus Proposition E.11
ist vollständig, also ist die Algebra einfach. Für $n=1,2$ wende man Satz E.13 an. $\square$

**Korollar E.25.** Für $n\ge4$ ist $\operatorname{Herm}_n(\mathbb O)$ keine Jordan-Algebra.

*Beweis.* Wäre sie eine, so wäre Lemma E.10(2) anwendbar: Sein Beweis benutzt nur die Jordan-Identität und die
Produkte aus Lemma E.20 (mit $c_k=E_{kk}$). Dann wäre
$\frac14((ab)c)[14]=(a[12]\circ b[23])\circ c[34]=a[12]\circ(b[23]\circ c[34])=\frac14(a(bc))[14]$ für alle
$a,b,c\in\mathbb O$, im Widerspruch zu Lemma E.15(2). $\square$

**Satz E.26 (Albert).** Sei $D$ eine Kompositionsalgebra. Gibt es eine assoziative reelle Algebra $A$ und eine
injektive lineare Abbildung $\varphi:\operatorname{Herm}_3(D)\to A$ mit
$\varphi(x\circ y)=\frac12(\varphi(x)\varphi(y)+\varphi(y)\varphi(x))$, so ist $D$ assoziativ. Insbesondere ist die
Albert-Algebra **exzeptionell**: Sie ist zu keiner Jordan-Unteralgebra irgendeiner assoziativen Algebra, etwa
$M_N(\mathbb C)$, isomorph.

*Beweis.* Man setze $e_i=\varphi(E_{ii})$ und $e=e_1+e_2+e_3=\varphi(1)$. Dann ist $e_i^2=e_i$, und $e_ie_j+e_je_i=0$
für $i\ne j$; Multiplikation mit $e_i$ von links bzw. von rechts liefert $e_ie_j+e_ie_je_i=0=e_ie_je_i+e_je_i$, also
$e_ie_j=e_je_i=0$. Für $x=\varphi(X)$ liefert $ex+xe=2x$, dass $exe=ex=xe=x$, also $x=\sum_{p,q}x_{pq}$ mit
$x_{pq}=e_pxe_q$. Für $X=a[ij]$ besagen die Relationen $e_i\circ x=e_j\circ x=\frac12x$ und $e_k\circ x=0$
($k\ne i,j$), dass $\frac12([p=i]+[q=i])x_{pq}=\frac12x_{pq}$ usw., was nur $x_{ij}$ und $x_{ji}$ übrig lässt. Man
definiere für geordnete Paare $\alpha_{ij}(a)=e_i\varphi(a[ij])e_j$. Für verschiedene $i,j,k$ ist
$\varphi(a[ij])\varphi(b[jk])=\alpha_{ij}(a)\alpha_{jk}(b)$ der einzige Term in $e_iAe_k$ von
$\varphi(a[ij])\varphi(b[jk])+\varphi(b[jk])\varphi(a[ij])=\varphi((ab)[ik])$, also
$$
\alpha_{ij}(a)\,\alpha_{jk}(b)=\alpha_{ik}(ab).
$$
Aus $\varphi(1[13])^2=\varphi(E_{11}+E_{33})$ erhalten wir $\alpha_{13}(1)\alpha_{31}(1)=e_1$,
$\alpha_{31}(1)\alpha_{13}(1)=e_3$; aus $\varphi(a[13])^2=N(a)(e_1+e_3)$ folgt $\alpha_{13}(a)\alpha_{31}(\bar a)=N(a)e_1$,
also ist $\alpha_{13}$ injektiv. Nun ist
$$
\alpha_{13}(ay)=\alpha_{12}(a)\alpha_{23}(y)=\alpha_{13}(a)\alpha_{32}(1)\alpha_{21}(1)\alpha_{13}(y)=\alpha_{13}(a)\alpha_{31}(1)\alpha_{13}(y),
$$
und daher, mit $\pi=\alpha_{31}(1)$ und der Assoziativität in $A$,
$\alpha_{13}((ab)c)=\alpha_{13}(a)\pi\alpha_{13}(b)\pi\alpha_{13}(c)=\alpha_{13}(a(bc))$. Also ist $(ab)c=a(bc)$. Nach
Lemma E.15(2) ist $\mathbb O$ nicht assoziativ. $\square$

Alberts ursprüngliches Argument (A. A. Albert, Ann. of Math. 35 (1934) 65–73) ist von dieser Art. Der Code prüft an
zufälligen Elementen die Jordan-Identität von $\operatorname{Herm}_3(\mathbb O)$ sowie die Assoziativität und Positivität
ihrer Spurform, zeigt numerisch, dass die Jordan-Identität in $\operatorname{Herm}_4(\mathbb O)$ verletzt ist
(Korollar E.25), prüft den Automorphismus aus Lemma E.23 und führt die drei Schritte (a)–(c) des Beweises von Satz
E.24(2) an einem zufälligen Element aus:

```python
n3, prod3 = herm(3, 8)
X, Y, Z = torch.randn(3, 200, 27, dtype=dt)
X2 = prod3(X, X)
assert (prod3(prod3(X2, Y), X) - prod3(X2, prod3(Y, X))).abs().max() < 1e-10    # Jordan identity, Theorem E.24(2)
n4, prod4 = herm(4, 8)
X4, Y4 = torch.randn(2, 200, n4, dtype=dt)
X42 = prod4(X4, X4)
defect4 = (prod4(prod4(X42, Y4), X4) - prod4(X42, prod4(Y4, X4))).abs().max().item()
print("Herm_4(O): max |Jordan defect| =", defect4)
assert defect4 > 1                                                               # Corollary E.25
retr = lambda v: v[..., :3].sum(-1)               # Re tr
assert torch.allclose(retr(prod3(prod3(X, Y), Z)), retr(prod3(X, prod3(Y, Z))))
sq = (X[:, :3] ** 2).sum(-1) + 2 * (X[:, 3:] ** 2).sum(-1)
assert torch.allclose(retr(prod3(X, X)), sq)                                    # (X|X) = sum of N(X_ij)

def moufang_auto(v, uu):                          # (L_u, R_u, B_u) on the entries a_12, a_23, a_13
    a12, a13, a23 = v[..., 3:11], v[..., 11:19], v[..., 19:27]   # pairs order (0,1),(0,2),(1,2)
    uu = uu.expand_as(a12)
    return torch.cat([v[..., :3], cd(uu, a12), cd(cd(uu, a13), uu), cd(a23, uu)], -1)
uu = torch.randn(8, dtype=dt); uu = uu / uu.norm()
assert torch.allclose(moufang_auto(prod3(X, Y), uu), prod3(moufang_auto(X, uu), moufang_auto(Y, uu)))

x = torch.randn(27, dtype=dt)                     # steps (a)-(c) of the proof of Theorem E.24(2)
s = x[11:19]
im = torch.cat([torch.zeros(1, dtype=dt), s[1:]])
v, th = im / im.norm(), torch.atan2(im.norm(), s[0]).item()
x1 = moufang_auto(x, math.cos(th / 2) * e[0] - math.sin(th / 2) * v)            # (a): (1,3) entry -> |s|
imq = torch.cat([torch.zeros(1, dtype=dt), x1[20:27]])
g = torch.cat([torch.zeros(1, dtype=dt), torch.randn(7, dtype=dt)])
g = g - (g @ imq) / (imq @ imq) * imq                                             # u perpendicular to 1 and q
x2 = moufang_auto(x1, g / g.norm())                                               # (b): (2,3) entry imaginary
x3 = moufang_auto(x2, x2[19:27] / x2[19:27].norm())                              # (c): (2,3) entry real
print("entries (1,3), (2,3) after (a)-(c):", x3[11:19].numpy().round(3), x3[19:27].numpy().round(3))
assert x3[12:19].abs().max() < 1e-12 and x3[20:27].abs().max() < 1e-12
assert abs(abs(x3[11].item()) - s.norm().item()) < 1e-12
```

## E.9 Die Klassifikation

**Satz E.27 (Jordan–von Neumann–Wigner).**

1. Jede endlichdimensionale formal reelle (äquivalent: euklidische) Jordan-Algebra ist eine direkte Summe einfacher
   Ideale, und jede einfache ist isomorph zu $\operatorname{Herm}_n(\mathbb R)$, $\operatorname{Herm}_n(\mathbb C)$ oder
   $\operatorname{Herm}_n(\mathbb H)$ für ein $n\ge1$, zu einem Spinfaktor $V_m$ mit $m\ge2$ oder zur Albert-Algebra
   $\operatorname{Herm}_3(\mathbb O)$.
2. Die Algebren $\mathbb R$; $V_m$ ($m\ge2$); $\operatorname{Herm}_n(D)$ ($D=\mathbb R,\mathbb C,\mathbb H$, $n\ge3$);
   $\operatorname{Herm}_3(\mathbb O)$ sind einfach, formal reell und paarweise nicht isomorph, sodass jede einfache
   formal reelle Jordan-Algebra zu genau einer von ihnen isomorph ist. Außerdem ist $\operatorname{Herm}_1(D)\cong\mathbb R$ und
   $\operatorname{Herm}_2(D)\cong V_{\dim D+1}$ für $D=\mathbb R,\mathbb C,\mathbb H,\mathbb O$.
3. In einer einfachen formal reellen Jordan-Algebra haben alle Jordan-Rahmen dieselbe Anzahl $r$ von Elementen und alle
   ihre Peirce-Räume $J_{ij}$ dieselbe Dimension $d$: $r=1$ für $\mathbb R$, $(r,d)=(2,m-1)$ für $V_m$ und
   $(n,\dim D)$ für $\operatorname{Herm}_n(D)$, $n\ge2$.

*Beweis.* (1) Nach Satz E.9 ist $J$ eine direkte Summe einfacher formal reeller Ideale. Ein solches Ideal besitzt einen
Jordan-Rahmen einer gewissen Größe $r$ (Lemma E.6(2)); es ist $\cong\mathbb R=\operatorname{Herm}_1(\mathbb R)$ für
$r=1$, $\cong V_{d+1}$ für $r=2$ (Satz E.13) und $\cong\operatorname{Herm}_r(D)$ mit $D\ne\mathbb O$ für $r\ge4$ und
$D\in\{\mathbb R,\mathbb C,\mathbb H,\mathbb O\}$ für $r=3$ (Satz E.22).

(2) Für $\operatorname{Herm}_n(D)$ siehe Satz E.24. $V_m$ ist eine Jordan-Algebra nach Lemma 13.3 und euklidisch für
$\langle(a,u),(b,v)\rangle=ab+u\cdot v$, da $\langle x\circ y,z\rangle=abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w$
symmetrisch in $x=(a,u)$, $y=(b,v)$, $z=(c,w)$ ist. Für einen Einheitsvektor $w$ ist $c_{\pm w}=\frac12(1,\pm w)$ ein
Jordan-Rahmen (Satz 13.4(1); $c_{\pm w}$ ist primitiv, weil eine Summe $c_{w'}+c_{-w'}$ zweier orthogonaler, von null
verschiedener Idempotente gleich $1$ ist) mit $J_{12}=\{(0,v):v\perp w\}$, der Dimension $m-1\ge1$ (die Dimensionen
summieren sich zu $m+1$); also ist $V_m$ nach Proposition E.11 einfach. Um die Algebren zu unterscheiden, sei $s(J)$ die
größte Anzahl paarweise orthogonaler, von null verschiedener Idempotente (endlich nach Lemma E.6(2)); $s$ und $\dim$
sind Isomorphieinvarianten. Offenbar ist $s(\mathbb R)=1$. In $V_m$ sind die von null verschiedenen Idempotente $1$ und
die $c_w$, mit $1\circ c_w=c_w$ und $c_w\circ c_{w'}=0$ nur für $w'=-w$ (Satz 13.4(1)); also $s(V_m)=2$. Für
$D\in\{\mathbb R,\mathbb C,\mathbb H\}$ ist $s(\operatorname{Herm}_n(D))=n$: Die $E_{ii}$ liefern $s\ge n$. Umgekehrt
erfüllen orthogonale Idempotente $p_k\ne0$ die Gleichung $p_kp_l=0$ für $k\ne l$ (wie im Beweis von Satz E.26), also
sind die $D$-Rechtsuntermoduln $V_k=p_kD^n$ von $D^n$ von null verschieden und unabhängig: Ist $\sum_kv_k=0$ mit
$v_k\in V_k$, so ist $v_l=p_l\sum_kv_k=0$. Ein von null verschiedenes $V_k$ enthält ein $v\ne0$ und damit $vD$, von
reeller Dimension $d$ ($D$ hat keine Nullteiler); also $sd\le nd$. Für $\operatorname{Herm}_3(\mathbb O)$ liefern die
$E_{ii}$, dass $s\ge3$; und $s$ paarweise orthogonale, von null verschiedene Idempotente lassen sich zusammen mit $1$
minus ihrer Summe wie im Beweis von Satz E.7 zu einem Jordan-Rahmen aus $r'\ge s$ Elementen verfeinern. Nach
Proposition E.11(2) und Satz E.22 ist $27=r'+d'\frac{r'(r'-1)}2$ mit $d'\ge1$, und $d'\in\{1,2,4\}$ für $r'\ge4$. Für
$r'=4,5,6$ liegt die rechte Seite in $\{10,16,28\}$, $\{15,25,45\}$, $\{21,36,66\}$, und für $r'\ge7$ übersteigt sie
$27$; also $s\le r'=3$. Die Paare $(s,\dim)$ sind $(1,1)$, $(2,m+1)$, $(n,\frac{n(n+1)}2)$, $(n,n^2)$, $(n,n(2n-1))$ für
$n\ge3$ und $(3,27)$, alle verschieden, da $\frac{n(n+1)}2<n^2<n(2n-1)$ und $6,9,15\ne27$. Die Isomorphismen in Rang
$1$ und $2$ sind Satz E.24(3).

(3) Ist $F$ ein Jordan-Rahmen aus $r_F$ Elementen eines einfachen $J$, so liefert der Beweis von (1), angewandt mit
$F$, einen Isomorphismus von $J$ mit einem Modell vom Rang $r_F$ aus (2), und das $s$ dieses Modells ist $r_F$. Also
ist $r_F=s(J)$ für jeden Rahmen, und dann $d=(\dim J-r)/\binom r2$ für $r\ge2$ (Proposition E.11(2)). $\square$

**Korollar E.28 (Rang und Spur).** In jeder formal reellen Jordan-Algebra haben alle Jordan-Rahmen dieselbe Anzahl von
Elementen, den **Rang**. Ist $J$ einfach vom Rang $r$, so erfüllt $\operatorname{tr}x=\frac r{\dim J}\operatorname{Tr}L(x)$
die Gleichung $\operatorname{tr}c=1$ für jedes primitive Idempotent $c$, folglich
$\operatorname{tr}(\sum_i\lambda_ic_i)=\sum_i\lambda_i$ für jede Zerlegung aus Satz E.7, und $\operatorname{tr}(x\circ y)$
ist ein assoziatives Skalarprodukt.

*Beweis.* Sei $J=\bigoplus_pJ^{(p)}$ wie in Satz E.9 und $c$ primitiv. Seine Komponenten $c^{(p)}\in J^{(p)}$ sind
paarweise orthogonale Idempotente ($c^2=\sum_p(c^{(p)})^2$), also ist genau eine von null verschieden. Daher bilden die
in $J^{(p)}$ liegenden Mitglieder eines Jordan-Rahmens von $J$ einen Jordan-Rahmen von $J^{(p)}$ (sie summieren sich
zur Komponente $1^{(p)}$ von $1$), und der Rang von $J$ ist die Summe der Ränge der $J^{(p)}$ (Satz E.27(3)). In einem
einfachen $J$ gehört ein primitives $c$ zu einem Jordan-Rahmen (man verfeinere $1-c$ wie im Beweis von Satz E.7), also
ist $\operatorname{Tr}L(c)=t=\dim J/r$ nach Proposition E.11(2) (für $r=1$ ist $J=\mathbb Rc$ und
$\operatorname{Tr}L(c)=1$). Schließlich ist $\operatorname{tr}(x\circ y)=\frac r{\dim J}\tau(x,y)$ (Satz E.5). $\square$

Der Code berechnet für jeden Typ den Jordan-Rahmen $(E_{ii})$ bzw. $(c_w,c_{-w})$, die Dimensionen der Peirce-Räume als
gemeinsame Eigenräume der $L(c_i)$ (Satz E.8), die Dimensionsformel aus Proposition E.11(2) und die Positivität der
Spurform $\tau$ (Satz E.5):

```python
def peirce(C, frame):
    """dimensions of the joint Peirce spaces J_ij = {x : c_i o x = c_j o x = x/2}"""
    I = torch.eye(C.shape[0], dtype=dt)
    Ls = [Lop(C, f) for f in frame]
    return [C.shape[0] - torch.linalg.matrix_rank(torch.cat([Ls[i] - I / 2, Ls[j] - I / 2])).item()
            for i in range(len(frame)) for j in range(i + 1, len(frame))]

cases = [(f"Herm_{n}({'RCHO'[k]})", herm(n, 2 ** k), n, 2 ** k) for n, k in ((2, 0), (3, 0), (4, 0), (3, 1), (4, 1), (3, 2), (2, 3), (3, 3))]
cases += [(f"V_{m}", spin(m), 2, m - 1) for m in (2, 3, 5)]
for name, (dim, prod), r, d in cases:
    C = structure(dim, prod)
    if name.startswith('V'):
        frame = [0.5 * torch.eye(dim, dtype=dt)[0] + s * 0.5 * torch.eye(dim, dtype=dt)[1] for s in (1, -1)]
    else:
        frame = list(torch.eye(dim, dtype=dt)[:r])
    for i, f in enumerate(frame):
        assert torch.allclose(mul(C, f, f), f) and all(torch.allclose(mul(C, f, g), torch.zeros(dim, dtype=dt)) for g in frame[i + 1:])
    dims = peirce(C, frame)
    t = torch.einsum('cbb->c', Lop(C, torch.eye(dim, dtype=dt)))
    gram = torch.einsum('abc,c->ab', C, t)
    print(f"{name:10s} dim {dim:2d}  rank {r}  Peirce dimensions {dims}  min eig of trace form {torch.linalg.eigvalsh(gram).min().item():.2f}")
    assert dims == [d] * (r * (r - 1) // 2) and dim == r + d * r * (r - 1) // 2
    assert torch.linalg.eigvalsh(gram).min() > 0
```

## E.10 Was dieser Anhang für Kapitel 13 beweist

| Aussage in Kapitel 13 | hier bewiesen in |
|---|---|
| formal reell $\Leftrightarrow$ euklidisch (Abschnitt 13.4) | Satz E.5 |
| Satz 13.7: Zerlegung in einfache Ideale, die Liste, „genau eine“ | Sätze E.9, E.27 |
| $\operatorname{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}$ für $\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O$ | Satz E.24(3) |
| $\operatorname{Herm}_r(\mathbb H)$ bettet sich in $\operatorname{Herm}_{2r}(\mathbb C)$ ein | Satz E.24(1) |
| die Albert-Algebra ist eine Jordan-Algebra und exzeptionell | Sätze E.24(2), E.26 |
| $\operatorname{Herm}_r(\mathbb O)$ ist nur für $r\le3$ eine Jordan-Algebra | Satz E.24(2), Korollar E.25 |
| Spektralsatz mit Jordan-Rahmen, Rang, $\operatorname{tr}z=\sum_i\mu_i$, Assoziativität von $\operatorname{tr}(x\circ y)$ (Abschnitt 13.6) | Satz E.7, Korollar E.28 |

**Anmerkungen.** Die Klassifikation stammt von P. Jordan, J. von Neumann und E. Wigner, *On an algebraic generalization
of the quantum mechanical formalism*, Ann. of Math. 35 (1934) 29–64, und die Exzeptionalität von
$\operatorname{Herm}_3(\mathbb O)$ von A. A. Albert (gleicher Band, 65–73). Der obige Beweis folgt dem Aufbau von
Faraut–Korányi, *Analysis on Symmetric Cones* (Oxford 1994), Kapitel II–V; der Verdopplungsbeweis des Satzes von
Hurwitz (Nachr. Ges. Wiss. Göttingen 1898, 309–316) folgt J. H. Conway, D. A. Smith, *On Quaternions and Octonions*
(2003) und T. A. Springer, F. D. Veldkamp, *Octonions, Jordan Algebras and Exceptional Groups* (2000), und der Beweis
des Satzes von Frobenius (1878) ist der von R. S. Palais, Amer. Math. Monthly 75 (1968) 366–368. Weiterführende
Literatur (im Kurs nirgends verwendet): C. M. Glennie, Pacific J. Math. 16 (1966) 47–59, für einen anderen Zugang zur
Exzeptionalität über polynomiale Identitäten, und McCrimmons *A Taste of Jordan Algebras* (2004) für die Geschichte und
die allgemeine Theorie.
