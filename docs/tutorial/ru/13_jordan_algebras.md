[English](../13_jordan_algebras.md) · [Español](../es/13_jordan_algebras.md) · [Français](../fr/13_jordan_algebras.md) · [Deutsch](../de/13_jordan_algebras.md) · [中文](../zh/13_jordan_algebras.md) · [日本語](../ja/13_jordan_algebras.md) · **Русский**

# 13. Йордановы алгебры: почему квантовая механика использует *-алгебры

Произведение $xy$ двух наблюдаемых не является наблюдаемой, если $x$ и $y$ не коммутируют, но симметризованное произведение
$x\circ y=\frac12(xy+yx)$ является ею всегда. В 1933 году П. Йордан предложил взять наблюдаемые с этим произведением в качестве
основной структуры квантовой механики, а Йордан, фон Нейман и Вигнер классифицировали возникающие при этом конечномерные алгебры.
В этой главе развивается конечномерная теория — с полными доказательствами там, где они коротки, — и ставится вопрос о том,
что упускает йорданова структура. Ответ — знак коммутатора: йорданова алгебра одной системы определяет коммутатор лишь с
точностью до знака (раздел 13.5), но составным системам этот знак нужен, и именно он решает, является ли
состояние запутанным.

**Вы узнаете** йорданово произведение и квадратичное представление $P(y)x=yxy$, которое, в отличие от $\circ$,
сохраняет положительность; спин-факторы, их спектральную теорию и конус Лоренца; йордановы подалгебры эрмитовых
матриц, не являющиеся самосопряжёнными частями $*$-алгебр; классификацию Йордана–фон Неймана–Вигнера; почему
транспонирование — йорданов автоморфизм, но не $*$-автоморфизм, и что оно делает на составной системе; а также следовое
неравенство из главы 11, записанное на йордановом языке и доказанное в спин-факторе.

Всюду далее $\operatorname{Herm}_n(\mathbb K)$ — вещественное векторное пространство эрмитовых матриц $n\times n$ над
$\mathbb K=\mathbb R$ или $\mathbb C$, $M_{sa}$ — самосопряжённая часть $\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$,
$X,Y,Z$ — матрицы Паули и $\sigma=(X,Y,Z)$, так что $u\cdot\sigma=u_1X+u_2Y+u_3Z$ для $u\in\mathbb R^3$.
$[a,b]=ab-ba$, а $|u|$ — евклидова норма $u\in\mathbb R^m$.

## 13.1 Наблюдаемые и йорданово произведение

**Определение.** **Йордановой алгеброй** называется вещественное векторное пространство $J$ с билинейным коммутативным произведением $x\circ y$,
удовлетворяющим **йорданову тождеству** $(x^2\circ y)\circ x=x^2\circ(y\circ x)$, где $x^2=x\circ x$. Она
**специальна**, если изоморфна подпространству некоторой ассоциативной алгебры, замкнутому относительно $\frac12(xy+yx)$. На $J$ положим
$L(y)x=y\circ x$; **квадратичное представление** — это $P(y)=2L(y)^2-L(y^2)$.

**Лемма 13.1.** Пусть $\mathcal A$ — ассоциативная вещественная или комплексная алгебра и $x\circ y=\frac12(xy+yx)$. Тогда
$\circ$ удовлетворяет йорданову тождеству, а йордановы степени $x^{\circ1}=x$, $x^{\circ(k+1)}=x\circ x^{\circ k}$ совпадают
с ассоциативными степенями $x^k$. В частности, $M_{sa}$, $\operatorname{Herm}_n(\mathbb C)$ и
$\operatorname{Herm}_n(\mathbb R)$ — специальные йордановы алгебры.

*Доказательство.* Раскрывая скобки, получаем $4(x^2\circ y)\circ x=x^2yx+yx^3+x^3y+xyx^2$ и
$4x^2\circ(y\circ x)=x^2yx+x^3y+yx^3+xyx^2$ — те же четыре слагаемых. По индукции
$x\circ x^k=\frac12(x^{k+1}+x^{k+1})=x^{k+1}$. Если $x,y$ эрмитовы, то $(xy+yx)^*=yx+xy$, так что все три пространства замкнуты
относительно $\circ$; для вещественных симметричных $x,y$ произведение вещественно. $\square$

**Предложение 13.2.** Пусть $x,y,z\in M_{sa}$.

1. $P(y)x=yxy$, и выполняется **фундаментальная формула** $P(P(y)x)=P(y)P(x)P(y)$.
2. $\operatorname{Tr}P(y)x=\operatorname{Tr}y^2x$.
3. Если $x\ge0$, то $P(y)x\ge0$. Если $y$ обратим, то $P(y)$ — биекция положительного конуса на себя с
   обратным отображением $P(y^{-1})$.
4. Йорданово произведение положительных элементов не обязано быть положительным: для $p=\begin{pmatrix}1&0\\0&0\end{pmatrix}$ и
   $q=\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}$ имеем $\det(p\circ q)=-\frac1{16}$.

*Доказательство.* (1) $2y\circ(y\circ x)=\frac12(y^2x+2yxy+xy^2)$ и $y^2\circ x=\frac12(y^2x+xy^2)$; вычтем. Далее
$P(P(y)x)z=yxy\,z\,yxy=P(y)P(x)P(y)z$ по ассоциативности. (2) Цикличность следа. (3)
$yxy=(x^{1/2}y)^*(x^{1/2}y)\ge0$ (следствие 1.7) и $P(y)P(y^{-1})=P(y^{-1})P(y)=\mathrm{id}$. (4)
$p\circ q=\frac14\begin{pmatrix}2&1\\1&0\end{pmatrix}$. $\square$

Если записать $P(y)$ как $2L(y)^2-L(y^2)$, каждое тождество в (1) выражается только через $\circ$, поэтому оно выполняется в любой специальной йордановой алгебре; фундаментальная
формула выполняется в любой йордановой алгебре (K. McCrimmon, *A Taste of Jordan Algebras*, Springer 2004; это доказательство мы не
воспроизводим). Модуль `torch_vn_algebra.jordan` реализует эти отображения на операторах `TypeIAlgebra` — лениво
и пакетно по всем секторам:

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

> **Физический смысл.** $x\circ y$ — симметризованное произведение наблюдаемых: единственная комбинация
> $\alpha xy+\beta yx$, симметричная по $x,y$ и сводящаяся к $x^2$ при $y=x$, и она снова является наблюдаемой. Квадратичное
> представление — это **правило Людерса**: измерение вопроса «да/нет» $p$ с сохранением исхода «да» переводит
> плотность $\rho$ в $P(p)\rho/\operatorname{Tr}P(p)\rho=p\rho p/\operatorname{Tr}(p\rho)$. Предложение 13.2(4)
> говорит, что само $\circ$ не имеет такого операционального смысла для положительных элементов.

## 13.2 Спин-факторы и конус Лоренца

**Определение.** Для $m\ge1$ **спин-фактор** $V_m=\mathbb R\oplus\mathbb R^m$ состоит из элементов $x=(a,u)$ с
произведением $(a,u)\circ(b,v)=(ab+u\cdot v,\;av+bu)$ и единицей $1=(1,0)$. Его **след** равен $\operatorname{tr}x=2a$.

Произведение билинейно и коммутативно; йорданово тождество следует из матричного представления.

**Лемма 13.3 (антикоммутирующие симметрии).** Пусть $p=\lfloor m/2\rfloor$ и $D=2^p$ ($D=2$ при $m=1$). На
$(\mathbb C^2)^{\otimes p}$ положим $s_{2j-1}=Z^{\otimes(j-1)}\otimes X\otimes1$, $s_{2j}=Z^{\otimes(j-1)}\otimes Y\otimes1$
при $j\le p$ и $s_m=Z^{\otimes p}$, если $m$ нечётно ($s_1=Z$ при $m=1$). Тогда $s_i=s_i^*$ и
$s_is_j+s_js_i=2\delta_{ij}1$. Отображение $\phi(a,u)=a1+\sum_iu_is_i$ инъективно и удовлетворяет
$\phi(x\circ y)=\phi(x)\circ\phi(y)$; следовательно, $V_m$ — специальная йорданова алгебра. При $m\ge2$ элементы $s_i$ порождают
$M_D(\mathbb C)$ как $*$-алгебру.

*Доказательство.* Квадраты матриц Паули равны $1$, и они попарно антикоммутируют; $Z$ антикоммутирует с $X$ и $Y$. Две различные строки
$s_i,s_j$ отличаются ровно в одном узле, причём антикоммутирующим образом: для $s_{2j-1},s_{2j}$ — в узле $j$ ($X$ против
$Y$); для строк с индексами $j<j'$ — в узле $j$ ($X$ или $Y$ против $Z$), тогда как в остальных узлах они несут одинаковые
множители или множитель $1$; то же верно для $s_m=Z^{\otimes p}$. Поэтому $s_is_j=-s_js_i$ при $i\ne j$ и
$s_i^2=1$. Значит, $(u\cdot s)\circ(v\cdot s)=\sum_{ij}u_iv_j\,s_i\circ s_j=(u\cdot v)1$, откуда
$\phi(x)\circ\phi(y)=\phi(x\circ y)$. Каждое $s_i$ имеет нулевой след (тензорное произведение с бесследовым множителем), и
$\operatorname{Tr}(s_is_j)=\operatorname{Tr}(s_i\circ s_j)=D\delta_{ij}$, так что $1,s_1,\dots,s_m$ ортогональны относительно
$\operatorname{Tr}(ab)$ и $\phi$ инъективно. Йорданово тождество выполняется в $\operatorname{Herm}_D$ (лемма 13.1) и
переносится обратно с помощью $\phi$. Наконец, $s_{2j-1}s_{2j}=1\otimes\cdots\otimes XY\otimes\cdots\otimes1=iZ_j$ ($Z$ в узле
$j$), так что порождённая $*$-алгебра содержит $Z_j$, затем $X_j=Z_1\cdots Z_{j-1}s_{2j-1}$ и $Y_j$, а значит, и все
тензорные произведения матриц Паули, образующие базис $M_D(\mathbb C)$. $\square$

С точностью до знаков $s_{2j-1}$ и $s_{2j}$ — это майорановские операторы $c_j+c_j^*$ и $i(c_j-c_j^*)$ фермионов
Йордана–Вигнера $c_j$ из раздела 9.4.

**Теорема 13.4 (спектральная теория $V_m$).** Пусть $x=(a,u)\in V_m$ и $\lambda_\pm=a\pm|u|$.

1. Идемпотенты $V_m$ — это $0$, $1$ и $c_w=\frac12(1,w)$ с $|w|=1$; при этом $c_w\circ c_{w'}=0$ тогда и только тогда, когда $w'=-w$.
2. Если $u\ne0$, то $x=\lambda_+e_++\lambda_-e_-$, где $e_\pm=c_{\pm u/|u|}$; если $u=0$, то $x=a1$. Если $x=\mu_1c_w+\mu_2c_{-w}$,
   то $\{\mu_1,\mu_2\}=\{\lambda_+,\lambda_-\}$. Числа $\lambda_\pm$ мы называем **собственными значениями** $x$.
3. Для $f:\{\lambda_+,\lambda_-\}\to\mathbb R$ положим $f(x)=f(\lambda_+)e_++f(\lambda_-)e_-$ ($f(a)1$ при $u=0$). Тогда
   $x^{\circ k}=\lambda_+^ke_++\lambda_-^ke_-$, $\phi(f(x))=f(\phi(x))$, собственные значения $\phi(x)$ равны $\lambda_\pm$,
   каждое кратности $D/2$, и $\operatorname{Tr}\phi(x)=\frac D2\operatorname{tr}x$. Кроме того,
   $\operatorname{tr}x=\lambda_++\lambda_-$ и $\det x:=\lambda_+\lambda_-=a^2-|u|^2$.
4. Конус квадратов $Q=\{x^2:x\in V_m\}$ — это **конус Лоренца** $\{(a,u):a\ge|u|\}=\{x:\lambda_-\ge0\}=\{x:\phi(x)\ge0\}$.
   Он **самодвойствен**: $y\in Q$ тогда и только тогда, когда $\operatorname{tr}(x\circ y)\ge0$ для всех $x\in Q$.
5. $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$; следовательно, $P(y)Q\subseteq Q$ для всякого $y\in V_m$.

*Доказательство.* (1) $(a,u)^2=(a^2+|u|^2,2au)$ равно $(a,u)$ тогда и только тогда, когда $u=0$, $a\in\{0,1\}$, или $a=\frac12$, $|u|=\frac12$.
$c_w\circ c_{w'}=\frac14(1+w\cdot w',\,w+w')$ обращается в нуль тогда и только тогда, когда $w'=-w$. (2)
$\lambda_+e_++\lambda_-e_-=(a,0)+|u|(0,u/|u|)=x$, и $\mu_1c_w+\mu_2c_{-w}=(\frac{\mu_1+\mu_2}2,\frac{\mu_1-\mu_2}2w)$, поэтому
$\mu_1+\mu_2=2a$ и $|\mu_1-\mu_2|=2|u|$. (3) В силу (1) $(\alpha e_++\beta e_-)\circ(\gamma e_++\delta e_-)=\alpha\gamma
e_++\beta\delta e_-$, что по индукции даёт степени. При $w=u/|u|$ имеем $(w\cdot s)^2=1$, так что
$\phi(e_\pm)=\frac12(1\pm w\cdot s)$ — ортогональные проекторы с суммой $1$ и произведением $0$, каждый со следом $D/2$; значит,
$\phi(x)=\lambda_+\phi(e_+)+\lambda_-\phi(e_-)$ — спектральное разложение, и
$f(\phi(x))=f(\lambda_+)\phi(e_+)+f(\lambda_-)\phi(e_-)=\phi(f(x))$ (теорема 1.5). (4) У $x^2$ собственные значения
$\lambda_\pm^2\ge0$; обратно, если $\lambda_-\ge0$, то $x=(\sqrt x)^2$. $\lambda_-\ge0$ тогда и только тогда, когда $a\ge|u|$, и, в силу (3), тогда и только тогда, когда
$\phi(x)\ge0$. Для самодвойственности: $\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$. Если $a\ge|u|$ и $b\ge|v|$, то
$ab\ge|u||v|\ge-u\cdot v$. Если $b<|v|$, то $x=(|v|,-v)\in Q$ даёт $2(|v|b-|v|^2)<0$ (при $v=0$, $b<0$ возьмём $x=1$).
(5) $P(y)x=2y\circ(y\circ x)-y^2\circ x$, и $\phi$ сохраняет $\circ$, поэтому $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ по
предложению 13.2(1); остаётся применить (4) и предложение 13.2(3). $\square$

**Предложение 13.5.** (a) При $m=3$ отображение $\phi(a,u)=a1+u\cdot\sigma$ — йорданов изоморфизм
$V_3\cong\operatorname{Herm}_2(\mathbb C)$. (b) Йордановы автоморфизмы $V_m$ — это отображения $(a,u)\mapsto(a,Ou)$,
$O\in O(m)$. (c) При $m\ge2$ размерность $\phi(V_m)$ равна $m+1$, а размерность самосопряжённой части порождённой им
$*$-алгебры равна $D^2=4^{\lfloor m/2\rfloor}$; они совпадают только при $m=3$.

*Доказательство.* (a) При $m=3$ имеем $(s_1,s_2,s_3)=(X,Y,Z)$, и $\phi$ — инъекция между пространствами размерности $4$. (b) Отображения
$(a,u)\mapsto(a,Ou)$ сохраняют $u\cdot v$, а значит, и $\circ$. Обратно, пусть $\alpha$ — йорданов автоморфизм. Тогда
$\alpha(1)$ — единица $V_m$, так что $\alpha(1)=1$. По формуле для $x^2$ равенство $x^2=1$ выполняется тогда и только тогда, когда $x=\pm1$ или $x=(0,w)$ с
$|w|=1$. При $|w|=1$ имеем $\alpha(0,w)^2=1$ и $\alpha(0,w)\ne\pm1=\alpha(\pm1)$, так что $\alpha(0,w)=(0,w')$ с $|w'|=1$.
Записав $\alpha(0,u)=(\ell(u),Ou)$ с линейными $\ell,O$, получаем $\ell=0$ и $|Ou|=|u|$, т. е. $O\in O(m)$.
(c) Лемма 13.3; $4^{\lfloor m/2\rfloor}>m+1$ при $m=2$ и $m\ge4$. $\square$

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

При $m=1$ элемент $s_1=Z$ порождает диагональные матрицы, размерность которых $2=m+1$, а $V_1\cong\mathbb R\oplus\mathbb R$ не
прост. Напечатанные размерности: $3,4$ при $m=2$; $5,16$ при $m=4$; $7,64$ при $m=6$.

> **Физический смысл.** Состояние на $V_m$ (положительный нормированный линейный функционал) в силу самодвойственности имеет вид
> $x\mapsto\operatorname{tr}(\rho\circ x)$ с $\rho=\frac12(1,r)\in Q$, $|r|\le1$: пространство состояний — шар
> размерности $m$. При $m=3$ это шар Блоха кубита (предложение 13.5(a)). Конус Лоренца — это световой конус
> будущего специальной теории относительности в $1+m$ измерениях, а предложение 13.5(b) говорит, что йордановы автоморфизмы —
> это ортогональные преобразования пространства, оставляющие на месте ось времени.

Конусы $Q$ задают ограничения в **коническом программировании второго порядка** (SOCP): линейной оптимизации на
произведениях конусов Лоренца. Методы внутренней точки для SOCP используют в точности теорему 13.4: собственные значения, $\det$ и
$P(y)$ для масштабирования (F. Alizadeh, D. Goldfarb, *Second-order cone programming*, Math. Program. 95 (2003) 3–51).

## 13.3 Йордановы подалгебры, не являющиеся *-подалгебрами

Йорданова подалгебра $J\subseteq\operatorname{Herm}_n(\mathbb C)$ порождает комплексную $*$-алгебру $\mathcal A$, причём
$J\subseteq\mathcal A_{sa}$. Предложение 13.5(c) дало примеры с $J\ne\mathcal A_{sa}$; вот ещё два.

**Предложение 13.6.** (a) $\operatorname{Herm}_n(\mathbb R)$ — йорданова подалгебра $\operatorname{Herm}_n(\mathbb C)$
размерности $n(n+1)/2$; порождённая ею $*$-алгебра есть $M_n(\mathbb C)$, самосопряжённая часть которой имеет размерность $n^2$.
(b) Для $x,y\in\operatorname{Herm}_2(\mathbb C)$ существует унитарная $U$, для которой $UxU^*$ и $UyU^*$ вещественны и
симметричны. Йорданова алгебра, порождённая $1,x,y$, есть $\operatorname{span}\{1,x,y\}$; если $x,y,1$ линейно
независимы, то это $U^*\operatorname{Herm}_2(\mathbb R)U\cong V_2$ размерности $3$, тогда как $*$-алгебра, порождённая
$x,y$, есть $M_2(\mathbb C)$ размерности $4$.

*Доказательство.* (a) Лемма 13.1. При $i\ne j$ матрица $E_{ij}=E_{ii}(E_{ij}+E_{ji})$ — произведение вещественных симметричных матриц, так что
порождённая алгебра содержит все матричные единицы. (b) Запишем $x=a1+u\cdot\sigma$, $y=b1+v\cdot\sigma$ и выберем единичный вектор
$n\perp u,v$. Матрицы $n\cdot\sigma$ и $Y$ эрмитовы, бесследовы и в квадрате дают $1$, поэтому обе имеют
собственные значения $1,-1$, и $U(n\cdot\sigma)U^*=Y$ для некоторой унитарной $U$ (теорема 1.3). Сопряжение с помощью $U$ сохраняет
$\operatorname{Tr}$ и $\circ$, поэтому переводит $u\cdot\sigma$ в $(Ru)\cdot\sigma$ с линейным $R$ и
$Ru\cdot Rv=\frac12\operatorname{Tr}\big((Ru\cdot\sigma)(Rv\cdot\sigma)\big)=u\cdot v$: $R\in O(3)$, и $Rn=e_2$.
Следовательно, $Ru,Rv\perp e_2$, и $UxU^*=a1+(Ru)_1X+(Ru)_3Z$ вещественна, поскольку таковы $X,Z$; аналогично для $UyU^*$. Далее,
$x\circ y=(ab+u\cdot v)1+av\cdot\sigma+bu\cdot\sigma\in\operatorname{span}\{1,x,y\}$, и так же для $x^2,y^2$, так что линейная оболочка
замкнута. Если её размерность равна $3$, то её сопряжение с помощью $U$ есть всё $\operatorname{Herm}_2(\mathbb R)=\operatorname{span}\{1,X,Z\}$,
которое $\cong V_2$ по лемме 13.3 ($X,Z$ антикоммутируют) и порождает $M_2(\mathbb C)$ по (a). $\square$

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

Три эрмитовы матрицы $2\times2$ общего положения порождают всё $\operatorname{Herm}_2(\mathbb C)$: размерность $4$ в обоих
столбцах. Поэтому в блоке $M_2$ алгебры $\mathcal M$ две наблюдаемые никогда не выявляют комплексную структуру; третья
может.

## 13.4 Классификация

Йорданова алгебра называется **формально вещественной**, если из $\sum_ix_i^2=0$ следует, что все $x_i=0$; в конечной размерности это
равносильно тому, что она **евклидова**, т. е. обладает скалярным произведением с $\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$
[FK94]. $M_{sa}$ (с $\operatorname{Tr}(xy)$) и $V_m$ (с $\operatorname{tr}(x\circ y)$)
евклидовы.

**Теорема 13.7 (Йордан–фон Нейман–Вигнер; без доказательства).** Всякая конечномерная формально вещественная йорданова
алгебра есть прямая сумма простых, а всякая простая изоморфна ровно одной из следующих: $\mathbb R$;
$\operatorname{Herm}_r(\mathbb K)$ при $\mathbb K=\mathbb R,\mathbb C,\mathbb H$ и $r\ge3$; спин-факторы
$V_m$, $m\ge2$; **алгебра Альберта** $\operatorname{Herm}_3(\mathbb O)$ эрмитовых матриц $3\times3$ над
октонионами (размерности $27$). В ранге $2$ имеем $\operatorname{Herm}_2(\mathbb K)\cong V_{\dim\mathbb K+1}$ для
$\mathbb K=\mathbb R,\mathbb C,\mathbb H,\mathbb O$. (P. Jordan, J. von Neumann, E. Wigner, *On an algebraic
generalization of the quantum mechanical formalism*, Ann. of Math. 35 (1934) 29–64; современное доказательство см. в [FK94,
Chapter V].)

Все они, кроме алгебры Альберта, специальны: $\operatorname{Herm}_r(\mathbb H)$ вкладывается в
$\operatorname{Herm}_{2r}(\mathbb C)$, если записать кватернионы как комплексные матрицы $2\times2$, а $V_m$ — по лемме 13.3.
Алгебра Альберта **исключительна**: она не изоморфна йордановой подалгебре никакой ассоциативной алгебры
(A. A. Albert, *On a certain algebra of quantum mechanics*, Ann. of Math. 35 (1934) 65–73). Причина, которую можно
изложить кратко: К. М. Гленни нашёл полиномиальные тождества, которые выполняются в любой специальной йордановой алгебре, но нарушаются
в $\operatorname{Herm}_3(\mathbb O)$ (*Some identities valid in special Jordan algebras but not valid in all Jordan
algebras*, Pacific J. Math. 16 (1966) 47–59; см. McCrimmon, *A Taste of Jordan Algebras*). Октонионы не
ассоциативны, и $\operatorname{Herm}_r(\mathbb O)$ удовлетворяет йорданову тождеству лишь при $r\le3$, поэтому
список на этом и обрывается. Итак, $\operatorname{Herm}_3(\mathbb O)$ нельзя точно представить вещественными или комплексными матрицами
(или внутри какой-либо ассоциативной алгебры) с произведением $\frac12(xy+yx)$; программам приходится работать с октонионными элементами.

**Теорема 13.8 (Кёхер–Винберг; без доказательства).** Внутренности конусов квадратов евклидовых йордановых
алгебр — это в точности **симметрические конусы**: открытые выпуклые конусы в евклидовом пространстве, которые самодвойственны и
однородны (их линейные автоморфизмы действуют транзитивно). [FK94, Chapter III.]

Теорема 13.4(4) и упражнение 1 доказывают оба свойства для конуса Лоренца. Таким образом, конусы полуопределённого программирования
(произведения $\operatorname{Herm}_n(\mathbb R)^+$ или $\operatorname{Herm}_n(\mathbb C)^+$), SOCP и линейного
программирования ($\mathbb R_+^n$ — конус $\mathbb R^n=\mathbb R\oplus\dots\oplus\mathbb R$) все имеют этот вид, и
методы внутренней точки можно написать один раз для всех евклидовых йордановых алгебр.

[FK94]: J. Faraut, A. Korányi, *Analysis on Symmetric Cones*, Oxford University Press, 1994.

## 13.5 Йордановы отображения и *-отображения

**Йордановым $*$-гомоморфизмом** алгебр фон Неймана называется комплексно-линейное отображение $\Phi$ с $\Phi(x^*)=\Phi(x)^*$ и
$\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$; на самосопряжённых частях это йорданов гомоморфизм, и всякий йорданов
гомоморфизм самосопряжённых частей продолжается до йорданова $*$-гомоморфизма по формуле $\Phi(x+iy)=\Phi(x)+i\Phi(y)$.

**Предложение 13.9 (транспонирование).** Пусть $T(A)=\bigoplus_cA_c^T$ на $\mathcal M$. (a) $T$ — унитальный йорданов
$*$-автоморфизм с $T(AB)=T(B)T(A)$; он сохраняет спектры, поэтому $T$ и $T^{-1}=T$ положительны. (b) На $M_k$,
$k\ge2$, матрица Чоя $\sum_{ij}E_{ij}\otimes T(E_{ij})$ есть оператор перестановки $F(e_a\otimes e_b)=e_b\otimes e_a$ с
собственными значениями $+1$ (кратности $k(k+1)/2$) и $-1$ (кратности $k(k-1)/2$). Поэтому $T$ не является вполне положительным.

*Доказательство.* (a) $(A^T)^*=(A^*)^T$, $(AB)^T=B^TA^T$, следовательно, $(A\circ B)^T=B^T\circ A^T=A^T\circ B^T$;
$\det(A^T-\lambda)=\det(A-\lambda)$. (b) $\sum_{ij}E_{ij}e_a\otimes E_{ji}e_b=e_b\otimes e_a$. $F^2=1$, и $F$ равен $+1$
на симметричных и $-1$ на антисимметричных тензорах, которые вместе порождают всё пространство. По теореме 6.2 (и предложению 6.1) $T$ не является
CP. $\square$

Итак, йорданов автоморфизм не обязан быть $*$-автоморфизмом. С точностью до этого явления положительность определяет йорданову
структуру:

**Теорема 13.10 (Кадисон; без доказательства).** Пусть $\Phi$ — линейная биекция между унитальными C\*-алгебрами с
$\Phi(1)=1$, причём $\Phi$ и $\Phi^{-1}$ положительны. Тогда $\Phi$ — йорданов $*$-изоморфизм. (R. V. Kadison,
*A generalized Schwarz inequality and algebraic invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503,
Corollary 5. В работе *Isometries of operator algebras*, Ann. of Math. 54 (1951) 325–338, он показал, что линейная
изометрия $\Phi$ одной унитальной C\*-алгебры на другую есть произведение унитарного элемента $\Phi(1)$ на йорданов $*$-изоморфизм.)

Обратное очевидно: йорданов $*$-изоморфизм переводит квадраты $h^2=h\circ h$ самосопряжённых элементов в квадраты,
а значит, положительный конус — на положительный конус. Для факторов $M_n$ йордановы автоморфизмы известны точно:
всякий йорданов $*$-автоморфизм $M_n(\mathbb C)$ имеет вид $A\mapsto UAU^*$ или $A\mapsto UA^TU^*$ с унитарной $U$
(N. Jacobson, C. E. Rickart, *Jordan homomorphisms of rings*, Trans. Amer. Math. Soc. 69 (1950) 479–502; I. N.
Herstein, *Jordan homomorphisms*, Trans. Amer. Math. Soc. 81 (1956) 331–341). Мы докажем случай $n=2$.

**Предложение 13.11.** Всякий йорданов автоморфизм $\operatorname{Herm}_2(\mathbb C)$ имеет вид $A\mapsto UAU^*$ или
$A\mapsto UA^TU^*$ с $U\in U(2)$.

*Доказательство.* По предложению 13.5(a,b) йорданов автоморфизм действует как $a1+u\cdot\sigma\mapsto a1+(Ou)\cdot\sigma$ с
$O\in O(3)$. Обозначим через $R_U\in O(3)$ матрицу отображения $A\mapsto UAU^*$ в этом смысле; $R_{UV}=R_UR_V$. (i) $\det R_U=1$:
$U=W\operatorname{diag}(e^{i\theta_1},e^{i\theta_2})W^*$ соединяется с $1$ путём
$U_t=W\operatorname{diag}(e^{it\theta_1},e^{it\theta_2})W^*$, а $\det R_{U_t}\in\{\pm1\}$ непрерывно зависит от $t$.
(ii) Всякое $R\in SO(3)$ есть некоторое $R_U$: как в предложении 13.6, существует $U_1$ с
$U_1ZU_1^*=(Re_3)\cdot\sigma$, т. е. $R_{U_1}e_3=Re_3$. Тогда $R_{U_1}^{-1}R\in SO(3)$ оставляет $e_3$ на месте, так что это поворот
на некоторый угол $\beta$ вокруг $e_3$. Прямое вычисление даёт для $U_2=\operatorname{diag}(e^{-i\beta/2},e^{i\beta/2})$
$U_2XU_2^*=\cos\beta\,X+\sin\beta\,Y$ и $U_2YU_2^*=-\sin\beta\,X+\cos\beta\,Y$, так что $R=R_{U_1U_2}$. (iii)
Транспонирование оставляет на месте $X,Z$ и переводит $Y\mapsto-Y$: $O_T=\operatorname{diag}(1,-1,1)$, $\det O_T=-1$. Следовательно, $O=R_U$, если
$\det O=1$, и $OO_T=R_U$, т. е. автоморфизм имеет вид $A\mapsto UA^TU^*$, если $\det O=-1$. $\square$

**Составные системы.** Йорданова структура фактора не позволяет отличить $T$ от тождественного отображения. Для двух систем эта
неоднозначность становится видимой:

**Предложение 13.12.** Пусть $x,x'\in\operatorname{Herm}_n(\mathbb C)$, $y,y'\in\operatorname{Herm}_m(\mathbb C)$. В
$\operatorname{Herm}_{nm}(\mathbb C)$, действующем на $\mathbb C^n\otimes\mathbb C^m$,
$$
(x\otimes y)\circ(x'\otimes y')=(x\circ x')\otimes(y\circ y')+\tfrac14[x,x']\otimes[y,y'].
$$
Пусть $\theta=\mathrm{id}\otimes T$ (**частичное транспонирование**) и $A\circ_\theta B=\theta(\theta(A)\circ\theta(B))$. Тогда
$\circ_\theta$ — йорданово произведение на $\operatorname{Herm}_{nm}(\mathbb C)$, изоморфное $\circ$ посредством $\theta$, которое
совпадает с $\circ$ на всех парах локальных наблюдаемых $x\otimes1$, $1\otimes y$, и
$(x\otimes y)\circ_\theta(x'\otimes y')=(x\circ x')\otimes(y\circ y')-\frac14[x,x']\otimes[y,y']$. Его конус квадратов есть
$\theta(\operatorname{Herm}_{nm}^+)$. При $n=m=2$ и $\Omega=e_1\otimes e_1+e_2\otimes e_2$ состояние
$\omega=\frac12|\Omega\rangle\langle\Omega|\ge0$ не лежит в этом конусе: $\theta(\omega)=\frac12F$ имеет собственное значение
$-\frac12$.

*Доказательство.* Положим $j=x\circ x'$, $c=\frac12[x,x']$, $k=y\circ y'$, $d=\frac12[y,y']$; тогда $xx'=j+c$, $x'x=j-c$, и аналогично для
$y$, так что $\frac12(xx'\otimes yy'+x'x\otimes y'y)=\frac12\big((j+c)\otimes(k+d)+(j-c)\otimes(k-d)\big)=j\otimes k+c\otimes d$.
$\theta$ — вещественно-линейная инволюция $\operatorname{Herm}_{nm}$: она переводит $x\otimes y$ в $x\otimes y^T$, а
элементы $x\otimes y$ порождают $\operatorname{Herm}_{nm}$ над $\mathbb R$ (произведения $h_\alpha\otimes g_\beta$ вещественных базисов
$\operatorname{Herm}_n,\operatorname{Herm}_m$ образуют комплексный базис $M_{nm}$, и эрмитова матрица имеет в нём вещественные
коэффициенты). Поэтому $\circ_\theta$ — йорданово произведение, а $\theta$ — изоморфизм. Поскольку $(y^T\circ y'^T)^T=y\circ y'$ и
$[y^T,y'^T]^T=-[y,y']$, применение $\theta$ к формуле для $(x\otimes y^T)\circ(x'\otimes y'^T)$ даёт формулу для
$\circ_\theta$. На локальных наблюдаемых коммутаторный член исчезает ($[1,1]=0$, $[x,1]=0$). Квадраты суть
$A\circ_\theta A=\theta(\theta(A)^2)$, и $\omega$ лежит в $\theta(\operatorname{Herm}_{nm}^+)$ тогда и только тогда, когда $\theta(\omega)\ge0$,
поскольку $\theta^2=\mathrm{id}$. Наконец, $\theta(\omega)=\frac12\sum_{ij}E_{ij}\otimes E_{ji}=\frac12F$
(предложение 13.9(b); это вычисление из доказательства предложения 6.1). $\square$

Оба произведения задают одинаково хорошие йордановы алгебры, они согласуются во всём, что можно измерить на одной подсистеме
за раз, и расходятся в том, какие совместные операторы положительны, т. е. в том, каковы состояния
составной системы. Различает их знак коммутаторного члена, и этот знак фиксируется
$*$-алгеброй $M_n\otimes M_m$, а не йордановыми алгебрами сомножителей. ($\circ_\theta$ — йорданово произведение
ассоциативного произведения $A\star B=\theta(\theta(A)\theta(B))$, для которого $(x\otimes y)\star(x'\otimes y')=xx'\otimes y'y$:
это алгебра $M_n\otimes M_m^{\rm op}$.) Именно конус $\theta(\operatorname{Herm}_{nm}^+)$
даёт **критерий PPT** (A. Peres, Phys. Rev. Lett. 77 (1996) 1413): сепарабельное состояние
$\rho=\sum_kp_ka_k\otimes b_k$, $p_k\ge0$, $a_k,b_k\ge0$, лежит в обоих конусах, поскольку
$\theta(\rho)=\sum_kp_ka_k\otimes b_k^T$ с $b_k^T\ge0$ (предложение 13.9(a)) и
$a\otimes b=(a^{1/2}\otimes b^{1/2})^2\ge0$. Поэтому состояние $\omega$, не лежащее в $\theta(\operatorname{Herm}_4^+)$,
запутано. О том, какие составные системы йордановых алгебр существуют, см. H. Hanche-Olsen, *On the structure and tensor
products of JC-algebras*, Canad. J. Math. 35 (1983) 1059–1074, и H. Barnum, M. Graydon, A. Wilce,
*Composites and categories of Euclidean Jordan algebras*, Quantum 4 (2020) 359.

**Динамика.** Можно было бы надеяться, что эволюция во времени вынуждает ассоциативное произведение. Это не так:

**Предложение 13.13.** Пусть $H\in\operatorname{Herm}_n(\mathbb C)$ и $\delta(z)=i[H,z]$. (a) $\delta$ переводит
$\operatorname{Herm}_n$ в себя и является дифференцированием относительно $\circ$; $e^{t\delta}(z)=e^{itH}ze^{-itH}$ — йорданов
автоморфизм. (b) $[L(x),L(y)]z=\frac14[[x,y],z]$ для $x,y,z\in\operatorname{Herm}_n$. (c) $\delta=4\sum_k[L(x_k),L(y_k)]$
для некоторых $x_k,y_k\in\operatorname{Herm}_n$: гамильтонова динамика порождается одним лишь йордановым произведением.

*Доказательство.* (a) $(i[H,z])^*=-i[z,H]=i[H,z]$; $[H,\cdot]$ — дифференцирование ассоциативного произведения, а значит, и $\circ$.
(b) $4x\circ(y\circ z)=xyz+xzy+yzx+zyx$; вычтем то же выражение с переставленными $x,y$. (c) В силу (b) достаточно записать
$iH-\frac{\operatorname{Tr}H}n i1$ (с теми же коммутаторами, что и у $iH$) в виде суммы $[x_k,y_k]$. Для $i\ne j$ и
$a=E_{ij}+E_{ji}$, $b=i(E_{ij}-E_{ji})$: $[a,E_{ii}]=E_{ji}-E_{ij}$, $[b,E_{ii}]=-i(E_{ij}+E_{ji})$ и
$[a,b]=2i(E_{jj}-E_{ii})$; они порождают над $\mathbb R$ бесследовые антиэрмитовы матрицы. $\square$

Чего йорданова алгебра не содержит, так это **соответствия** $H\mapsto i[H,\cdot]$ между наблюдаемой и
порождаемой ею динамикой. В $\operatorname{Herm}_n(\mathbb R)$, например, дифференцирования
$[L(x),L(y)]=\frac14[[x,y],\cdot]$ происходят от вещественных антисимметричных $[x,y]$, которые не являются произведением $i$ на вещественную симметричную
матрицу. Э. Альфсен и Ф. Шульц показали, что JBW-алгебра (бесконечномерный вариант формально
вещественной йордановой алгебры; в конечной размерности эти понятия совпадают) является самосопряжённой частью алгебры фон Неймана тогда и только тогда, когда
она допускает такое «динамическое соответствие», и геометрически охарактеризовали пространства состояний C\*-алгебр
(*Geometry of State Spaces of Operator Algebras*, Birkhäuser, 2003; без доказательства). Соответствие не
единственно: на $\operatorname{Herm}_n(\mathbb C)$ подходят и $H\mapsto i[H,\cdot]$, и $H\mapsto-i[H,\cdot]$;
транспонирование переставляет их, поскольку $T(i[H,T(z)])=-i[H^T,z]$, и они относятся к $M_n$ и к противоположной алгебре — это
два знака из предложения 13.12.

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

> **Физический смысл.** Для бесспиновой частицы транспонирование (комплексное сопряжение в вещественном базисе) есть **обращение
> времени**. Это симметрия наблюдаемых одной системы, но она антиунитарна, и, применённая лишь к одной
> подсистеме, не является физической операцией: предложение 13.12 показывает, что она превращает состояние Белла в
> оператор с отрицательным собственным значением. Квантовая механика использует комплексные $*$-алгебры потому, что составным системам
> нужен коммутатор со своим знаком, а не только симметризованное произведение.

## 13.6 Снова о следовых неравенствах

Предложение 11.9 (эксп. 3) утверждает, что $\operatorname{Tr}|YXY|\le\operatorname{Tr}(Y|X|Y)$ для $X=X^*$ и $Y\ge0$ (здесь
$X,Y$ — операторы из $\mathcal M$, а не матрицы Паули). Его доказательство нигде не использует $Y\ge0$: равенство $YX_\pm Y=T^*T$ с
$T=X_\pm^{1/2}Y$ требует лишь $Y=Y^*$, а лемма 11.8 применяется к $YX_+Y$ и $YX_-Y$. Поэтому оно верно для всех $Y=Y^*$
и на йордановом языке записывается как
$$
\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|,
$$
что имеет смысл в любой евклидовой йордановой алгебре. Мы докажем это в $V_m$ рассуждением, использующим лишь три
свойства.

**Предложение 13.14 (эксп. 3 в спин-факторе).** Для $x,y\in V_m$
$\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|$, причём равенство имеет место, если $x\in Q$ или $-x\in Q$.

*Доказательство.* Пусть $x_\pm=f_\pm(x)$ с $f_+(t)=\max(t,0)$, $f_-(t)=\max(-t,0)$ (теорема 13.4(3)): $x_\pm\in Q$, $x=x_+-x_-$ и
$|x|=x_++x_-$. Положим $A=P(y)x_+$ и $B=P(y)x_-$; они лежат в $Q$ по теореме 13.4(5), и $P(y)x=A-B$. Пусть
$P(y)x=\mu_1c_1+\mu_2c_2$ — спектральное разложение (теорема 13.4(2); $c_1=c_w$, $c_2=c_{-w}$, так что $c_1+c_2=1$ и
$c_1\circ c_2=0$), $c$ — сумма тех $c_i$, для которых $\mu_i>0$, и
$s=2c-1=\sum_i\operatorname{sgn}'(\mu_i)c_i$, где $\operatorname{sgn}'(t)=1$ при $t>0$ и $-1$ в остальных случаях. Тогда
$s\circ P(y)x=\sum_i|\mu_i|c_i=|P(y)x|$ (в силу теоремы 13.4(2,3) это не зависит от выбора $w$ при
$\mu_1=\mu_2$), и $c,1-c\in Q$ (идемпотенты — квадраты). По лёгкой половине самодвойственности (теорема 13.4(4):
$\operatorname{tr}(a\circ b)\ge0$ для $a,b\in Q$), применённой к $\operatorname{tr}(A\circ(1-c))$ и $\operatorname{tr}(B\circ c)$,
$$
\operatorname{tr}|P(y)x|=\operatorname{tr}(A\circ c)-\operatorname{tr}(A\circ(1-c))-\operatorname{tr}(B\circ c)+\operatorname{tr}(B\circ(1-c))
\le\operatorname{tr}A+\operatorname{tr}B=\operatorname{tr}P(y)|x|.
$$
Если $x\in Q$, то $x_-=0$, $B=0$ и $|P(y)x|=A$, так что обе части равны $\operatorname{tr}A$; аналогично при $-x\in Q$.
$\square$

Доказательство использует лишь следующее: (i) $P(y)Q\subseteq Q$ для всякого $y$; (ii) спектральную теорему: всякое $z$ есть
$\sum_i\mu_ic_i$ с ортогональными идемпотентами $c_i$, в сумме дающими $1$, и $f(z)=\sum_if(\mu_i)c_i$ (так что $x_\pm$ и
идемпотенты — квадраты); (iii) $\operatorname{tr}$ линеен, $\operatorname{tr}(\sum_i\mu_ic_i)=\sum_i\mu_i\operatorname{tr}c_i$
с $\operatorname{tr}c_i\ge0$, и $\operatorname{tr}(a\circ b)\ge0$ для $a,b\in Q$. В евклидовой йордановой алгебре
ранга $r$ свойство (ii) выполняется с $r$ примитивными идемпотентами $c_i$ (йордановым репером), след равен
$\operatorname{tr}z=\sum_i\mu_i$, форма $\operatorname{tr}(x\circ y)$ ассоциативна,
$\operatorname{tr}((x\circ y)\circ z)=\operatorname{tr}(x\circ(y\circ z))$, и выполняется (i) [FK94, Chapters II–III];
эти факты мы не доказываем. Тогда последняя часть (iii) следует из (i): раскрывая $P(a^{1/2})=2L(a^{1/2})^2-L(a)$
и используя ассоциативность, получаем $\operatorname{tr}(a\circ b)=\operatorname{tr}P(a^{1/2})b$ — след элемента $Q$,
а значит, $\ge0$. Таким образом, неравенство выполняется в любой евклидовой йордановой алгебре — включая
$\operatorname{Herm}_r(\mathbb H)$ и алгебру Альберта, которые не изоморфны как йордановы алгебры самосопряжённым частям
комплексных $*$-алгебр. В $M_{sa}$ то же рассуждение с $\operatorname{Tr}$ и $P(y)x=yxy$ заново доказывает
предложение 11.9. В $V_m$ неравенство следует также из предложения 11.9
через $\phi$, поскольку $\phi(P(y)x)=\phi(y)\phi(x)\phi(y)$ и $\phi(|x|)=|\phi(x)|$; код проверяет и то, и другое:

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

Разность $\operatorname{tr}P(y)|x|-\operatorname{tr}|P(y)x|$ обращается в нуль для знакоопределённых $x$ и положительна для типичных
знаконеопределённых $x$ (напечатанные максимумы для этих гауссовых выборок имеют порядок $10^1$–$10^2$).

## 13.7 С чем сравнивать

Насколько нам известно, распространённой численной библиотеки для йордановых алгебр не существует; операции этой главы
обычно пишут вручную в терминах матриц. Алгебраическую сторону покрывает компьютерная алгебра: `JordanAlgebra` в SageMath строит
специальные йордановы алгебры, йордановы алгебры симметричных билинейных форм (к которым относятся спин-факторы) и, по
алгебре октонионов, 27-мерную исключительную йорданову алгебру; в Magma также есть конструкторы йордановых алгебр,
включая исключительную. Это точные алгебраические инструменты, а не пакетные численные. Конусы Лоренца стандартны в
оптимизации: решатели SOCP (ECOS, SCS, Clarabel, MOSEK, доступные через ограничение `SOC` в cvxpy) используют спектральное
разложение $V_m$ (для проекции на конус) или его квадратичное представление (для масштабирования Нестерова–Тодда в
методах внутренней точки), а решатели полуопределённых задач делают то же для $\operatorname{Herm}_n$. Библиотека добавляет йордановы произведения, квадратичное
представление и транспонирование как пакетные дифференцируемые операторы на $\bigoplus_cM_{k_c}$; спин-факторы с их
функциональным исчислением и матричным представлением; и рандомизированное вычисление йордановой алгебры, порождённой данными
матрицами, в сравнении с порождённой ими $*$-алгеброй.

## Упражнения

1. (**Однородность конуса Лоренца.**) Пусть $y\in V_m$ и $\det y\ne0$. Покажите, что $y^{-1}=(a,-u)/\det y$ удовлетворяет
   $y\circ y^{-1}=1$, что $P(y)P(y^{-1})=\mathrm{id}$, что $P(y)$ отображает внутренность $Q$ на себя и что
   $P(x^{1/2})1=x$ для $x$ из внутренности. Выведите отсюда, что внутренность $Q$ однородна.
2. (**Изотропные состояния.**) Для $p\in[0,1]$ пусть $\rho_p=p\,\omega+(1-p)\frac14$ на $\mathbb C^2\otimes\mathbb C^2$
   (предложение 13.12). Вычислите спектр $\theta(\rho_p)$ и покажите, что $\rho_p$ запутано при $p>\frac13$.
3. (**Произведения граничных точек.**) Пусть $|w|=|w'|=1$ в $\mathbb R^m$. Покажите, что $(1,w)\circ(1,w')\in Q$ тогда и только тогда, когда
   $w'=w$ или $w'=-w$: йорданово произведение двух граничных лучей конуса Лоренца почти никогда не лежит в конусе.
4. (**Вещественная квантовая механика.**) Покажите, что для $x,y\in\operatorname{Herm}_n(\mathbb R)$ дифференцирование
   $[L(x),L(y)]$ алгебры $\operatorname{Herm}_n(\mathbb R)$ есть $z\mapsto[K,z]$ с вещественной антисимметричной $K$ и что при
   $n\ge2$ такое отображение никогда не имеет вида $z\mapsto i[H,z]$ с $H\in\operatorname{Herm}_n(\mathbb R)$, если только оно не равно $0$.

<details>
<summary>Решения</summary>

1. $(a,u)\circ(a,-u)=(a^2-|u|^2,0)$, откуда $y\circ y^{-1}=1$. $\phi(y)$ и $\phi(y^{-1})$ коммутируют (оба лежат в
   $\operatorname{span}\{1,u\cdot s\}$), поэтому $\phi(y)\phi(y^{-1})=\phi(y)\circ\phi(y^{-1})=1$, и по теореме 13.4(5)
   $\phi(P(y)P(y^{-1})x)=\phi(x)$; $\phi$ инъективно. $P(y)$ отображает $Q$ в $Q$, как и обратное к нему $P(y^{-1})$, так что
   $P(y)Q=Q$; будучи линейным гомеоморфизмом, оно отображает внутренность на внутренность. Наконец,
   $P(y)1=2y\circ y-y^2=y^2$, так что $P(x^{1/2})1=x$, и $P(x^{1/2})$ обратимо, когда $x$ — внутренняя точка
   ($\det x^{1/2}=\sqrt{\lambda_+\lambda_-}>0$). Любые две внутренние точки переводятся друг в друга отображением
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

2. $\theta(\rho_p)=\frac p2F+\frac{1-p}4$ по предложению 13.12. $F$ имеет собственные значения $1$ (трижды) и $-1$ (один раз), так что
   спектр состоит из $\frac{1+p}4$ (трижды) и $\frac{1-3p}4$, которое отрицательно тогда и только тогда, когда $p>\frac13$; в этом случае $\rho_p$ не
   сепарабельно по критерию PPT. (Для систем $2\times2$ верно и обратное: M., P. and R. Horodecki, Phys.
   Lett. A 223 (1996) 1; мы это не доказываем.)

```python
for pp in (0.2, 1 / 3, 0.5):
    rho_p = pp * omega + (1 - pp) * torch.eye(4, dtype=dt) / 4
    lam_min = torch.linalg.eigvalsh(ptranspose(rho_p, 2, 2))[0].item()
    assert abs(lam_min - (1 - 3 * pp) / 4) < 1e-12
```

3. $(1,w)\circ(1,w')=(1+t,\,w+w')$, где $t=w\cdot w'\in[-1,1]$ и $|w+w'|=\sqrt{2+2t}$. Это лежит в $Q$ тогда и только тогда, когда
   $1+t\ge\sqrt{2(1+t)}$, т. е. когда $t=-1$ (произведение равно $0$) или $1+t\ge2$, т. е. $t=1$.

```python
V3c = SpinFactor(3)
wv = torch.nn.functional.normalize(torch.randn(2, 1000, 3, dtype=torch.float64), dim=-1)
prod = V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[1]))
assert not V3c.in_cone(prod).any()
assert V3c.in_cone(V3c.product(V3c.element(1.0, wv[0]), V3c.element(1.0, wv[0]))).all()
```

4. По предложению 13.13(b) $[L(x),L(y)]z=[K,z]$ с $K=\frac14[x,y]$, и $K^T=\frac14[y^T,x^T]=-K$ вещественна. Для вещественных $H$
   и $z$ выражение $i[H,z]$ чисто мнимо, тогда как $[K,z]$ вещественно, поэтому они совпадают, лишь если оба равны нулю. Таким образом, в вещественной
   квантовой механике генераторы динамики $K$ не являются наблюдаемыми.

</details>

**Что дальше.** Это последняя глава курса. Стандартный источник по конечномерной
теории — [FK94]; K. McCrimmon, *A Taste of Jordan Algebras* (Springer 2004), излагает алгебру и историю, а
H. Hanche-Olsen and E. Størmer, *Jordan Operator Algebras* (Pitman 1984), и Alfsen–Shultz (2003) — бесконечномерную
теорию и её связь с C\*-алгебрами. О роли составных систем в реконструкциях
квантовой теории из операциональных аксиом см. цитированные выше работы Barnum и Wilce. О следовых неравенствах за пределами
главы 11 см. E. Carlen, *Trace inequalities and quantum entropy* (Contemp. Math. 529, 2010). Квадратичное
представление и конус Лоренца снова появляются в конической оптимизации: Alizadeh–Goldfarb (2003) и обзор йордановых
методов внутренней точки в ней.
