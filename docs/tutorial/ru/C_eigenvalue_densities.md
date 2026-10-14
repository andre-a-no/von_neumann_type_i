[English](../C_eigenvalue_densities.md) · [Español](../es/C_eigenvalue_densities.md) · [Français](../fr/C_eigenvalue_densities.md) · [Deutsch](../de/C_eigenvalue_densities.md) · [中文](../zh/C_eigenvalue_densities.md) · [日本語](../ja/C_eigenvalue_densities.md) · **Русский**

# Приложение C. Плотности собственных значений и формула Пейджа

В главе 5 используется совместная плотность собственных значений гауссовых ансамблей, а в главе 7 — формула
Пейджа для средней энтропии запутанности случайного чистого состояния. В этом приложении доказывается и то, и другое.
Путь таков: формула интегрирования Вейля (якобиан «диагонализации»), гауссовы
ансамбли, одноточечная функция детерминантной плотности через ортогональные многочлены, ансамбль
Уишарта–Лагерра и индуцированная мера на редуцированных состояниях и, наконец, точное
вычисление с многочленами Лагерра.

**Вы узнаете**, почему собственные значения инвариантной случайной матрицы отталкиваются по закону
$\prod_{i<j}|\lambda_i-\lambda_j|^\beta$, почему собственные векторы при этом распределены по Хаару и не зависят
от собственных значений, как интегрируется определитель из ортонормированных многочленов и как
$\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}$ получается из нескольких интегралов.

Помимо предварительных сведений курса (теория меры, формула замены переменных в
$\mathbb R^N$ для инъективных $C^1$-отображений с обратимой производной, теоремы об обратной и о неявной
функции), мы используем из главы 5 лемму 5.1 (единственность и правоинвариантность меры Хаара),
лемму 5.2 и теорему 5.3 (существование меры Хаара на $U(n)$ через гауссово QR-разложение), из приложения B —
теорему B.17 (то же для $O(n)$), из приложения A — предложение A.16(a) (унитарная инвариантность комплексных
гауссовых векторов) и теорему A.20(b) (распределение $|z|^2$), а из главы 7 — представление
$\rho_A=GG^*/\operatorname{Tr}GG^*$ из доказательства теоремы 7.11.

**Обозначения.** $\beta\in\{1,2\}$ и соответственно $\mathbb F=\mathbb R$ или $\mathbb C$. $\mathcal H$ —
пространство вещественных симметричных ($\beta=1$) или комплексных эрмитовых ($\beta=2$) матриц размера $n\times n$ с
мерой Лебега $dH$ в координатах $h_{ii}$, $\operatorname{Re}h_{ij}$ и (при $\beta=2$)
$\operatorname{Im}h_{ij}$, $i<j$; его размерность равна $n+\beta n(n-1)/2$. $K=O(n)$ или $U(n)$, $T$ —
группа диагональных элементов $K$ (знаков или фаз), $\mathfrak k=\{X\in M_n(\mathbb F):X^*=-X\}$,
$\mathfrak t$ — её диагональная часть, а $\mathfrak m$ — внедиагональная часть (с нулевой диагональю) с
координатами $\operatorname{Re}X_{ij}$, $\operatorname{Im}X_{ij}$, $i<j$. Для $\lambda\in\mathbb R^n$
$\Lambda=\operatorname{diag}(\lambda)$ и $\Delta(\lambda)=\prod_{i<j}(\lambda_j-\lambda_i)$; $C=\{\lambda:\lambda_1<\dots<\lambda_n\}$,
и на этом множестве $\Delta>0$. Для $H\in\mathcal H$ через $\lambda(H)$ обозначается вектор её собственных значений в порядке
возрастания. $\|X\|_2^2=\operatorname{Tr}X^*X$. Для $V\in K$ отображение $H\mapsto VHV^*$ линейно на $\mathcal H$
и сохраняет $\|\cdot\|_2$; у линейной изометрии $|\det|=1$, поэтому оно **сохраняет $dH$**.

## C.1 Простые спектры и собственные реперы

**Лемма C.1.** (a) Множество $\Omega\subseteq\mathcal H$ матриц с $n$ различными собственными значениями
открыто, а $\mathcal H\setminus\Omega$ — множество лебеговой меры нуль. (b) На $\Omega$ отображение $H\mapsto\lambda(H)$ принадлежит классу
$C^\infty$, а спектральные проекторы $P_j(H)=\prod_{k\neq j}\frac{H-\lambda_k(H)}{\lambda_j(H)-\lambda_k(H)}$ непрерывны.

*Доказательство.* (a) Пусть $A_{ij}=\lambda_j^{i-1}$, где $\lambda=\lambda(H)$. Тогда $\det A=\Delta(\lambda)$ (определитель Вандермонда) и
$(AA^T)_{ij}=\sum_k\lambda_k^{i+j-2}=\operatorname{Tr}H^{i+j-2}$, поэтому
$\operatorname{disc}(H)=\det\big(\operatorname{Tr}H^{i+j-2}\big)_{i,j=1}^n=\Delta(\lambda(H))^2$ — многочлен
от координат $H$, отличный от нуля в точке $\operatorname{diag}(1,\dots,n)$, и множество $\Omega=\{\operatorname{disc}\neq0\}$ открыто.
Множество нулей ненулевого многочлена $p$ на $\mathbb R^N$ имеет меру нуль; это доказывается индукцией по $N$: при $N=1$ оно
конечно; в противном случае запишем $p=\sum_jp_j(x')x_N^j$, где некоторый $p_{j_0}\neq0$; для $x'$ вне множества нулей
$p_{j_0}$ (меры нуль по предположению индукции) функция $x_N\mapsto p(x',x_N)$ имеет конечное число нулей, и теорема Тонелли, применённая
к индикатору множества нулей, даёт меру $0$.
(b) $\chi(H,t)=\det(t-H)=\prod_i(t-\lambda_i)$ — вещественный многочлен от $(H,t)$. В точке $H_0\in\Omega$
$\partial_t\chi(H_0,\lambda_j(H_0))=\prod_{k\neq j}(\lambda_j-\lambda_k)\neq0$, поэтому теорема о неявной функции даёт
$C^\infty$-функции $r_j$ в окрестности $H_0$ с $\chi(H,r_j(H))=0$ и $r_j(H_0)=\lambda_j(H_0)$. По непрерывности
$r_1<\dots<r_n$ вблизи $H_0$: это $n$ различных корней многочлена степени $n$, а значит, все его корни, причём по
порядку, и $r_j=\lambda_j$. Если $H=\sum_k\lambda_kQ_k$ — спектральное разложение, то интерполяция Лагранжа даёт
$P_j(H)=Q_j$; непрерывность следует из непрерывности $\lambda$. $\square$

Для $H\in\Omega$ **собственный репер** — это такой $V\in K$, что $H=V\Lambda V^*$, $\Lambda=\operatorname{diag}\lambda(H)$. Если $V,V'$
— собственные реперы, то $M=V'^*V$ коммутирует с $\Lambda$, так что $(\lambda_i-\lambda_j)M_{ij}=0$ и $M\in T$; обратно, $VD$
— собственный репер для любого $D\in T$. **Борелевский собственный репер** $F(H)$ получается, если взять в качестве $j$-го столбца
$P_j(H)e_m/\|P_j(H)e_m\|$, где $m$ — наименьший индекс, для которого $P_j(H)e_m\neq0$: он непрерывен на каждом из
борелевских множеств, где эти индексы фиксированы, его столбцы — ортонормированные собственные векторы, и он вещественен, если вещественна $H$.

**Лемма C.2 (меры Хаара).** У $U(n)$, $O(n)$ и $T$ есть меры Хаара, и $\int_Kf(VW)\,dV=\int_Kf\,dV$
для всех $W\in K$ и ограниченных борелевских $f$.

*Доказательство.* $U(n)$: теорема 5.3. $O(n)$: теорема B.17 (вещественная гауссова конструкция через QR: вещественная гауссова матрица
почти наверное обратима, и ортогональный множитель её QR-разложения распределён по Хаару). $T$: равномерная мера на $\{\pm1\}^n$ или произведение нормированных мер длины дуги на
окружностях инвариантны относительно умножения на фиксированный элемент. Правоинвариантность — это лемма 5.1. $\square$

**Лемма C.3 (собственные реперы инвариантных ансамблей).** Пусть $H$ — случайный элемент $\mathcal H$ с
$P(H\in\Omega)=1$, распределение которого инвариантно относительно $H\mapsto VHV^*$ для любого $V\in K$. Пусть $D$ равномерно распределён на $T$ и
не зависит от $H$, и $\tilde U=F(H)D$. Тогда $\tilde U$ — собственный репер $H$, он распределён по Хаару на
$K$ и не зависит от $\lambda(H)$.

*Доказательство.* Зафиксируем $V\in K$. И $F(VHV^*)$, и $VF(H)$ — собственные реперы $VHV^*$, поэтому
$F(VHV^*)=VF(H)E(H)$, где $E(H)=(VF(H))^*F(VHV^*)\in T$ — борелевская функция от $H$. Для ограниченных борелевских
$\varphi$ на $\mathbb R^n$ и $\psi$ на $K$, поскольку $VHV^*$ распределена так же, как $H$, и $\lambda(VHV^*)=\lambda(H)$,
$$
\mathbb E\big[\varphi(\lambda(H))\psi(\tilde U)\big]=\mathbb E\big[\varphi(\lambda(H))\psi(VF(H)E(H)D)\big]
=\mathbb E\big[\varphi(\lambda(H))\psi(V\tilde U)\big],
$$
где на втором шаге сначала интегрируем по $D$ (Фубини): при фиксированном $H$ элемент $E(H)D$ равномерно распределён на $T$. Пусть теперь $V$
распределён по Хаару и не зависит от всего остального; по теореме Фубини и лемме C.2
$\mathbb E[\varphi(\lambda(H))\psi(\tilde U)]=\mathbb E\big[\varphi(\lambda(H))\int_K\psi(V\tilde U)\,dV\big]
=\mathbb E[\varphi(\lambda(H))]\int_K\psi\,dV$. При $\varphi=1$ получаем, что $\tilde U$ распределён по Хаару; форма произведения означает
независимость. $\square$

## C.2 Формула интегрирования Вейля

Для $\varepsilon>0$ положим $\mathfrak m_\varepsilon=\{X\in\mathfrak m:\|X\|_2<\varepsilon\}$ и
$$
\Psi(X,\lambda)=e^X\Lambda e^{-X},\qquad (X,\lambda)\in\mathfrak m_\varepsilon\times C .
$$
Для $X\in\mathfrak k$ имеем $(e^X)^*=e^{-X}=(e^X)^{-1}$, поэтому $e^X\in K$.

**Лемма C.4 (карта).** Существуют $\varepsilon>0$ и непрерывная функция $j:\overline{\mathfrak m_\varepsilon}\to(0,\infty)$
с $j(0)=1$, такие что: (a) $|\det D\Psi(X,\lambda)|=\Delta(\lambda)^\beta j(X)$; (b) $\Psi$ инъективно; (c) для
$H\in\Omega$ включение $H\in\Psi(\mathfrak m_\varepsilon\times C)$ равносильно $F(H)\in S=\{e^XD:X\in\mathfrak m_\varepsilon,D\in T\}$, причём
$S$ — борелевское множество и $ST=S$.

*Доказательство.* (a) $\Psi$ принадлежит классу $C^\infty$ (степенные ряды). Для $Y\in\mathfrak m$ положим $A_X(Y)=e^{-X}\frac{d}{dt}\big|_0e^{X+tY}$.
Дифференцируя $W(t)^*W(t)=1$ для $W(t)=e^{X+tY}\in K$, получаем $A_X(Y)\in\mathfrak k$, и с учётом
$\frac{d}{dt}W^{-1}=-W^{-1}W'W^{-1}$ для $g=e^X$
$$
D\Psi(X,\lambda)(Y,\mu)=g\big([A_X(Y),\Lambda]+\operatorname{diag}(\mu)\big)g^{-1}.
$$
Диагональные матрицы коммутируют с $\Lambda$, поэтому $[A_X(Y),\Lambda]=[B_XY,\Lambda]$, где $B_X=\pi\circ A_X:\mathfrak m\to\mathfrak m$, а
$\pi$ — проекция $\mathfrak k$ на $\mathfrak m$ вдоль $\mathfrak t$. Таким образом, $D\Psi=\operatorname{Ad}_g\circ L_\lambda\circ(B_X\times1)$,
где $L_\lambda(Z,\mu)=[Z,\Lambda]+\operatorname{diag}(\mu)$. Поскольку $([Z,\Lambda])_{ij}=Z_{ij}(\lambda_j-\lambda_i)$, в координатах
$L_\lambda$ умножает $\operatorname{Re}Z_{ij}$ (и $\operatorname{Im}Z_{ij}$) на $\lambda_j-\lambda_i$ и отображает $\mu$ на диагональ:
$|\det L_\lambda|=\Delta(\lambda)^\beta$. $\operatorname{Ad}_g$ сохраняет $dH$. Положим $j(X)=|\det B_X|$: эта функция непрерывна, и
$B_0=1$, так как $\frac d{dt}e^{tY}|_0=Y$; выберем $\varepsilon$ так, чтобы $j>0$ на $\overline{\mathfrak m_\varepsilon}$.
(b) Если $\Psi(X,\lambda)=\Psi(X',\lambda')=H$, то $\lambda=\lambda'=\lambda(H)$, а $e^X,e^{X'}$ — собственные реперы $H\in\Omega$,
поэтому $e^X=e^{X'}D$ с $D\in T$. Пусть $G(X,t)=e^Xe^t$ на $\mathfrak m\times\mathfrak t$ (при $\beta=1$ $\mathfrak t=\{0\}$); $DG(0,0)(Y,s)=Y+s$ инъективно,
и так как $Y$ внедиагональна, а $s$ диагональна, $\|Y+s\|_2^2=\|Y\|_2^2+\|s\|_2^2$, так что $\|DG(0,0)v\|=\|v\|$ для нормы произведения; положим $m=1$. Выберем $\delta>0$ так, чтобы $\|DG(z)-DG(0)\|\le m/2$ при $\|z\|<\delta$; для $u,v$ из этого
шара равенство $G(u)-G(v)=\int_0^1DG(v+s(u-v))(u-v)\,ds$ даёт $\|G(u)-G(v)\|\ge\frac m2\|u-v\|$, так что $G$ там инъективно.
($\delta$ зависит только от $n,\beta$; $\varepsilon$ выбирается после него.)
Далее, $\|D-1\|_2=\|e^{X'}(D-1)\|_2=\|e^X-e^{X'}\|_2\le\|e^X-1\|_2+\|e^{X'}-1\|_2\le2(e^\varepsilon-1)$ в силу степенного ряда и
$\|AB\|_2\le\|A\|_2\|B\|_2$. При $\beta=1$ $D=\operatorname{diag}(\pm1)$, и элемент $-1$ дал бы
$\|D-1\|_2\ge2$; при $\varepsilon<\log2$ имеем $D=1=e^t$ с $t=0$. При $\beta=2$ $D=e^t$ с $t=i\operatorname{diag}(\theta)$,
$\theta_j\in(-\pi,\pi]$, и $|\theta_j|\le\frac\pi2|e^{i\theta_j}-1|$ (так как $|e^{i\theta}-1|=2|\sin\frac\theta2|$ и $\sin x\ge2x/\pi$ на $[0,\pi/2]$), поэтому
$\|t\|_2\le\frac\pi2\|D-1\|_2\le\pi(e^\varepsilon-1)$. При малом $\varepsilon$ ($\varepsilon^2+\pi^2(e^\varepsilon-1)^2<\delta^2$) точки $(X,0)$ и $(X',t)$ лежат в $\delta$-шаре и
$G(X,0)=e^X=e^{X'}D=G(X',t)$, поэтому $X=X'$.
(c) Если $H=\Psi(X,\lambda')$, то $e^X$ — собственный репер, поэтому $F(H)\in e^XT\subseteq S$. Если $F(H)=e^XD$, то
$H=F(H)\Lambda F(H)^*=\Psi(X,\lambda(H))$. $S=G(\mathfrak m_\varepsilon\times\mathfrak t)$ (или $e^{\mathfrak m_\varepsilon}T$, $T$ конечно) —
непрерывный образ $\sigma$-компактного множества, а значит, $\sigma$-компактно и борелевское; $ST=S$ очевидно. $\square$

**Теорема C.5 (формула интегрирования Вейля).** Существует константа $c_{n,\beta}>0$, такая что для любой борелевской
$f:\mathcal H\to[0,\infty]$
$$
\int_{\mathcal H}f(H)\,dH=c_{n,\beta}\int_C\Delta(\lambda)^\beta\Big(\int_Kf(V\Lambda V^*)\,dV\Big)d\lambda .
$$
Равносильно, $\Phi(V,\lambda)=V\Lambda V^*$ переводит меру $c_{n,\beta}\,dV\otimes\Delta(\lambda)^\beta d\lambda$ на $K\times C$
в меру Лебега на $\mathcal H$. Множитель $\Delta(\lambda)^\beta$ — это якобиан карты $\Psi$ при $X=0$
(лемма C.4(a)); константа $c_{n,\beta}$ отражает нормировку меры Хаара, и её явное значение ниже нигде
не нужно, поскольку каждая получаемая из формулы плотность нормируется своей полной массой.

*Доказательство.* *Эталонный ансамбль.* Пусть $H_0$ имеет плотность $p_0(H)=Z_0^{-1}e^{-\operatorname{Tr}H^2/2}$; поскольку
$\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$, $Z_0$ — конечное произведение гауссовых интегралов. $p_0$ и $dH$
инвариантны относительно $H\mapsto VHV^*$, а значит, инвариантно и распределение $H_0$; $P(H_0\in\Omega)=1$ по лемме C.1. По лемме C.3
$\tilde U$ распределён по Хаару, не зависит от $\lambda=\lambda(H_0)$, и $H_0=\tilde U\Lambda\tilde U^*$.
*Распределение $\lambda$.* Пусть $A\subseteq C$ — борелевское множество. По лемме C.4(c) и в силу $ST=S$
$\{H_0\in\Psi(\mathfrak m_\varepsilon\times A)\}=\{\lambda\in A,\ \tilde U\in S\}$ с точностью до события нулевой вероятности, и вероятность этого события равна
$P(\lambda\in A)\,h$, где $h=\int_K1_S\,dV$. С другой стороны, $\Psi$ инъективно и принадлежит классу $C^1$ с ненулевым якобианом на
открытом множестве $\mathfrak m_\varepsilon\times C$, и $p_0(\Psi(X,\lambda))=p_0(\Lambda)$; формула замены переменных и лемма C.4(a) дают
$$
P\big(H_0\in\Psi(\mathfrak m_\varepsilon\times A)\big)=\int_{\mathfrak m_\varepsilon}\!\int_Ap_0(\Lambda)\Delta(\lambda)^\beta j(X)\,d\lambda\,dX
=\kappa\int_Ap_0(\Lambda)\Delta(\lambda)^\beta\,d\lambda,\qquad \kappa=\int_{\mathfrak m_\varepsilon}j>0 .
$$
При $A=C$ получаем $h=\kappa\int_Cp_0\Delta^\beta>0$. Следовательно, $\lambda$ имеет на $C$ плотность $c\,p_0(\Lambda)\Delta(\lambda)^\beta$, где $c=\kappa/h$.
*Общая $f$.* Так как $p_0>0$, $\int f\,dH=\mathbb E\big[(f/p_0)(\tilde U\Lambda\tilde U^*)\big]$. Совместное распределение $(\tilde U,\lambda)$ равно
$dV\otimes c\,p_0(\Lambda)\Delta^\beta d\lambda$, и $p_0(V\Lambda V^*)=p_0(\Lambda)$; теорема Тонелли даёт формулу с $c_{n,\beta}=c$.
Константа определяется самой формулой (возьмите $f=p_0$). $\square$

**Следствие C.6 (инвариантные плотности).** Пусть $H$ имеет плотность $p$ на $\mathcal H$, причём $p(VHV^*)=p(H)$ для всех
$V\in K$. Тогда $\lambda(H)$ имеет на $C$ плотность $c_{n,\beta}\,p(\Lambda)\Delta(\lambda)^\beta$; собственный репер $\tilde U$ из леммы C.3
распределён по Хаару и не зависит от $\lambda(H)$; и для симметричных борелевских $g\ge0$ на $\mathbb R^n$
$$
\mathbb E\,g(\lambda(H))=\frac{c_{n,\beta}}{n!}\int_{\mathbb R^n}g(x)\,p(\operatorname{diag}x)\,|\Delta(x)|^\beta\,dx .
$$

*Доказательство.* $P(H\in\Omega)=1$ по лемме C.1. Для борелевского $A\subseteq C$ применим теорему C.5 к $f=1_A(\lambda(\cdot))\,p$, для которой
$f(V\Lambda V^*)=1_A(\lambda)p(\Lambda)$. Распределение $H$ инвариантно, поэтому применима лемма C.3. С точностью до множества меры нуль, где две
координаты совпадают, $\mathbb R^n$ — дизъюнктное объединение $n!$ образов $C$ при перестановках координат, а
$g$, $|\Delta|$ и $x\mapsto p(\operatorname{diag}x)$ инвариантны относительно перестановок (матрица перестановки лежит в $K$). $\square$

Проверка ниже численно дифференцирует $\Psi$ для $n=2,3$ и обоих $\beta$: при $X=0$ якобиан равен $\Delta(\lambda)^\beta$, а
при фиксированном $X\neq0$ его отношение к $\Delta(\lambda)^\beta$ не зависит от $\lambda$ (лемма C.4(a)).

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

## C.3 Гауссовы ансамбли

**Теорема C.7 (GOE и GUE).** Пусть $H=(G+G^T)/2$, где $G$ — вещественная матрица с независимыми элементами $N(0,1)$ ($\beta=1$, GOE), или
$H=(G+G^*)/\sqrt2$, где $G$ — матрица Жинибра ($\beta=2$, GUE), как в разделе 5.3. Тогда $H$ имеет на $\mathcal H$ плотность
$(2\pi)^{-n/2}\pi^{-\beta n(n-1)/4}e^{-\operatorname{Tr}H^2/2}$, её собственные значения имеют плотность
$$
Z_{n,\beta}^{-1}\,\Delta(\lambda)^\beta\,e^{-\sum_i\lambda_i^2/2}\quad\text{on }C,\qquad Z_{n,\beta}=\int_C\Delta^\beta e^{-|\lambda|^2/2}\,d\lambda\in(0,\infty),
$$
а собственный репер $\tilde U$ распределён по Хаару на $K$ и не зависит от собственных значений.

*Доказательство.* GOE: $h_{ii}=G_{ii}\sim N(0,1)$ и $h_{ij}=(G_{ij}+G_{ji})/2\sim N(0,\frac12)$ при $i<j$. GUE: $h_{ii}=\sqrt2\operatorname{Re}G_{ii}\sim N(0,1)$,
а $\operatorname{Re}h_{ij}=(\operatorname{Re}G_{ij}+\operatorname{Re}G_{ji})/\sqrt2$ и $\operatorname{Im}h_{ij}=(\operatorname{Im}G_{ij}-\operatorname{Im}G_{ji})/\sqrt2$ распределены по $N(0,\frac12)$.
Различные координаты являются функциями непересекающихся наборов независимых величин, поэтому они независимы, и
плотность равна $\prod_i(2\pi)^{-1/2}e^{-h_{ii}^2/2}\prod_{i<j}\pi^{-\beta/2}e^{-|h_{ij}|^2}$, что совпадает с указанной, так как
$\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$. Она инвариантна, и применимо следствие C.6; константа определяется
нормировкой. $Z_{n,\beta}<\infty$, поскольку $|\lambda_j-\lambda_i|\le(1+|\lambda_i|)(1+|\lambda_j|)$, так что
$|\Delta(\lambda)|\le\prod_i(1+|\lambda_i|)^{n-1}$. $\square$

**Следствие C.8 (плотность уровней).** Для GOE или GUE размера $n$ положим $q(x)=(n!Z_{n,\beta})^{-1}|\Delta(x)|^\beta e^{-|x|^2/2}$ на
$\mathbb R^n$ и $\rho(x)=n\int_{\mathbb R^{n-1}}q(x,y)\,dy$. Тогда $\mathbb E\sum_i\varphi(\lambda_i(H))=\int\varphi\rho$ для ограниченных борелевских
$\varphi$, $\int\rho=n$, и $\rho$ непрерывна, причём $0\le\rho(x)\le C_n(1+|x|)^{\beta(n-1)}e^{-x^2/2}$. В частности,
предположение теоремы 5.15 выполнено для независимых блоков GOE или GUE.

*Доказательство.* По следствию C.6 с $g(x)=\sum_i\varphi(x_i)$ (для $\varphi\ge0$; в общем случае применяем его к $\varphi^\pm=\max(\pm\varphi,0)$) и в силу симметрии $q$
$\mathbb E\sum_i\varphi(\lambda_i)=n\int\varphi(x_1)q(x)\,dx=\int\varphi\rho$ (Фубини). По оценке для $|\Delta|$
$q(x,y)\le(n!Z_{n,\beta})^{-1}(1+|x|)^{\beta(n-1)}e^{-x^2/2}\prod_{i}(1+|y_i|)^{\beta(n-1)}e^{-y_i^2/2}$, что после
интегрирования по $y$ даёт требуемую оценку, а поскольку $\sup_x(1+|x|)^{\beta(n-1)}e^{-x^2/2}<\infty$, — и интегрируемую мажоранту, не зависящую от $x$: непрерывность
следует из теоремы о мажорируемой сходимости. $\square$

**Предложение C.9 (наименьший интервал).** Для GOE или GUE размера $n\ge2$ пусть $\delta=\min_i(\lambda_{i+1}-\lambda_i)$. Существуют
константы $0<c\le C$ (зависящие от $n,\beta$), такие что $c\,\varepsilon^{\beta+1}\le P(\delta<\varepsilon)\le C\varepsilon^{\beta+1}$ при $0<\varepsilon\le1$.

*Доказательство.* *Верхняя оценка.* $P(\delta<\varepsilon)\le\sum_iP(\lambda_{i+1}-\lambda_i<\varepsilon)$. На $\{0<u=x_{i+1}-x_i<\varepsilon\}$ оценим множитель
пары $(i,i+1)$ через $u^\beta$, остальные — как выше, и воспользуемся $1+|x_{i+1}|\le2(1+|x_i|)$ и $e^{-x_{i+1}^2/2}\le1$; подставляя
$x_{i+1}=x_i+u$, получаем
$$
P(\lambda_{i+1}-\lambda_i<\varepsilon)\le\frac{2^{\beta(n-1)}}{Z_{n,\beta}}\int_0^\varepsilon u^\beta du\int(1+|x|)^{2\beta(n-1)}e^{-\frac{x^2}2}dx
\Big(\int(1+|x|)^{\beta(n-1)}e^{-\frac{x^2}2}dx\Big)^{n-2}.
$$
*Нижняя оценка.* Пусть $R_\varepsilon=\{x:0<x_1<1,\ 0<x_2-x_1<\varepsilon,\ 3k<x_k<3k+1\ (3\le k\le n)\}\subseteq C$. На $R_\varepsilon$
$\delta\le x_2-x_1$, все пары, кроме $(1,2)$, находятся на расстоянии $\ge1$, и $|x|^2\le M$ с $M$, не зависящим от $\varepsilon\le1$; поэтому
$P(\delta<\varepsilon)\ge P(\lambda\in R_\varepsilon)\ge Z_{n,\beta}^{-1}e^{-M/2}\int_0^\varepsilon u^\beta du$. $\square$

Для GOE это и есть порядок $\varepsilon^2$, упомянутый в конце раздела 5.4. Проверим теорему C.7 при $n=3$: плотность уровней GOE
вычисляется из совместной плотности интегрированием по двум другим собственным значениям на сетке и сравнивается
со средним числом собственных значений в каждом интервале гистограммы; проверим также, что собственные векторы не коррелируют с
собственными значениями ($|\tilde U_{11}|^2$ не зависит от фаз $D$, поэтому его можно взять из любого собственного репера, возвращаемого
`eigh`). Для распределённой по Хаару $U\in O(3)$ $\mathbb E|U_{11}|^2=\frac13$: умножение слева на матрицу перестановки переставляет
элементы первого столбца и сохраняет меру Хаара, поэтому три величины $\mathbb E|U_{i1}|^2$ равны, а их сумма равна $1$.

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

## C.4 Ортогональные многочлены и одноточечная функция

**Лемма C.10 (одноточечная функция).** Пусть $I\subseteq\mathbb R$ — интервал, а $w\ge0$ — борелевская функция на $I$, положительная на
множестве положительной меры, причём $\int_I|x|^mw<\infty$ для всех $m$. Пусть $p_0,\dots,p_{n-1}$ — вещественные многочлены, $\deg p_k=k$,
со старшим коэффициентом $\kappa_k\neq0$, ортонормированные в $L^2(w)$ (они существуют: применим процесс Грама–Шмидта к $1,x,\dots,x^{n-1}$, которые
линейно независимы в $L^2(w)$, поскольку ненулевой многочлен имеет конечное число нулей; знак каждого $p_k$ можно выбрать любым). Пусть $q=Z_n^{-1}\Delta(x)^2\prod_iw(x_i)$ на
$I^n$, где $Z_n=\int_{I^n}\Delta^2\prod_iw(x_i)\,dx$. Тогда $Z_n=n!\prod_k\kappa_k^{-2}$ и
$$
n\int_{I^{n-1}}q(x,y)\,dy=K_n(x)\,w(x),\qquad K_n(x)=\sum_{k=0}^{n-1}p_k(x)^2=v(x)^TM^{-1}v(x),
$$
где $v(x)=(1,x,\dots,x^{n-1})$ и $M_{ij}=\int_Ix^{i+j}w$, $0\le i,j<n$.

*Доказательство.* $\Delta(x)=\det(x_j^{i-1})_{i,j}$. Прибавляя к каждой строке линейные комбинации строк над ней, заменяем $x_j^{i-1}$ на
$p_{i-1}(x_j)/\kappa_{i-1}$, так что $\Delta(x)=(\prod_k\kappa_k)^{-1}\det(p_{i-1}(x_j))$ и
$$
\Delta(x)^2=\Big(\prod_k\kappa_k\Big)^{-2}\sum_{\sigma,\tau\in S_n}\operatorname{sgn}\sigma\operatorname{sgn}\tau\prod_{j=1}^np_{\sigma(j)-1}(x_j)\,p_{\tau(j)-1}(x_j).
$$
При интегрировании по $x_2,\dots,x_n$ с весом $\prod_{j\ge2}w(x_j)$ множитель с номером $j$ даёт $\delta_{\sigma(j)\tau(j)}$, поэтому остаются только слагаемые с $\sigma=\tau$,
и результат равен $(\prod_k\kappa_k)^{-2}\sum_\sigma p_{\sigma(1)-1}(x_1)^2=(\prod_k\kappa_k)^{-2}(n-1)!\,K_n(x_1)$.
Интегрирование с весом $w(x_1)$ даёт $Z_n=(\prod_k\kappa_k)^{-2}n!$, а деление даёт формулу. Наконец, если
$p_k(x)=\sum_iC_{ki}x^i$, то ортонормированность записывается как $CMC^T=1$ с обратимой $C$, поэтому $M^{-1}=C^TC$ и
$v^TM^{-1}v=|Cv|^2=\sum_kp_k^2$. $\square$

По следствию C.6 неупорядоченные собственные значения GUE имеют плотность $q$ из леммы C.10 с $w(x)=e^{-x^2/2}$ на
$\mathbb R$, поэтому плотность уровней GUE равна $K_n(x)e^{-x^2/2}$. Проверим это при $n=3$ с помощью моментов
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

## C.5 Матрицы Уишарта и индуцированная мера

Пусть $a\le b$, $G$ — матрица Жинибра размера $a\times b$ (независимые элементы с независимыми вещественной и мнимой частями $N(0,\frac12)$,
плотность $\pi^{-ab}e^{-\operatorname{Tr}GG^*}$) и $W=GG^*$. Пусть $T_+$ — множество нижнетреугольных матриц размера $a\times a$ с положительной
диагональю — открытое подмножество $\mathbb R^{a^2}$ в координатах $L_{ii}$, $\operatorname{Re}L_{ij}$, $\operatorname{Im}L_{ij}$ ($i>j$), а
$P_a$ — множество положительно определённых матриц, открытое в пространстве эрмитовых матриц размера $a\times a$.

**Лемма C.11 (разложение Бартлетта).** Почти наверное строки $r_1,\dots,r_a$ матрицы $G$ линейно независимы, и
процесс Грама–Шмидта даёт $G=LV$ с $L\in T_+$ и $V$ с ортонормированными строками, так что $W=LL^*$. Элементы $L_{ij}$, $i\ge j$,
независимы, $L_{ij}$ ($i>j$) — стандартные комплексные гауссовы величины, а $L_{ii}^2$ имеет плотность $x^{b-i}e^{-x}/(b-i)!$ на $(0,\infty)$.
Следовательно, $L$ имеет на $T_+$ плотность $c_L\prod_iL_{ii}^{2(b-i)+1}e^{-\operatorname{Tr}LL^*}$.

*Доказательство.* Пусть $\mathcal F_{i-1}=\sigma(r_1,\dots,r_{i-1})$, и предположим (индукция), что $r_1,\dots,r_{i-1}$ почти наверное линейно независимы, с
векторами Грама–Шмидта $v_1,\dots,v_{i-1}$. Пусть $u_1,\dots,u_b$ — ортонормированный базис, полученный процессом Грама–Шмидта из
$v_1,\dots,v_{i-1},e_1,\dots,e_b$ (с отбрасыванием линейно зависимых векторов); $u_k=v_k$ при $k<i$, и $u_k$ — борелевские функции от
$r_1,\dots,r_{i-1}$. Положим $c_k=\langle u_k,r_i\rangle$. Поскольку $r_i$ не зависит от $\mathcal F_{i-1}$, а его плотность $\pi^{-b}e^{-\|r\|^2}$
инвариантна относительно унитарных преобразований (предложение A.16(a)), условно относительно $\mathcal F_{i-1}$ вектор $c$ — стандартный комплексный гауссов (Фубини). Далее,
$L_{ij}=c_j$ при $j<i$ и $r_i-\sum_{j<i}c_jv_j=\sum_{k\ge i}c_ku_k$, поэтому почти наверное $L_{ii}^2=\sum_{k=i}^b|c_k|^2>0$. Таким образом, условное распределение
$i$-й строки $L$ относительно $\mathcal F_{i-1}$ — фиксированное распределение, так что $i$-я строка не зависит от предыдущих строк, а её координаты
$c_1,\dots,c_{i-1}$ и $L_{ii}$ независимы; по индукции все строки $L$ независимы. По теореме A.20(b) с $N=m=b-i+1$
$L_{ii}^2=\|(c_i,\dots,c_b)\|^2$ имеет плотность $s^{m-1}e^{-s}/(m-1)!$. Тогда плотность $L_{ii}=\sqrt s$ равна
$2t^{2(b-i)+1}e^{-t^2}/(b-i)!$, и $\operatorname{Tr}LL^*=\sum_{i\ge j}|L_{ij}|^2$. $\square$

**Лемма C.12 (отображение Холецкого).** $\chi(L)=LL^*$ — биекция класса $C^\infty$ множества $T_+$ на $P_a$, и
$|\det D\chi(L)|=2^a\prod_{i=1}^aL_{ii}^{2(a-i)+1}$.

*Доказательство.* *Сюръективность:* для $W\in P_a$ применим лемму 5.2 к $W^{1/2}=QR$; тогда $W=R^*R$ и $L=R^*\in T_+$. *Инъективность:* если
$LL^*=MM^*$, то $M^{-1}L=M^*(L^*)^{-1}$ одновременно нижне- и верхнетреугольна, а значит, диагональна с положительной диагональю, и
унитарна, а значит, равна $1$. *Якобиан:* при $i\ge j$ $W_{ij}=L_{ij}L_{jj}+\sum_{k<j}L_{ik}\overline{L_{jk}}$. Упорядочим координаты $L$ и
нижнего треугольника $W$ по столбцам, а внутри столбца $j$ поставим диагональный элемент первым. Тогда $W_{ij}$ зависит только от
$L_{ij}$ и от предшествующих координат, поэтому $D\chi$ треугольна, с диагональными элементами $2L_{jj}$ ($\partial W_{jj}/\partial L_{jj}$) и
для каждого $i>j$ дважды $L_{jj}$ (вещественная и мнимая части $W_{ij}$). Произведение равно $\prod_j2L_{jj}\cdot L_{jj}^{2(a-j)}$. $\square$

**Теорема C.13 (ансамбли Уишарта и Лагерра).** $W=GG^*$ имеет плотность $c_W\det(W)^{b-a}e^{-\operatorname{Tr}W}$ на $P_a$
(нулевую вне его), а её собственные значения имеют плотность
$$
Z_{a,b}^{-1}\,\Delta(\lambda)^2\prod_{i=1}^a\lambda_i^{\,b-a}e^{-\lambda_i}\quad\text{on }C_+=\{0<\lambda_1<\dots<\lambda_a\}.
$$

*Доказательство.* Для борелевской $g\ge0$ лемма C.11 и замена переменных $W=\chi(L)$ (лемма C.12) дают
$\mathbb E\,g(W)=\int_{P_a}g(W)\,c_L\prod_iL_{ii}^{2(b-i)+1}e^{-\operatorname{Tr}W}\big/\big(2^a\prod_iL_{ii}^{2(a-i)+1}\big)dW$ с
$L=\chi^{-1}(W)$, причём $\prod_iL_{ii}^{2(b-a)}=\det(W)^{b-a}$. Эта плотность инвариантна относительно $W\mapsto VWV^*$, $V\in U(a)$, и
$\Lambda\in P_a$ тогда и только тогда, когда $\lambda_1>0$; применим следствие C.6 с $\beta=2$. $\square$

**Теорема C.14 (индуцированная мера).** Пусть $\psi$ — случайное по Хаару чистое состояние $\mathbb C^a\otimes\mathbb C^b$ (раздел 7.5),
$\rho_A=\operatorname{Tr}_2|\psi\rangle\langle\psi|=W/s$, где $s=\operatorname{Tr}W$ (доказательство теоремы 7.11), и
$\mu_1<\dots<\mu_a$ — собственные значения $\rho_A$. Тогда $s$ и $\mu$ независимы, $s$ имеет плотность
$s^{ab-1}e^{-s}/(ab-1)!$, а $(\mu_1,\dots,\mu_{a-1})$ имеет плотность
$$
Z'^{-1}_{a,b}\,\Delta(\mu)^2\prod_{i=1}^a\mu_i^{\,b-a}\quad\text{on }D_a=\{0<\mu_1<\dots<\mu_{a-1}<\mu_a\},\qquad\mu_a=1-\sum_{i<a}\mu_i ,
$$
то есть плотность «$\propto\delta(1-\sum_i\mu_i)\,\Delta(\mu)^2\prod_i\mu_i^{b-a}$» физической литературы.

*Доказательство.* $\Theta(s,\mu_1,\dots,\mu_{a-1})=s\,(\mu_1,\dots,\mu_{a-1},1-\sum_{i<a}\mu_i)$ — биекция класса $C^\infty$ множества $(0,\infty)\times D_a$ на
$C_+$ с обратной $\lambda\mapsto(\sum\lambda_i,\lambda_1/\sum\lambda_i,\dots)$. Прибавление всех строк $D\Theta$ к последней превращает её в
$(1,0,\dots,0)$, поэтому $|\det D\Theta|=s^{a-1}$. Поскольку $\Delta$ однородна степени $a(a-1)/2$, плотность $(s,\mu)$ равна
$$
Z_{a,b}^{-1}\,s^{a(a-1)+a(b-a)+a-1}e^{-s}\cdot\Delta(\mu)^2\prod_i\mu_i^{b-a}=Z_{a,b}^{-1}\,s^{ab-1}e^{-s}\cdot\Delta(\mu)^2\prod_i\mu_i^{b-a},
$$
то есть является произведением, что доказывает независимость и после нормировки определяет оба множителя. $\square$

Проверим лемму C.11 на множителе Холецкого, вычисленном `torch.linalg.cholesky` (это единственная $L\in T_+$ по
лемме C.12), и теорему C.14 при $a=2$, когда $\mu_1\in(0,\frac12)$ имеет плотность $\propto(1-2t)^2(t(1-t))^{b-2}$:

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

## C.6 Многочлены Лагерра и формула Пейджа

**Лемма C.15 (многочлены Лагерра).** Для целых $\alpha,k\ge0$ положим
$L_k^{(\alpha)}(x)=\sum_{j=0}^k(-1)^j\binom{k+\alpha}{k-j}\frac{x^j}{j!}$ и $L_{-1}^{(\alpha)}=0$. Тогда:

(a) $L_k^{(\alpha)}(x)\,x^\alpha e^{-x}=\frac1{k!}\frac{d^k}{dx^k}\big(x^{k+\alpha}e^{-x}\big)$; $L_k^{(\alpha)}$ имеет степень $k$ и старший коэффициент $(-1)^k/k!$.

(b) $\frac{d}{dx}L_k^{(\alpha)}=-L_{k-1}^{(\alpha+1)}$, следовательно, $\frac{d^i}{dx^i}L_k^{(\alpha)}=(-1)^iL_{k-i}^{(\alpha+i)}$ при $i\le k$; и $L_k^{(\alpha)}=L_k^{(\alpha+1)}-L_{k-1}^{(\alpha+1)}$.

(c) $\int_0^\infty x^mL_k^{(\alpha)}x^\alpha e^{-x}dx=0$ при $m<k$, и $\int_0^\infty L_k^{(\alpha)}L_l^{(\alpha)}x^\alpha e^{-x}dx=\delta_{kl}\frac{(k+\alpha)!}{k!}$.

Следовательно, по лемме C.10 и теореме C.13 собственные значения $W$ имеют плотность уровней
($\mathbb E\sum_i\varphi(\lambda_i(W))=\int_0^\infty\varphi\,\rho_{a,b}$ для борелевских $\varphi\ge0$ и для всех $\varphi$ с $\int|\varphi|\rho_{a,b}<\infty$)
$$
\rho_{a,b}(x)=\sum_{k=0}^{a-1}\frac{k!}{(k+\alpha)!}L_k^{(\alpha)}(x)^2\,x^\alpha e^{-x},\qquad\alpha=b-a .
$$

*Доказательство.* (a) По формуле Лейбница $\frac{d^k}{dx^k}(x^{k+\alpha}e^{-x})=\sum_j\binom kj(-1)^j\frac{(k+\alpha)!}{(\alpha+j)!}x^{\alpha+j}e^{-x}$ (где $j$ производных
приходятся на $e^{-x}$), и $\frac1{k!}\binom kj\frac{(k+\alpha)!}{(\alpha+j)!}=\binom{k+\alpha}{k-j}\frac1{j!}$. (b) Почленное дифференцирование даёт
$-\sum_{i=0}^{k-1}(-1)^i\binom{k+\alpha}{k-1-i}\frac{x^i}{i!}=-L_{k-1}^{(\alpha+1)}$; второе тождество — правило Паскаля
$\binom{k+\alpha+1}{k-j}-\binom{k+\alpha}{k-1-j}=\binom{k+\alpha}{k-j}$. (c) Проинтегрируем (a) по частям $k$ раз против $x^m$:
граничные члены обращаются в нуль, так как $\frac{d^r}{dx^r}(x^{k+\alpha}e^{-x})$ содержит множитель $x^{k+\alpha-r}$, $r<k$, и убывает на $\infty$.
Получаем $\frac{(-1)^k}{k!}\int\frac{d^kx^m}{dx^k}x^{k+\alpha}e^{-x}dx$, что равно $0$ при $m<k$ и $(-1)^k(k+\alpha)!$ при $m=k$. При $l<k$
$L_l^{(\alpha)}$ — линейная комбинация $x^m$, $m<k$; при $l=k$ важен только старший коэффициент $(-1)^k/k!$. Поэтому
$p_k=\big(k!/(k+\alpha)!\big)^{1/2}L_k^{(\alpha)}$ ортонормированы для $w(x)=x^\alpha e^{-x}$ на $(0,\infty)$, и неупорядоченные
собственные значения $W$ имеют плотность $q$ из леммы C.10 (с $n=a$) по следствию C.6, так что, как в доказательстве следствия C.8,
$\mathbb E\sum_i\varphi(\lambda_i)=a\int\varphi(x_1)q(x)\,dx=\int\varphi\,K_aw=\int\varphi\,\rho_{a,b}$ для борелевских $\varphi\ge0$, а по линейности — и для
$\varphi$ с $\int|\varphi|\rho_{a,b}<\infty$. $\square$

**Лемма C.16 (логарифмические интегралы).** Пусть $I_m=\int_0^\infty x^m\log x\,e^{-x}dx$ и $H_m=\sum_{i=1}^m\frac1i$. Тогда
$I_m=m!\,(I_0+H_m)$. Для целых $\alpha\ge0$ и $0\le k\le m$
$$
A^{(\alpha)}_{k,m}:=\int_0^\infty\log x\,L_k^{(\alpha)}(x)L_m^{(\alpha)}(x)\,x^\alpha e^{-x}dx=
\begin{cases}I_{k+\alpha}/k!, & m=k,\\ -\dfrac{(k+\alpha)!}{k!\,(m-k)}, & m>k.\end{cases}
$$
($I_0$ сокращается в теореме C.17, поэтому его значение не нужно.)

*Доказательство.* При $m\ge1$ интегрирование по частям с $(-e^{-x})'=e^{-x}$ даёт $I_m=mI_{m-1}+(m-1)!$, т. е.
$I_m/m!=I_{m-1}/(m-1)!+1/m$. Для $A_{k,m}$ запишем $L_m^{(\alpha)}x^\alpha e^{-x}=\frac1{m!}g^{(m)}$ с $g=x^{m+\alpha}e^{-x}$
(лемма C.15(a)) и $f=\log x\,L_k^{(\alpha)}$ и проинтегрируем по частям $m$ раз. Граничные члены
$f^{(i)}g^{(m-1-i)}$ обращаются в нуль на $\infty$ и в $0$, где $f^{(i)}=O(|\log x|+x^{-i})$ и $g^{(m-1-i)}=O(x^{\alpha+1+i})$. Поэтому
$A_{k,m}=\frac{(-1)^m}{m!}\int_0^\infty f^{(m)}g\,dx$ (подынтегральное выражение есть $O(x^\alpha(1+|\log x|))$ в $0$), где
$f^{(m)}=\sum_{i\le k}\binom mi(\log x)^{(m-i)}\frac{d^i}{dx^i}L_k^{(\alpha)}$. При $m-i\ge1$ $(\log x)^{(m-i)}=(-1)^{m-i-1}(m-i-1)!\,x^{i-m}$,
и по лемме C.15(b,c)
$$
\int_0^\infty x^{i-m}\frac{d^iL_k^{(\alpha)}}{dx^i}\,x^{m+\alpha}e^{-x}dx=(-1)^i\int_0^\infty L_{k-i}^{(\alpha+i)}\,x^{\alpha+i}e^{-x}dx=(-1)^i(\alpha+i)!\,\delta_{ik}.
$$
Если $m>k$, вклад даёт только $i=k$: $A_{k,m}=\frac{(-1)^m}{m!}\binom mk(-1)^{m-k-1}(m-k-1)!(-1)^k(\alpha+k)!=-\frac{(k+\alpha)!}{k!(m-k)}$.
Если $m=k$, слагаемые с $i<k$ равны нулю, а слагаемое с $i=k$ равно $\log x\cdot(-1)^k$, что даёт $A_{k,k}=\frac{(-1)^k}{k!}(-1)^kI_{k+\alpha}$. $\square$

**Теорема C.17 (формула Пейджа).** Для случайного по Хаару чистого состояния $\mathbb C^a\otimes\mathbb C^b$ с $a\le b$
$$
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
$$

*Доказательство.* *Сведение.* В обозначениях теоремы C.14 $\lambda_i=s\mu_i$ и $\sum_i\mu_i=1$, поэтому
$\sum_i\lambda_i\log\lambda_i=s\log s-s\,S(\rho_A)$. Все слагаемые интегрируемы ($|x\log x|\le1+x^2$, и $\mathbb Es^2<\infty$,
$0\le S(\rho_A)\le a/e$, поскольку $0\le-x\log x\le1/e$ на $[0,1]$), $\mathbb E\sum_i\lambda_i\log\lambda_i=\int x\log x\,\rho_{a,b}(x)\,dx$ по лемме C.15, и $s$ не зависит от
$S(\rho_A)$, поэтому $\mathbb E\sum_i\lambda_i\log\lambda_i=\mathbb E[s\log s]-\mathbb E[s]\,\mathbb E\,S(\rho_A)$. Из плотности $s$
$\mathbb Es=ab$ и $\mathbb E[s\log s]=I_{ab}/(ab-1)!=ab\,(I_0+H_{ab})$ (лемма C.16). Следовательно,
$$
\mathbb E\,S(\rho_A)=I_0+H_{ab}-\frac1{ab}\int_0^\infty x\log x\,\rho_{a,b}(x)\,dx .
$$
*Интеграл.* По лемме C.15(b) $L_k^{(\alpha)}=L_k^{(\alpha+1)}-L_{k-1}^{(\alpha+1)}$, поэтому с учётом $x\cdot x^\alpha=x^{\alpha+1}$ и леммы C.16
для $\alpha+1$
$$
\frac{k!}{(k+\alpha)!}\int_0^\infty x\log x\,L_k^{(\alpha)}(x)^2x^\alpha e^{-x}dx=\frac{k!}{(k+\alpha)!}\Big(A^{(\alpha+1)}_{k,k}-2A^{(\alpha+1)}_{k-1,k}+A^{(\alpha+1)}_{k-1,k-1}\Big)
$$
$$
=(k+\alpha+1)(I_0+H_{k+\alpha+1})+2k+k\,(I_0+H_{k+\alpha}),
$$
причём и при $k=0$. Поскольку $\sum_{k=0}^{a-1}(2k+\alpha+1)=a(a-1)+a(b-a+1)=ab$,
$$
\int_0^\infty x\log x\,\rho_{a,b}\,dx=ab\,I_0+\Sigma,\qquad \Sigma=\sum_{k=0}^{a-1}\big[(k+\alpha+1)H_{k+\alpha+1}+kH_{k+\alpha}+2k\big].
$$
*Сумма.* Положим $m=\alpha+1=b-a+1$ и $j=k+m\in\{m,\dots,b\}$. Поскольку $(j-m)H_{j-1}=(j-m)H_j-1+\frac mj$,
$\Sigma=\sum_{j=m}^b(2j-m)H_j-a+m(H_b-H_{m-1})+a(a-1)$. Для $F(j)=j(j+1-m)$ имеем $F(j)-F(j-1)=2j-m$, $F(m-1)=0$,
$F(b)=ab$, и суммирование по частям даёт
$$
\sum_{j=m}^b(2j-m)H_j=abH_b-\sum_{j=m}^{b-1}\frac{F(j)}{j+1}=abH_b-\sum_{j=m}^{b-1}\Big(j-m+\frac m{j+1}\Big)=abH_b-\frac{(a-1)(a-2)}2-m(H_b-H_m).
$$
Так как $m(H_m-H_{m-1})=1$, $\Sigma=abH_b-\frac{(a-1)(a-2)}2+1-a+a(a-1)=abH_b+\frac{a(a-1)}2$. Следовательно,
$\mathbb E\,S(\rho_A)=H_{ab}-H_b-\frac{a-1}{2b}$. $\square$

**Следствие C.18.** При $a\le b$
$\log a-\frac{a-1}{ab+1}-\frac{a-1}{2b}\le\mathbb E\,S(\rho_A)\le\log a-\frac{a-1}{2b}$; в частности,
$\big|\mathbb E\,S(\rho_A)-\log a+\frac a{2b}\big|\le\frac1{2b}$. Для двух кубитов $\mathbb E\,S(\rho_A)=\frac13$.

*Доказательство.* Поскольку $1/x$ убывает, $\int_k^{k+1}\frac{dx}x\le\frac1k\le\int_{k-1}^k\frac{dx}x$; суммируя по $b<k\le ab$, получаем
$\log\frac{ab+1}{b+1}\le H_{ab}-H_b\le\log a$, и $\log\frac{ab+1}{b+1}=\log a-\log\big(1+\frac{a-1}{ab+1}\big)\ge\log a-\frac{a-1}{ab+1}$.
Наконец, $\frac{a-1}{ab+1}\le\frac1b$, а $a=b=2$ даёт $\frac13+\frac14-\frac14$. $\square$

Формула была высказана как гипотеза D. N. Page (Phys. Rev. Lett. 71, 1291, 1993) и доказана S. K. Foong и S. Kanno
(Phys. Rev. Lett. 72, 1148, 1994), J. Sánchez-Ruiz (Phys. Rev. E 52, 5653, 1995) и S. Sen (Phys. Rev. Lett. 77, 1, 1996);
приведённое выше вычисление следует сведению S. Sen к ансамблю Лагерра, но суммирует одноточечную плотность почленно
вместо использования формулы Кристоффеля–Дарбу. Индуцированная мера принадлежит K. Życzkowski и H.-J. Sommers
(J. Phys. A 34, 7111, 2001); разложение леммы C.11 принадлежит M. S. Bartlett (1933), а теорема C.5 — H. Weyl.

Мы проверяем плотность уровней Лагерра по гистограмме, значения $A^{(\alpha)}_{k,m}$ из леммы C.16 и промежуточный
результат $\int x\log x\,\rho_{a,b}=ab\,I_0+\Sigma$ численным интегрированием (причём сам $I_0$ тоже вычисляется численным интегрированием), тождество для
$\Sigma$ в точной арифметике, формулу Пейджа — сравнением с прямым интегрированием
индуцированной плотности при $a=2$ (точные рациональные значения) и с методом Монте-Карло, — и следствие C.18:

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

**Назад к:** [главе 5](05_random_matrices.md) (теорема C.7, следствие C.8 и предложение C.9),
[главе 7](07_tensor_products.md) (теоремы C.14 и C.17 и следствие C.18).
