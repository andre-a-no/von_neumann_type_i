[English](../12_equivariant_learning.md) · [Español](../es/12_equivariant_learning.md) · [Français](../fr/12_equivariant_learning.md) · [Deutsch](../de/12_equivariant_learning.md) · [中文](../zh/12_equivariant_learning.md) · [日本語](../ja/12_equivariant_learning.md) · **Русский**

# Руководство 12. Эквивариантное машинное обучение

Это руководство сопровождает [главу 12](../../lecture_notes/ru/12_equivariant_learning.md) конспекта лекций, где доказывается, что эквивариантные линейные
отображения между двумя представлениями конечной группы являются блочными матрицами ([теорема 12.4](../../lecture_notes/ru/12_equivariant_learning.md#r-12-4)), что
усреднение по группе есть ортогональный проектор на них ([теорема 12.5](../../lecture_notes/ru/12_equivariant_learning.md#r-12-5)), и подсчитываются инварианты.
Здесь мы строим и обучаем эквивариантные слои с помощью `torch_vn_algebra.decompose` для $S_4$, действующей на признаках вершин и
рёбер графов, и численно проверяем эти результаты.

В коде пакет векторов — это тензор формы `(batch, n)`, на который матрица $T$ действует как `x @ T.T`.

<a id="s-12-1"></a>

## 12.1 Эквивариантные линейные отображения

*Конспект:* [раздел 12.1](../../lecture_notes/ru/12_equivariant_learning.md#s-12-1): [лемма 12.1](../../lecture_notes/ru/12_equivariant_learning.md#r-12-1) (орбиты) и [следствие 12.2](../../lecture_notes/ru/12_equivariant_learning.md#r-12-2) (DeepSets, числа Белла).

Слой DeepSets (Zaheer et al. 2017) — это [конспект лекций, следствие 12.2](../../lecture_notes/ru/12_equivariant_learning.md#r-12-2)(a), базисы Maron et al. (2019) для графовых сетей — это (b).
Код подсчитывает орбиты перебором, сравнивает результат с численно вычисленным коммутантом ([лемма 12.1](../../lecture_notes/ru/12_equivariant_learning.md#r-12-1))
и раскладывает оба представления $S_4$:

```python
import itertools
import torch
import torch_vn_algebra.decompose as dc

torch.manual_seed(0)
dt = torch.complex128
n = 4
perms = list(itertools.permutations(range(n)))                   # S_4, |G| = 24

def perm_matrix(p):
    """pi(sigma) e_x = e_{sigma(x)} on C^n."""
    P = torch.zeros(n, n, dtype=dt)
    P[list(p), list(range(n))] = 1
    return P

G_node = torch.stack([perm_matrix(p) for p in perms])            # S_4 on C^4: vertex features
G_edge = torch.stack([torch.kron(P, P) for P in G_node])         # S_4 on C^16: X -> P X P^T, row-major

def n_orbits(tuples):
    """Number of orbits of S_4 on tuples of vertices: the canonical form is the smallest image."""
    return len({min(tuple(p[i] for i in t) for p in perms) for t in tuples})

bell = {r: n_orbits(itertools.product(range(n), repeat=r)) for r in (2, 3, 4)}
print(bell)
assert bell == {2: 2, 3: 5, 4: 15}                               # Corollary 12.2(b)
for G, r in ((G_node, 2), (G_edge, 4)):                          # Lemma 12.1: dim pi(G)' = number of orbits
    assert dc.commutant_basis(G, group=True).shape[0] == dc.commutant_dimension(G) == bell[r]
for X in dc.commutant_basis(G_node, group=True):                 # Corollary 12.2(a): X = a 1 + b (J - 1)
    off = X[~torch.eye(n, dtype=torch.bool)]
    assert (X.diagonal() - X[0, 0]).abs().max() < 1e-12 and (off - off[0]).abs().max() < 1e-12

D_node, D_edge = dc.from_group(G_node), dc.from_group(G_edge)
print(D_node)
print(D_edge)
assert (D_node.k_factors, D_node.multiplicities) == ((1, 3), (1, 1))
assert (D_edge.k_factors, D_edge.multiplicities) == ((1, 2, 3, 3), (2, 1, 1, 3))
```

Итак, в разложении, вычисленном `dc.from_group` (его результат проверяют утверждения assert), $\mathbb C^{16}$ содержит
представление размерности $1$ дважды — это тривиальное представление, кратность которого по предложению 12.7(a) и
лемме 12.1 равна числу $2$ орбит на $[4]^2$ (код раздела 12.4 определяет этот сектор), — представление размерности $2$
один раз и два представления размерности $3$ — один раз и три раза: $16=1\cdot2+2\cdot1+3\cdot1+3\cdot3$, а
$15=2^2+1^2+1^2+3^2$ согласуется с подсчётом орбит по лемме 12.1 и с $\dim\pi(G)'=\sum_cm_c^2$ (теорема 2.6). Для двух
представлений всё сводится к одному коммутанту:

Согласно [конспекту лекций, лемма 12.3](../../lecture_notes/ru/12_equivariant_learning.md#r-12-3), $\mathrm{Hom}_G(\pi_1,\pi_2)=Q_2^*\,\pi(G)'\,Q_1$ для $\pi=\pi_1\oplus\pi_2$, а по
[теореме 12.4](../../lecture_notes/ru/12_equivariant_learning.md#r-12-4) в подходящих базисах это $\bigoplus_cM_{b_c\times a_c}(\mathbb C)$, размерности $\sum_\rho m_\rho(\pi_1)\,m_\rho(\pi_2)$.

<a id="s-12-2"></a>

## 12.2 Эквивариантный слой средствами библиотеки

*Конспект:* [раздел 12.2](../../lecture_notes/ru/12_equivariant_learning.md#s-12-2): [теорема 12.5](../../lecture_notes/ru/12_equivariant_learning.md#r-12-5) (твирлинг); слои следуют доказательствам [леммы 12.3](../../lecture_notes/ru/12_equivariant_learning.md#r-12-3) и [теоремы 12.4](../../lecture_notes/ru/12_equivariant_learning.md#r-12-4).

Доказательство теоремы 12.4 является алгоритмом, коль скоро известно разложение из теоремы 2.6; `dc.from_group`
вычисляет его рандомизированным методом из раздела 4.3 (в точной арифметике он верен с вероятностью единица; здесь
вычисленный результат проверяют утверждения assert). При $\pi_1=\pi_2$ это просто теорема 2.6: веса — это блоки
$w_c\in M_{m_c}$, т. е. оператор абстрактного коммутанта `D.commutant_algebra()` $=\bigoplus_cM_{m_c}$, а
слой — это $x\mapsto V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*x$, вычисляемый `D.from_commutant`. Модуль
ниже — код главы, а не часть библиотеки.

```python
class EquivariantLinear(torch.nn.Module):
    """x -> T x, T = V ((+)_c 1_{k_c} (x) w_c) V^* in pi(G)'; one weight block w_c of size m_c x m_c per sector."""
    def __init__(self, D):
        super().__init__()
        self.D, self.CA = D, D.commutant_algebra()
        self.w = torch.nn.ParameterList(torch.nn.Parameter(torch.randn(m, m, dtype=dt) / m ** 0.5)
                                        for m in D.multiplicities)

    def matrix(self):
        return self.D.from_commutant(self.CA.from_blocks(list(self.w)))[0]

    def forward(self, x):                                        # x: (batch, n)
        return x @ self.matrix().T

def equivariance_error(layer, G_in, G_out):
    """max_g |layer(pi_in(g) x) - pi_out(g) layer(x)| on a random batch."""
    x = torch.randn(8, G_in.shape[-1], dtype=dt)
    with torch.no_grad():
        return max((layer(x @ g.T) - layer(x) @ h.T).abs().max().item() for g, h in zip(G_in, G_out))

edge_layer = EquivariantLinear(D_edge)
n_params = sum(p.numel() for p in edge_layer.parameters())
err = equivariance_error(edge_layer, G_edge, G_edge)
print(f"edge -> edge: {n_params} parameters, equivariance error {err:.1e}")
assert n_params == D_edge.commutant_dim == bell[4] and err < 1e-12
```

При $\pi_1\ne\pi_2$ доказательство воспроизводится строка за строкой: разложить $\pi_1\oplus\pi_2$, вычислить блоки $p_c$
проектора $P_1$ с помощью `to_commutant`, диагонализовать их и оставить левые нижние углы. След $p_c$ равен
$a_c=m_c(\pi_1)$.

Для следствия 12.2(b) при $n=4$ в вычисленном разложении (оно печатается и проверяется ниже): признаки вершин $\pi_1$
имеют $a=(1,0,0,1)$ в четырёх секторах $\pi_1\oplus\pi_2$ (размерностей $1,2,3,3$; первый — тривиальный, поскольку по
следствию 12.2(a) единственная одномерная составляющая $\pi_1$ тривиальна), признаки рёбер $\pi_2$ имеют
$b=(2,1,1,3)$, так что $\dim\mathrm{Hom}_G=1\cdot2+1\cdot3=5=B_3$, как и доказывает следствие 12.2(b).

```python
class EquivariantMap(torch.nn.Module):
    """x -> T x with T in Hom_G(pi_1, pi_2) (Theorem 12.4): one weight block y_c of size b_c x a_c per sector."""
    def __init__(self, G1, G2):
        super().__init__()
        self.n1 = n1 = G1.shape[-1]
        self.D = dc.from_group(torch.stack([torch.block_diag(g1, g2) for g1, g2 in zip(G1, G2)]))
        P1 = torch.block_diag(torch.eye(n1, dtype=dt), torch.zeros(G2.shape[-1], G2.shape[-1], dtype=dt))
        p = self.D.to_commutant(P1).matrix[0]                    # blocks p_c of P_1, padded
        self.u, self.a = [], []
        for c, m in enumerate(self.D.multiplicities):
            lam, u = torch.linalg.eigh(p[c, :m, :m])             # eigenvalues 0, ..., 0, 1, ..., 1
            self.a.append(round(lam.sum().item()))               # a_c = rank p_c = m_c(pi_1)
            self.u.append(u.flip(-1))                            # u_c^* p_c u_c = 1_{a_c} (+) 0_{b_c}
        self.y = torch.nn.ParameterList(torch.nn.Parameter(torch.randn(m - a, a, dtype=dt))
                                        for m, a in zip(self.D.multiplicities, self.a))

    def matrix(self):
        blocks = []
        for u, a, y in zip(self.u, self.a, self.y):
            x = torch.zeros(u.shape, dtype=dt)
            x[a:, :a] = y                                        # kappa_c y_c iota_c^*
            blocks.append(u @ x @ u.mH)
        T = self.D.from_commutant(self.D.commutant_algebra().from_blocks(blocks))[0]
        return T[self.n1:, :self.n1]                             # Q_2^* (...) Q_1

    def forward(self, x):
        return x @ self.matrix().T

node_to_edge = EquivariantMap(G_node, G_edge)
n_params = sum(p.numel() for p in node_to_edge.parameters())
err = equivariance_error(node_to_edge, G_node, G_edge)
print(f"node -> edge: sectors k = {node_to_edge.D.k_factors}, a = {node_to_edge.a}, "
      f"{n_params} parameters, equivariance error {err:.1e}")
assert n_params == bell[3] and err < 1e-12
b = tuple(m - a for m, a in zip(node_to_edge.D.multiplicities, node_to_edge.a))
assert node_to_edge.D.k_factors == (1, 2, 3, 3) and node_to_edge.a == [1, 0, 0, 1] and b == (2, 1, 1, 3)
corners = node_to_edge.D.basis('commutant')[:, n:, :n].reshape(-1, 16 * n)   # Lemma 12.3: Q_2^* pi(G)' Q_1
assert torch.linalg.matrix_rank(corners).item() == bell[3]
```

Второй способ сделать слой эквивариантным — усреднить его по группе. **Твирлинг** — это
$$
\mathcal P(W)=\frac1{|G|}\sum_{g\in G}\pi_2(g)\,W\,\pi_1(g)^*,\qquad W\in M_{n_2\times n_1}(\mathbb C).
$$

Согласно [конспекту лекций, теорема 12.5](../../lecture_notes/ru/12_equivariant_learning.md#r-12-5), $\mathcal P$ — ортогональный относительно скалярного произведения Гильберта–Шмидта проектор на $\mathrm{Hom}_G(\pi_1,\pi_2)$, и
$\dim\mathrm{Hom}_G(\pi_1,\pi_2)=|G|^{-1}\sum_g\overline{\chi_{\pi_1}(g)}\chi_{\pi_2}(g)$; при $\pi_1=\pi_2$ эту формулу через характеры
вычисляет `dc.commutant_dimension`. Параметризация $W\mapsto\mathcal P(W)$ тратит впустую на ядро $n_1n_2-\dim\mathrm{Hom}_G$ параметров ($256-15$ для признаков рёбер), а одно
вычисление стоит $O(|G|n_1n_2(n_1+n_2))$ для плотных $\pi(g)$, причём $|G|=n!$ для $S_n$. Библиотека вычисляет тот же
проектор по разложению, `D.conditional_expectation(W, 'commutant')`, за четыре произведения матриц $n\times n$.
Код проверяет, что `dc.twirl` и `D.conditional_expectation` совпадают, что $\mathcal P$ тождествен на $\pi(G)'$
и самосопряжён, а также формулу размерности от признаков вершин к признакам рёбер.

```python
W = torch.randn(16, 16, dtype=dt)
PW = dc.twirl(G_edge, W)
S = torch.randn(16, 16, dtype=dt)
hs = lambda A, B: (A.conj() * B).sum()                           # Hilbert-Schmidt inner product
T_layer = edge_layer.matrix().detach()
assert (PW - D_edge.conditional_expectation(W, 'commutant')[0]).abs().max() < 1e-12
assert (dc.twirl(G_edge, T_layer) - T_layer).abs().max() < 1e-12           # P = id on pi(G)'
assert abs(hs(S, PW) - hs(dc.twirl(G_edge, S), W)) < 1e-10                  # P is self-adjoint
assert abs(hs(T_layer, W - PW)) < 1e-10                                     # W - P(W) is orthogonal to pi(G)'
chi_node, chi_edge = (G.diagonal(dim1=-2, dim2=-1).sum(-1) for G in (G_node, G_edge))
assert abs((chi_node.conj() * chi_edge).mean() - bell[3]) < 1e-12           # Theorem 12.5: 5
```

<a id="s-12-3"></a>

## 12.3 Слои с ограничениями и спектральные слои

*Конспект:* [раздел 12.3](../../lecture_notes/ru/12_equivariant_learning.md#s-12-3): [предложение 12.6](../../lecture_notes/ru/12_equivariant_learning.md#r-12-6).

Веса $w=(w_c)_c$ слоя `EquivariantLinear` живут в $*$-алгебре, поэтому к ним применима любая конструкция глав 1 и 11.
Положим $\Phi(w)=V\big(\bigoplus_c1_{k_c}\otimes w_c\big)V^*$ для $w\in\bigoplus_cM_{m_c}$.

Согласно [конспекту лекций, предложение 12.6](../../lecture_notes/ru/12_equivariant_learning.md#r-12-6), (a) $\Phi$ — $*$-изоморфизм $\bigoplus_cM_{m_c}$ на $\pi(G)'$, и $\Phi(w)$
самосопряжён, унитарен или положителен тогда и только тогда, когда таков каждый $w_c$; (b) $f(\Phi(w))=\Phi\big((f(w_c))_c\big)$ для нормального $w$; (c)
эквивариантные унитарные операторы — это в точности $\Phi(u_0e^{s-s^*})$.

Итак, ортогональный (унитарный) эквивариантный слой — это `optimize.UnitaryParam` на `D.commutant_algebra()` с последующим
$\Phi$; по (c) эта параметризация достигает всякого эквивариантного унитарного оператора, и в точной арифметике каждая
итерация градиентного метода унитарна и эквивариантна (код проверяет и то и другое с точностью $10^{-12}$ после обучения).
Аналогично `PositiveParam` даёт положительные эквивариантные слои, нормированные как
в главе 11 (спектр в $[\ell,1]$ с максимумом $1$; по (a) общий положительный слой имеет вид $\Phi(t^*t)$), а
спектральный слой $f(L)$ самосопряжённого эквивариантного генератора $L$, например диффузия $e^{-tL}$, требует
спектральных разложений блоков $w_c$ стоимостью порядка $\sum_cm_c^3$ (здесь $37$) вместо $n^3=4096$; библиотека
дифференцирует через `apply_function` (о её обратном проходе мы ничего не доказываем; код ниже проверяет одну
производную). Код подгоняет эквивариантный унитарный оператор под $32$ пары вход–выход с помощью Adam — это
*численное наблюдение* данного запуска, а не теорема: функция потерь уменьшается более чем в $10^3$ раз, что проверяет
assert, и никаких утверждений о сходимости не делается, — и сверяет (b) и градиент по $t$ с плотными матрицами и с
формулой $\frac d{dt}e^{-tL}=-Le^{-tL}$ (для $L=U\mathrm{diag}(\lambda)U^*$ по теореме 1.5
$e^{-tL}=U\mathrm{diag}(e^{-t\lambda_i})U^*$, и $U$ не зависит от $t$).

```python
from torch_vn_algebra import optimize as opt

CA = D_edge.commutant_algebra()                                  # M_2 (+) M_1 (+) M_1 (+) M_3
U_target = D_edge.from_commutant(CA.random_unitary_operator(1))[0]
X = torch.randn(32, 16, dtype=dt)
Y = X @ U_target.T

U = opt.UnitaryParam(CA, 1)                                      # u_0 Haar-random, s = 0
adam = torch.optim.Adam(U.parameters(), lr=0.05)
for step in range(150):
    U_full = D_edge.from_commutant(U.operator())[0]
    loss = (X @ U_full.T - Y).abs().pow(2).mean()
    if step == 0:
        loss0 = loss.item()
    adam.zero_grad()
    loss.backward()
    adam.step()
U_full = D_edge.from_commutant(U.operator())[0].detach()
unit_err = (U_full.mH @ U_full - torch.eye(16, dtype=dt)).abs().max().item()
equi_err = max((g @ U_full - U_full @ g).abs().max().item() for g in G_edge)
print(f"loss {loss0:.3f} -> {loss.item():.1e};  |U*U - 1| = {unit_err:.1e}, |[U, pi(g)]| = {equi_err:.1e}")
assert loss.item() < 1e-3 * loss0 and unit_err < 1e-12 and equi_err < 1e-12

h = torch.randn(1, CA.C, CA.k_max, CA.k_max, dtype=dt)
h = CA.operator(h + h.mH, is_self_adjoint=True)                  # a self-adjoint generator h = (h_c)
L = D_edge.from_commutant(h)[0]
t = torch.tensor(0.7, dtype=torch.float64, requires_grad=True)
E = D_edge.from_commutant(h.apply_function(lambda w: torch.exp(-t * w)))[0]   # e^{-tL}, in the small algebra
Z = torch.randn(16, 16, dtype=dt)
(Z.conj() * E).sum().real.backward()                             # d/dt Re <Z, e^{-tL}>
dE = -L @ E.detach()
assert (E - torch.linalg.matrix_exp(-0.7 * L)).abs().max() < 1e-12          # Proposition 12.6(b)
assert abs(t.grad.item() - (Z.conj() * dE).sum().real.item()) < 1e-10
```

<a id="s-12-4"></a>

## 12.4 Инварианты

*Конспект:* [раздел 12.4](../../lecture_notes/ru/12_equivariant_learning.md#s-12-4): [предложение 12.7](../../lecture_notes/ru/12_equivariant_learning.md#r-12-7) (инвариантные векторы и квадратичные формы) и [предложение 12.8](../../lecture_notes/ru/12_equivariant_learning.md#r-12-8) (формула Молина).

Для графов инвариантные признаки рёбер порождаются индикаторами диагонали и её дополнения
($B_2=2$), а инвариантные линейные считывания — это суммы диагональных и внедиагональных элементов. Код
находит тривиальный сектор по блокам $\rho_c(g)$ (`to_algebra`) и проверяет (a), (b) и оба подсчёта из (c).

```python
rho = D_edge.to_algebra(G_edge).matrix                           # blocks rho_c(g), shape (24, C, 3, 3)
c0 = [c for c, k in enumerate(D_edge.k_factors) if k == 1 and (rho[:, c, 0, 0] - 1).abs().max() < 1e-12]
P0 = G_edge.mean(0)
assert c0 == [0] and (P0 - D_edge.central_projections()[0]).abs().max() < 1e-12      # (a)
assert abs(chi_edge.mean() - D_edge.multiplicities[0]) < 1e-12
diag = torch.eye(n, dtype=dt).reshape(16)
assert (P0 @ diag - diag).abs().max() < 1e-12 and (P0 @ (1 - diag) - (1 - diag)).abs().max() < 1e-12

w = [torch.randn(m, m, dtype=dt) for m in D_edge.multiplicities]
A = D_edge.from_commutant(CA.from_blocks([b + b.mH for b in w]))[0]                   # (b)
x = torch.randn(16, dtype=dt)
q = lambda v: (v.conj() @ A @ v).real
assert (A - A.mH).abs().max() < 1e-12
assert max(abs(q(g @ x) - q(x)).item() for g in G_edge) < 1e-12

def degree2(G):                                                  # (c): characters
    chi, chi_sq = (H.diagonal(dim1=-2, dim2=-1).sum(-1) for H in (G, G @ G))
    return ((chi ** 2 + chi_sq) / 2).mean().real.item()

pairs = list(itertools.product(range(n), repeat=2))
multisets = [(a, b) for i, a in enumerate(pairs) for b in pairs[i:]]
orbits2 = len({min(tuple(sorted(tuple(p[v] for v in e) for e in (a, b))) for p in perms) for a, b in multisets})
print(f"invariant quadratic polynomials: vertices {degree2(G_node):.4f}, edges {degree2(G_edge):.4f} "
      f"(orbits of S_4 on monomials: {orbits2})")
assert abs(degree2(G_node) - 2) < 1e-12 and abs(degree2(G_edge) - orbits2) < 1e-12 and orbits2 == 11
```

На $\mathbb C^4$ это $\sum_ix_i^2$ и $\sum_{i<j}x_ix_j$. Доказательство (c) работает в любой степени:

Согласно [конспекту лекций, предложение 12.8](../../lecture_notes/ru/12_equivariant_learning.md#r-12-8) (формула Молина),
$\sum_{d\ge0}\dim\mathbb C[x]^G_d\,t^d=|G|^{-1}\sum_g\det\big(1-t\pi(g)\big)^{-1}$, а для перестановочного представления
$\dim\mathbb C[x]^G_d$ равна числу орбит $G$ на мономах степени $d$.

(Это теорема Молина 1897 года; для дальнейшего чтения см. B. Sturmfels, *Algorithms in Invariant Theory*, Springer 2008,
раздел 2.2.) Код раскладывает правую часть до $t^4$, умножая на один геометрический ряд для каждого
собственного значения, и сравнивает результат с орбитами мономов; для $S_4$ на $\mathbb C^4$ числа $1,1,2,3,5$ — это
числа разбиений $d$ не более чем на четыре слагаемых (орбиты мультимножеств определяются их упорядоченными кратностями).

```python
def molien(G, d_max):
    """Coefficients of |G|^{-1} sum_g det(1 - t pi(g))^{-1} up to t^d_max: products of geometric series."""
    lam = torch.linalg.eigvals(G)                                # (|G|, n), |lambda| = 1
    c = torch.zeros(G.shape[0], d_max + 1, dtype=dt)
    c[:, 0] = 1
    for i in range(lam.shape[1]):                                # multiply by 1 / (1 - t lambda_i)
        for d in range(1, d_max + 1):
            c[:, d] = c[:, d] + lam[:, i] * c[:, d - 1]
    return c.mean(0)

def monomial_orbits(G, d):
    """Orbits of G on the monomials of degree d, i.e. on the multisets of d points of Omega."""
    act = [g.real.argmax(0).tolist() for g in G]                 # pi(g) e_x = e_{act[x]}
    mons = itertools.combinations_with_replacement(range(G.shape[-1]), d)
    return len({min(tuple(sorted(a[x] for x in m)) for a in act) for m in mons})

for G in (G_node, G_edge):
    series, counts = molien(G, 4), [monomial_orbits(G, d) for d in range(5)]
    print("Molien series:", [f"{v.real:.4f}" for v in series], " orbits of monomials:", counts)
    assert (series - torch.tensor(counts, dtype=dt)).abs().max() < 1e-9          # Proposition 12.8
assert [monomial_orbits(G_node, d) for d in range(5)] == [1, 1, 2, 3, 5] and monomial_orbits(G_edge, 2) == orbits2
```

<a id="s-12-5"></a>

## 12.5 С чем сравнивать

Эквивариантные слои для конкретных групп предоставляются зрелыми библиотеками. **escnn** (Cesa, Lang и Weiler 2022,
преемник e2cnn) строит управляемые (steerable) свёртки для подгрупп $O(2)$ и $O(3)$ из неприводимых
представлений и гармонических базисов ядер; **e3nn** (Geiger и Smidt 2022) строит $E(3)$-эквивариантные сети
из сферических гармоник и тензорных произведений Клебша–Гордана. **EMLP** (Finzi, Welling и Wilson 2021) работает с любой
матричной группой: он налагает $\pi_2(g)T=T\pi_1(g)$ для дискретных образующих и $d\pi_2(X)T=Td\pi_1(X)$ для образующих
алгебры Ли и вычисляет базис решений как численное ядро — это линейная задача
`dc.commutant_basis`. Групповые свёртки (Cohen и Welling 2016) — это случай $\Omega_1=\Omega_2=G$ леммы 12.1,
где $G$ действует на себе левыми сдвигами: $(y,x)$ и $(y',x')$ лежат в одной орбите тогда и только тогда, когда
$x^{-1}y=x'^{-1}y'$ (в этом случае $g=x'x^{-1}$ переводит одну пару в другую), так что эквивариантные отображения — это
$(Tv)(y)=\sum_xf(x^{-1}y)v(x)$, свёртки с функциями $f$ на $G$. Kondor и Trivedi (2018) изучают аналогичный вопрос
для однородных пространств компактных групп, который в этом курсе не рассматривается. Для перестановок DeepSets (Zaheer et al. 2017) и базисы
Maron et al. (2019) — это следствие 12.2.

Библиотека не предоставляет свёрток на решётках, управляемых ядер на непрерывных пространствах, разложений над
$\mathbb R$ (вещественный, комплексный и кватернионный типы; упражнение 2 рассматривает вещественные матрицы) и не рассчитана на большие $n$: $V$ — плотная
матрица $n\times n$, и разложение стоит $O(|G|n^3)$ на одну случайную выборку для группы, заданной всеми своими элементами
(образующие алгебры Ли обрабатываются через `dc.from_generators`; раздел 4.3). Что она добавляет, так это
**структуру** пространства решений: базис коммутанта — это лишь список матриц, тогда как разложение
отождествляет коммутант с $*$-алгеброй $\bigoplus_cM_{m_c}$. Произведения, сопряжения, унитарность, положительность и
функциональное исчисление тогда действуют поблочно на маленьких матрицах (предложение 12.6), библиотека
дифференцирует через них, а
`to_commutant`/`from_commutant` переводят из одной картины в другую.

<a id="x-12"></a>

## Упражнения

[Упражнения](../../lecture_notes/ru/12_equivariant_learning.md#x-12) 1, 3 и 4 конспекта лекций, численно:

- *Упражнение 1* (DeepSets с каналами): проверьте с помощью `EquivariantMap`, что эквивариантные отображения из
  $\mathbb C^n\otimes\mathbb C^{d_1}$ в $\mathbb C^n\otimes\mathbb C^{d_2}$ — это в точности $1\otimes A+J\otimes B$.
- *Упражнение 3* (нелинейности): покоординатный ReLU коммутирует с перестановочным представлением, но не
  с поворотом $\mathbb R^2$ на $\pi/2$; $x\mapsto\varphi(\|x\|)x$ коммутирует с любым унитарным представлением.
- *Упражнение 4*: для представления на рёбрах одна гауссова выборка оставляет неопределённым $2$-мерное пространство
  эквивариантных слоёв, две выборки определяют слой.

<details>
<summary>Решения</summary>

*Упражнение 1.*

```python
G_in = torch.stack([torch.kron(P, torch.eye(2, dtype=dt)) for P in G_node])     # d_1 = 2
G_out = torch.stack([torch.kron(P, torch.eye(3, dtype=dt)) for P in G_node])    # d_2 = 3
deepsets = EquivariantMap(G_in, G_out)
assert sum(p.numel() for p in deepsets.parameters()) == 2 * 2 * 3
assert equivariance_error(deepsets, G_in, G_out) < 1e-12
T = deepsets.matrix().detach().reshape(n, 3, n, 2)                # blocks T[i, :, j, :] = delta_ij A + B
A_, B_ = T[0, :, 0] - T[0, :, 1], T[0, :, 1]
assert all((T[i, :, j] - (A_ if i == j else 0) - B_).abs().max() < 1e-12 for i in range(n) for j in range(n))
```

*Упражнение 3.*

```python
relu = lambda v: v.clamp(min=0)
R = torch.tensor([[0.0, -1.0], [1.0, 0.0]], dtype=torch.float64)
v = torch.tensor([1.0, -1.0], dtype=torch.float64)
assert (relu(R @ v) - R @ relu(v)).abs().max() == 1                             # (b)
xr = torch.randn(16, dtype=torch.float64)
assert all(torch.equal(relu(g.real @ xr), g.real @ relu(xr)) for g in G_edge)   # (a): exact
gate = lambda v: torch.tanh(v.norm()) * v
assert max((gate(g @ x) - g @ gate(x)).abs().max().item() for g in G_edge) < 1e-12   # (c)
```

*Упражнение 4.*

```python
basis = D_edge.basis('commutant')                                # (15, 16, 16), Hilbert-Schmidt orthonormal
for N, lost in ((1, 2), (2, 0)):
    Xs = torch.randn(16, N, dtype=dt)
    M = torch.stack([(B @ Xs).reshape(-1) for B in basis], dim=1)               # w -> (T x_1, ..., T x_N)
    assert torch.linalg.matrix_rank(M).item() == 15 - lost
```

</details>

**Далее:** [Руководство 13. Йордановы алгебры: почему квантовая механика использует *-алгебры](13_jordan_algebras.md).
