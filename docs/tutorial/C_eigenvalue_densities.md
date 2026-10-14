**English** · [Español](es/C_eigenvalue_densities.md) · [Français](fr/C_eigenvalue_densities.md) · [Deutsch](de/C_eigenvalue_densities.md) · [中文](zh/C_eigenvalue_densities.md) · [日本語](ja/C_eigenvalue_densities.md) · [Русский](ru/C_eigenvalue_densities.md)

# Appendix C. Eigenvalue densities and Page's formula

Chapter 5 uses the joint density of the eigenvalues of the Gaussian ensembles, and Chapter 7 uses
Page's formula for the mean entanglement entropy of a random pure state. This appendix proves both.
The route is: the Weyl integration formula (the Jacobian of "diagonalisation"), the Gaussian
ensembles, the one-point function of a determinantal density by orthogonal polynomials, the
Wishart–Laguerre ensemble and the induced measure on reduced states, and finally an exact
computation with Laguerre polynomials.

**You will learn** why the eigenvalues of an invariant random matrix repel like
$\prod_{i<j}|\lambda_i-\lambda_j|^\beta$, why the eigenvectors are then Haar distributed and independent
of the eigenvalues, how a determinant of orthonormal polynomials integrates out, and how
$\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}$ follows from a handful of integrals.

Besides the prerequisites of the course (measure theory, the change-of-variables formula in
$\mathbb R^N$ for injective $C^1$ maps with invertible derivative, the inverse and implicit function
theorems), we use from Chapter 5 Lemma 5.1 (uniqueness and right invariance of Haar measure),
Lemma 5.2 and Theorem 5.3 (existence of Haar measure on $U(n)$ by Gaussian QR), from Appendix B
Theorem B.17 (the same for $O(n)$), from Appendix A Proposition A.16(a) (unitary invariance of complex
Gaussian vectors) and Theorem A.20(b) (the law of $|z|^2$), and from Chapter 7 the
representation $\rho_A=GG^*/\operatorname{Tr}GG^*$ of the proof of Theorem 7.11.

**Notation.** $\beta\in\{1,2\}$ and $\mathbb F=\mathbb R$ or $\mathbb C$ accordingly. $\mathcal H$ is the
space of real symmetric ($\beta=1$) or complex Hermitian ($\beta=2$) $n\times n$ matrices, with
Lebesgue measure $dH$ in the coordinates $h_{ii}$, $\operatorname{Re}h_{ij}$ and (for $\beta=2$)
$\operatorname{Im}h_{ij}$, $i<j$; its dimension is $n+\beta n(n-1)/2$. $K=O(n)$ or $U(n)$, $T$ is the
group of diagonal elements of $K$ (signs, or phases), $\mathfrak k=\{X\in M_n(\mathbb F):X^*=-X\}$,
$\mathfrak t$ its diagonal part and $\mathfrak m$ its off-diagonal part (zero diagonal), with
coordinates $\operatorname{Re}X_{ij}$, $\operatorname{Im}X_{ij}$, $i<j$. For $\lambda\in\mathbb R^n$,
$\Lambda=\operatorname{diag}(\lambda)$ and $\Delta(\lambda)=\prod_{i<j}(\lambda_j-\lambda_i)$; $C=\{\lambda:\lambda_1<\dots<\lambda_n\}$,
on which $\Delta>0$. For $H\in\mathcal H$, $\lambda(H)$ is the vector of its eigenvalues in ascending
order. $\|X\|_2^2=\operatorname{Tr}X^*X$. For $V\in K$ the map $H\mapsto VHV^*$ is linear on $\mathcal H$
and preserves $\|\cdot\|_2$; a linear isometry has $|\det|=1$, so it **preserves $dH$**.

## C.1 Simple spectra and eigenframes

**Lemma C.1.** (a) The set $\Omega\subseteq\mathcal H$ of matrices with $n$ distinct eigenvalues is
open, and $\mathcal H\setminus\Omega$ is a Lebesgue null set. (b) On $\Omega$ the map $H\mapsto\lambda(H)$ is
$C^\infty$, and the spectral projections $P_j(H)=\prod_{k\neq j}\frac{H-\lambda_k(H)}{\lambda_j(H)-\lambda_k(H)}$ are continuous.

*Proof.* (a) Let $A_{ij}=\lambda_j^{i-1}$ with $\lambda=\lambda(H)$. Then $\det A=\Delta(\lambda)$ (Vandermonde) and
$(AA^T)_{ij}=\sum_k\lambda_k^{i+j-2}=\operatorname{Tr}H^{i+j-2}$, so
$\operatorname{disc}(H)=\det\big(\operatorname{Tr}H^{i+j-2}\big)_{i,j=1}^n=\Delta(\lambda(H))^2$ is a polynomial
in the coordinates of $H$, nonzero at $\operatorname{diag}(1,\dots,n)$, and $\Omega=\{\operatorname{disc}\neq0\}$ is open.
The zero set of a nonzero polynomial $p$ on $\mathbb R^N$ is null, by induction on $N$: for $N=1$ it is
finite; otherwise write $p=\sum_jp_j(x')x_N^j$ with some $p_{j_0}\neq0$; for $x'$ outside the zero set of
$p_{j_0}$ (null by induction) $x_N\mapsto p(x',x_N)$ has finitely many zeros, and Tonelli's theorem applied
to the indicator of the zero set gives measure $0$.
(b) $\chi(H,t)=\det(t-H)=\prod_i(t-\lambda_i)$ is a real polynomial in $(H,t)$. At $H_0\in\Omega$,
$\partial_t\chi(H_0,\lambda_j(H_0))=\prod_{k\neq j}(\lambda_j-\lambda_k)\neq0$, so the implicit function theorem gives
$C^\infty$ functions $r_j$ near $H_0$ with $\chi(H,r_j(H))=0$ and $r_j(H_0)=\lambda_j(H_0)$. By continuity
$r_1<\dots<r_n$ near $H_0$: these are $n$ distinct roots of a polynomial of degree $n$, hence all of them, in
order, and $r_j=\lambda_j$. If $H=\sum_k\lambda_kQ_k$ is the spectral decomposition, Lagrange interpolation gives
$P_j(H)=Q_j$; continuity follows from that of $\lambda$. $\square$

For $H\in\Omega$, an **eigenframe** is a $V\in K$ with $H=V\Lambda V^*$, $\Lambda=\operatorname{diag}\lambda(H)$. If $V,V'$
are eigenframes, $M=V'^*V$ commutes with $\Lambda$, so $(\lambda_i-\lambda_j)M_{ij}=0$ and $M\in T$; conversely $VD$
is an eigenframe for every $D\in T$. A **Borel eigenframe** $F(H)$ is obtained by taking as $j$-th column
$P_j(H)e_m/\|P_j(H)e_m\|$ with $m$ the least index such that $P_j(H)e_m\neq0$: it is continuous on each of the
Borel sets where these indices are fixed, its columns are orthonormal eigenvectors, and it is real when $H$ is.

**Lemma C.2 (Haar measures).** $U(n)$, $O(n)$ and $T$ have Haar measures, and $\int_Kf(VW)\,dV=\int_Kf\,dV$
for all $W\in K$ and bounded Borel $f$.

*Proof.* $U(n)$: Theorem 5.3. $O(n)$: Theorem B.17 (the real Gaussian QR construction: a real Gaussian matrix is
invertible almost surely, and the orthogonal factor of its QR decomposition is Haar distributed). $T$: the uniform measure on $\{\pm1\}^n$, or the product of normalised arc length on the
circles, is invariant under multiplication by a fixed element. Right invariance is Lemma 5.1. $\square$

**Lemma C.3 (eigenframes of invariant ensembles).** Let $H$ be a random element of $\mathcal H$ with
$P(H\in\Omega)=1$ whose law is invariant under $H\mapsto VHV^*$ for every $V\in K$. Let $D$ be uniform on $T$ and
independent of $H$, and $\tilde U=F(H)D$. Then $\tilde U$ is an eigenframe of $H$, it is Haar distributed on
$K$, and it is independent of $\lambda(H)$.

*Proof.* Fix $V\in K$. Both $F(VHV^*)$ and $VF(H)$ are eigenframes of $VHV^*$, so
$F(VHV^*)=VF(H)E(H)$ with $E(H)=(VF(H))^*F(VHV^*)\in T$ a Borel function of $H$. For bounded Borel
$\varphi$ on $\mathbb R^n$ and $\psi$ on $K$, since $VHV^*$ has the law of $H$ and $\lambda(VHV^*)=\lambda(H)$,
$$
\mathbb E\big[\varphi(\lambda(H))\psi(\tilde U)\big]=\mathbb E\big[\varphi(\lambda(H))\psi(VF(H)E(H)D)\big]
=\mathbb E\big[\varphi(\lambda(H))\psi(V\tilde U)\big],
$$
where the second step integrates first over $D$ (Fubini): for fixed $H$, $E(H)D$ is uniform on $T$. Now let $V$
be Haar distributed and independent of everything; by Fubini and Lemma C.2,
$\mathbb E[\varphi(\lambda(H))\psi(\tilde U)]=\mathbb E\big[\varphi(\lambda(H))\int_K\psi(V\tilde U)\,dV\big]
=\mathbb E[\varphi(\lambda(H))]\int_K\psi\,dV$. With $\varphi=1$, $\tilde U$ is Haar; the product form is the
independence. $\square$

## C.2 The Weyl integration formula

For $\varepsilon>0$ let $\mathfrak m_\varepsilon=\{X\in\mathfrak m:\|X\|_2<\varepsilon\}$ and
$$
\Psi(X,\lambda)=e^X\Lambda e^{-X},\qquad (X,\lambda)\in\mathfrak m_\varepsilon\times C .
$$
For $X\in\mathfrak k$, $(e^X)^*=e^{-X}=(e^X)^{-1}$, so $e^X\in K$.

**Lemma C.4 (a chart).** There are $\varepsilon>0$ and a continuous $j:\overline{\mathfrak m_\varepsilon}\to(0,\infty)$
with $j(0)=1$ such that: (a) $|\det D\Psi(X,\lambda)|=\Delta(\lambda)^\beta j(X)$; (b) $\Psi$ is injective; (c) for
$H\in\Omega$, $H\in\Psi(\mathfrak m_\varepsilon\times C)$ iff $F(H)\in S=\{e^XD:X\in\mathfrak m_\varepsilon,D\in T\}$, and
$S$ is a Borel set with $ST=S$.

*Proof.* (a) $\Psi$ is $C^\infty$ (power series). For $Y\in\mathfrak m$ put $A_X(Y)=e^{-X}\frac{d}{dt}\big|_0e^{X+tY}$.
Differentiating $W(t)^*W(t)=1$ for $W(t)=e^{X+tY}\in K$ gives $A_X(Y)\in\mathfrak k$, and with
$\frac{d}{dt}W^{-1}=-W^{-1}W'W^{-1}$, for $g=e^X$,
$$
D\Psi(X,\lambda)(Y,\mu)=g\big([A_X(Y),\Lambda]+\operatorname{diag}(\mu)\big)g^{-1}.
$$
Diagonal matrices commute with $\Lambda$, so $[A_X(Y),\Lambda]=[B_XY,\Lambda]$ with $B_X=\pi\circ A_X:\mathfrak m\to\mathfrak m$ and
$\pi$ the projection of $\mathfrak k$ onto $\mathfrak m$ along $\mathfrak t$. Thus $D\Psi=\operatorname{Ad}_g\circ L_\lambda\circ(B_X\times1)$
with $L_\lambda(Z,\mu)=[Z,\Lambda]+\operatorname{diag}(\mu)$. Since $([Z,\Lambda])_{ij}=Z_{ij}(\lambda_j-\lambda_i)$, in the coordinates
$L_\lambda$ multiplies $\operatorname{Re}Z_{ij}$ (and $\operatorname{Im}Z_{ij}$) by $\lambda_j-\lambda_i$ and maps $\mu$ to the diagonal:
$|\det L_\lambda|=\Delta(\lambda)^\beta$. $\operatorname{Ad}_g$ preserves $dH$. Put $j(X)=|\det B_X|$: it is continuous, and
$B_0=1$ because $\frac d{dt}e^{tY}|_0=Y$; choose $\varepsilon$ with $j>0$ on $\overline{\mathfrak m_\varepsilon}$.
(b) If $\Psi(X,\lambda)=\Psi(X',\lambda')=H$, then $\lambda=\lambda'=\lambda(H)$, and $e^X,e^{X'}$ are eigenframes of $H\in\Omega$,
so $e^X=e^{X'}D$ with $D\in T$. Let $G(X,t)=e^Xe^t$ on $\mathfrak m\times\mathfrak t$ (for $\beta=1$, $\mathfrak t=\{0\}$); $DG(0,0)(Y,s)=Y+s$ is injective,
and as $Y$ is off-diagonal and $s$ diagonal, $\|Y+s\|_2^2=\|Y\|_2^2+\|s\|_2^2$, so $\|DG(0,0)v\|=\|v\|$ for the product norm; put $m=1$. Choose $\delta>0$ with $\|DG(z)-DG(0)\|\le m/2$ for $\|z\|<\delta$; for $u,v$ in this
ball, $G(u)-G(v)=\int_0^1DG(v+s(u-v))(u-v)\,ds$ gives $\|G(u)-G(v)\|\ge\frac m2\|u-v\|$, so $G$ is injective there.
($\delta$ depends only on $n,\beta$; $\varepsilon$ is chosen after it.)
Now $\|D-1\|_2=\|e^{X'}(D-1)\|_2=\|e^X-e^{X'}\|_2\le\|e^X-1\|_2+\|e^{X'}-1\|_2\le2(e^\varepsilon-1)$, by the power series and
$\|AB\|_2\le\|A\|_2\|B\|_2$. For $\beta=1$, $D=\operatorname{diag}(\pm1)$ and an entry $-1$ would give
$\|D-1\|_2\ge2$; for $\varepsilon<\log2$, $D=1=e^t$ with $t=0$. For $\beta=2$, $D=e^t$ with $t=i\operatorname{diag}(\theta)$,
$\theta_j\in(-\pi,\pi]$, and $|\theta_j|\le\frac\pi2|e^{i\theta_j}-1|$ (as $|e^{i\theta}-1|=2|\sin\frac\theta2|$ and $\sin x\ge2x/\pi$ on $[0,\pi/2]$), so
$\|t\|_2\le\frac\pi2\|D-1\|_2\le\pi(e^\varepsilon-1)$. For $\varepsilon$ small ($\varepsilon^2+\pi^2(e^\varepsilon-1)^2<\delta^2$), $(X,0)$ and $(X',t)$ lie in the $\delta$-ball and
$G(X,0)=e^X=e^{X'}D=G(X',t)$, so $X=X'$.
(c) If $H=\Psi(X,\lambda')$, then $e^X$ is an eigenframe, so $F(H)\in e^XT\subseteq S$. If $F(H)=e^XD$, then
$H=F(H)\Lambda F(H)^*=\Psi(X,\lambda(H))$. $S=G(\mathfrak m_\varepsilon\times\mathfrak t)$ (or $e^{\mathfrak m_\varepsilon}T$, $T$ finite) is a
continuous image of a $\sigma$-compact set, hence $\sigma$-compact and Borel; $ST=S$ is clear. $\square$

**Theorem C.5 (Weyl integration formula).** There is a constant $c_{n,\beta}>0$ such that for every Borel
$f:\mathcal H\to[0,\infty]$,
$$
\int_{\mathcal H}f(H)\,dH=c_{n,\beta}\int_C\Delta(\lambda)^\beta\Big(\int_Kf(V\Lambda V^*)\,dV\Big)d\lambda .
$$
Equivalently, $\Phi(V,\lambda)=V\Lambda V^*$ maps the measure $c_{n,\beta}\,dV\otimes\Delta(\lambda)^\beta d\lambda$ on $K\times C$
to Lebesgue measure on $\mathcal H$. The factor $\Delta(\lambda)^\beta$ is the Jacobian of the chart $\Psi$ at $X=0$
(Lemma C.4(a)); the constant $c_{n,\beta}$ reflects the normalisation of the Haar measure and is never needed
explicitly below, since every density obtained from the formula is normalised by its total mass.

*Proof.* *Reference ensemble.* Let $H_0$ have the density $p_0(H)=Z_0^{-1}e^{-\operatorname{Tr}H^2/2}$; since
$\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$, $Z_0$ is a finite product of Gaussian integrals. $p_0$ and $dH$
are invariant under $H\mapsto VHV^*$, hence so is the law of $H_0$; $P(H_0\in\Omega)=1$ by Lemma C.1. By Lemma C.3,
$\tilde U$ is Haar, independent of $\lambda=\lambda(H_0)$, and $H_0=\tilde U\Lambda\tilde U^*$.
*The law of $\lambda$.* Let $A\subseteq C$ be Borel. By Lemma C.4(c) and $ST=S$,
$\{H_0\in\Psi(\mathfrak m_\varepsilon\times A)\}=\{\lambda\in A,\ \tilde U\in S\}$ up to a null event, whose probability is
$P(\lambda\in A)\,h$ with $h=\int_K1_S\,dV$. On the other hand $\Psi$ is injective and $C^1$ with nonzero Jacobian on the
open set $\mathfrak m_\varepsilon\times C$, and $p_0(\Psi(X,\lambda))=p_0(\Lambda)$; the change-of-variables formula and Lemma C.4(a) give
$$
P\big(H_0\in\Psi(\mathfrak m_\varepsilon\times A)\big)=\int_{\mathfrak m_\varepsilon}\!\int_Ap_0(\Lambda)\Delta(\lambda)^\beta j(X)\,d\lambda\,dX
=\kappa\int_Ap_0(\Lambda)\Delta(\lambda)^\beta\,d\lambda,\qquad \kappa=\int_{\mathfrak m_\varepsilon}j>0 .
$$
With $A=C$, $h=\kappa\int_Cp_0\Delta^\beta>0$. Hence $\lambda$ has the density $c\,p_0(\Lambda)\Delta(\lambda)^\beta$ on $C$, $c=\kappa/h$.
*General $f$.* As $p_0>0$, $\int f\,dH=\mathbb E\big[(f/p_0)(\tilde U\Lambda\tilde U^*)\big]$. The joint law of $(\tilde U,\lambda)$ is
$dV\otimes c\,p_0(\Lambda)\Delta^\beta d\lambda$, and $p_0(V\Lambda V^*)=p_0(\Lambda)$; Tonelli's theorem gives the formula with $c_{n,\beta}=c$.
The constant is determined by the formula itself (take $f=p_0$). $\square$

**Corollary C.6 (invariant densities).** Let $H$ have a density $p$ on $\mathcal H$ with $p(VHV^*)=p(H)$ for all
$V\in K$. Then $\lambda(H)$ has the density $c_{n,\beta}\,p(\Lambda)\Delta(\lambda)^\beta$ on $C$; the eigenframe $\tilde U$ of Lemma C.3
is Haar distributed and independent of $\lambda(H)$; and for symmetric Borel $g\ge0$ on $\mathbb R^n$,
$$
\mathbb E\,g(\lambda(H))=\frac{c_{n,\beta}}{n!}\int_{\mathbb R^n}g(x)\,p(\operatorname{diag}x)\,|\Delta(x)|^\beta\,dx .
$$

*Proof.* $P(H\in\Omega)=1$ by Lemma C.1. For Borel $A\subseteq C$ apply Theorem C.5 to $f=1_A(\lambda(\cdot))\,p$, for which
$f(V\Lambda V^*)=1_A(\lambda)p(\Lambda)$. The law of $H$ is invariant, so Lemma C.3 applies. Up to the null set where two
coordinates coincide, $\mathbb R^n$ is the disjoint union of the $n!$ images of $C$ under coordinate permutations, and
$g$, $|\Delta|$ and $x\mapsto p(\operatorname{diag}x)$ are permutation invariant (a permutation matrix lies in $K$). $\square$

The check below differentiates $\Psi$ numerically for $n=2,3$ and both $\beta$: at $X=0$ the Jacobian is $\Delta(\lambda)^\beta$, and
at a fixed $X\neq0$ its ratio to $\Delta(\lambda)^\beta$ does not depend on $\lambda$ (Lemma C.4(a)).

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

## C.3 The Gaussian ensembles

**Theorem C.7 (GOE and GUE).** Let $H=(G+G^T)/2$ with $G$ real with independent $N(0,1)$ entries ($\beta=1$, GOE), or
$H=(G+G^*)/\sqrt2$ with $G$ Ginibre ($\beta=2$, GUE), as in Section 5.3. Then $H$ has the density
$(2\pi)^{-n/2}\pi^{-\beta n(n-1)/4}e^{-\operatorname{Tr}H^2/2}$ on $\mathcal H$, its eigenvalues have the density
$$
Z_{n,\beta}^{-1}\,\Delta(\lambda)^\beta\,e^{-\sum_i\lambda_i^2/2}\quad\text{on }C,\qquad Z_{n,\beta}=\int_C\Delta^\beta e^{-|\lambda|^2/2}\,d\lambda\in(0,\infty),
$$
and the eigenframe $\tilde U$ is Haar distributed on $K$ and independent of the eigenvalues.

*Proof.* GOE: $h_{ii}=G_{ii}\sim N(0,1)$ and $h_{ij}=(G_{ij}+G_{ji})/2\sim N(0,\frac12)$ for $i<j$. GUE: $h_{ii}=\sqrt2\operatorname{Re}G_{ii}\sim N(0,1)$,
$\operatorname{Re}h_{ij}=(\operatorname{Re}G_{ij}+\operatorname{Re}G_{ji})/\sqrt2$ and $\operatorname{Im}h_{ij}=(\operatorname{Im}G_{ij}-\operatorname{Im}G_{ji})/\sqrt2$ are $N(0,\frac12)$.
Distinct coordinates are functions of disjoint sets of independent variables, so they are independent, and the
density is $\prod_i(2\pi)^{-1/2}e^{-h_{ii}^2/2}\prod_{i<j}\pi^{-\beta/2}e^{-|h_{ij}|^2}$, which is the stated one because
$\operatorname{Tr}H^2=\sum_ih_{ii}^2+2\sum_{i<j}|h_{ij}|^2$. It is invariant, and Corollary C.6 applies; the constant is fixed by
normalisation. $Z_{n,\beta}<\infty$ because $|\lambda_j-\lambda_i|\le(1+|\lambda_i|)(1+|\lambda_j|)$, so
$|\Delta(\lambda)|\le\prod_i(1+|\lambda_i|)^{n-1}$. $\square$

**Corollary C.8 (level density).** For the GOE or GUE of size $n$, let $q(x)=(n!Z_{n,\beta})^{-1}|\Delta(x)|^\beta e^{-|x|^2/2}$ on
$\mathbb R^n$ and $\rho(x)=n\int_{\mathbb R^{n-1}}q(x,y)\,dy$. Then $\mathbb E\sum_i\varphi(\lambda_i(H))=\int\varphi\rho$ for bounded Borel
$\varphi$, $\int\rho=n$, and $\rho$ is continuous with $0\le\rho(x)\le C_n(1+|x|)^{\beta(n-1)}e^{-x^2/2}$. In particular the
hypothesis of Theorem 5.15 holds for independent GOE or GUE blocks.

*Proof.* By Corollary C.6 with $g(x)=\sum_i\varphi(x_i)$ (for $\varphi\ge0$; in general apply it to $\varphi^\pm=\max(\pm\varphi,0)$) and the symmetry of $q$,
$\mathbb E\sum_i\varphi(\lambda_i)=n\int\varphi(x_1)q(x)\,dx=\int\varphi\rho$ (Fubini). By the bound on $|\Delta|$,
$q(x,y)\le(n!Z_{n,\beta})^{-1}(1+|x|)^{\beta(n-1)}e^{-x^2/2}\prod_{i}(1+|y_i|)^{\beta(n-1)}e^{-y_i^2/2}$, which gives the bound after
integrating over $y$, and, since $\sup_x(1+|x|)^{\beta(n-1)}e^{-x^2/2}<\infty$, an integrable majorant independent of $x$: continuity
follows by dominated convergence. $\square$

**Proposition C.9 (smallest spacing).** For the GOE or GUE of size $n\ge2$ let $\delta=\min_i(\lambda_{i+1}-\lambda_i)$. There are
constants $0<c\le C$ (depending on $n,\beta$) with $c\,\varepsilon^{\beta+1}\le P(\delta<\varepsilon)\le C\varepsilon^{\beta+1}$ for $0<\varepsilon\le1$.

*Proof.* *Upper bound.* $P(\delta<\varepsilon)\le\sum_iP(\lambda_{i+1}-\lambda_i<\varepsilon)$. On $\{0<u=x_{i+1}-x_i<\varepsilon\}$ bound the factor of
the pair $(i,i+1)$ by $u^\beta$, the others as above, $1+|x_{i+1}|\le2(1+|x_i|)$ and $e^{-x_{i+1}^2/2}\le1$; substituting
$x_{i+1}=x_i+u$,
$$
P(\lambda_{i+1}-\lambda_i<\varepsilon)\le\frac{2^{\beta(n-1)}}{Z_{n,\beta}}\int_0^\varepsilon u^\beta du\int(1+|x|)^{2\beta(n-1)}e^{-\frac{x^2}2}dx
\Big(\int(1+|x|)^{\beta(n-1)}e^{-\frac{x^2}2}dx\Big)^{n-2}.
$$
*Lower bound.* Let $R_\varepsilon=\{x:0<x_1<1,\ 0<x_2-x_1<\varepsilon,\ 3k<x_k<3k+1\ (3\le k\le n)\}\subseteq C$. On $R_\varepsilon$,
$\delta\le x_2-x_1$, every pair other than $(1,2)$ is at distance $\ge1$, and $|x|^2\le M$ with $M$ independent of $\varepsilon\le1$; so
$P(\delta<\varepsilon)\ge P(\lambda\in R_\varepsilon)\ge Z_{n,\beta}^{-1}e^{-M/2}\int_0^\varepsilon u^\beta du$. $\square$

For the GOE this is the order $\varepsilon^2$ quoted at the end of Section 5.4. We check Theorem C.7 for $n=3$: the GOE
level density is computed from the joint density by integrating over the other two eigenvalues on a grid, and compared
with the mean number of eigenvalues in each bin; we also check that the eigenvectors are uncorrelated with the
eigenvalues ($|\tilde U_{11}|^2$ does not depend on the phases $D$, so it can be read off any eigenframe returned by
`eigh`). For Haar $U\in O(3)$, $\mathbb E|U_{11}|^2=\frac13$: left multiplication by a permutation matrix permutes the
entries of the first column and preserves the Haar measure, so the three $\mathbb E|U_{i1}|^2$ are equal, and they sum to $1$.

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

## C.4 Orthogonal polynomials and the one-point function

**Lemma C.10 (one-point function).** Let $I\subseteq\mathbb R$ be an interval and $w\ge0$ a Borel function on $I$, positive on a
set of positive measure, with $\int_I|x|^mw<\infty$ for all $m$. Let $p_0,\dots,p_{n-1}$ be real polynomials, $\deg p_k=k$
with leading coefficient $\kappa_k\neq0$, orthonormal in $L^2(w)$ (they exist: apply Gram–Schmidt to $1,x,\dots,x^{n-1}$, which are
independent in $L^2(w)$ since a nonzero polynomial has finitely many zeros; any sign may be chosen for each $p_k$). Let $q=Z_n^{-1}\Delta(x)^2\prod_iw(x_i)$ on
$I^n$ with $Z_n=\int_{I^n}\Delta^2\prod_iw(x_i)\,dx$. Then $Z_n=n!\prod_k\kappa_k^{-2}$ and
$$
n\int_{I^{n-1}}q(x,y)\,dy=K_n(x)\,w(x),\qquad K_n(x)=\sum_{k=0}^{n-1}p_k(x)^2=v(x)^TM^{-1}v(x),
$$
where $v(x)=(1,x,\dots,x^{n-1})$ and $M_{ij}=\int_Ix^{i+j}w$, $0\le i,j<n$.

*Proof.* $\Delta(x)=\det(x_j^{i-1})_{i,j}$. Adding to each row multiples of the rows above it replaces $x_j^{i-1}$ by
$p_{i-1}(x_j)/\kappa_{i-1}$, so $\Delta(x)=(\prod_k\kappa_k)^{-1}\det(p_{i-1}(x_j))$ and
$$
\Delta(x)^2=\Big(\prod_k\kappa_k\Big)^{-2}\sum_{\sigma,\tau\in S_n}\operatorname{sgn}\sigma\operatorname{sgn}\tau\prod_{j=1}^np_{\sigma(j)-1}(x_j)\,p_{\tau(j)-1}(x_j).
$$
Integrating over $x_2,\dots,x_n$ against $\prod_{j\ge2}w(x_j)$, the factor $j$ gives $\delta_{\sigma(j)\tau(j)}$, so only $\sigma=\tau$
survives, and the result is $(\prod_k\kappa_k)^{-2}\sum_\sigma p_{\sigma(1)-1}(x_1)^2=(\prod_k\kappa_k)^{-2}(n-1)!\,K_n(x_1)$.
Integrating against $w(x_1)$ gives $Z_n=(\prod_k\kappa_k)^{-2}n!$, and dividing gives the formula. Finally, if
$p_k(x)=\sum_iC_{ki}x^i$, orthonormality reads $CMC^T=1$ with $C$ invertible, so $M^{-1}=C^TC$ and
$v^TM^{-1}v=|Cv|^2=\sum_kp_k^2$. $\square$

By Corollary C.6, the unordered eigenvalues of the GUE have the density $q$ of Lemma C.10 with $w(x)=e^{-x^2/2}$ on
$\mathbb R$, so the GUE level density is $K_n(x)e^{-x^2/2}$. We check this for $n=3$ with the moments
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

## C.5 Wishart matrices and the induced measure

Let $a\le b$, $G$ an $a\times b$ Ginibre matrix (independent entries with independent $N(0,\frac12)$ real and imaginary parts,
density $\pi^{-ab}e^{-\operatorname{Tr}GG^*}$) and $W=GG^*$. Let $T_+$ be the set of lower triangular $a\times a$ matrices with positive
diagonal, an open subset of $\mathbb R^{a^2}$ in the coordinates $L_{ii}$, $\operatorname{Re}L_{ij}$, $\operatorname{Im}L_{ij}$ ($i>j$), and
$P_a$ the set of positive definite matrices, open in the Hermitian $a\times a$ matrices.

**Lemma C.11 (Bartlett decomposition).** Almost surely the rows $r_1,\dots,r_a$ of $G$ are linearly independent, and
Gram–Schmidt gives $G=LV$ with $L\in T_+$ and $V$ with orthonormal rows, so $W=LL^*$. The entries $L_{ij}$, $i\ge j$, are
independent, $L_{ij}$ ($i>j$) is standard complex Gaussian and $L_{ii}^2$ has the density $x^{b-i}e^{-x}/(b-i)!$ on $(0,\infty)$.
Hence $L$ has the density $c_L\prod_iL_{ii}^{2(b-i)+1}e^{-\operatorname{Tr}LL^*}$ on $T_+$.

*Proof.* Let $\mathcal F_{i-1}=\sigma(r_1,\dots,r_{i-1})$ and assume (induction) that $r_1,\dots,r_{i-1}$ are a.s. independent, with
Gram–Schmidt vectors $v_1,\dots,v_{i-1}$. Let $u_1,\dots,u_b$ be the orthonormal basis obtained by Gram–Schmidt from
$v_1,\dots,v_{i-1},e_1,\dots,e_b$ (dropping dependent vectors); $u_k=v_k$ for $k<i$, and the $u_k$ are Borel functions of
$r_1,\dots,r_{i-1}$. Put $c_k=\langle u_k,r_i\rangle$. Since $r_i$ is independent of $\mathcal F_{i-1}$ and its density $\pi^{-b}e^{-\|r\|^2}$ is
invariant under unitaries (Proposition A.16(a)), conditionally on $\mathcal F_{i-1}$ the vector $c$ is standard complex Gaussian (Fubini). Now
$L_{ij}=c_j$ for $j<i$ and $r_i-\sum_{j<i}c_jv_j=\sum_{k\ge i}c_ku_k$, so $L_{ii}^2=\sum_{k=i}^b|c_k|^2>0$ a.s. The conditional law of
row $i$ of $L$ given $\mathcal F_{i-1}$ is thus a fixed law, so row $i$ is independent of the earlier rows, with the coordinates
$c_1,\dots,c_{i-1}$ and $L_{ii}$ independent; by induction all rows of $L$ are independent. By Theorem A.20(b) with $N=m=b-i+1$,
$L_{ii}^2=\|(c_i,\dots,c_b)\|^2$ has the density $s^{m-1}e^{-s}/(m-1)!$. The density of $L_{ii}=\sqrt s$ is then
$2t^{2(b-i)+1}e^{-t^2}/(b-i)!$, and $\operatorname{Tr}LL^*=\sum_{i\ge j}|L_{ij}|^2$. $\square$

**Lemma C.12 (Cholesky map).** $\chi(L)=LL^*$ is a $C^\infty$ bijection of $T_+$ onto $P_a$ with
$|\det D\chi(L)|=2^a\prod_{i=1}^aL_{ii}^{2(a-i)+1}$.

*Proof.* *Onto:* for $W\in P_a$ apply Lemma 5.2 to $W^{1/2}=QR$; then $W=R^*R$ and $L=R^*\in T_+$. *Injective:* if
$LL^*=MM^*$, then $M^{-1}L=M^*(L^*)^{-1}$ is lower and upper triangular, hence diagonal with positive diagonal, and
unitary, hence $1$. *Jacobian:* for $i\ge j$, $W_{ij}=L_{ij}L_{jj}+\sum_{k<j}L_{ik}\overline{L_{jk}}$. Order the coordinates of $L$ and of
the lower triangle of $W$ column by column, and within column $j$ the diagonal first. Then $W_{ij}$ depends only on
$L_{ij}$ and on earlier coordinates, so $D\chi$ is triangular, with diagonal entries $2L_{jj}$ ($\partial W_{jj}/\partial L_{jj}$) and,
for each $i>j$, $L_{jj}$ twice (real and imaginary parts of $W_{ij}$). The product is $\prod_j2L_{jj}\cdot L_{jj}^{2(a-j)}$. $\square$

**Theorem C.13 (Wishart and Laguerre ensembles).** $W=GG^*$ has the density $c_W\det(W)^{b-a}e^{-\operatorname{Tr}W}$ on $P_a$
(zero outside), and its eigenvalues have the density
$$
Z_{a,b}^{-1}\,\Delta(\lambda)^2\prod_{i=1}^a\lambda_i^{\,b-a}e^{-\lambda_i}\quad\text{on }C_+=\{0<\lambda_1<\dots<\lambda_a\}.
$$

*Proof.* For Borel $g\ge0$, Lemma C.11 and the change of variables $W=\chi(L)$ (Lemma C.12) give
$\mathbb E\,g(W)=\int_{P_a}g(W)\,c_L\prod_iL_{ii}^{2(b-i)+1}e^{-\operatorname{Tr}W}\big/\big(2^a\prod_iL_{ii}^{2(a-i)+1}\big)dW$ with
$L=\chi^{-1}(W)$, and $\prod_iL_{ii}^{2(b-a)}=\det(W)^{b-a}$. This density is invariant under $W\mapsto VWV^*$, $V\in U(a)$, and
$\Lambda\in P_a$ iff $\lambda_1>0$; apply Corollary C.6 with $\beta=2$. $\square$

**Theorem C.14 (induced measure).** Let $\psi$ be a Haar-random pure state of $\mathbb C^a\otimes\mathbb C^b$ (Section 7.5),
$\rho_A=\operatorname{Tr}_2|\psi\rangle\langle\psi|=W/s$ with $s=\operatorname{Tr}W$ (proof of Theorem 7.11), and
$\mu_1<\dots<\mu_a$ the eigenvalues of $\rho_A$. Then $s$ and $\mu$ are independent, $s$ has the density
$s^{ab-1}e^{-s}/(ab-1)!$, and $(\mu_1,\dots,\mu_{a-1})$ has the density
$$
Z'^{-1}_{a,b}\,\Delta(\mu)^2\prod_{i=1}^a\mu_i^{\,b-a}\quad\text{on }D_a=\{0<\mu_1<\dots<\mu_{a-1}<\mu_a\},\qquad\mu_a=1-\sum_{i<a}\mu_i ,
$$
which is the density "$\propto\delta(1-\sum_i\mu_i)\,\Delta(\mu)^2\prod_i\mu_i^{b-a}$" of the physics literature.

*Proof.* $\Theta(s,\mu_1,\dots,\mu_{a-1})=s\,(\mu_1,\dots,\mu_{a-1},1-\sum_{i<a}\mu_i)$ is a $C^\infty$ bijection of $(0,\infty)\times D_a$ onto
$C_+$, with inverse $\lambda\mapsto(\sum\lambda_i,\lambda_1/\sum\lambda_i,\dots)$. Adding all rows of $D\Theta$ to the last one turns it into
$(1,0,\dots,0)$, so $|\det D\Theta|=s^{a-1}$. As $\Delta$ is homogeneous of degree $a(a-1)/2$, the density of $(s,\mu)$ is
$$
Z_{a,b}^{-1}\,s^{a(a-1)+a(b-a)+a-1}e^{-s}\cdot\Delta(\mu)^2\prod_i\mu_i^{b-a}=Z_{a,b}^{-1}\,s^{ab-1}e^{-s}\cdot\Delta(\mu)^2\prod_i\mu_i^{b-a},
$$
a product, which proves independence and identifies both factors after normalisation. $\square$

We check Lemma C.11 on the Cholesky factor computed by `torch.linalg.cholesky` (which is the unique $L\in T_+$ by
Lemma C.12), and Theorem C.14 for $a=2$, where $\mu_1\in(0,\frac12)$ has the density $\propto(1-2t)^2(t(1-t))^{b-2}$:

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

## C.6 Laguerre polynomials and Page's formula

**Lemma C.15 (Laguerre polynomials).** For integers $\alpha,k\ge0$ let
$L_k^{(\alpha)}(x)=\sum_{j=0}^k(-1)^j\binom{k+\alpha}{k-j}\frac{x^j}{j!}$, and $L_{-1}^{(\alpha)}=0$. Then:

(a) $L_k^{(\alpha)}(x)\,x^\alpha e^{-x}=\frac1{k!}\frac{d^k}{dx^k}\big(x^{k+\alpha}e^{-x}\big)$; $L_k^{(\alpha)}$ has degree $k$ and leading coefficient $(-1)^k/k!$.

(b) $\frac{d}{dx}L_k^{(\alpha)}=-L_{k-1}^{(\alpha+1)}$, hence $\frac{d^i}{dx^i}L_k^{(\alpha)}=(-1)^iL_{k-i}^{(\alpha+i)}$ for $i\le k$; and $L_k^{(\alpha)}=L_k^{(\alpha+1)}-L_{k-1}^{(\alpha+1)}$.

(c) $\int_0^\infty x^mL_k^{(\alpha)}x^\alpha e^{-x}dx=0$ for $m<k$, and $\int_0^\infty L_k^{(\alpha)}L_l^{(\alpha)}x^\alpha e^{-x}dx=\delta_{kl}\frac{(k+\alpha)!}{k!}$.

Consequently, by Lemma C.10 and Theorem C.13, the eigenvalues of $W$ have the level density
($\mathbb E\sum_i\varphi(\lambda_i(W))=\int_0^\infty\varphi\,\rho_{a,b}$ for Borel $\varphi\ge0$, and for every $\varphi$ with $\int|\varphi|\rho_{a,b}<\infty$)
$$
\rho_{a,b}(x)=\sum_{k=0}^{a-1}\frac{k!}{(k+\alpha)!}L_k^{(\alpha)}(x)^2\,x^\alpha e^{-x},\qquad\alpha=b-a .
$$

*Proof.* (a) By Leibniz, $\frac{d^k}{dx^k}(x^{k+\alpha}e^{-x})=\sum_j\binom kj(-1)^j\frac{(k+\alpha)!}{(\alpha+j)!}x^{\alpha+j}e^{-x}$ (with $j$ derivatives
on $e^{-x}$), and $\frac1{k!}\binom kj\frac{(k+\alpha)!}{(\alpha+j)!}=\binom{k+\alpha}{k-j}\frac1{j!}$. (b) Differentiating term by term,
$-\sum_{i=0}^{k-1}(-1)^i\binom{k+\alpha}{k-1-i}\frac{x^i}{i!}=-L_{k-1}^{(\alpha+1)}$; the second identity is Pascal's rule
$\binom{k+\alpha+1}{k-j}-\binom{k+\alpha}{k-1-j}=\binom{k+\alpha}{k-j}$. (c) Integrate (a) by parts $k$ times against $x^m$: the
boundary terms vanish, because $\frac{d^r}{dx^r}(x^{k+\alpha}e^{-x})$ has the factor $x^{k+\alpha-r}$, $r<k$, and decays at $\infty$. The
result is $\frac{(-1)^k}{k!}\int\frac{d^kx^m}{dx^k}x^{k+\alpha}e^{-x}dx$, which is $0$ for $m<k$ and $(-1)^k(k+\alpha)!$ for $m=k$. For $l<k$,
$L_l^{(\alpha)}$ is a combination of $x^m$, $m<k$; for $l=k$ only the leading coefficient $(-1)^k/k!$ counts. So
$p_k=\big(k!/(k+\alpha)!\big)^{1/2}L_k^{(\alpha)}$ are orthonormal for $w(x)=x^\alpha e^{-x}$ on $(0,\infty)$, and the unordered
eigenvalues of $W$ have the density $q$ of Lemma C.10 (with $n=a$) by Corollary C.6, so that, as in the proof of Corollary C.8,
$\mathbb E\sum_i\varphi(\lambda_i)=a\int\varphi(x_1)q(x)\,dx=\int\varphi\,K_aw=\int\varphi\,\rho_{a,b}$ for Borel $\varphi\ge0$, and by linearity for
$\varphi$ with $\int|\varphi|\rho_{a,b}<\infty$. $\square$

**Lemma C.16 (logarithmic integrals).** Let $I_m=\int_0^\infty x^m\log x\,e^{-x}dx$ and $H_m=\sum_{i=1}^m\frac1i$. Then
$I_m=m!\,(I_0+H_m)$. For integers $\alpha\ge0$ and $0\le k\le m$,
$$
A^{(\alpha)}_{k,m}:=\int_0^\infty\log x\,L_k^{(\alpha)}(x)L_m^{(\alpha)}(x)\,x^\alpha e^{-x}dx=
\begin{cases}I_{k+\alpha}/k!, & m=k,\\ -\dfrac{(k+\alpha)!}{k!\,(m-k)}, & m>k.\end{cases}
$$
($I_0$ cancels in Theorem C.17, so its value is not needed.)

*Proof.* For $m\ge1$, integration by parts with $(-e^{-x})'=e^{-x}$ gives $I_m=mI_{m-1}+(m-1)!$, i.e.
$I_m/m!=I_{m-1}/(m-1)!+1/m$. For $A_{k,m}$, write $L_m^{(\alpha)}x^\alpha e^{-x}=\frac1{m!}g^{(m)}$ with $g=x^{m+\alpha}e^{-x}$
(Lemma C.15(a)) and $f=\log x\,L_k^{(\alpha)}$, and integrate by parts $m$ times. The boundary terms
$f^{(i)}g^{(m-1-i)}$ vanish at $\infty$ and at $0$, where $f^{(i)}=O(|\log x|+x^{-i})$ and $g^{(m-1-i)}=O(x^{\alpha+1+i})$. So
$A_{k,m}=\frac{(-1)^m}{m!}\int_0^\infty f^{(m)}g\,dx$ (the integrand is $O(x^\alpha(1+|\log x|))$ at $0$) with
$f^{(m)}=\sum_{i\le k}\binom mi(\log x)^{(m-i)}\frac{d^i}{dx^i}L_k^{(\alpha)}$. For $m-i\ge1$, $(\log x)^{(m-i)}=(-1)^{m-i-1}(m-i-1)!\,x^{i-m}$,
and by Lemma C.15(b,c)
$$
\int_0^\infty x^{i-m}\frac{d^iL_k^{(\alpha)}}{dx^i}\,x^{m+\alpha}e^{-x}dx=(-1)^i\int_0^\infty L_{k-i}^{(\alpha+i)}\,x^{\alpha+i}e^{-x}dx=(-1)^i(\alpha+i)!\,\delta_{ik}.
$$
If $m>k$, only $i=k$ contributes: $A_{k,m}=\frac{(-1)^m}{m!}\binom mk(-1)^{m-k-1}(m-k-1)!(-1)^k(\alpha+k)!=-\frac{(k+\alpha)!}{k!(m-k)}$.
If $m=k$, the terms $i<k$ vanish and the term $i=k$ is $\log x\cdot(-1)^k$, giving $A_{k,k}=\frac{(-1)^k}{k!}(-1)^kI_{k+\alpha}$. $\square$

**Theorem C.17 (Page's formula).** For a Haar-random pure state of $\mathbb C^a\otimes\mathbb C^b$ with $a\le b$,
$$
\mathbb E\,S(\rho_A)=\sum_{k=b+1}^{ab}\frac1k-\frac{a-1}{2b}.
$$

*Proof.* *Reduction.* With the notation of Theorem C.14, $\lambda_i=s\mu_i$ and $\sum_i\mu_i=1$, so
$\sum_i\lambda_i\log\lambda_i=s\log s-s\,S(\rho_A)$. All terms are integrable ($|x\log x|\le1+x^2$, and $\mathbb Es^2<\infty$,
$0\le S(\rho_A)\le a/e$ as $0\le-x\log x\le1/e$ on $[0,1]$), $\mathbb E\sum_i\lambda_i\log\lambda_i=\int x\log x\,\rho_{a,b}(x)\,dx$ by Lemma C.15, and $s$ is independent of
$S(\rho_A)$, so $\mathbb E\sum_i\lambda_i\log\lambda_i=\mathbb E[s\log s]-\mathbb E[s]\,\mathbb E\,S(\rho_A)$. From the density of $s$,
$\mathbb Es=ab$ and $\mathbb E[s\log s]=I_{ab}/(ab-1)!=ab\,(I_0+H_{ab})$ (Lemma C.16). Hence
$$
\mathbb E\,S(\rho_A)=I_0+H_{ab}-\frac1{ab}\int_0^\infty x\log x\,\rho_{a,b}(x)\,dx .
$$
*The integral.* By Lemma C.15(b), $L_k^{(\alpha)}=L_k^{(\alpha+1)}-L_{k-1}^{(\alpha+1)}$, so with $x\cdot x^\alpha=x^{\alpha+1}$ and Lemma C.16
for $\alpha+1$,
$$
\frac{k!}{(k+\alpha)!}\int_0^\infty x\log x\,L_k^{(\alpha)}(x)^2x^\alpha e^{-x}dx=\frac{k!}{(k+\alpha)!}\Big(A^{(\alpha+1)}_{k,k}-2A^{(\alpha+1)}_{k-1,k}+A^{(\alpha+1)}_{k-1,k-1}\Big)
$$
$$
=(k+\alpha+1)(I_0+H_{k+\alpha+1})+2k+k\,(I_0+H_{k+\alpha}),
$$
also for $k=0$. Since $\sum_{k=0}^{a-1}(2k+\alpha+1)=a(a-1)+a(b-a+1)=ab$,
$$
\int_0^\infty x\log x\,\rho_{a,b}\,dx=ab\,I_0+\Sigma,\qquad \Sigma=\sum_{k=0}^{a-1}\big[(k+\alpha+1)H_{k+\alpha+1}+kH_{k+\alpha}+2k\big].
$$
*The sum.* Put $m=\alpha+1=b-a+1$ and $j=k+m\in\{m,\dots,b\}$. Since $(j-m)H_{j-1}=(j-m)H_j-1+\frac mj$,
$\Sigma=\sum_{j=m}^b(2j-m)H_j-a+m(H_b-H_{m-1})+a(a-1)$. With $F(j)=j(j+1-m)$ we have $F(j)-F(j-1)=2j-m$, $F(m-1)=0$,
$F(b)=ab$, and summation by parts gives
$$
\sum_{j=m}^b(2j-m)H_j=abH_b-\sum_{j=m}^{b-1}\frac{F(j)}{j+1}=abH_b-\sum_{j=m}^{b-1}\Big(j-m+\frac m{j+1}\Big)=abH_b-\frac{(a-1)(a-2)}2-m(H_b-H_m).
$$
As $m(H_m-H_{m-1})=1$, $\Sigma=abH_b-\frac{(a-1)(a-2)}2+1-a+a(a-1)=abH_b+\frac{a(a-1)}2$. Therefore
$\mathbb E\,S(\rho_A)=H_{ab}-H_b-\frac{a-1}{2b}$. $\square$

**Corollary C.18.** For $a\le b$,
$\log a-\frac{a-1}{ab+1}-\frac{a-1}{2b}\le\mathbb E\,S(\rho_A)\le\log a-\frac{a-1}{2b}$; in particular
$\big|\mathbb E\,S(\rho_A)-\log a+\frac a{2b}\big|\le\frac1{2b}$. For two qubits $\mathbb E\,S(\rho_A)=\frac13$.

*Proof.* Since $1/x$ decreases, $\int_k^{k+1}\frac{dx}x\le\frac1k\le\int_{k-1}^k\frac{dx}x$; summing over $b<k\le ab$,
$\log\frac{ab+1}{b+1}\le H_{ab}-H_b\le\log a$, and $\log\frac{ab+1}{b+1}=\log a-\log\big(1+\frac{a-1}{ab+1}\big)\ge\log a-\frac{a-1}{ab+1}$.
Finally $\frac{a-1}{ab+1}\le\frac1b$, and $a=b=2$ gives $\frac13+\frac14-\frac14$. $\square$

The formula was conjectured by D. N. Page (Phys. Rev. Lett. 71, 1291, 1993) and proved by S. K. Foong and S. Kanno
(Phys. Rev. Lett. 72, 1148, 1994), J. Sánchez-Ruiz (Phys. Rev. E 52, 5653, 1995) and S. Sen (Phys. Rev. Lett. 77, 1, 1996);
the computation above follows Sen's reduction to the Laguerre ensemble but sums the one-point density term by term
instead of using the Christoffel–Darboux formula. The induced measure is that of K. Życzkowski and H.-J. Sommers
(J. Phys. A 34, 7111, 2001); the decomposition of Lemma C.11 is due to M. S. Bartlett (1933), and Theorem C.5 to H. Weyl.

We check the Laguerre level density against a histogram, the values $A^{(\alpha)}_{k,m}$ of Lemma C.16 and the intermediate
result $\int x\log x\,\rho_{a,b}=ab\,I_0+\Sigma$ by quadrature (with $I_0$ itself computed by quadrature), the identity for
$\Sigma$ in exact arithmetic, Page's formula against a direct integration of the
induced density for $a=2$ (exact rational values) and against Monte Carlo, and Corollary C.18:

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

**Back:** [Chapter 5](05_random_matrices.md) uses Theorem C.7, Corollary C.8 and Proposition C.9;
[Chapter 7](07_tensor_products.md) uses Theorems C.14 and C.17 and Corollary C.18.
