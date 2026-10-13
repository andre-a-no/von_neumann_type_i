**English** · [Español](es/D_gksl.md) · [Français](fr/D_gksl.md) · [Deutsch](de/D_gksl.md) · [中文](zh/D_gksl.md) · [日本語](ja/D_gksl.md) · [Русский](ru/D_gksl.md)

# Appendix D. Generators of quantum dynamical semigroups

Chapter 8 proved that a generator of Gorini–Kossakowski–Sudarshan–Lindblad (GKSL) form produces a
quantum dynamical semigroup (Theorem 8.6) and used the converse (Theorem 8.7) and its version for
$\mathcal M=\bigoplus_cM_{k_c}$. This appendix proves both in finite dimensions, together with the
freedom in the GKSL form and a criterion for the uniqueness of the stationary state. The tool is the
Choi matrix of Chapter 6: a generator is not completely positive, but its Choi matrix is positive on
the orthogonal complement of the maximally entangled vector. The same kind of argument — positivity
of a family of operators, then a decomposition — proves in the last section the theorem on operator
algebra error correction stated in Chapter 10 (Theorem 10.7), with an explicit recovery channel.

**You will learn** why the generator of a completely positive semigroup is conditionally completely
positive, how to read off $H$ and the jump operators from the Choi matrix of a generator, in what sense
they are unique, why every generator on an algebra with sectors has the form used in Section 8.3, when
the stationary state is unique, and which algebras of observables a recovery channel can protect.

Throughout, $M_n=M_n(\mathbb C)$, $E_{ij}=|i\rangle\langle j|$ are its matrix units, a
**superoperator** is a linear map on $M_n$ (or on $\mathcal M$), and
$$
\Omega=\sum_{i=1}^ne_i\otimes e_i\in\mathbb C^n\otimes\mathbb C^n,\qquad
J(\Psi)=\sum_{i,j=1}^nE_{ij}\otimes\Psi(E_{ij})=[\Psi(E_{ij})]_{i,j=1}^n
$$
are the maximally entangled vector and the **Choi matrix** of a superoperator $\Psi$ (Theorem 6.2;
here $\Psi$ need not be positive). For vectors $v,w$, $|v\rangle\langle w|$ is the operator
$x\mapsto\langle w,x\rangle v$. A superoperator $\Psi$ is **Hermiticity preserving** if
$\Psi(A^*)=\Psi(A)^*$ for all $A$. Quantum dynamical semigroups and their generators are those of
Section 8.2 (Lemma 8.5): $T_t=e^{t\mathcal L}$ with channels $T_t$.

## D.1 Choi matrices of the building blocks

**Lemma D.1.** Let $\Psi$ be a superoperator on $M_n$.

(a) $\Psi\mapsto J(\Psi)$ is a linear bijection from the superoperators on $M_n$ onto $M_n(M_n)=M_{n^2}$.
(b) Every $v\in\mathbb C^n\otimes\mathbb C^n$ is $v=(1\otimes K)\Omega$ for exactly one $K\in M_n$, and
$\langle(1\otimes K)\Omega,(1\otimes L)\Omega\rangle=\operatorname{Tr}(K^*L)$. In particular
$\langle\Omega,(1\otimes K)\Omega\rangle=\operatorname{Tr}K$.
(c) For $K,L\in M_n$ the superoperator $A\mapsto KAL^*$ has the Choi matrix
$|(1\otimes K)\Omega\rangle\langle(1\otimes L)\Omega|$. In particular $J(\mathrm{id})=|\Omega\rangle\langle\Omega|$.
(d) $\Psi$ is Hermiticity preserving iff $J(\Psi)$ is self-adjoint. Every CP map is Hermiticity preserving.

*Proof.* (a) The $(i,j)$ block of $J(\Psi)$ is $\Psi(E_{ij})$ and the $E_{ij}$ form a basis of $M_n$, so
$J$ is linear and injective, and every block matrix $[X_{ij}]$ is the Choi matrix of the superoperator
with $E_{ij}\mapsto X_{ij}$. (b) $(1\otimes K)\Omega=\sum_ie_i\otimes Ke_i$ has the $i$-th block $Ke_i$, the
$i$-th column of $K$, which determines $K$, and every $v$ arises in this way. Then
$\langle\sum_ie_i\otimes Ke_i,\sum_je_j\otimes Le_j\rangle=\sum_i\langle Ke_i,Le_i\rangle=\operatorname{Tr}(K^*L)$;
for the last claim take $1$ in place of $K$. (c) The $(i,j)$ block of the right-hand side is
$Ke_i(Le_j)^*=KE_{ij}L^*$. (d) The $(i,j)$ block of $J(\Psi)^*$ is $\Psi(E_{ji})^*$. If $\Psi$ is
Hermiticity preserving, $\Psi(E_{ji})^*=\Psi(E_{ji}^*)=\Psi(E_{ij})$. Conversely, if
$\Psi(E_{ji})^*=\Psi(E_{ij})$ for all $i,j$ and $A=\sum a_{ij}E_{ij}$, then
$\Psi(A)^*=\sum\bar a_{ij}\Psi(E_{ij})^*=\sum\bar a_{ij}\Psi(E_{ji})=\Psi(A^*)$. A CP map has a Kraus form
(Theorem 6.2), and $(\sum_lK_lAK_l^*)^*=\sum_lK_lA^*K_l^*$. $\square$

**Lemma D.2 (decomposition).** Let $q_1,\dots,q_C$ be orthonormal vectors in
$\mathbb C^n\otimes\mathbb C^n$, written $q_c=(1\otimes Q_c)\Omega$, let $W=\operatorname{span}\{q_c\}$, and let
$\Psi$ be a Hermiticity preserving superoperator with $\langle v,J(\Psi)v\rangle\ge0$ for all $v\perp W$.
Put $J=J(\Psi)$ and $\xi_c=Jq_c-\frac12\sum_d\langle q_d,Jq_c\rangle q_d$, and let $K_c\in M_n$ be defined by
$(1\otimes K_c)\Omega=\xi_c$. Then there are $L_1,\dots,L_r\in M_n$, $r=\operatorname{rank}(P_W^\perp JP_W^\perp)$
with $P_W^\perp$ the projection onto $W^\perp$, such that
$$
\Psi(A)=\sum_{j=1}^rL_jAL_j^*+\sum_{c=1}^C\big(K_cAQ_c^*+Q_cAK_c^*\big)\qquad(A\in M_n),
$$
$\operatorname{Tr}(Q_c^*L_j)=0$ for all $c,j$, and $\operatorname{Tr}(L_j^*L_k)=0$ for $j\neq k$.

*Proof.* $J$ is self-adjoint by Lemma D.1(d). Let $P=\sum_c|q_c\rangle\langle q_c|$ be the projection onto $W$
and $P^\perp=P_W^\perp=1-P$. Expanding $J=(P+P^\perp)J(P+P^\perp)$ and using $PJP^\perp=PJ-PJP$ and
$P^\perp JP=JP-PJP$,
$$
J=P^\perp JP^\perp+\big(PJ+JP-PJP\big).
$$
The first term is self-adjoint, and $\langle v,P^\perp JP^\perp v\rangle=\langle P^\perp v,JP^\perp v\rangle\ge0$
because $P^\perp v\perp W$; so $P^\perp JP^\perp\ge0$. By Theorem 1.3,
$P^\perp JP^\perp=\sum_{j=1}^r|l_j\rangle\langle l_j|$, where the $l_j$ are the eigenvectors of the nonzero
eigenvalues $\mu_j$ scaled by $\sqrt{\mu_j}$: they are mutually orthogonal and lie in the range of
$P^\perp JP^\perp$, which is orthogonal to $W$. Write $l_j=(1\otimes L_j)\Omega$; by Lemma D.1(b),
$\operatorname{Tr}(Q_c^*L_j)=\langle q_c,l_j\rangle=0$ and $\operatorname{Tr}(L_j^*L_k)=\langle l_j,l_k\rangle=0$ for $j\ne k$,
and by Lemma D.1(c), $P^\perp JP^\perp$ is the Choi matrix of $A\mapsto\sum_jL_jAL_j^*$. For the second term,
$JP=\sum_c|Jq_c\rangle\langle q_c|$ and $PJP=\sum_{c,d}\langle q_d,Jq_c\rangle|q_d\rangle\langle q_c|$, so
$$
\sum_c|\xi_c\rangle\langle q_c|=JP-\tfrac12PJP,\qquad\sum_c|q_c\rangle\langle\xi_c|=PJ-\tfrac12PJP
$$
(the second is the adjoint of the first; $PJP$ is self-adjoint), and their sum is $PJ+JP-PJP$. By
Lemma D.1(c), $|\xi_c\rangle\langle q_c|$ and $|q_c\rangle\langle\xi_c|$ are the Choi matrices of
$A\mapsto K_cAQ_c^*$ and $A\mapsto Q_cAK_c^*$. So both sides of the claimed identity have the Choi matrix $J$,
and Lemma D.1(a) gives the identity. $\square$

## D.2 Conditional complete positivity

A superoperator $\Psi$ on $M_n$ is **conditionally completely positive** (CCP) if
$\langle v,J(\Psi)v\rangle\ge0$ for every $v\perp\Omega$. By Lemma D.1(b) these $v$ are the vectors
$(1\otimes K)\Omega$ with $\operatorname{Tr}K=0$.

**Theorem D.3.** For a superoperator $\mathcal L$ on $M_n$ the following are equivalent:

1. $e^{t\mathcal L}$ is CP for every $t\ge0$;
2. $\mathcal L$ is Hermiticity preserving and conditionally completely positive;
3. there are a CP map $\Phi$ on $M_n$ and $G\in M_n$ with $\mathcal L(A)=\Phi(A)+GA+AG^*$ for all $A$.

In 3 one can take $\Phi(A)=\sum_{j=1}^rL_jAL_j^*$ with $r\le n^2-1$, $\operatorname{Tr}L_j=0$ and
$\operatorname{Tr}(L_j^*L_k)=0$ for $j\neq k$, and $G$ with $\operatorname{Tr}G\in\mathbb R$.

*Proof.* (1)$\Rightarrow$(2). By Lemma 8.1, $(e^{t\mathcal L}-\mathrm{id})/t\to\mathcal L$ as $t\to0$, and $J$ is
linear on a finite-dimensional space, hence continuous, so
$J(e^{t\mathcal L})=J(\mathrm{id})+t\,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)$ with the last Choi matrix tending to
$J(\mathcal L)$. For $v\perp\Omega$, $\langle v,J(\mathrm{id})v\rangle=|\langle\Omega,v\rangle|^2=0$ (Lemma D.1(c)), so
for $t>0$
$$
0\le\frac1t\langle v,J(e^{t\mathcal L})v\rangle=\big\langle v,J\big((e^{t\mathcal L}-\mathrm{id})/t\big)v\big\rangle
\;\xrightarrow[t\to0]{}\;\langle v,J(\mathcal L)v\rangle .
$$
Each $e^{t\mathcal L}$ is CP, hence Hermiticity preserving (Lemma D.1(d)), so
$\mathcal L(A^*)=\lim_t(e^{t\mathcal L}(A)^*-A^*)/t=\mathcal L(A)^*$.

(2)$\Rightarrow$(3). Apply Lemma D.2 with $C=1$ and $q_1=\Omega/\sqrt n$, i.e. $Q_1=1/\sqrt n$:
$\mathcal L(A)=\sum_jL_jAL_j^*+GA+AG^*$ with $G=K_1/\sqrt n$, and $\Phi=\sum_jL_j\cdot L_j^*$ is CP
(Theorem 6.2). Further $\operatorname{Tr}L_j=\sqrt n\operatorname{Tr}(Q_1^*L_j)=0$, $r\le\operatorname{rank}P_W^\perp=n^2-1$,
and by Lemma D.1(b)
$\operatorname{Tr}G=\langle\Omega,\xi_1\rangle/\sqrt n=\langle q_1,\xi_1\rangle=\frac12\langle q_1,J(\mathcal L)q_1\rangle\in\mathbb R$,
as $J(\mathcal L)$ is self-adjoint.

(3)$\Rightarrow$(1). Steps 1–4 of the proof of Theorem 8.6 use only that the generator is
$\mathcal L_0+\Phi$ with $\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)$ for some matrix
$H_{\mathrm{eff}}$ and $\Phi$ CP. With $H_{\mathrm{eff}}=iG$, $\mathcal L_0(\rho)=G\rho+\rho G^*$, and these
steps show that $e^{t\mathcal L}=\sum_nS_n(t)$ is CP. $\square$

Positivity alone would not do. Let $\theta$ be the transpose and $\mathcal L=\theta-\mathrm{id}$. As $\theta^2=\mathrm{id}$,
$$
e^{t\mathcal L}=e^{-t}\big(\cosh t\,\mathrm{id}+\sinh t\,\theta\big)=\tfrac{1+e^{-2t}}2\,\mathrm{id}+\tfrac{1-e^{-2t}}2\,\theta ,
$$
a convex combination of positive trace preserving maps (Proposition 6.1). But $J(\theta)=[E_{ji}]$ maps
$v=e_1\otimes e_2-e_2\otimes e_1$, which is orthogonal to $\Omega$, to $-v$ (proof of Proposition 6.1), so
$\langle v,J(\mathcal L)v\rangle=-2$: $\mathcal L$ is not CCP, and indeed
$\langle v,J(e^{t\mathcal L})v\rangle=-(1-e^{-2t})<0$ for every $t>0$. A semigroup of positive trace preserving maps
need not consist of channels.

We check Theorem D.3 on a GKSL generator with arbitrary jumps on $M_3$, built by `lindblad_superoperator`
(the matrix $S$ of $\mathcal L$ on row-major $\operatorname{vec}\rho$, so that $J_{(i,a),(j,b)}=S_{(a,b),(i,j)}$ as in
Section 6.1), and on $\theta-\mathrm{id}$:

```python
import math
import numpy as np
import scipy.linalg as sla
import torch
from torch_vn_algebra import TypeIAlgebra, channels, dynamics

torch.manual_seed(0)
rng = np.random.default_rng(0)
f64, c128 = torch.float64, torch.complex128

def choi_of(S, n):
    """J_{(i,a),(j,b)} = S_{(a,b),(i,j)}: the Choi matrix of the map with row-major matrix S (Section 6.1)."""
    return S.reshape(n, n, n, n).transpose(2, 0, 3, 1).reshape(n * n, n * n)

def superop(f, n):
    """Row-major matrix of a linear map f on n x n matrices."""
    E = np.eye(n * n).reshape(n * n, n, n)
    return np.stack([f(X).reshape(-1) for X in E], axis=1)

vec_of = lambda K: K.T.reshape(-1)                        # (1 (x) K) Omega = sum_i e_i (x) K e_i
op_of = lambda v, n: v.reshape(n, n).T                    # its inverse (Lemma D.1(b))
cplx = lambda *s: rng.standard_normal(s) + 1j * rng.standard_normal(s)

n = 3
M3 = TypeIAlgebra([n], [n], precision='double')
H = cplx(n, n)
H = (H + H.conj().T) / 2
Ls = [cplx(n, n) / 2 for _ in range(2)]                   # arbitrary jumps, not traceless
to_op = lambda X: M3.operator(torch.from_numpy(X).reshape(1, 1, n, n))
S = dynamics.lindblad_superoperator(M3, to_op(H), [to_op(L) for L in Ls])[0][0].numpy()
J = choi_of(S, n)
q = vec_of(np.eye(n)) / math.sqrt(n)                      # the unit vector Omega / sqrt(n)
Pperp = np.eye(n * n) - np.outer(q, q)
assert np.abs(J - J.conj().T).max() < 1e-12               # Hermiticity preserving (Lemma D.1(d))
print(f"min eigenvalue of J(L): {np.linalg.eigvalsh(J).min():.2f}; "
      f"on the complement of Omega: {np.linalg.eigvalsh(Pperp @ J @ Pperp).min():.1e}")
assert np.linalg.eigvalsh(J).min() < -0.1                 # L itself is not CP ...
assert np.linalg.eigvalsh(Pperp @ J @ Pperp).min() > -1e-12        # ... but CCP (Theorem D.3)
for t in (0.01, 0.1, 1.0):
    assert np.linalg.eigvalsh(choi_of(sla.expm(t * S), n)).min() > -1e-12

L_theta = superop(lambda X: X.T, n) - np.eye(n * n)     # L = transpose - id: Hermiticity preserving, not CCP
anti = vec_of(np.eye(n)[:, [1]] @ np.eye(n)[[0], :] - np.eye(n)[:, [0]] @ np.eye(n)[[1], :])   # e1(x)e2 - e2(x)e1
assert abs(anti @ q) < 1e-15 and abs(anti @ choi_of(L_theta, n) @ anti + 2) < 1e-12
X = cplx(n, n)
for t in (0.01, 0.1, 1.0):
    T_t = sla.expm(t * L_theta)
    b = (1 - math.exp(-2 * t)) / 2
    assert np.abs(T_t - ((1 - b) * np.eye(n * n) + b * superop(lambda Y: Y.T, n))).max() < 1e-12
    assert abs(np.trace((T_t @ X.reshape(-1)).reshape(n, n)) - np.trace(X)) < 1e-12   # trace preserving
    assert np.linalg.eigvalsh(choi_of(T_t, n)).min() < -0.9 * b                     # but not CP
```

## D.3 The GKSL theorem

**Theorem D.4 (Gorini–Kossakowski–Sudarshan; Lindblad; Theorem 8.7).** A superoperator $\mathcal L$ on
$M_n$ is the generator of a quantum dynamical semigroup on $M_n$ iff there are $H=H^*$ and
$L_1,\dots,L_r$ in $M_n$ with
$$
\mathcal L(\rho)=-i[H,\rho]+\sum_{j=1}^r\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
$$
One can choose $r\le n^2-1$, $\operatorname{Tr}H=0$, $\operatorname{Tr}L_j=0$, and $\operatorname{Tr}(L_j^*L_k)=0$ for $j\neq k$.

*Proof.* "If" is Theorem 8.6. Conversely, let $(T_t)$ be a quantum dynamical semigroup with generator
$\mathcal L$, so $T_t=e^{t\mathcal L}$ (Lemma 8.5). As every $T_t$ is CP, Theorem D.3 gives
$\mathcal L(\rho)=\sum_jL_j\rho L_j^*+G\rho+\rho G^*$ with $r\le n^2-1$ traceless, mutually orthogonal $L_j$ and
$\operatorname{Tr}G\in\mathbb R$. Put $H=\frac i2(G-G^*)$ and $R=\frac12(G+G^*)$, both self-adjoint; then
$G=-iH+R$, $G^*=iH+R$, and $G\rho+\rho G^*=-i[H,\rho]+\{R,\rho\}$. Moreover
$\operatorname{Tr}H=\frac i2(\operatorname{Tr}G-\overline{\operatorname{Tr}G})=0$. Since every $T_t$ is TP,
$\operatorname{Tr}\mathcal L(\rho)=\frac d{dt}\operatorname{Tr}T_t(\rho)\big|_{t=0}=0$ (Lemma 8.1), and by cyclicity
$$
0=\operatorname{Tr}\mathcal L(\rho)=\operatorname{Tr}\Big(\rho\,\Big(\sum_jL_j^*L_j+2R\Big)\Big)\qquad\text{for all }\rho .
$$
With $\rho=X=\sum_jL_j^*L_j+2R$, which is self-adjoint, $\operatorname{Tr}(X^*X)=0$, so $X=0$ and
$R=-\frac12\sum_jL_j^*L_j$, which is the GKSL form. $\square$

**Theorem D.5 (freedom in the GKSL form).** Let $(H;L_1,\dots,L_r)$ and $(H';L'_1,\dots,L'_s)$ be GKSL
data ($H,H'$ self-adjoint) with generators $\mathcal L$, $\mathcal L'$ as in Theorem D.4.

(a) $\mathcal L'=\mathcal L$ in each of the following cases:
(i) $L'_j=L_j$ and $H'=H+c1$, $c\in\mathbb R$;
(ii) $H'=H$ and, after padding the shorter family with zeros, $L'_j=\sum_ku_{jk}L_k$ with a unitary $u$;
(iii) $L'_j=L_j+a_j1$ with $a_j\in\mathbb C$, and $H'=H-\frac i2\sum_j(\bar a_jL_j-a_jL_j^*)$, which is self-adjoint.

(b) If $\operatorname{Tr}H=\operatorname{Tr}H'=0$, all $L_j,L'_j$ are traceless and $\mathcal L=\mathcal L'$, then $H=H'$ and,
after padding, $L'_j=\sum_ku_{jk}L_k$ with a unitary $u$.

(c) Hence two GKSL data have the same generator iff one is obtained from the other by finitely many
steps (i)–(iii). In particular the operator $H$ of the normalised form of Theorem D.4 and its
dissipative part $\Phi(\rho)=\sum_jL_j\rho L_j^*$ are determined by $\mathcal L$.

*Proof.* (a) (i) $[c1,\rho]=0$. (ii) By Theorem 6.7, $\Phi'=\sum_jL'_j\cdot L_j'^*$ equals
$\Phi=\sum_jL_j\cdot L_j^*$; the dual does not depend on the Kraus operators (Section 6.2), so
$\sum_jL_j'^*L'_j=\Phi^*(1)=\sum_jL_j^*L_j$. (iii) For one jump $L$ and $a\in\mathbb C$, expanding gives
$$
(L+a)\rho(L+a)^*-\tfrac12\{(L+a)^*(L+a),\rho\}-\Big(L\rho L^*-\tfrac12\{L^*L,\rho\}\Big)=\tfrac12[B,\rho],
\qquad B=\bar aL-aL^* ,
$$
since the terms $\bar aL\rho+a\rho L^*$ from the first product combine with
$-\frac12(aL^*\rho+a\rho L^*+\bar aL\rho+\bar a\rho L)$ and the terms with $|a|^2$ cancel. As $B^*=-B$, the
matrix $\frac i2B$ is self-adjoint, so $H'=H-\frac i2\sum_jB_j$ is self-adjoint, and
$-i[H'-H,\rho]=-\frac12\sum_j[B_j,\rho]$ cancels the change of the dissipator.

(b) Let $\Phi=\sum_jL_j\cdot L_j^*$ and $G=-iH-\frac12\Phi^*(1)$, so that $\mathcal L(\rho)=\Phi(\rho)+G\rho+\rho G^*$.
By Lemma D.1(c),
$J(\mathcal L)=\sum_j|l_j\rangle\langle l_j|+|g\rangle\langle\Omega|+|\Omega\rangle\langle g|$ with $l_j=(1\otimes L_j)\Omega$ and
$g=(1\otimes G)\Omega$, and $l_j\perp\Omega$ because $\operatorname{Tr}L_j=0$ (Lemma D.1(b)). With $P^\perp$ the projection
onto $\Omega^\perp$, therefore $P^\perp J(\mathcal L)P^\perp=\sum_j|l_j\rangle\langle l_j|=J(\Phi)$, and likewise
$P^\perp J(\mathcal L')P^\perp=J(\Phi')$. So $J(\Phi)=J(\Phi')$, $\Phi=\Phi'$ (Lemma D.1(a)), and Theorem 6.7 gives the
unitary $u$. Then $\Phi^*(1)=\Phi'^*(1)$, hence $[H-H',\rho]=0$ for all $\rho$. A matrix $X$ that commutes with all
$E_{ij}$ is a multiple of $1$ (compare the entries of $XE_{ij}$ and $E_{ij}X$), so $H-H'=\lambda1$, and
$\lambda=0$ by taking traces.

(c) By (a), the steps (i)–(iii) do not change the generator, and the inverse of each step is a step of the
same kind ((iii) with $-a_j$, since $\bar a_jL'_j-a_jL_j'^*=\bar a_jL_j-a_jL_j^*$). Any data are brought to
normalised form by (iii) with $a_j=-\operatorname{Tr}L_j/n$ followed by (i) with $c=-\operatorname{Tr}H''/n$, $H''$ the
Hamiltonian after (iii); two normalised data with the same generator are related by (ii), by (b). $\square$

The proof of Theorem D.4 is an algorithm: $\Phi$ from the compression of $J(\mathcal L)$ to $\Omega^\perp$, then
$G$ from $\xi_1$, then $H$. We run it on the generator above, rebuild the generator with the library, and
compare with the normalised form of the original data computed by Theorem D.5(iii) and (i):

```python
def gksl_from_generator(S, n, tol=1e-10):
    """(H, [L_j]) with Tr H = Tr L_j = 0 from the matrix of a GKSL generator: the proof of Theorem D.4."""
    J = choi_of(S, n)
    q = vec_of(np.eye(n)) / math.sqrt(n)
    Pperp = np.eye(n * n) - np.outer(q, q)
    w, U = np.linalg.eigh(Pperp @ J @ Pperp)              # = J(Phi) >= 0
    jumps = [math.sqrt(x) * op_of(U[:, k], n) for k, x in enumerate(w) if x > tol]
    xi = J @ q - 0.5 * (q.conj() @ J @ q) * q
    G = op_of(xi, n) / math.sqrt(n)
    return 0.5j * (G - G.conj().T), jumps

H_rec, L_rec = gksl_from_generator(S, n)
S_rec = dynamics.lindblad_superoperator(M3, to_op(H_rec), [to_op(L) for L in L_rec])[0][0].numpy()
assert np.abs(S_rec - S).max() < 1e-12                    # the same generator
assert abs(np.trace(H_rec)) < 1e-12 and all(abs(np.trace(L)) < 1e-12 for L in L_rec)

a = [-np.trace(L) / n for L in Ls]                        # Theorem D.5(iii): make the jumps traceless
L_norm = [L + aj * np.eye(n) for L, aj in zip(Ls, a)]
H_norm = H - 0.5j * sum(np.conj(aj) * L - aj * L.conj().T for L, aj in zip(Ls, a))
H_norm = H_norm - np.trace(H_norm) / n * np.eye(n)        # Theorem D.5(i)
assert np.abs(H_norm - H_norm.conj().T).max() < 1e-12
assert np.abs(H_rec - H_norm).max() < 1e-12               # H is unique (Theorem D.5(b))
jump_choi = lambda Ks: sum(np.outer(vec_of(K), vec_of(K).conj()) for K in Ks)
assert np.abs(jump_choi(L_rec) - jump_choi(L_norm)).max() < 1e-12      # so is sum_j L_j . L_j^*
Wr, Wn = np.stack([vec_of(K) for K in L_rec], 1), np.stack([vec_of(K) for K in L_norm], 1)
u = np.linalg.lstsq(Wr, Wn, rcond=None)[0]                # L_norm_j = sum_k u_kj L_rec_k
assert np.abs(Wr @ u - Wn).max() < 1e-12 and np.abs(u.conj().T @ u - np.eye(2)).max() < 1e-12   # unitary
print("eigenvalues of the normalised H:", np.round(np.linalg.eigvalsh(H_rec), 4), " jumps:", len(L_rec))
```

**Theorem D.6 (uniqueness of the stationary state).** Let $\mathcal L$ be the GKSL generator of Theorem D.4
with data $(H;L_1,\dots,L_r)$. If $\mathcal L$ has a faithful stationary state $\sigma$ ($\sigma>0$) and the only
matrices commuting with $H$ and all $L_j$, $L_j^*$ are the multiples of $1$, then $\ker\mathcal L=\mathbb C\sigma$:
$\sigma$ is the only stationary state. Moreover $\ker\mathcal L^*=\mathbb C1$.

*Proof.* The dual for the trace pairing is
$\mathcal L^*(A)=i[H,A]+\sum_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})$ (Chapter 8, Exercise 1, with its solution).
A direct expansion, in which the Hamiltonian terms cancel because $A\mapsto i[H,A]$ is a derivation, gives
for every $X$
$$
\mathcal L^*(X^*X)-X^*\mathcal L^*(X)-\mathcal L^*(X^*)X=\sum_j[L_j,X]^*[L_j,X]:
$$
for one jump, the left side is $L^*X^*XL-X^*L^*XL-L^*X^*LX+X^*L^*LX$ (the anticommutator terms cancel in
pairs), which is $(XL-LX)^*(XL-LX)$. Let $X\in\ker\mathcal L^*$. Then $\mathcal L^*(X^*)=\mathcal L^*(X)^*=0$, as $\mathcal L^*$
is Hermiticity preserving, so $\mathcal L^*(X^*X)=\sum_j[L_j,X]^*[L_j,X]=:Y\ge0$, and
$\operatorname{Tr}(\sigma Y)=\operatorname{Tr}(\mathcal L(\sigma)X^*X)=0$. Since $\sigma>0$,
$\operatorname{Tr}(\sigma Y)=\operatorname{Tr}(\sigma^{1/2}Y\sigma^{1/2})$ with $\sigma^{1/2}Y\sigma^{1/2}\ge0$ vanishes only if
$\sigma^{1/2}Y\sigma^{1/2}=0$, i.e. $Y=0$; so $[L_j,X]=0$ for all $j$. The same argument for $X^*\in\ker\mathcal L^*$
gives $[L_j,X^*]=0$, i.e. $[L_j^*,X]=0$. Then $X$ commutes with $L_j^*L_j$, the dissipative part of $\mathcal L^*(X)$
vanishes, and $0=\mathcal L^*(X)=i[H,X]$. By hypothesis $X\in\mathbb C1$; conversely $\mathcal L^*(1)=0$. So
$\ker\mathcal L^*=\mathbb C1$. In the bases $(E_{ij})$ and $(E_{ji})$ of $M_n$, which are dual for the pairing
$(\rho,A)\mapsto\operatorname{Tr}(\rho A)$, the matrix of $\mathcal L^*$ is the transpose of that of $\mathcal L$; so both have the
same rank, and $\dim\ker\mathcal L=\dim\ker\mathcal L^*=1$. A state in $\mathbb C\sigma$ with trace $1$ is $\sigma$. $\square$

This is the finite-dimensional case of a criterion of Frigerio [*Commun. Math. Phys.* 63, 269 (1978)],
quoted in Section 8.5. Both hypotheses are needed: in Exercise 4 of Chapter 8 the state $\frac12$ is faithful
and stationary, but $H$ and $\sigma_z$ commute with all diagonal matrices, and the stationary states form a
two-dimensional family. For the random generator above the code checks the hypotheses and the conclusion
(the matrix of $\mathcal L^*$ on row-major vectors of $A^T$ is $S^T$):

```python
def kernel_dim(S, tol=1e-9):
    return int((np.linalg.svd(S, compute_uv=False) < tol).sum())

def commutant_dim(ops, n):                    # dim {X : XA = AX for all A in ops}
    return kernel_dim(np.concatenate([np.kron(A, np.eye(n)) - np.kron(np.eye(n), A.T) for A in ops]))

w, V = np.linalg.eig(S)
sigma = V[:, np.argmin(abs(w))].reshape(n, n)             # row-major: vec(X)_{(a,b)} = X_ab
sigma = sigma / np.trace(sigma)
assert np.abs(S @ sigma.reshape(-1)).max() < 1e-12        # stationary
assert np.abs(sigma - sigma.conj().T).max() < 1e-10 and np.linalg.eigvalsh(sigma).min() > 1e-3   # faithful
assert commutant_dim([H] + Ls + [L.conj().T for L in Ls], n) == 1
assert kernel_dim(S) == 1 and kernel_dim(S.T) == 1        # Theorem D.6
assert np.abs(S.T @ np.eye(n).reshape(-1)).max() < 1e-12  # ker L^* = C 1
```

## D.4 Algebras with sectors

Let $\mathcal M=\bigoplus_{c=1}^CM_{k_c}\subseteq M_D$ be the block-diagonal matrices, $D=\sum_ck_c$, with
minimal central projections $P_c$, and let $E(X)=\sum_cP_cXP_c$ be the conditional expectation of $M_D$ onto
$\mathcal M$. As in the proof of Proposition 8.8, $E(B)=B$, $E(XB)=E(X)B$ and $E(BX)=BE(X)$ for $B\in\mathcal M$,
and $E(X^*)=E(X)^*$. For a CP map $\Phi:\mathcal M\to\mathcal M$, $\Phi^*$ is its blunt dual (Proposition 6.8).

**Theorem D.7.** Let $\mathcal L$ be the generator of a quantum dynamical semigroup $(T_t)$ on $\mathcal M$. There
are a CP map $\Phi:\mathcal M\to\mathcal M$ and $H=H^*\in\mathcal M$ with
$$
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]\qquad(\rho\in\mathcal M).
$$
By Proposition 6.4, $\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l$ with Kraus operators between the sectors;
conversely every such $\mathcal L$ generates a quantum dynamical semigroup on $\mathcal M$ (Proposition 8.8(1)).
So the generators of Section 8.3 (where $\Phi=\sum_jg_jL_j\cdot L_j^*+\sum_mg_m\Phi_m$) are exactly the
generators of quantum dynamical semigroups on $\mathcal M$.

*Proof.* *Step 1: a CP map on $M_D$.* By Proposition 6.4,
$T_t(\rho)_d=\sum_{c,l}K^{dc}_l(t)\rho_cK^{dc}_l(t)^*$. Let $N^{dc}_l(t)\in M_D$ have the block $K^{dc}_l(t)$ in
position $(d,c)$ and zeros elsewhere. For $X\in M_D$, $N^{dc}_lXN^{dc*}_l$ is the block
$K^{dc}_lX_{cc}K^{dc*}_l$ in position $(d,d)$, where $X_{cc}$ is the $c$-th diagonal block, i.e. the $c$-th block
of $E(X)$. Hence $\hat T_t=T_t\circ E=\sum_{d,c,l}N^{dc}_l\cdot N^{dc*}_l$ is CP on $M_D$ (Theorem 6.2), with
$\hat T_0=E$. By Lemma 8.5, $(\hat T_t-E)/t=\big((T_t-\mathrm{id})/t\big)\circ E\to\mathcal L\circ E$.

*Step 2: positivity off the sector vectors.* Let $p_c=(1\otimes P_c)\Omega\in\mathbb C^D\otimes\mathbb C^D$ and
$W=\operatorname{span}\{p_c\}$. By Lemma D.1(c), $J(E)=\sum_c|p_c\rangle\langle p_c|$, so for $v\perp W$ and $t>0$,
$0\le\frac1t\langle v,J(\hat T_t)v\rangle=\langle v,J((\hat T_t-E)/t)v\rangle\to\langle v,J(\mathcal L\circ E)v\rangle$.
Every $T_t$ is Hermiticity preserving (its Kraus form), hence so is $\mathcal L$ (as in Theorem D.3), and so are $E$
and $\mathcal L\circ E$.

*Step 3: decomposition.* By Lemma D.1(b), $\langle p_c,p_d\rangle=\operatorname{Tr}(P_cP_d)=\delta_{cd}k_c$, so
$q_c=p_c/\sqrt{k_c}=(1\otimes Q_c)\Omega$ with $Q_c=P_c/\sqrt{k_c}$ are orthonormal. Lemma D.2 for
$\Psi=\mathcal L\circ E$ gives a CP map $\hat\Phi=\sum_j\hat L_j\cdot\hat L_j^*$ on $M_D$ and $K_c\in M_D$ with
$$
\mathcal L(E(A))=\hat\Phi(A)+\sum_c\big(K_cAP_c+P_cAK_c^*\big)/\sqrt{k_c}\qquad(A\in M_D).
$$
For $A\in\mathcal M$, $E(A)=A$ and $AP_c=P_cA$, so the sum is $\tilde GA+A\tilde G^*$ with
$\tilde G=\sum_cK_cP_c/\sqrt{k_c}$.

*Step 4: compression.* For $A\in\mathcal M$, $\mathcal L(A)\in\mathcal M$, so
$\mathcal L(A)=E(\mathcal L(A))=E(\hat\Phi(A))+GA+AG^*$ with $G=E(\tilde G)\in\mathcal M$. The map
$\Phi=E\circ\hat\Phi$ on $\mathcal M$ is
$\Phi(A)_d=\sum_{c,j}(P_d\hat L_jP_c)A_c(P_d\hat L_jP_c)^*$: it is CP by Proposition 6.4, with Kraus operators the
blocks of the $\hat L_j$.

*Step 5: trace preservation.* As in Theorem D.4, $\operatorname{Tr}\mathcal L(A)=0$ for $A\in\mathcal M$, and
$\operatorname{Tr}\Phi(A)=\operatorname{Tr}(A\Phi^*(1))$, so $\operatorname{Tr}(AX)=0$ for all $A\in\mathcal M$, where
$X=\Phi^*(1)+G+G^*\in\mathcal M$ is self-adjoint. With $A=X$, $X=0$ (the blunt trace is faithful, Theorem 3.1).
With $H=\frac i2(G-G^*)\in\mathcal M$, $G=-iH-\frac12\Phi^*(1)$ and $G\rho+\rho G^*=-i[H,\rho]-\frac12\{\Phi^*(1),\rho\}$. $\square$

The vectors $p_c$ replace $\Omega=\sum_cp_c$ because $\hat T_0=E$, not $\mathrm{id}$: $J(\mathcal L\circ E)$ need not be
positive on $\Omega^\perp$, and in the example below it is not.

**Corollary D.8 (Heisenberg picture).** Let $(S_t)_{t\ge0}$ be a continuous family of unital CP maps on
$\mathcal M$ with $S_0=\mathrm{id}$ and $S_{s+t}=S_sS_t$. Its generator is $A\mapsto\Psi(A)+K^*A+AK$ with
$\Psi$ CP on $\mathcal M$, $K\in\mathcal M$ and $\Psi(1)+K+K^*=0$.

*Proof.* Let $T_t=S_t^*$ (blunt duals). They are CP (Proposition 6.8) and TP (Corollary 6.9(b)),
$T_0=\mathrm{id}$, and $T_{s+t}=(S_tS_s)^*=S_s^*S_t^*=T_sT_t$, because
$\operatorname{Tr}(S_tS_s(x)y)=\operatorname{Tr}(xS_s^*S_t^*(y))$ and $S_tS_s=S_{s+t}$. As $\Psi\mapsto\Psi^*$ is linear on a
finite-dimensional space, $t\mapsto T_t$ is continuous and $(T_t-\mathrm{id})/t\to\mathcal L_S^*$, $\mathcal L_S$ the
generator of $(S_t)$. Theorem D.7 gives $\mathcal L_S^*(\rho)=\Phi(\rho)+G\rho+\rho G^*$ with $G+G^*=-\Phi^*(1)$, and by
cyclicity, $\operatorname{Tr}((G\rho+\rho G^*)A)=\operatorname{Tr}(\rho(G^*A+AG))$, so by Corollary 6.9(c)
$\mathcal L_S(A)=\Phi^*(A)+G^*A+AG$. Take $\Psi=\Phi^*$ (CP by Proposition 6.8) and $K=G$. $\square$

This is the finite-dimensional case of the theorem of Christensen and Evans quoted after Theorem 8.7. We
check Theorem D.7 on $M_1\oplus M_2\subseteq M_3$ for a generator of `lindblad_rhs` with a jump inside the
sectors and a random jump between them, following the five steps of the proof:

```python
alg = TypeIAlgebra([1, 2], [1, 2], precision='double')    # M = M_1 (+) M_2 in M_3
D, ks = 3, alg.k_factors
P = [np.diag([1.0, 0, 0]), np.diag([0, 1.0, 1.0])]
E = lambda X: sum(Pc @ X @ Pc for Pc in P)               # the conditional expectation onto M
def blocks(X):                                            # dense block-diagonal 3 x 3 -> library layout
    out = torch.zeros(1, 2, 2, 2, dtype=c128)
    out[0, 0, :1, :1] = torch.from_numpy(X[:1, :1])
    out[0, 1] = torch.from_numpy(X[1:, 1:])
    return out
dense = lambda Y: sla.block_diag(Y[0, 0, :1, :1].numpy(), Y[0, 1].numpy())
Hs = alg.operator_from_eigenvalues(lambda d: torch.randn(d, dtype=f64), force_self_adjoint=True)
inside = alg.operator(torch.randn(1, 2, 2, 2, dtype=c128))       # a jump inside the sectors
between = channels.random_inter_sector_channel(alg, alg, kraus_rank=2)
rhs = dynamics.lindblad_rhs(Hs, [inside, between], rates=[0.7, 1.3], algebra=alg)
L_M = lambda X: dense(rhs(0.0, blocks(X)))               # the generator on M (X block diagonal)
J = choi_of(superop(lambda X: L_M(E(X)), D), D)          # J(L o E) on M_3
p = [vec_of(Pc) for Pc in P]
Q = np.stack([pc / np.linalg.norm(pc) for pc in p], 1)    # orthonormal q_c = p_c / sqrt(k_c)
Pw = np.eye(D * D) - Q @ Q.T                              # projection onto the complement of W
assert np.linalg.eigvalsh(Pw @ J @ Pw).min() > -1e-12     # Step 2
w, U = np.linalg.eigh(Pw @ J @ Pw)                        # Step 3 (Lemma D.2)
jumps = [math.sqrt(x) * op_of(U[:, k], D) for k, x in enumerate(w) if x > 1e-10]
xis = [J @ Q[:, c] - 0.5 * sum((Q[:, d] @ J @ Q[:, c]) * Q[:, d] for d in range(2)) for c in range(2)]
G = E(sum(op_of(xi, D) @ Pc / math.sqrt(k) for xi, Pc, k in zip(xis, P, ks)))   # Step 4
Phi = lambda X: E(sum(K @ X @ K.conj().T for K in jumps))
Phi_dual_one = E(sum(K.conj().T @ K for K in jumps))      # Phi^*(1)
H_M = 0.5j * (G - G.conj().T)
assert np.abs(G + G.conj().T + Phi_dual_one).max() < 1e-12      # Step 5
for _ in range(3):
    X = E(cplx(D, D))
    gksl = Phi(X) - 0.5 * (Phi_dual_one @ X + X @ Phi_dual_one) - 1j * (H_M @ X - X @ H_M)
    assert np.abs(L_M(X) - gksl).max() < 1e-12            # Theorem D.7
    dual = lambda A: E(sum(K.conj().T @ A @ K for K in jumps)) + G.conj().T @ A + A @ G   # Corollary D.8
    assert abs(np.trace(L_M(X) @ (Y := E(cplx(D, D)))) - np.trace(X @ dual(Y))) < 1e-12
Pomega = np.eye(D * D) - np.outer(sum(p), sum(p)) / D     # projection onto the complement of Omega
print(f"numerical observation: min eigenvalue of J(L o E) on the complement of Omega: "
      f"{np.linalg.eigvalsh(Pomega @ J @ Pomega).min():.2f}")
assert np.linalg.eigvalsh(Pomega @ J @ Pomega).min() < -0.1
```

## D.5 Operator algebra error correction

We use the setting of Section 10.3. $\Phi(\rho)=\sum_{i=1}^rE_i\rho E_i^*$ is a channel on $M_D$, $\mathcal C=P\mathbb C^D$
a code, $B(\mathcal C)=\{X\in M_D:X=PXP\}$ the operators on $\mathcal C$, $\mathfrak A\subseteq B(\mathcal C)$ a $*$-algebra
with unit $P$ (a von Neumann algebra on $\mathcal C$ by Theorem 2.2), and $\mathfrak A'$ its commutant in
$B(\mathcal C)$. $\mathfrak A$ is **correctable** for $\Phi$ if there is a channel $\mathcal R$ on $M_D$ with
$P(\mathcal R\circ\Phi)^*(A)P=A$ for all $A\in\mathfrak A$. If $\mathcal R$ has the Kraus operators $R_m$, then
$\mathcal R\circ\Phi$ has the Kraus operators $R_mE_i$, and $\sum_mR_m^*R_m=1$.

**Theorem D.9 (Bény–Kempf–Kribs; Theorem 10.7).** $\mathfrak A$ is correctable for $\Phi$ iff
$PE_i^*E_jP\in\mathfrak A'$ for all $i,j$. In that case the channel $\mathcal R$ constructed in the proof corrects $\mathfrak A$.

*Proof.* *Necessity.* Let $\mathcal R$ correct $\mathfrak A$ and put $W_{mi}=R_mE_iP$, so that
$\sum_{m,i}W_{mi}^*XW_{mi}=P(\mathcal R\circ\Phi)^*(X)P$ for all $X\in M_D$. This equals $X$ for $X\in\mathfrak A$, and $P$ for
$X=1$, because $(\mathcal R\circ\Phi)^*$ is unital (Corollary 6.9(b)). Let $p\in\mathfrak A$ be a projection; then $p\le P$,
and $\sum W_{mi}^*(1-p)W_{mi}=P-p$. For $x=px$,
$\sum_{m,i}\|(1-p)W_{mi}x\|^2=\langle x,(P-p)x\rangle=0$, and for $x=(P-p)x$,
$\sum_{m,i}\|pW_{mi}x\|^2=\langle x,px\rangle=0$. Hence $(1-p)W_{mi}p=0$ and $pW_{mi}(P-p)=0$, and since
$W_{mi}=W_{mi}P$,
$$
W_{mi}p=pW_{mi}p=pW_{mi}P=pW_{mi}.
$$
Every self-adjoint $A\in\mathfrak A$ is a linear combination of projections in $\mathfrak A$: by Theorem 1.3 on $\mathcal C$,
$A=\sum_k\lambda_kp_k$ with distinct $\lambda_k$ and orthogonal projections $p_k$, $\sum_kp_k=P$, and
$p_k=\prod_{l\neq k}(A-\lambda_lP)/(\lambda_k-\lambda_l)\in\mathfrak A$ (apply both sides to the eigenvectors in $\mathcal C$ and to $\mathcal C^\perp$; the empty product is $P$). Every element
of $\mathfrak A$ is $A_1+iA_2$ with $A_1,A_2\in\mathfrak A$ self-adjoint. So $W_{mi}A=AW_{mi}$ for all $A\in\mathfrak A$, and,
taking adjoints, $W_{mi}^*A=AW_{mi}^*$. As $\sum_mW_{mi}^*W_{mj}=PE_i^*(\sum_mR_m^*R_m)E_jP=PE_i^*E_jP$,
$$
A\,PE_i^*E_jP=\sum_mW_{mi}^*AW_{mj}=PE_i^*E_jP\,A\qquad(A\in\mathfrak A),
$$
i.e. $PE_i^*E_jP\in\mathfrak A'$.

*Sufficiency.* By Theorem 2.6 for $\mathfrak A$ on $\mathcal C$ there is an orthonormal basis $|c;a,b\rangle$
($c=1,\dots,C$, $a=1,\dots,k_c$, $b=1,\dots,m_c$) of $\mathcal C$ — the images of $e_a\otimes f_b$ in the $c$-th summand under
the unitary $V$ — in which $\mathfrak A$ consists of the operators $\sum_{c,a,a',b}(A_c)_{aa'}|c;a,b\rangle\langle c;a',b|$,
$A_c\in M_{k_c}$, and $\mathfrak A'$ of the operators $\sum_{c,a,b,b'}(Y_c)_{bb'}|c;a,b\rangle\langle c;a,b'|$, $Y_c\in M_{m_c}$.
The hypothesis therefore gives matrices $\Lambda^c_{ij}\in M_{m_c}$ with
$$
\langle c;a,b|E_i^*E_j|d;a',b'\rangle=\delta_{cd}\,\delta_{aa'}\,(\Lambda^c_{ij})_{bb'} .
$$
*Step 1.* For each $c$ and each pair $\alpha=(i,b)$ let $F^c_\alpha:\mathbb C^{k_c}\to\mathbb C^D$,
$F^c_\alpha e_a=E_i|c;a,b\rangle$, and $\Lambda^c_{\alpha\beta}=(\Lambda^c_{ij})_{bb'}$ for $\beta=(j,b')$. The display says
$F^{c*}_\alpha F^d_\beta=\delta_{cd}\Lambda^c_{\alpha\beta}1_{k_c}$: a Knill–Laflamme condition for the "errors" $F^c_\alpha$ on the
code $\mathbb C^{k_c}$, with in addition orthogonality between the sectors.

*Step 2.* As in the proof of Theorem 10.5, $\Lambda^c$ is self-adjoint and positive, since
$\sum\bar z_\alpha z_\beta\Lambda^c_{\alpha\beta}1=(\sum z_\alpha F^c_\alpha)^*(\sum z_\beta F^c_\beta)\ge0$. Write
$\Lambda^c=U^c\operatorname{diag}(d^c_\kappa)U^{c*}$ with $U^c$ unitary and $d^c_\kappa\ge0$, and
$G^c_\kappa=\sum_\alpha u^c_{\alpha\kappa}F^c_\alpha$. Then $G^{c*}_\kappa G^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}d^c_\kappa1$, so
$G^c_\kappa=0$ if $d^c_\kappa=0$, and otherwise $V^c_\kappa=G^c_\kappa/\sqrt{d^c_\kappa}$ is an isometry with
$V^{c*}_\kappa V^d_\lambda=\delta_{cd}\delta_{\kappa\lambda}1$. Inverting, $F^c_\alpha=\sum_\kappa\bar u^c_{\alpha\kappa}G^c_\kappa
=\sum_{\kappa:d^c_\kappa>0}\bar u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,V^c_\kappa$, hence
$$
V^{c'*}_\kappa F^d_\beta=\delta_{c'd}\,\bar u^d_{\beta\kappa}\sqrt{d^d_\kappa}\,1,\qquad
F^{c*}_\alpha V^{c'}_\kappa=\delta_{cc'}\,u^c_{\alpha\kappa}\sqrt{d^c_\kappa}\,1 .
$$
*Step 3: the recovery.* The $V^c_\kappa V^{c*}_\kappa$ ($d^c_\kappa>0$) are mutually orthogonal projections; let $Q$ be their
sum. Let $\iota_c:\mathbb C^{k_c}\to\mathbb C^D$, $\iota_ce_a=|c;a,1\rangle$, an isometry, and
$$
\mathcal R(\sigma)=\sum_{c}\sum_{\kappa:d^c_\kappa>0}\iota_cV^{c*}_\kappa\,\sigma\,V^c_\kappa\iota_c^*+(1-Q)\sigma(1-Q).
$$
It is a channel: $\sum_{c,\kappa}V^c_\kappa\iota_c^*\iota_cV^{c*}_\kappa+(1-Q)^2=Q+1-Q=1$. It measures the sector $c$ and
the syndrome $\kappa$, undoes the error on $\mathbb C^{k_c}$ and writes the result into the copy $b=1$.

*Step 4: correctness.* $\mathcal R^*(X)=\sum_{c,\kappa}V^c_\kappa\iota_c^*X\iota_cV^{c*}_\kappa+(1-Q)X(1-Q)$, and for
$A\in\mathfrak A$ with blocks $A_c$, $\iota_c^*A\iota_c=A_c$. Since every $F^d_\beta$ is a combination of the $V^d_\kappa$,
$(1-Q)F^d_\beta=0$. With Step 2, for $\alpha=(i,b)$, $\beta=(j,b')$,
$$
\langle c;a,b|E_i^*\mathcal R^*(A)E_j|d;a',b'\rangle=\sum_{c',\kappa}\langle e_a,F^{c*}_\alpha V^{c'}_\kappa A_{c'}V^{c'*}_\kappa F^d_\beta e_{a'}\rangle
=\delta_{cd}\sum_\kappa u^c_{\alpha\kappa}d^c_\kappa\bar u^c_{\beta\kappa}(A_c)_{aa'}=\delta_{cd}(\Lambda^c_{ij})_{bb'}(A_c)_{aa'} .
$$
Taking $j=i$ and summing over $i$: by the display with $\sum_iE_i^*E_i=1$, $\sum_i\Lambda^c_{ii}=1_{m_c}$, so
$$
\langle c;a,b|(\mathcal R\circ\Phi)^*(A)|d;a',b'\rangle=\delta_{cd}\,\delta_{bb'}(A_c)_{aa'}=\langle c;a,b|A|d;a',b'\rangle .
$$
Both $P(\mathcal R\circ\Phi)^*(A)P$ and $A$ lie in $B(\mathcal C)$ and have the same matrix in the basis $|c;a,b\rangle$, so they are
equal. $\square$

**Corollary D.10.** (a) For $\mathfrak A=B(\mathcal C)$, $\mathfrak A'=\mathbb CP$, and Theorem D.9 is the Knill–Laflamme
theorem 10.5 ($C=1$, $m_1=1$; Step 3 is the recovery of Theorem 10.5). (b) For $\mathfrak A=\operatorname{span}\{P_1,\dots,P_C\}$ with
orthogonal projections $P_c\ne0$, $\sum_cP_c=P$, it is the equivalence of 2 and 3 in Proposition 10.6. (c) (*Noiseless
subsystems.*) In the notation of Proposition 10.1(a), let $\mathcal C=V(\mathbb C^{k_c}\otimes\mathbb C^{m_c})$ be the $c$-th
summand for the interaction algebra $\mathcal A$ of $\Phi$, and $\mathfrak A=V(1_{k_c}\otimes M_{m_c})V^*$ on $\mathcal C$, an
algebra of the form $M_{m}\otimes1_{k}$ after exchanging the two factors. Then $\mathfrak A$ is corrected by
$\mathcal R=\mathrm{id}$, and $PE_i^*E_jP=V(E_{i,c}^*E_{j,c}\otimes1_{m_c})V^*\in\mathfrak A'$.

*Proof.* (a) $\mathfrak A'=\mathbb CP$ is the case $C=1$, $m_1=1$ of the description of $\mathfrak A'$ above, and
$PE_i^*E_jP\in\mathbb CP$ is the Knill–Laflamme condition. (b) $T\in B(\mathcal C)$ commutes with all $P_c$ iff
$P_cTP_d=0$ for $c\ne d$, which is condition 2 of Proposition 10.6. (c) Every $A\in\mathfrak A$, extended by $0$ to the
other summands, lies in $\mathcal A'=V(\bigoplus_c1_{k_c}\otimes M_{m_c})V^*$ (Theorem 2.6), so $\Phi^*(A)=A$ by Proposition
10.1(a) and $PAP=A$. The formula for $PE_i^*E_jP$ follows from $V^*E_iV=\bigoplus_cE_{i,c}\otimes1_{m_c}$, and
$E_{i,c}^*E_{j,c}\otimes1_{m_c}$ commutes with $1_{k_c}\otimes M_{m_c}$. $\square$

By Theorem 2.6 every $\mathfrak A$ is $*$-isomorphic to some $\bigoplus_cM_{k_c}$. When $C\ge2$ and some $k_c\ge2$, a
correctable $\mathfrak A$ protects the sector label $c$ (classical information) together with the quantum information of
each $M_{k_c}$; such codes are called **hybrid** codes. As an example we take $\mathfrak A=(M_2\otimes1_2)\oplus M_1$ on a random
five-dimensional code in $\mathbb C^8$ and noise that acts on the gauge factors $\mathbb C^{m_c}$ through
$\sum_iB_{c,i}^*B_{c,i}=1_{m_c}$, maps the sectors to orthogonal subspaces, scrambles $\mathcal C^\perp$, and is then mixed by a
random unitary (Theorem 6.7). The Knill–Laflamme condition fails, the condition of Theorem D.9 holds, and the
recovery of the proof restores every $A\in\mathfrak A$ but not a generic operator on $\mathcal C$:

```python
D, dims = 8, [(2, 2), (1, 1)]                             # A = (M_2 (x) 1_2) (+) (M_1 (x) 1_1), dim C = 5
U8 = np.linalg.qr(cplx(D, D))[0]
Vc = np.linalg.qr(cplx(D, 5))[0]                          # orthonormal basis |c; a, b> of the code
W = [Vc[:, :4], Vc[:, 4:]]                                # sector c: columns indexed by (a, b), row-major
Pcode = sum(Wc @ Wc.conj().T for Wc in W)
Y = [U8[:, :4], U8[:, 4:6]]                               # orthogonal images C^{k_c} (x) C^{s_c}, s_c = 2
r = 3
B = [np.linalg.qr(cplx(2 * r, m))[0].reshape(r, 2, m) for _, m in dims]   # sum_i B_i^* B_i = 1_{m_c}
E_code = [sum(Yc @ np.kron(np.eye(k), Bc[i]) @ Wc.conj().T for Yc, Bc, Wc, (k, m) in zip(Y, B, W, dims))
          for i in range(r)]
kraus = E_code + [np.linalg.qr(cplx(D, D))[0] @ (np.eye(D) - Pcode)]
mix = np.linalg.qr(cplx(r + 1, r + 1))[0]                 # Theorem 6.7: an equivalent family
kraus = [sum(mix[i, l] * K for l, K in enumerate(kraus)) for i in range(r + 1)]
assert np.abs(sum(K.conj().T @ K for K in kraus) - np.eye(D)).max() < 1e-12      # a channel

units = []                                                 # matrix units of A
for Wc, (k, m) in zip(W, dims):
    for a1 in range(k):
        for a2 in range(k):
            units.append(Wc @ np.kron(np.eye(k)[:, [a1]] @ np.eye(k)[[a2], :], np.eye(m)) @ Wc.conj().T)
comm = lambda X: max(np.abs(X @ A - A @ X).max() for A in units)
assert max(comm(Pcode @ Ki.conj().T @ Kj @ Pcode) for Ki in kraus for Kj in kraus) < 1e-12   # Theorem D.9
KL_dev = max(np.abs(Pcode @ Ki.conj().T @ Kj @ Pcode - np.trace(Pcode @ Ki.conj().T @ Kj) / 5 * Pcode).max()
             for Ki in kraus for Kj in kraus)
print(f"Knill-Laflamme condition violated by {KL_dev:.2f}")
assert KL_dev > 0.05

R_kraus, Qsum = [], np.zeros((D, D), complex)             # Steps 1-3 of the proof
for Wc, (k, m) in zip(W, dims):
    F = [Ki @ Wc[:, b::m] for Ki in kraus for b in range(m)]   # F_(i,b) e_a = E_i |c; a, b>
    Lam = np.array([[(Fa.conj().T @ Fb)[0, 0] for Fb in F] for Fa in F])
    assert max(np.abs(Fa.conj().T @ Fb - Lam[x, y] * np.eye(k)).max()
               for x, Fa in enumerate(F) for y, Fb in enumerate(F)) < 1e-12      # Step 1
    d, u = np.linalg.eigh(Lam)
    iota = Wc[:, 0::m]                                    # e_a -> |c; a, 1>
    for kk in range(len(d)):
        if d[kk] > 1e-12:
            Vk = sum(u[x, kk] * Fx for x, Fx in enumerate(F)) / math.sqrt(d[kk])
            R_kraus.append(iota @ Vk.conj().T)
            Qsum += Vk @ Vk.conj().T
R_kraus.append(np.eye(D) - Qsum)
assert np.abs(sum(R.conj().T @ R for R in R_kraus) - np.eye(D)).max() < 1e-10    # R is a channel
dual = lambda X: sum(Ki.conj().T @ R.conj().T @ X @ R @ Ki for R in R_kraus for Ki in kraus)
for _ in range(3):
    A = sum(complex(z) * Ue for z, Ue in zip(cplx(len(units)), units))
    assert np.abs(Pcode @ dual(A) @ Pcode - A).max() < 1e-10                       # A is corrected
B_gen = Pcode @ cplx(D, D) @ Pcode
print(f"a generic operator on the code: error {np.abs(Pcode @ dual(B_gen) @ Pcode - B_gen).max():.2f}")
assert np.abs(Pcode @ dual(B_gen) @ Pcode - B_gen).max() > 0.1
assert max(comm(R @ Ki @ Pcode) for R in R_kraus for Ki in kraus) < 1e-10   # proof of necessity: W_mi in A'
```

The results of this appendix are due to V. Gorini, A. Kossakowski and E. C. G. Sudarshan
[*J. Math. Phys.* 17, 821 (1976)] and G. Lindblad [*Commun. Math. Phys.* 48, 119 (1976)] (Theorems D.3–D.5),
E. Christensen and D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)] (Corollary D.8 in general von Neumann
algebras), A. Frigerio (Theorem D.6), and C. Bény, A. Kempf and D. W. Kribs [*Phys. Rev. Lett.* 98, 100502
(2007); *Phys. Rev. A* 76, 042303 (2007)] (Theorem D.9); the proofs above are complete for the finite-dimensional
statements used in the course.
