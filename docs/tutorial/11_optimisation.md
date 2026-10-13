**English** · [Español](es/11_optimisation.md) · [Français](fr/11_optimisation.md) · [Deutsch](de/11_optimisation.md) · [中文](zh/11_optimisation.md) · [日本語](ja/11_optimisation.md) · [Русский](ru/11_optimisation.md)

# 11. Optimisation and trace inequalities

A trace inequality is a statement about the extreme values of a functional of operators. Section 3.7 showed one
that is plausible and false, and a sampling test that misses the violation. This chapter explains why sampling
fails in high dimension and replaces it by optimisation over exactly parametrised constraint sets; where the
extremes are known in closed form, we prove them, so that the optimiser is checked against a theorem.

**You will learn** how the library parametrises unitaries and spectra with a prescribed contrast, and why every
admissible operator is reached; von Neumann's trace inequality and the exact extremes it gives; a trace inequality
with its equality condition; why Haar samples stay far from the extremes; and what a multistart search proves.

Throughout, $\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)$ with $D=\sum_ck_c\ge2$, $\operatorname{Tr}$ is the
blunt trace (Section 3.1) and $\|A\|_1=\sum_c\|A_c\|_1$, the trace norm of $A$ as a block-diagonal matrix in
$M_D(\mathbb C)$. The singular values of $A\in M_k(\mathbb C)$ are $s_1(A)\ge\dots\ge s_k(A)$; $x^\downarrow$ and
$x^\uparrow$ are the decreasing and increasing rearrangements of a real vector $x$. $\Delta$ is the Michelson
contrast (Section 3.6), and for $\delta\in[0,1)$
$$
\ell(\delta)=\frac{1-\delta}{1+\delta},\qquad 1-\ell(\delta)=\frac{2\delta}{1+\delta};
$$
a positive operator with $\lambda_{\max}=1$ has contrast $\ell(\lambda_{\min})$, so it has contrast $\delta$ iff
$\lambda_{\min}=\ell(\delta)$, since $\ell$ is a strictly decreasing involution of $[0,1]$. The examples are the three functionals of the companion paper
(section *From sampling to optimisation*):

| | functional | variables |
|---|---|---|
| Exp. 1 | $z=\lvert\operatorname{Tr}(XUY)\rvert-\operatorname{Tr}(XY)$ | $X,Y\ge0$, $U\in\mathcal M$ unitary |
| Exp. 2 | $z=\operatorname{Tr}\lvert XY\rvert-\operatorname{Tr}(\lvert X\rvert\,\lvert Y\rvert)$ | $X=UX_0$, $Y=VY_0$, $X_0,Y_0\ge0$, $U,V$ unitary |
| Exp. 3 | $z=\operatorname{Tr}(Y\lvert X\rvert Y)-\operatorname{Tr}\lvert YXY\rvert$ | $X=X^*$, $Y\ge0$ |

with $|X|$, $|Y|$ normalised to $\lambda_{\max}=1$ and their contrasts as parameters. The **envelopes** are
$\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}$ and $\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}$, all other
variables being free.

## 11.1 Constraint sets and their parametrisations

**Definition.** For $\delta\in[0,1)$ and $\ell=\ell(\delta)$, the **contrast shell** is
$$
\mathcal P_\delta=\{X\in\mathcal M:X\ge0,\ \lambda_{\max}(X)=1,\ \Delta(X)=\delta\}
=\{X=X^*:\operatorname{spec}X\subseteq[\ell,1],\ 1,\ell\in\operatorname{spec}X\}.
$$
The **slots** of $\mathcal M$ are the $D$ pairs $(c,i)$, $i\le k_c$; for $\lambda\in\mathbb R^D$,
$\operatorname{diag}(\lambda)=\bigoplus_c\operatorname{diag}(\lambda_{c,1},\dots,\lambda_{c,k_c})$.

`optimize.PositiveParam` with prescribed contrast $\delta$ draws, per start, two distinct slots
$s_\pm=(c_\pm,i_\pm)$ and a Haar-random unitary $U_0\in\mathcal M$, and maps $\theta\in\mathbb R^D$ and
$G=\bigoplus_cG_c$ to
$$
X(\theta,G)=U\operatorname{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell,\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
$$
`UnitaryParam` is the map $G\mapsto U_0e^{G-G^*}$ alone. Unlike a sigmoid, $(1+\sin\theta)/2$ reaches $0$ and $1$
at finite $\theta$, which matters because extremal spectra often sit at $\ell$ and $1$. A free contrast in $[a,b]$
is one more parameter, $\delta=a+(b-a)(1+\sin\eta)/2$.

**Lemma 11.1.** For every unitary $U_0\in\mathcal M$, $G\mapsto U_0e^{G-G^*}$ maps onto the unitary group of
$\mathcal M$.

*Proof.* Let $U$ be unitary. By Theorem 1.3 (blockwise) and Corollary 1.4(2),
$U_0^*U=\bigoplus_cW_c\operatorname{diag}(e^{i\varphi_{c,j}})_jW_c^*$ with $W_c$ unitary, $\varphi_{c,j}\in\mathbb R$. Then
$S=\bigoplus_cW_c\operatorname{diag}(i\varphi_{c,j})_jW_c^*$ satisfies $S^*=-S$ and $e^S=U_0^*U$ (the series
$\sum_nS^n/n!$ is computed term by term in the eigenbasis), and $G=S/2$ gives $G-G^*=S$. $\square$

**Proposition 11.2 (the parametrisation is onto).** For fixed $\delta$, $s_\pm$ and $U_0$, the image of
$(\theta,G)\mapsto X(\theta,G)$ is $\{X\in\mathcal P_\delta:1\in\operatorname{spec}X_{c_+},\ \ell\in\operatorname{spec}X_{c_-}\}$,
and $\mathcal P_\delta$ is the union of these images over all pairs of distinct slots.

*Proof.* The image lies in this set: $X$ has eigenvalues $\lambda_s\in[\ell,1]$, with $\lambda_{s_\pm}$ in blocks
$c_\pm$. Conversely, diagonalise $X_c=W_c\operatorname{diag}(\mu_c)W_c^*$ and permute each $\mu_c$ (replace $W_c$ by
$W_c\Pi_c$, $\Pi_c$ a permutation matrix) so that $\mu_{s_+}=1$ and $\mu_{s_-}=\ell$. This is possible: if $c_+=c_-$
and $\delta>0$, then $1\ne\ell$ are two eigenvalues of the same block and occupy different positions; if $\delta=0$,
all eigenvalues are $1$. If $\ell<1$, put $\theta_s=\arcsin\big(2(\mu_s-\ell)/(1-\ell)-1\big)$ for the other slots,
so that $\lambda=\mu$ (if $\ell=1$, any $\theta$). Lemma 11.1 gives $G$ with $U_0e^{G-G^*}=\bigoplus_cW_c$.
*Union:* for $X\in\mathcal P_\delta$ choose blocks $c_\pm$ containing the eigenvalues $1$, $\ell$ and distinct slots
in them; if $c_+=c_-$ and $\delta>0$ the block has two different eigenvalues, so $k_{c_+}\ge2$; if $\delta=0$, $X=1$
and any two slots will do. $\square$

So random slots let a multistart search cover the whole shell, and the constraint holds exactly at every iterate.
`SelfAdjointParam` uses $X=U\operatorname{diag}(\sigma\lambda)U^*$ with signs $\sigma\in\{\pm1\}^D$ drawn per start
(and changed only by the sign flips of Section 11.5), so $|X|=U\operatorname{diag}(\lambda)U^*$, and (by the same proof)
the union over signs and slots is $\{X=X^*:|X|\in\mathcal P_\delta\}$. Only $S=G-G^*$ matters. The map
$S\mapsto U_0e^S$ on skew-adjoint $S$ has differential $S\mapsto U_0S$ at $S=0$, a bijection onto the tangent space
$U_0\cdot\{S^*=-S\}$ of the unitary group, so it is a local diffeomorphism near $S=0$; far from $0$ it is not injective,
and its differential is singular where two eigenvalues of $S$ differ by a nonzero multiple of $2\pi i$. `extremize`
therefore **re-centres** after every round ($U_0\leftarrow U_0e^{G-G^*}$, $G\leftarrow0$). The Cayley chart
$S\mapsto U_0C(-S/2)$ (Exercise 1) and the polar retraction $S\mapsto$ unitary factor of $U_0(1+S)$ (Theorem 1.8) have
the same differential at $0$.

The code checks the constraint, and Proposition 11.2 constructively: the parameters of the proof reproduce a target.

```python
import math
import scipy.linalg
import torch
from torch_vn_algebra import TypeIAlgebra, optimize as opt

torch.manual_seed(0)
alg = TypeIAlgebra([3, 3], [3, 3], precision='double')          # M_3 (+) M_3, D = 6
ell = lambda delta: (1 - delta) / (1 + delta)                     # lambda_min for contrast delta
d = torch.tensor([0.0, 0.25, 0.5, 0.9], dtype=torch.float64)
x = opt.PositiveParam(alg, 4, delta=d).operator()                 # four starts, one contrast each
assert torch.allclose(x.michelson_contrast, d, atol=1e-12)        # the constraint holds exactly
assert torch.allclose(x.lambda_max, torch.ones(4, dtype=torch.float64), atol=1e-12)

P = opt.PositiveParam(alg, 1, delta=torch.tensor([0.5], dtype=torch.float64))
lo = ell(0.5)
mu = lo + (1 - lo) * torch.rand(1, 2, 3, dtype=torch.float64)    # target spectrum in [lo, 1] ...
(cp, ip), (cm, im) = P.spectrum.pin_max[0].tolist(), P.spectrum.pin_min[0].tolist()
mu[0, cp, ip], mu[0, cm, im] = 1.0, lo                            # ... with 1 and lo in the drawn slots
W = alg.random_unitary_operator(1).matrix                         # target eigenvectors
target = W @ torch.diag_embed(mu.to(W.dtype)) @ W.mH
with torch.no_grad():                                             # the parameters of the proof
    P.spectrum.theta.copy_(torch.asin((2 * (mu - lo) / (1 - lo) - 1).clamp(-1, 1)))
    V = P.unitary.base.mH @ W                                     # U0^* W = exp(S), G = S / 2
    S = torch.stack([torch.from_numpy(scipy.linalg.logm(V[0, c].numpy())) for c in range(2)])
    P.unitary.G.copy_(S[None] / 2)
err = (P.operator().matrix - target).abs().max().item()
print(f"slots of 1 and lo: {(cp, ip)}, {(cm, im)};  reconstruction error {err:.1e}")
assert err < 1e-10
```

> **Physical meaning.** The Michelson contrast is the visibility $(I_{\max}-I_{\min})/(I_{\max}+I_{\min})$ of
> interference fringes. For $0\le X\le1$ with $\lambda_{\max}=1$ — a filter, an attenuator, a detector efficiency —
> it fixes the dynamic range $\lambda_{\min}/\lambda_{\max}$, the inverse of a condition number.

## 11.2 Von Neumann's trace inequality

A matrix $(c_{ij})$ with $c_{ij}\ge0$ and all row and column sums $1$ is **doubly stochastic**; for a unitary $G$,
$(|G_{ij}|^2)$ and $(|G_{ji}|^2)$ are doubly stochastic, as rows and columns of $G$ are unit vectors.

**Lemma 11.3 (rearrangement).** Let $p_1\ge\dots\ge p_k$, $q_1\ge\dots\ge q_k$ be real and $(c_{ij})$ doubly
stochastic. Then $\sum_ip_iq_{k+1-i}\le\sum_{i,j}p_iq_jc_{ij}\le\sum_ip_iq_i$.

*Proof.* *Upper bound.* For $m\le k$ let $r_j=\sum_{i\le m}c_{ij}\in[0,1]$, so $\sum_jr_j=m$. Since $q_j\ge q_m$ for
$j\le m$ and $q_j\le q_m$ for $j>m$,
$$
\sum_{j\le m}q_j-\sum_jq_jr_j=\sum_{j\le m}q_j(1-r_j)-\sum_{j>m}q_jr_j\ge q_m\Big(m-\sum_jr_j\Big)=0 .
$$
With $x_i=\sum_jq_jc_{ij}$ this says $\sum_{i\le m}x_i\le\sum_{j\le m}q_j$, with equality for $m=k$. By Abel
summation and $p_m-p_{m+1}\ge0$,
$$
\sum_ip_ix_i=p_k\sum_{i\le k}x_i+\sum_{m<k}(p_m-p_{m+1})\sum_{i\le m}x_i
\le p_k\sum_{j\le k}q_j+\sum_{m<k}(p_m-p_{m+1})\sum_{j\le m}q_j=\sum_ip_iq_i .
$$
*Lower bound.* Apply the upper bound to $p$, the decreasing $q'_j=-q_{k+1-j}$ and the doubly stochastic
$c'_{ij}=c_{i,k+1-j}$: $-\sum_{ij}p_iq_jc_{ij}=\sum_{ij}p_iq'_jc'_{ij}\le\sum_ip_iq'_i=-\sum_ip_iq_{k+1-i}$. $\square$

Every $A\in M_k(\mathbb C)$ has a **singular value decomposition** $A=W\Sigma V^*$, $W,V$ unitary,
$\Sigma=\operatorname{diag}(s_i(A))$: write $A=W'|A|$ (Theorem 1.8), $|A|=V\Sigma V^*$ (Theorem 1.3), $W=W'V$.

**Theorem 11.4 (von Neumann's trace inequality).** For $A,B\in M_k(\mathbb C)$,
$|\operatorname{Tr}(AB)|\le\sum_is_i(A)\,s_i(B)$.

*Proof.* Let $A=W_A\Sigma_AV_A^*$, $B=W_B\Sigma_BV_B^*$, $a_i=s_i(A)$, $b_j=s_j(B)$, and $G=V_A^*W_B$, $H=V_B^*W_A$.
By cyclicity $\operatorname{Tr}(AB)=\operatorname{Tr}(\Sigma_AG\Sigma_BH)=\sum_{i,j}a_iG_{ij}b_jH_{ji}$, so, by
$|g||h|\le\frac12(|g|^2+|h|^2)$ and $a_ib_j\ge0$,
$$
|\operatorname{Tr}(AB)|\le\sum_{i,j}a_ib_j\,d_{ij},\qquad d_{ij}=\tfrac12\big(|G_{ij}|^2+|H_{ji}|^2\big).
$$
$(d_{ij})$ is doubly stochastic, and Lemma 11.3 gives $\sum_{ij}a_ib_jd_{ij}\le\sum_ia_ib_i$. $\square$

**Corollary 11.5.** (a) For self-adjoint $P,Q\in M_k(\mathbb C)$ with eigenvalues $p,q$:
$\sum_ip^\downarrow_iq^\uparrow_i\le\operatorname{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i$.
(b) For $P,Q\ge0$: $\operatorname{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i$.

*Proof.* (a) As in the proof of Theorem 3.5, $\operatorname{Tr}(PQ)=\sum_{ij}p_iq_jc_{ij}$ with the doubly stochastic
$c_{ij}=|\langle u_i,v_j\rangle|^2$ for eigenbases ordered decreasingly; apply Lemma 11.3 (its lower bound is
$\sum_ip^\downarrow_iq^\downarrow_{k+1-i}=\sum_ip^\downarrow_iq^\uparrow_i$). (b) The first inequality is
Theorem 1.9. By Theorem 1.9 and Theorem 11.4, $\|PQ\|_1=\max_U|\operatorname{Tr}(P\,QU)|\le\sum_is_i(P)s_i(QU)$, and
$s_i(QU)=q^\downarrow_i$ because $|QU|=(U^*Q^2U)^{1/2}=U^*QU$ (Theorem 1.6). $\square$

So the violation in Proposition 3.9 is at most $\sum_ix^\downarrow_i(y^\downarrow_i-y^\uparrow_i)$. For Exp. 2 this
number is the exact answer:

**Proposition 11.6 (Exp. 2 for fixed spectra).** In Exp. 2, $|X|=X_0$, $|Y|=Y_0$ and
$z=\|X_0VY_0\|_1-\operatorname{Tr}(X_0Y_0)$. If the spectra $x_c$ of $(X_0)_c$ and $y_c$ of $(Y_0)_c$ are fixed and the
eigenbases and $V$ vary, then
$$
\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}\big(y^\downarrow_{c,i}-y^\uparrow_{c,i}\big),
$$
and both are attained.

*Proof.* $X^*X=X_0U^*UX_0=X_0^2$, so $|X|=X_0$ (Theorem 1.6); likewise $|Y|=Y_0$, and $|UA|=|A|$ gives
$\operatorname{Tr}|XY|=\|X_0VY_0\|_1$. All variables act blockwise and $z$ is a sum over blocks, so let $C=1$.
Put $Y'=VY_0V^*\ge0$ (spectrum $y$); since $|AV|=V^*|A|V$, $\|X_0VY_0\|_1=\|X_0Y'V\|_1=\|X_0Y'\|_1$. Corollary 11.5
for the pairs $(X_0,Y')$ and $(X_0,Y_0)$, with $\|X_0Y'\|_1\ge\operatorname{Tr}(X_0Y')$, gives
$$
\sum_ix^\downarrow_iy^\uparrow_i-\sum_ix^\downarrow_iy^\downarrow_i\le\|X_0Y'\|_1-\operatorname{Tr}(X_0Y_0)\le\sum_ix^\downarrow_iy^\downarrow_i-\sum_ix^\downarrow_iy^\uparrow_i .
$$
Let $J$ be the reversal permutation, $Je_i=e_{k+1-i}$; then $J\operatorname{diag}(w)=\operatorname{diag}(w_k,\dots,w_1)J$, and
the trace norm of $\operatorname{diag}(w)J$ is $\sum_i|w_i|$. With $X_0=\operatorname{diag}(x^\downarrow)$ and $V=J$:
for $Y_0=\operatorname{diag}(y^\uparrow)$, $X_0VY_0=\operatorname{diag}(x^\downarrow_iy^\downarrow_i)J$ and
$\operatorname{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\uparrow_i$ (the upper bound); for $Y_0=\operatorname{diag}(y^\downarrow)$,
$X_0VY_0=\operatorname{diag}(x^\downarrow_iy^\uparrow_i)J$ and $\operatorname{Tr}(X_0Y_0)=\sum_ix^\downarrow_iy^\downarrow_i$
(the lower bound). $\square$

**Corollary 11.7 (envelope of Exp. 2).** Let $\ell_X,\ell_Y\in(0,1]$ and some $k_c\ge2$. Over all $X_0,Y_0\ge0$
with $\lambda_{\max}=1$, $\lambda_{\min}(X_0)=\ell_X$, $\lambda_{\min}(Y_0)=\ell_Y$, and all unitaries,
$$
\sup z=-\inf z=\sum_c\Big\lfloor\frac{k_c}2\Big\rfloor(1-\ell_X)(1-\ell_Y).
$$

*Proof.* In one block write $x_i=x^\downarrow_{c,i}$, $y_i=y^\downarrow_{c,i}$, $k=k_c$, so $y^\uparrow_{c,i}=y_{k+1-i}$.
Pairing $i$ with $k+1-i$, $x_i(y_i-y_{k+1-i})+x_{k+1-i}(y_{k+1-i}-y_i)=(x_i-x_{k+1-i})(y_i-y_{k+1-i})$, and for odd
$k$ the middle index contributes $x_i(y_i-y_i)=0$. Hence
$$
\sum_ix_i(y_i-y_{k+1-i})=\sum_{i\le\lfloor k/2\rfloor}(x_i-x_{k+1-i})(y_i-y_{k+1-i})\le\Big\lfloor\frac k2\Big\rfloor(1-\ell_X)(1-\ell_Y),
$$
as $x_i\in[\ell_X,1]$, $y_i\in[\ell_Y,1]$. Equality holds if in every block the $\lfloor k_c/2\rfloor$ largest
eigenvalues are $1$ and the $\lfloor k_c/2\rfloor$ smallest are $\ell_X$ (resp. $\ell_Y$), the middle one being
arbitrary; these spectra have global extremes $1$ and $\ell_X$ ($\ell_Y$) because some $k_c\ge2$. Proposition 11.6
concludes. $\square$

Here $1-\ell_X=2\delta/(1+\delta)$ for $\Delta(|X|)=\delta$; for odd $k_c$ the middle eigenvalue is irrelevant. The code
checks Theorem 11.4, Corollary 11.5 and the extremal configurations of Proposition 11.6 for $k=3$, $C=2$.

```python
M5 = TypeIAlgebra([5], [5], precision='double')
A = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))         # Ginibre, no structure
B = M5.operator(torch.randn(2000, 1, 5, 5, dtype=torch.complex128))
sA, sB = torch.linalg.svdvals(A.matrix[:, 0]), torch.linalg.svdvals(B.matrix[:, 0])   # decreasing
assert ((A @ B).trace.abs() <= (sA * sB).sum(-1) + 1e-12).all()             # Theorem 11.4
Ph, Qh = (M5.operator_from_eigenvalues(lambda k: torch.randn(2000, k), batch_size=2000, force_self_adjoint=True)
          for _ in range(2))
p, q = Ph.eigenvalues()[0], Qh.eigenvalues()[0]                             # increasing order
t = (Ph @ Qh).trace.real
assert (t <= (p * q).sum(-1) + 1e-12).all() and (t >= (p * q.flip(-1)).sum(-1) - 1e-12).all()   # Cor. 11.5(a)
Pp, Qp = Ph @ Ph, Qh @ Qh                                                   # positive, spectra p^2, q^2
nuc = (Pp @ Qp).trace_norm()
assert ((Pp @ Qp).trace.real <= nuc + 1e-12).all()                          # Cor. 11.5(b)
assert (nuc <= ((p ** 2).sort(-1)[0] * (q ** 2).sort(-1)[0]).sum(-1) + 1e-10).all()

alg3 = TypeIAlgebra([3, 3], [3, 3], precision='double')
x_dn = torch.tensor([1.0, 0.6, ell(0.5)], dtype=torch.float64)              # Delta = 0.5, middle value 0.6
y_dn = torch.tensor([1.0, 0.7, ell(0.8)], dtype=torch.float64)              # Delta = 0.8, middle value 0.7
J = torch.eye(3, dtype=torch.float64).flip(0)                               # reversal permutation
blocks = lambda m: alg3.from_blocks([m, m])                                 # the same block in both sectors

def z_exp2(X0, V, Y0):
    return ((X0 @ V @ Y0).trace_norm() - (X0 @ Y0).trace.real).item()

X0 = blocks(torch.diag(x_dn))
z_max = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn.flip(0))))
z_min = z_exp2(X0, blocks(J), blocks(torch.diag(y_dn)))
exact = 2 * (x_dn * (y_dn - y_dn.flip(0))).sum().item()                     # Proposition 11.6
envelope = 2 * (1 - ell(0.5)) * (1 - ell(0.8))                              # Corollary 11.7, floor(3/2) = 1
print(f"Exp. 2, k = 3, C = 2: z_max = {z_max:.6f}, z_min = {z_min:.6f}, formula {envelope:.6f}")
assert abs(z_max - exact) < 1e-12 and abs(z_min + exact) < 1e-12 and abs(exact - envelope) < 1e-12
```

## 11.3 A trace inequality with an equality condition

For $X=X^*$ let $X_\pm=f_\pm(X)$ with $f_+(t)=\max(t,0)$, $f_-(t)=\max(-t,0)$. Then $X_\pm\ge0$, $X=X_+-X_-$,
$|X|=X_++X_-$ and $X_+X_-=0$ (Theorem 1.5(c)).

**Lemma 11.8.** For $A,B\ge0$ in $M_n(\mathbb C)$, $\|A-B\|_1\le\operatorname{Tr}A+\operatorname{Tr}B$, with equality iff
$AB=0$.

*Proof.* Let $H=A-B$ and $Q=g(H)$ with $g=1$ on $(0,\infty)$, $g=0$ elsewhere. As $|t|=t(2g(t)-1)$, $|H|=H(2Q-1)$
and
$$
\|H\|_1=\operatorname{Tr}\big(A(2Q-1)\big)-\operatorname{Tr}\big(B(2Q-1)\big)
=\operatorname{Tr}A+\operatorname{Tr}B-2\operatorname{Tr}\big((1-Q)A(1-Q)\big)-2\operatorname{Tr}(QBQ).
$$
The subtracted traces are $\ge0$, and $\operatorname{Tr}(RAR)=\|A^{1/2}R\|_2^2$ vanishes iff $AR=0$; so equality holds
iff $A(1-Q)=0$ and $BQ=0$, and then $AB=AQB=A(BQ)^*=0$. Conversely, if $AB=0$ then $BA=0$,
$(A-B)^2=A^2+B^2=(A+B)^2$, and Theorem 1.6 gives $|A-B|=A+B$. $\square$

**Proposition 11.9 (Exp. 3).** For $X=X^*$ and $Y\ge0$ in $\mathcal M$,
$\operatorname{Tr}|YXY|\le\operatorname{Tr}(Y|X|Y)$, with equality iff $YX_+Y^2X_-Y=0$, in particular if $XY=YX$.

*Proof.* In $M_D(\mathbb C)$, $A=YX_+Y$ and $B=YX_-Y$ are positive (Corollary 1.7 with $T=X_\pm^{1/2}Y$), $YXY=A-B$
and $\operatorname{Tr}(Y|X|Y)=\operatorname{Tr}A+\operatorname{Tr}B$; apply Lemma 11.8. If $XY=YX$, then $Y$ commutes with
$X_\pm$ (Theorem 1.5(d)), and $AB=Y^4X_+X_-=0$. $\square$

So $z\ge0$ in Exp. 3, and $\underline z=0$ (take $Y=1$). The equality condition is weaker than commutation: for
$X=\operatorname{diag}(1,-1,0)$ and $Y$ with rows $(1,0,a)$, $(0,1,0)$, $(a,0,1)$, $|a|<1$,
$AB=(Ye_1)\langle Ye_1,Ye_2\rangle(Ye_2)^*=0$, while $XY\ne YX$ for $a\ne0$.

```python
def z_exp3(X, Y):
    return (Y @ X.abs() @ Y).trace.real - (Y @ X @ Y).trace_norm()

Xs = alg3.operator_from_eigenvalues(lambda k: torch.randn(500, k), batch_size=500, force_self_adjoint=True)
Yp = alg3.operator_from_eigenvalues(lambda k: torch.rand(500, k), batch_size=500, force_positive=True)
assert (z_exp3(Xs, Yp) >= -1e-12).all()                                     # Proposition 11.9

M3 = TypeIAlgebra([3], [3], precision='double')
Xe = M3.operator(torch.diag(torch.tensor([1.0, -1.0, 0.0]))[None, None])
Ye = M3.operator(torch.tensor([[1.0, 0, 0.5], [0, 1, 0], [0.5, 0, 1]])[None, None])
print(f"equality without commutation: z = {z_exp3(Xe, Ye).item():.1e}, "
      f"||XY - YX|| = {(Xe @ Ye - Ye @ Xe).operator_norm().item():.2f}")
assert abs(z_exp3(Xe, Ye).item()) < 1e-12 and (Xe @ Ye - Ye @ Xe).operator_norm().item() > 0.1
```

## 11.4 Sampling versus optimisation

How large is $|\operatorname{Tr}(AU)|$ typically, compared with its maximum $\|A\|_1$? The library samples $U$ as
the factor $Q$ of the QR decomposition $Z=QR$, $R_{ii}>0$, of a **Ginibre matrix** $Z$ (independent entries with
density $\pi^{-1}e^{-|z|^2}$; Mezzadri's recipe), independently in every block.

**Lemma 11.10.** $Q$ is Haar distributed on $U(k)$; in particular, for every fixed unitary $V\in M_k(\mathbb C)$, $VQ$ has
the law of $Q$. The unitary $\bigoplus_cQ_c$ with independent blocks is Haar distributed on the unitary group of
$\mathcal M$.

*Proof.* This is Theorem 5.3 (the left invariance is the equivariance step of its proof) and Corollary 5.4. $\square$

(For concentration of measure and more on these groups see E. Meckes, *The Random Matrix Theory of the Classical
Compact Groups*, Cambridge 2019.)

**Proposition 11.11.** For Haar $U\in U(k)$ and $A\in M_k(\mathbb C)$: $\mathbb E\operatorname{Tr}(AU)=0$ and
$\mathbb E|\operatorname{Tr}(AU)|^2=\|A\|_2^2/k$.

*Proof.* $\mathbb EU_{ij}=0$ by Corollary 5.5. Proposition 5.7 (with $C=1$) for $X=|e_j\rangle\langle e_m|$ gives the
Weingarten formula of order one,
$\mathbb E[U_{ij}\bar U_{lm}]=\mathbb E(UXU^*)_{il}=\delta_{il}\delta_{jm}/k$; for $k\ge2$ it also follows from the
order-two formula, Theorem 5.10, by setting $i_2=i'_2$, $j_2=j'_2$ and summing over both, since
$\sum_{i_2,j_2}|U_{i_2j_2}|^2=k$. With $\operatorname{Tr}(AU)=\sum_{ij}A_{ji}U_{ij}$,
$\mathbb E|\operatorname{Tr}(AU)|^2=\sum_{ijlm}A_{ji}\bar A_{ml}\mathbb E[U_{ij}\bar U_{lm}]=\frac1k\sum_{ij}|A_{ji}|^2$. $\square$

**Corollary 11.12.** Let $X,Y\ge0$ in $M_k(\mathbb C)$ have spectra in $[a,b]$, $a>0$, and $U$ be Haar. Then
$$
\frac{\big(\mathbb E|\operatorname{Tr}(XUY)|^2\big)^{1/2}}{\max_V|\operatorname{Tr}(XVY)|}=\frac{\|YX\|_2}{\sqrt k\,\|YX\|_1}\le\min\Big(\frac1{\sqrt k},\frac{b^2}{a^2k}\Big),
\qquad
\mathbb P\big(|\operatorname{Tr}(XUY)|\ge\operatorname{Tr}(XY)\big)\le\frac{b^4}{a^4k^2}.
$$

*Proof.* The equality is Proposition 11.11 for $A=YX$ and Proposition 3.9. $\sum_is_i^2\le(\sum_is_i)^2$ gives
$1/\sqrt k$. As $s_1(A)\le\|Y\|\|X\|\le b^2$ and $s_k(A)=\min_{\|v\|=1}\|YXv\|\ge a^2$, $\|A\|_2^2\le kb^4$ and
$\|A\|_1\ge ka^2$. Finally $\operatorname{Tr}(XY)=\operatorname{Tr}(X^{1/2}YX^{1/2})\ge a\operatorname{Tr}X\ge ka^2$, and
Markov's inequality for $|\operatorname{Tr}(XUY)|^2$ bounds the probability by $(kb^4/k)/(ka^2)^2$. $\square$

So $N$ Haar samples find a violation of $|\operatorname{Tr}(XUY)|\le\operatorname{Tr}(XY)$ with probability at most
$Nb^4/(a^4k^2)$, although one exists for every non-commuting pair. The truth is more drastic: $f(U)=|\operatorname{Tr}(AU)|$
is Lipschitz with constant $L=\|A\|_2$ for $\|U-V\|_2$, so $\mathbb P(f\ge\mathbb Ef+t)\le e^{-ckt^2/L^2}$ with an
absolute $c>0$ (concentration of measure, stated without proof; Meckes, *loc. cit.*, Ch. 5): an exponent of order
$k^2a^4/b^4$ here. The code checks Proposition 11.11 and prints, for spectra uniform in $[0.2,1.2]$ ($b/a=6$), the
ratio of Corollary 11.12, $\sup z$ and the best $z$ of $2000$ samples, both relative to $\operatorname{Tr}(XY)$.

```python
def random_pair(alg_k):
    pos = lambda k: 0.2 + torch.rand(1, k)                                  # spectrum in [a, b] = [0.2, 1.2]
    return tuple(alg_k.operator_from_eigenvalues(pos, force_positive=True) for _ in range(2))

M6 = TypeIAlgebra([6], [6], precision='double')
Xf, Yf = random_pair(M6)
N = 20_000
t2 = (Yf @ Xf @ M6.random_unitary_operator(N)).trace.abs() ** 2
m2, err, exact2 = t2.mean().item(), t2.std().item() / math.sqrt(N), (Yf @ Xf).frobenius_norm().item() ** 2 / 6
print(f"E|Tr(AU)|^2 = {m2:.4f} +- {err:.4f}   (exact ||A||_2^2 / k = {exact2:.4f})")
assert abs(m2 - exact2) < 5 * err                                           # Proposition 11.11

print(" k   rms/max   sup z/TrXY   best of 2000 samples")
for k in (2, 4, 8, 16, 32):
    Mk = TypeIAlgebra([k], [k], precision='double')
    Xk, Yk = random_pair(Mk)
    trXY, nuc = (Xk @ Yk).trace.real.item(), (Yk @ Xk).trace_norm().item()
    ratio = (Yk @ Xk).frobenius_norm().item() / (math.sqrt(k) * nuc)
    best = (Xk @ Mk.random_unitary_operator(2000) @ Yk).trace.abs().max().item()
    print(f"{k:2d}   {ratio:.4f}    {(nuc - trXY) / trXY:.2e}     {(best - trXY) / trXY:+.3f}")
    assert ratio <= min(1 / math.sqrt(k), 36 / k) + 1e-12                   # Corollary 11.12
    assert nuc > trXY and best < nuc and (k < 4 or best < trXY)            # from k = 4 on: no violation found
```

The ratio falls like $1/k$, and from $k=4$ on the best sample does not even reach $z=0$. The extremes are sets of
tiny Haar measure, which a gradient method reaches directly. For the pair `Xf, Yf` both extremes over $U$ are known:

**Lemma 11.13.** If $k\ge2$, then $\min_U|\operatorname{Tr}(AU)|=0$ for every $A\in M_k(\mathbb C)$.

*Proof.* Let $A=W\Sigma V^*$ and $P$ the cyclic shift, $Pe_i=e_{(i\bmod k)+1}$, whose diagonal vanishes for $k\ge2$. For
$U=VPW^*$, $\operatorname{Tr}(AU)=\operatorname{Tr}(\Sigma P)=0$. $\square$

So $\inf_Uz=-\operatorname{Tr}(XY)$ and $\sup_Uz=\|YX\|_1-\operatorname{Tr}(XY)$. One call of `extremize` finds both, with
$12$ starts maximising $z$ and $12$ maximising $-z$. Sampling comes close to the infimum, which lies next to the
typical values ($|\operatorname{Tr}(XUY)|$ is typically of order $\|YX\|_2/\sqrt k$), but not to the supremum:

```python
trXY = (Xf @ Yf).trace.real
sup_exact, inf_exact = ((Yf @ Xf).trace_norm() - trXY).item(), -trXY.item()
U6 = opt.UnitaryParam(M6, 24)                                               # 24 starts
pm = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(12)   # 12 maximise z, 12 maximise -z
best = opt.extremize(lambda: pm * ((Xf @ U6.operator() @ Yf).trace.abs() - trXY), [U6],
                     maximize=True, rounds=6, steps=30).reshape(2, 12).max(1)[0]
best_up, best_down = best[0].item(), -best[1].item()
z_mc = (Xf @ M6.random_unitary_operator(N) @ Yf).trace.abs() - trXY
print(f"sup z: exact {sup_exact:.6e}, ascent {best_up:.6e}, best sample {z_mc.max().item():+.4f}")
print(f"inf z: exact {inf_exact:.6f}, descent {best_down:.6f}, worst sample {z_mc.min().item():+.6f}")
assert best_up <= sup_exact + 1e-12 and best_down >= inf_exact - 1e-12    # attained values: one-sided
assert sup_exact - best_up < 1e-3 * sup_exact and best_down - inf_exact < 1e-3 * abs(inf_exact)
assert z_mc.min().item() - inf_exact < 1e-2 * abs(inf_exact)               # sampling finds the typical side ...
assert z_mc.max().item() < 0                                               # ... but misses the violation
```

> **Physical meaning.** Concentration of measure is the mathematical form of **typicality** in statistical
> mechanics: in high dimension almost every state or evolution gives the same value of a macroscopic quantity. A
> sharp inequality is decided by fine-tuned configurations — aligned eigenbases, two-level spectra — that a random
> unitary (a random circuit, a random scattering matrix) essentially never produces.

## 11.5 Multistart search for envelopes

`optimize.extremize(objective, params, maximize)` runs Adam in `rounds` rounds of `steps` steps (learning rate
`lr * decay**round`), re-centres the unitary charts after every round, and returns **per start** the best value
attained. Every batch element is an independent start with its own slots, signs, $U_0$ and contrast. The signs of
`SelfAdjointParam` are discrete: between rounds, `discrete_step` tries eight single sign flips per start, one after
the other, and keeps each flip that improves the objective.

Since every iterate is admissible (Proposition 11.2), every returned value is attained: a maximisation gives a
**lower bound** for the supremum, a minimisation an upper bound for the infimum, and nothing more — a local search
can stop at a local extremum. The other direction needs a proof, like Corollary 11.7 and Proposition 11.9, or a
convex relaxation (Section 11.6). A computed envelope is checked (i) against closed forms; (ii) against
configurations built by hand, such as direct sums of smaller solutions, which it must not fall below; (iii) against
structured samples (two-level spectra, permutation bases) that come close to the extremes in every dimension,
unlike Haar samples (`scripts/bounds_checks.py`).

**Proposition 11.14 (lower envelope of Exp. 1).** Let all $k_c\ge2$. Over $X\in\mathcal P_\delta$, $Y\ge0$ with
$\lambda_{\max}(Y)=1$ and unitaries $U\in\mathcal M$, $\ \min z=-\big(D-2\delta/(1+\delta)\big)$.

*Proof.* $Y\le1$ gives $X^{1/2}YX^{1/2}\le X$, so $\operatorname{Tr}(XY)\le\operatorname{Tr}X\le(D-1)+\ell(\delta)$, as one
eigenvalue of $X$ is $\ell(\delta)$ and the others are $\le1$; hence $z\ge-\operatorname{Tr}(XY)\ge-(D-2\delta/(1+\delta))$.
Equality holds for $X$ with eigenvalues $\ell(\delta)$ once and $1$ otherwise, $Y=1$, and $U=\bigoplus_cU_c$ with
$\operatorname{Tr}(X_cU_c)=0$ (Lemma 11.13). $\square$

The code checks Proposition 11.14 on $M_2\oplus M_2$, Corollary 11.7 on $M_3\oplus M_3$ (both contrasts prescribed,
$k=3$ odd) and Proposition 11.9; suprema and infima share one call, as above. The keyword `is_self_adjoint` selects
the algorithm for the trace norm (singular values, or eigenvalues for $YXY$) without a test at every step. The last
check is of type (ii): $z$ of Exp. 3 is additive over direct sums, and the contrast of a direct sum is the largest
contrast of its summands (Exercise 3), so
$$
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\max_{\delta'\le\delta}\overline z_{M_2}(\delta').
$$

```python
torch.manual_seed(1)
S = 6                                                                       # starts per grid point
grid = torch.tensor([0.0, 0.5, 0.9], dtype=torch.float64)
alg22 = TypeIAlgebra([2, 2], [2, 2], precision='double')                    # D = 4

# Exp. 1, lower envelope along Delta(X); Delta(Y) free in [0, 0.999]
X1, Y1, U1 = (opt.PositiveParam(alg22, 3 * S, delta=grid.repeat_interleave(S)),
              opt.PositiveParam(alg22, 3 * S), opt.UnitaryParam(alg22, 3 * S))
def z1():
    x, y = X1.operator(), Y1.operator()
    return (x @ U1.operator() @ y).trace.abs() - (x @ y).trace.real
inf1 = opt.extremize(z1, [X1, Y1, U1], maximize=False, rounds=4, steps=30).reshape(3, S).min(1)[0]
exact1 = -(4 - 2 * grid / (1 + grid))                                       # Proposition 11.14
print("Exp. 1  inf:", [f"{v:.4f}" for v in inf1.tolist()], " exact:", [f"{v:.4f}" for v in exact1.tolist()])
assert (inf1 >= exact1 - 1e-12).all() and (inf1 - exact1).max() < 4e-3

# Exp. 2 on M_3 (+) M_3, Delta(X0) = 0.5, Delta(Y0) = 0.8: sup and inf in one call
sign = torch.tensor([1.0, -1.0], dtype=torch.float64).repeat_interleave(S)
X2, Y2, V2 = (opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.5, dtype=torch.float64)),
              opt.PositiveParam(alg3, 2 * S, delta=torch.full((2 * S,), 0.8, dtype=torch.float64)),
              opt.UnitaryParam(alg3, 2 * S))
def z2():
    x0, y0 = X2.operator(), Y2.operator()
    M = alg3.operator((x0 @ V2.operator() @ y0).matrix, is_self_adjoint=False)
    return sign * (M.trace_norm() - (x0 @ y0).trace.real)
best2 = opt.extremize(z2, [X2, Y2, V2], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]
sup2, inf2 = best2[0].item(), -best2[1].item()
print(f"Exp. 2  sup {sup2:.6f}, inf {inf2:.6f}, exact +-{envelope:.6f}")
assert sup2 <= envelope + 1e-12 and inf2 >= -envelope - 1e-12               # Corollary 11.7: one-sided ...
assert envelope - sup2 < 1e-3 * envelope and inf2 + envelope < 1e-3 * envelope   # ... and reached

# Exp. 3: maximise sign * z at the contrasts `deltas` of |X|; Delta(Y) free in [0, 0.999]
def exp3(alg_, deltas, sign):
    X3, Y3 = opt.SelfAdjointParam(alg_, len(deltas), delta=deltas), opt.PositiveParam(alg_, len(deltas))
    def z3():
        y = Y3.operator()
        M = alg_.operator((y @ X3.operator() @ y).matrix, is_self_adjoint=True)
        return sign * ((y @ X3.abs_operator() @ y).trace.real - M.trace_norm())
    return opt.extremize(z3, [X3, Y3], maximize=True, rounds=4, steps=30).reshape(2, S).max(1)[0]

best3 = exp3(alg22, torch.full((2 * S,), 0.5, dtype=torch.float64), sign)
sup_M22, inf3 = best3[0].item(), -best3[1].item()
s_M2 = exp3(TypeIAlgebra([2], [2], precision='double'),
            torch.tensor([0.0, 0.5], dtype=torch.float64).repeat_interleave(S), 1.0)
direct_sum = (s_M2[1] + s_M2.max()).item()                                 # delta = 0.5, best delta' <= 0.5
conj = [(1 - dl) * (2 * 0.999 / 1.999) ** 2 for dl in (0.0, 0.5)]          # unproven closed form, d = 0.999
print(f"Exp. 3  inf {inf3:.1e};  sup for M_2 at delta = 0, 0.5: {s_M2[0]:.5f}, {s_M2[1]:.5f}"
      f"  (closed form {conj[0]:.5f}, {conj[1]:.5f})")
print(f"        sup for M_2 (+) M_2 at 0.5: search {sup_M22:.5f}, direct sum {direct_sum:.5f}")
assert inf3 >= -1e-12 and inf3 < 1e-6                                       # Proposition 11.9
assert sup_M22 >= direct_sum - 5e-3 * direct_sum                           # (ii), up to 0.5 %
```

The search on $M_2\oplus M_2$ falls short of the direct sum by about $0.2\,\%$: it has not fully converged, and the
best lower bound for $\overline z_{M_2\oplus M_2}(0.5)$ is the larger number — detecting this is the purpose of check
(ii). The values for $M_2$ agree to $10^{-4}$ with $(1-\delta)\big(2d/(1+d)\big)^2=(1-\delta)(1-\ell(d))^2$, where
$d=0.999$ is the largest admitted contrast of $Y$: a closed form that the companion paper observed numerically (in
real algebras) and could not prove; the code does not assert it. For larger blocks the upper envelopes of Exp. 1 and
Exp. 3 are, to our knowledge, not known in closed form; computed curves are then lower bounds and conjectures,
nothing more.

## 11.6 What to compare with

**Manopt**, **pymanopt** and **geoopt** provide matrix manifolds (Stiefel manifolds, orthogonal and unitary groups)
with retractions such as QR, polar or Cayley, and Riemannian solvers (conjugate gradients, trust regions, Riemannian
Adam) that are often more efficient than Adam in a re-centred chart. **cvxpy** with an SDP solver solves convex
problems in operator variables to global optimality, with a dual certificate up to the solver tolerance, and
semidefinite relaxations give rigorous bounds in the direction a local search cannot. **NCTSSOS**, **ncpol2sdpa** and **NCSOStools** bound eigenvalues and traces of noncommutative
polynomials in all dimensions at once (Section 0.3). The library adds exact parametrisations of the spectral and
contrast constraints of $\bigoplus_cM_{k_c}$ (the extremes may sit in any sector), unitaries that respect the
blocks, and one batched computation for many starts and constraint values, on a CPU or a GPU. It produces candidate
extremisers and lower bounds for suprema; proofs or relaxations must supply the rest.

## Exercises

1. (**Cayley transform.**) For $S^*=-S$ show that $C(S)=(1-S)(1+S)^{-1}$ is unitary, and that $S\mapsto C(S)$ is a
   bijection onto the unitaries without eigenvalue $-1$. So, unlike the exponential map (Lemma 11.1), a Cayley chart
   $U_0C(S)$ misses $-U_0$.
2. (**Real Haar measure.**) For $O$ Haar-distributed on $O(k)$ (QR of a real Ginibre matrix) and real $A$, show
   $\mathbb E\operatorname{Tr}(AO)^2=\|A\|_2^2/k$, and check it with a real algebra.
3. (**Direct sums.**) Show that $z$ of Exp. 3 is additive over direct sums and that
   $\Delta(A\oplus A')=\max(\Delta(A),\Delta(A'))$ for $A,A'\ge0$ with $\lambda_{\max}=1$. Deduce the direct-sum
   inequality of Section 11.5 (with $\Delta(Y)$ free in a fixed interval $[0,b]$).

<details>
<summary>Solutions</summary>

1. $S$ has imaginary eigenvalues, so $1\pm S$ are invertible; they commute, and $(1+S)^*=1-S$. Hence
   $C^*C=(1-S)^{-1}(1+S)(1-S)(1+S)^{-1}=1$, and $1+C=2(1+S)^{-1}$ is invertible, so $-1\notin\operatorname{spec}C$.
   Conversely, for a unitary $U$ with $1+U$ invertible, $S=(1-U)(1+U)^{-1}$ is skew-adjoint:
   $S^*=(1+U^*)^{-1}(1-U^*)=(U+1)^{-1}U\,U^*(U-1)=-S$. Since $C$ commutes with $S$, $C(1+S)=1-S$ is equivalent to
   $S(1+C)=1-C$, i.e. $S=(1-C)(1+C)^{-1}$: the two maps are inverse to each other.

```python
Sk = torch.randn(4, 4, dtype=torch.complex128)
Sk = Sk - Sk.mH                                                   # skew-adjoint
I4 = torch.eye(4, dtype=torch.complex128)
Cay = (I4 - Sk) @ torch.linalg.inv(I4 + Sk)
assert torch.allclose(Cay.mH @ Cay, I4, atol=1e-12)
assert torch.allclose((I4 - Cay) @ torch.linalg.inv(I4 + Cay), Sk, atol=1e-10)
assert torch.linalg.eigvals(Cay).add(1).abs().min() > 1e-6        # -1 is never an eigenvalue
```

2. Theorem 5.3 and its proof hold for a real Gaussian matrix (density $\propto e^{-\|Z\|_2^2/2}$), the orthogonal
   group and orthogonal $V$ (the remark after Theorem 5.3), so $VO$ has the law of $O$.
   Left multiplication by $\operatorname{diag}(\pm1)$ gives $\mathbb E[O_{ij}O_{lm}]=0$ for $i\ne l$, permutations show
   that $\mathbb E[O_{ij}O_{im}]$ does not depend on $i$, and $\sum_iO_{ij}O_{im}=\delta_{jm}$; so
   $\mathbb E[O_{ij}O_{lm}]=\delta_{il}\delta_{jm}/k$ and $\mathbb E\operatorname{Tr}(AO)^2=\sum_{ij}A_{ji}^2/k$.

```python
Mr = TypeIAlgebra([5], [5], complex_valued=False, precision='double')
Ar = Mr.operator(torch.randn(1, 1, 5, 5, dtype=torch.float64))
tr2 = (Ar @ Mr.random_unitary_operator(20_000)).trace ** 2
assert abs(tr2.mean().item() - Ar.frobenius_norm().item() ** 2 / 5) < 5 * tr2.std().item() / math.sqrt(20_000)
```

3. Traces and trace norms of block-diagonal operators are sums over blocks, and $|X\oplus X'|=|X|\oplus|X'|$;
   $\lambda_{\min}(A\oplus A')=\min(\lambda_{\min}(A),\lambda_{\min}(A'))$, and $\Delta$ decreases with $\lambda_{\min}$.
   So if $(X,Y)$ is admissible at $\delta$ and $(X',Y')$ at $\delta'\le\delta$, then $(X\oplus X',Y\oplus Y')$ is
   admissible at $\delta$, with $z=z(X,Y)+z(X',Y')$; take suprema over both pairs and over $\delta'\le\delta$.

</details>

**Next:** Chapter 12 (in preparation) treats equivariant machine learning: layers that commute with a symmetry are
elements of a commutant (Chapter 2), built from the same batched, differentiable operators.
