**English** · [Español](es/F_cones_and_order.md) · [Français](fr/F_cones_and_order.md) · [Deutsch](de/F_cones_and_order.md) · [中文](zh/F_cones_and_order.md) · [日本語](ja/F_cones_and_order.md) · [Русский](ru/F_cones_and_order.md)

# Appendix F. Symmetric cones, order isomorphisms and Jordan maps

Chapter 13 connected the order of observables with their Jordan product through three theorems: the
Koecher–Vinberg theorem (Theorem 13.8), Kadison's theorem (Theorem 13.10) and the description of the Jordan
automorphisms of $M_n(\mathbb C)$ (Proposition 13.11 proves it for $n=2$). This appendix proves the first in full
generality and the other two for the algebras $\bigoplus_cM_{k_c}(\mathbb C)$ of the course, for all block sizes. On
the way it proves the facts about Euclidean Jordan algebras that Section 13.6 needs: a spectral theorem, the
self-duality of the cone of squares, and $P(y)Q\subseteq Q$.

**You will learn** why every Euclidean Jordan algebra has power associativity and a spectral theorem; why its cone
of squares is self-dual and homogeneous; how a Jordan product is recovered from a self-dual homogeneous cone alone,
through the third derivative of the logarithm of its characteristic function; why a unital linear bijection that
preserves positivity in both directions preserves squares; and why every such map of $\bigoplus_cM_{k_c}$ permutes
equal blocks and acts on each block by $UxU^*$ or by $Ux^TU^*$.

Besides linear algebra and calculus, Section F.3 uses Lebesgue integration (Fatou's lemma, dominated convergence and
differentiation under the integral sign, the linear change of variables formula); nothing else is assumed. A
**Euclidean space** is a finite-dimensional real vector space $V\ne\{0\}$ with an inner product
$\langle\cdot,\cdot\rangle$; $|x|$ is the norm, $\|A\|$ the operator norm of a linear map $A$ of $V$, $A^T$ its adjoint
($\langle Ax,y\rangle=\langle x,A^Ty\rangle$), and $A\ge0$ ($A>0$) means that $A=A^T$ is positive semidefinite
(definite). $e^A=\exp A=\sum_{k\ge0}A^k/k!$.

## F.1 Euclidean Jordan algebras

As in Section 13.1, a **Jordan algebra** is a real vector space with a bilinear commutative product $\circ$ such that
$x\circ(x^2\circ y)=x^2\circ(x\circ y)$, i.e. $[L(x),L(x^2)]=0$, where $L(x)y=x\circ y$ and $x^2=x\circ x$; and
$P(x)=2L(x)^2-L(x^2)$.

**Definition.** A **Euclidean Jordan algebra** is a finite-dimensional Jordan algebra $V$ with a unit $e$
($e\circ x=x$) and an inner product for which every $L(x)$ is symmetric: $\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$.
Its **cone of squares** is $Q=\{x^2:x\in V\}$.

Examples: $M_{sa}$ and $\operatorname{Herm}_n(\mathbb R)$ with $\langle x,y\rangle=\operatorname{Tr}(xy)$ (Lemma 13.1;
$\operatorname{Tr}((x\circ y)z)=\frac12\operatorname{Tr}(xyz+yxz)$ is symmetric in $x,y,z$ by cyclicity and
$\operatorname{Tr}(a^T)=\operatorname{Tr}a$); the spin factors $V_m$ with $\langle x,y\rangle=\operatorname{tr}(x\circ y)$,
since $\operatorname{tr}((x\circ y)\circ z)=2(abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w)$ for $x=(a,u)$, $y=(b,v)$, $z=(c,w)$
is symmetric; and $\mathbb R^N$ with the coordinatewise product.

**Lemma F.1.** In every Jordan algebra, for all $a,b,d,x,y$:

1. $[L(a),L(b\circ d)]+[L(b),L(d\circ a)]+[L(d),L(a\circ b)]=0$;
2. $L(x^2\circ y)=L(x^2)L(y)-2L(x)L(y)L(x)+2L(x\circ y)L(x)$.

*Proof.* (1) For $s,t,r\in\mathbb R$ and $x=sa+tb+rd$, $[L(x),L(x^2)]=0$ is a polynomial identity in $(s,t,r)$, so each
coefficient vanishes; the coefficient of $str$ is twice the left side of (1). (2) Put $a=b=x$ in (1):
$2[L(x),L(x\circ d)]+[L(d),L(x^2)]=0$. Applied to $y$ this reads
$2x\circ((x\circ d)\circ y)-2(x\circ d)\circ(x\circ y)+d\circ(x^2\circ y)-x^2\circ(d\circ y)=0$, that is,
$\big(2L(x)L(y)L(x)-2L(x\circ y)L(x)+L(x^2\circ y)-L(x^2)L(y)\big)d=0$ for all $d$. $\square$

**Proposition F.2 (power associativity).** Let $V$ be a Jordan algebra with unit $e$, $x\in V$, $x^0=e$ and
$x^{k+1}=x\circ x^k$. Then $x^p\circ x^q=x^{p+q}$ for all $p,q\ge0$, and every $L(x^p)$ is a polynomial in the
commuting operators $L(x)$ and $L(x^2)$. Hence $p(x)\circ q(x)=(pq)(x)$ for real polynomials $p,q$, and the operators
$L(p(x))$ commute with each other.

*Proof.* $L(x)$ and $L(x^2)$ commute by the Jordan identity. We show by induction on $N\ge2$: (A$_N$)
$x^i\circ x^j=x^{i+j}$ for $i+j\le N$, and (B$_N$) $L(x^i)$ is a polynomial in $L(x),L(x^2)$ for $i\le N$. Both are
clear for $N=2$. Assume them for $N$, and let $i+j=N+1$ with $i,j\ge1$ (for $i=0$ or $j=0$ there is nothing to prove).
Since $i\le N$, $L(x^i)$ commutes with $L(x)$, so by (A$_N$)
$$
x^i\circ x^j=L(x^i)L(x)x^{j-1}=L(x)L(x^i)x^{j-1}=x\circ x^N=x^{N+1}.
$$
This is (A$_{N+1}$). Lemma F.1(2) with $y=x^{N-1}$ and $x^2\circ x^{N-1}=x^{N+1}$ gives
$L(x^{N+1})=L(x^2)L(x^{N-1})-2L(x)L(x^{N-1})L(x)+2L(x^N)L(x)$, which is (B$_{N+1}$). The last sentence follows by
bilinearity. $\square$

**Theorem F.3 (spectral theorem).** Let $V$ be a Euclidean Jordan algebra and $x\in V$. There are unique
$k\ge1$, real numbers $\lambda_1>\dots>\lambda_k$ and nonzero $c_1,\dots,c_k\in V$ with
$$
c_i\circ c_j=\delta_{ij}c_i,\qquad \textstyle\sum_ic_i=e,\qquad x=\sum_i\lambda_ic_i .
$$
The $c_i$ are polynomials in $x$, and $p(x)=\sum_ip(\lambda_i)c_i$ for every real polynomial $p$.

We call $(c_i)$ a **complete system of orthogonal idempotents** and the $\lambda_i$ the **eigenvalues** of $x$, and put
$f(x)=\sum_if(\lambda_i)c_i$ for $f:\{\lambda_1,\dots,\lambda_k\}\to\mathbb R$; for instance $|x|$, $x_\pm=\max(\pm x,0)$,
$\sqrt x$ (if all $\lambda_i\ge0$), $\exp x$, and $\log x$ (if all $\lambda_i>0$).

*Proof.* Existence. $\mathbb R[x]=\{p(x)\}$ is a subspace that contains $e$ and, by Proposition F.2, satisfies
$p(x)\circ q(x)=(pq)(x)$. So $M=L(x)|_{\mathbb R[x]}$ maps $\mathbb R[x]$ into itself, is symmetric, and satisfies
$p(M)q(x)=p(x)\circ q(x)$. A symmetric operator is diagonalisable with real eigenvalues; let $\lambda_1>\dots>\lambda_k$
be those of $M$ and $m(t)=\prod_i(t-\lambda_i)$ its minimal polynomial. Then $m(x)=m(M)e=0$. With the Lagrange
polynomials $\ell_i(t)=\prod_{j\ne i}\frac{t-\lambda_j}{\lambda_i-\lambda_j}$, the polynomials
$\ell_i\ell_j-\delta_{ij}\ell_i$, $\sum_i\ell_i-1$ and $\sum_i\lambda_i\ell_i-t$ vanish at every $\lambda_l$, so they are
multiples of $m$; evaluating at $x$ shows that $c_i=\ell_i(x)$ have the three properties. If $c_i=0$, then
$\ell_i(M)q(x)=\ell_i(x)\circ q(x)=0$ for all $q$, so $\ell_i(M)=0$ although $\deg\ell_i<\deg m$; hence $c_i\ne0$.
Finally $p(x)=\sum_ip(\lambda_i)c_i$: true for $p=1$ and $p=t$, and $x\circ\sum_i\lambda_i^rc_i=\sum_i\lambda_i^{r+1}c_i$.

Uniqueness. Let $x=\sum_j\mu_jd_j$ with distinct $\mu_j$ and a complete system of nonzero orthogonal idempotents $d_j$.
As above $p(x)=\sum_jp(\mu_j)d_j$. The $d_j$ are linearly independent (multiply $\sum_ja_jd_j=0$ by $d_i$), so
$p(x)=0$ iff $p(\mu_j)=0$ for all $j$. Thus $m(x)=0$ gives $\{\mu_j\}\subseteq\{\lambda_i\}$, and
$m'(t)=\prod_j(t-\mu_j)$ satisfies $m'(M)q(x)=m'(x)\circ q(x)=0$, so $m$ divides $m'$ and $\{\lambda_i\}\subseteq\{\mu_j\}$.
Then $d_j=\ell_j(x)=c_j$. $\square$

**Lemma F.4 (Peirce decomposition).** Let $V$ be a Euclidean Jordan algebra and $c$ an idempotent ($c^2=c$).

1. $L(c)(2L(c)-I)(L(c)-I)=0$. Hence $V$ is the orthogonal sum of the eigenspaces $V(c,\nu)=\{z:c\circ z=\nu z\}$,
   $\nu\in\{0,\frac12,1\}$, and $L(c)\ge0$.
2. $V(c,1)\circ V(c,0)=\{0\}$; $V(c,\nu)\circ V(c,\nu)\subseteq V(c,\nu)$ and $V(c,\nu)\circ V(c,\frac12)\subseteq V(c,\frac12)$
   for $\nu\in\{0,1\}$; $V(c,\frac12)\circ V(c,\frac12)\subseteq V(c,0)+V(c,1)$.
3. $s=2c-e$ satisfies $s^2=e$, and $P(s)$ is a Jordan automorphism of $V$: it is $+1$ on $V(c,0)+V(c,1)$, $-1$ on
   $V(c,\frac12)$, and $P(s)(z\circ w)=P(s)z\circ P(s)w$.
4. Let $x=\sum_{i=1}^k\lambda_ic_i$ be as in Theorem F.3. The $L(c_i)$ commute, and $V$ is the orthogonal sum of the
   subspaces $V_{ii}=\{z:c_i\circ z=z\}$ and $V_{ij}=\{z:c_i\circ z=c_j\circ z=\frac12z\}$ ($i<j$). For
   $a=\sum_la_lc_l$ ($a_l\in\mathbb R$), $L(a)$ acts on $V_{ij}$ ($i\le j$) as the scalar $\frac12(a_i+a_j)$, and $P(a)$
   as the scalar $a_ia_j$.

*Proof.* (1) Lemma F.1(2) with $x=y=c$ gives $L(c)=L(c)^2-2L(c)^3+2L(c)^2$, i.e. $2L(c)^3-3L(c)^2+L(c)=0$, which is the
stated factorisation. So the eigenvalues of the symmetric operator $L(c)$ lie in $\{0,\frac12,1\}$.

(2) Lemma F.1(1) with $a=b=c$, $d=z\in V(c,\nu)$ gives $2\nu[L(c),L(z)]+[L(z),L(c)]=0$, so $L(z)$ commutes with $L(c)$
if $\nu\in\{0,1\}$ and then preserves every $V(c,\mu)$. This gives the second and third inclusions, and
$V(c,1)\circ V(c,0)\subseteq V(c,0)\cap V(c,1)=\{0\}$. For $z,w\in V(c,\frac12)$, apply Lemma F.1(1) with $(a,b,d)=(z,w,c)$
to $c$: since $L(w\circ c)=\frac12L(w)$ and $L(c\circ z)=\frac12L(z)$, the first two terms cancel, and
$[L(c),L(u)]c=L(c)^2u-L(c)u=0$ for $u=z\circ w$. Writing $u=u_1+u_{1/2}+u_0$ by (1), this is $-\frac14u_{1/2}=0$.

(3) $s^2=4c-4c+e=e$, $L(s)=2L(c)-I$ is $1,0,-1$ on $V(c,1),V(c,\frac12),V(c,0)$, and $P(s)=2L(s)^2-I$. Write
$z=z'+z''$, $w=w'+w''$ with $z',w'\in V(c,0)+V(c,1)$ and $z'',w''\in V(c,\frac12)$. By (2), $z'\circ w'+z''\circ w''$ lies in
$V(c,0)+V(c,1)$ and $z'\circ w''+z''\circ w'$ in $V(c,\frac12)$, so
$P(s)(z\circ w)=(z'-z'')\circ(w'-w'')=P(s)z\circ P(s)w$.

(4) The $c_l$ are polynomials in $x$, so the $L(c_l)$ commute (Proposition F.2); being symmetric, they have a common
orthonormal eigenbasis. On a common eigenvector $L(c_l)$ has an eigenvalue $\varepsilon_l\in\{0,\frac12,1\}$ by (1), and
$\sum_l\varepsilon_l=1$ because $\sum_lL(c_l)=L(e)=I$. So either one $\varepsilon_i=1$ and the others vanish, or
$\varepsilon_i=\varepsilon_j=\frac12$ for one pair $i<j$ and the others vanish; the common eigenspaces are the $V_{ij}$. Since
$a^2=\sum_la_l^2c_l$, on $V_{ij}$ ($i<j$) $L(a)=\frac12(a_i+a_j)$, $L(a^2)=\frac12(a_i^2+a_j^2)$ and
$P(a)=\frac12(a_i+a_j)^2-\frac12(a_i^2+a_j^2)=a_ia_j$; on $V_{ii}$, $L(a)=a_i$ and $P(a)=2a_i^2-a_i^2$. $\square$

## F.2 The cone of squares is symmetric

Let $\Omega$ be a nonempty open **convex cone** in a Euclidean space $V$ ($tx\in\Omega$ for $x\in\Omega$, $t>0$). Its
(open) **dual cone** is $\Omega^*=\{y:\langle x,y\rangle>0\text{ for all }x\in\bar\Omega\setminus\{0\}\}$; $\Omega$ is
**self-dual** if $\Omega^*=\Omega$. $G(\Omega)=\{g\in GL(V):g\Omega=\Omega\}$, and $\Omega$ is **homogeneous** if
$G(\Omega)$ acts transitively on $\Omega$. A **symmetric cone** is a self-dual homogeneous one. If $\Omega_1\subseteq\Omega_2$
then $\Omega_2^*\subseteq\Omega_1^*$.

**Theorem F.5 (self-duality).** Let $V$ be a Euclidean Jordan algebra and $y=\sum_i\lambda_ic_i$ (Theorem F.3).

1. The following are equivalent: (i) $y\in Q$; (ii) all $\lambda_i\ge0$; (iii) $L(y)\ge0$; (iv) $\langle y,x\rangle\ge0$
   for all $x\in Q$. In particular $Q$ is a closed convex cone.
2. The following are equivalent: (i) $y\in\Omega:=\operatorname{int}Q$; (ii) all $\lambda_i>0$; (iii) $L(y)>0$; (iv)
   $\langle y,x\rangle>0$ for all $x\in Q\setminus\{0\}$. Moreover $e\in\Omega$, $\bar\Omega=Q$ and $\Omega^*=\Omega$:
   $\Omega$ is a self-dual open convex cone.
3. If $x,y\in Q$ and $\langle x,y\rangle=0$, then $x\circ y=0$.
4. $x^2=0$ implies $x=0$, and $Q\cap(-Q)=\{0\}$.

*Proof.* We use $\langle y,x^2\rangle=\langle y\circ x,x\rangle=\langle L(y)x,x\rangle$ and
$\langle x^2,e\rangle=\langle x,x\circ e\rangle=|x|^2$; the latter gives (4) together with (1)(ii), since the
eigenvalues of $-y$ are the $-\lambda_i$ by uniqueness.

(1) (i)$\Rightarrow$(iii): if $y=x^2$ with $x=\sum_i\mu_id_i$, then $y=\sum_i\mu_i^2d_i$ and
$L(y)=\sum_i\mu_i^2L(d_i)\ge0$ by Lemma F.4(1). (iii)$\Leftrightarrow$(iv) by the first identity.
(iii)$\Rightarrow$(ii): $0\le\langle L(y)c_i,c_i\rangle=\lambda_i|c_i|^2$. (ii)$\Rightarrow$(i):
$y=(\sum_i\sqrt{\lambda_i}c_i)^2$. By (iv), $Q$ is an intersection of closed half-spaces through $0$.

(2) (ii)$\Rightarrow$(iii): $L(y)=\sum_i\lambda_iL(c_i)\ge(\min_i\lambda_i)\sum_iL(c_i)=(\min_i\lambda_i)I$.
(iii)$\Rightarrow$(iv): a nonzero $x\in Q$ is $z^2$ with $z\ne0$, and $\langle y,x\rangle=\langle L(y)z,z\rangle>0$.
(iv)$\Rightarrow$(ii): $\langle y,c_i\rangle=\lambda_i|c_i|^2$ and $c_i=c_i^2\in Q\setminus\{0\}$. (iii)$\Rightarrow$(i): if
$L(y)\ge mI$ with $m>0$, then $L(y')\ge(m-\|L\|\,|y-y'|)I$, so a neighbourhood of $y$ consists of elements with
$L>0$, which lie in $Q$ by (1). (i)$\Rightarrow$(iii): $L(y)\ge0$ by (1); if $\langle L(y)z,z\rangle=0$ with $z\ne0$, then
$x=z^2\ne0$ and $\langle y,x\rangle=0$; for small $\varepsilon>0$, $y-\varepsilon x\in Q$ but
$\langle y-\varepsilon x,x\rangle=-\varepsilon|x|^2<0$, contradicting (1)(iv). Next, $L(e)=I>0$, and for $y\in Q$,
$y+\varepsilon e\in\Omega$ tends to $y$, so $\bar\Omega=Q$; then $\Omega^*=\{y:\langle y,x\rangle>0\ \forall x\in Q\setminus\{0\}\}=\Omega$
by (iv). $\Omega=\{L>0\}$ is convex.

(3) Let $x=\sum_i\lambda_ic_i$ with $\lambda_i\ge0$. Then $0=\sum_i\lambda_i\langle c_i,y\rangle$ with all terms $\ge0$ by (1).
If $\lambda_i>0$, then $0=\langle y,c_i\rangle=\langle L(y)c_i,c_i\rangle$ with $L(y)\ge0$, so $L(y)c_i=0$. Hence
$x\circ y=\sum_i\lambda_i\,c_i\circ y=0$. $\square$

**Theorem F.6 (homogeneity and the quadratic representation).** Let $V$ be a Euclidean Jordan algebra, $Q$ and
$\Omega$ as in Theorem F.5.

1. For every $b\in V$: $e^{L(b)}Q=Q$, $e^{L(b)}\Omega=\Omega$, and $e^{L(b)}e=\exp b$.
2. Every $y\in\Omega$ equals $\exp(\log y)$. Hence $\Omega$ is homogeneous, and a symmetric cone.
3. For $y\in\Omega$, $P(y)=e^{2L(\log y)}$; so $P(y)Q=Q$, $P(y)\Omega=\Omega$, and $P(y)e=y^2$.
4. $P(y)Q\subseteq Q$ for every $y\in V$.

*Proof.* (1) Let $0\le s<1/\|L(b)\|$ (any $s\ge0$ if $L(b)=0$). The symmetric operator $I-sL(b)$ has eigenvalues
$\ge1-s\|L(b)\|>0$. Let $x\in Q$ and $z=(I-sL(b))^{-1}x=\sum_j\nu_jd_j$ (Theorem F.3). Using $z\circ d_j=\nu_jd_j$ and
$\langle b,d_j\rangle=\langle b,d_j^2\rangle=\langle L(b)d_j,d_j\rangle\le\|L(b)\|\,|d_j|^2$,
$$
0\le\langle x,d_j\rangle=\langle z-s\,b\circ z,d_j\rangle=\nu_j|d_j|^2-s\langle b,z\circ d_j\rangle
=\nu_j\big(|d_j|^2-s\langle b,d_j\rangle\big),
$$
where the left inequality is Theorem F.5(1)(iv). The bracket is $\ge(1-s\|L(b)\|)|d_j|^2>0$, so all $\nu_j\ge0$ and
$z\in Q$. For $t\ge0$, $(I-\frac tnL(b))^{-n}\to e^{tL(b)}$ as $n\to\infty$ (in an orthonormal eigenbasis of $L(b)$ this
is $(1-t\beta/n)^{-n}\to e^{t\beta}$), and $Q$ is closed, so $e^{tL(b)}Q\subseteq Q$. With $-b$ in place of $b$,
$e^{-L(b)}Q\subseteq Q$, so $e^{L(b)}Q=Q$, and the homeomorphism $e^{L(b)}$ maps $\operatorname{int}Q$ onto itself.
Finally $L(b)^re=b^r$ (Jordan powers), so $e^{L(b)}e=\sum_rb^r/r!=\sum_ie^{\beta_i}d_i=\exp b$ for $b=\sum_i\beta_id_i$,
by Theorem F.3.

(2) If $y=\sum_i\lambda_ic_i$ with $\lambda_i>0$, then $\log y=\sum_i(\log\lambda_i)c_i$ and $\exp\log y=y$. So
$y=e^{L(\log y)}e$ with $e^{L(\log y)}\in G(\Omega)$ by (1); self-duality is Theorem F.5(2).

(3) By Lemma F.4(4) for the $c_i$ of $y$, both $P(y)$ and $e^{2L(\log y)}$ act on $V_{ij}$ as
$\lambda_i\lambda_j=e^{\log\lambda_i+\log\lambda_j}$. Use (1), and $P(y)e=2y\circ y-y^2=y^2$.

(4) Let $y=\sum_i\lambda_ic_i$, $c=\sum_{\lambda_i\ge0}c_i$ (an idempotent), $s=2c-e=\sum_i\varepsilon_ic_i$ with
$\varepsilon_i=\pm1$, and $|y|=\sum_i|\lambda_i|c_i$. By Lemma F.4(4), $P(s)$, $P(|y|)$ and $P(y)$ act on $V_{ij}$ as
$\varepsilon_i\varepsilon_j$, $|\lambda_i\lambda_j|$ and $\lambda_i\lambda_j$, so $P(y)=P(s)P(|y|)$. For $\varepsilon>0$,
$|y|+\varepsilon e\in\Omega$, so $P(|y|+\varepsilon e)Q\subseteq Q$ by (3); $P(\cdot)$ is a polynomial and $Q$ is closed, so
$P(|y|)Q\subseteq Q$ ($\varepsilon\to0$). $P(s)$ is a Jordan automorphism (Lemma F.4(3)), so $P(s)x^2=(P(s)x)^2\in Q$. $\square$

Section 13.6 proved $\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|$ in $V_m$ from three properties. They hold in
every Euclidean Jordan algebra for the linear form $\tau(x)=\langle x,e\rangle$, which is $\operatorname{Tr}$ on
$M_{sa}$ and $\operatorname{tr}$ on $V_m$ (with the inner products above):

**Corollary F.7 (Proposition 13.14 in every Euclidean Jordan algebra).** For $x,y$ in a Euclidean Jordan algebra,
$\tau(|P(y)x|)\le\tau(P(y)|x|)$, with equality if $x\in Q$ or $-x\in Q$.

*Proof.* $x_\pm\in Q$, $x=x_+-x_-$, $|x|=x_++x_-$ (Theorem F.3). $A=P(y)x_+$ and $B=P(y)x_-$ lie in $Q$ (Theorem F.6(4)),
and $z=P(y)x=A-B$. Let $z=\sum_i\mu_ic_i$, $c=\sum_{\mu_i>0}c_i$; then $|z|=z\circ(2c-e)$, and
$\tau(a\circ b)=\langle a,b\rangle$. Since $\langle Q,Q\rangle\ge0$ (Theorem F.5(1)) and $c,e-c\in Q$,
$$
\tau(|z|)=\langle A,c\rangle-\langle A,e-c\rangle-\langle B,c\rangle+\langle B,e-c\rangle\le\tau(A)+\tau(B)=\tau(P(y)|x|).
$$
If $x\in Q$, then $x_-=0$ and $z=A\in Q$, so $|z|=z=P(y)|x|$; similarly if $-x\in Q$. $\square$

The code checks Theorem F.5 for $\operatorname{Herm}_4(\mathbb C)$ and for the Lorentz cone $Q\subseteq V_4$ (the
negative pairing is the witness $c_i$ of (iv)$\Rightarrow$(ii)), and Theorem F.6 in $V_4$, where in the coordinates
$(a,u)$ the operators $L(b)$ and $P(y)$ are $5\times5$ symmetric matrices:

```python
import math
import torch
from torch_vn_algebra import SpinFactor

torch.manual_seed(0)
dt = torch.complex128
f64 = torch.float64
close = lambda A, B, tol=1e-10: (A - B).abs().max().item() < tol

G = torch.randn(2, 500, 4, 4, dtype=dt)
X, Y = G[0] @ G[0].mH, G[1] @ G[1].mH                             # squares in Herm_4(C)
assert ((X @ Y).diagonal(dim1=-2, dim2=-1).sum(-1).real > -1e-10).all()   # <Q, Q> >= 0
Z = torch.randn(500, 4, 4, dtype=dt)
Z = Z + Z.mH
w, U = torch.linalg.eigh(Z)
c1 = U[:, :, :1] @ U[:, :, :1].mH                                  # eigenprojection of the smallest eigenvalue
pair = (Z @ c1).diagonal(dim1=-2, dim2=-1).sum(-1).real            # <Z, c1> = lambda_min |c1|^2
assert close(pair, w[:, 0]) and (pair < 0).float().mean() > 0.9   # witnesses of Z not in Q

V = SpinFactor(4)                                                  # <x, y> = tr(x o y) = V.inner
x, y = V.random(2000, cone=True), V.random(2000, cone=True)
assert (V.inner(x, y) > -1e-12).all()
z = V.random(2000)
lam, idem = V.spectral_decomposition(z)                            # (lambda-, lambda+), (c-, c+); |c-|^2 = 1
assert close(V.inner(z, idem[:, 0]), lam[:, 0])
assert ((lam[:, 0] < 0) == ~V.in_cone(z, tol=0.0)).all()

I5 = torch.eye(5, dtype=f64)
Lmat = lambda b: V.product(b.expand(5, 5), I5).T                   # columns b o e_k
Pmat = lambda b: V.quadratic_representation(b.expand(5, 5), I5).T
for _ in range(20):
    b, yi, yg = V.random(), V.random(cone=True), V.random()
    E = torch.linalg.matrix_exp(Lmat(b))
    assert close(E @ V.unit(), V.exp(b))                           # Theorem F.6(1)
    assert V.in_cone((E @ x.T).T, tol=1e-9).all()
    assert close(Pmat(yi), torch.linalg.matrix_exp(2 * Lmat(V.log(yi))))      # Theorem F.6(3)
    lam_g, idem_g = V.spectral_decomposition(yg)
    s = (torch.where(lam_g >= 0, 1.0, -1.0)[:, None] * idem_g).sum(0)
    assert close(Pmat(yg), Pmat(s) @ Pmat(V.apply_function(yg, torch.abs)))  # P(y) = P(s) P(|y|)
    assert V.in_cone((Pmat(yg) @ x.T).T, tol=1e-9).all()                     # Theorem F.6(4)

xs = torch.randn(300, 5, 5, dtype=f64)
xs, ys = xs + xs.mT, torch.randn(300, 5, 5, dtype=f64)
ys = ys + ys.mT                                                    # Corollary F.7 in Herm_5(R), tau = Tr
absm = lambda A: (lambda w, U: U @ torch.diag_embed(w.abs()) @ U.mT)(*torch.linalg.eigh(A))
lhs = absm(ys @ xs @ ys).diagonal(dim1=-2, dim2=-1).sum(-1)
rhs = (ys @ absm(xs) @ ys).diagonal(dim1=-2, dim2=-1).sum(-1)
assert (rhs - lhs > -1e-9 * rhs).all()
```

## F.3 The Koecher–Vinberg theorem

Section F.2 proved that the interior of the cone of squares of a Euclidean Jordan algebra is symmetric. The converse
reconstructs the Jordan product from the cone. Throughout this section $\Omega$ is a self-dual open convex cone in a
Euclidean space $V$ of dimension $N$, and $dy$ is Lebesgue measure in orthonormal coordinates (orthogonal changes
of coordinates have $|\det|=1$).

**Lemma F.8.** Let $\Omega$ be self-dual.

1. $0\notin\Omega$; $\langle x,y\rangle\ge0$ for $x,y\in\bar\Omega$; for each $x\in\Omega$ there is $\delta>0$ with
   $\langle x,y\rangle\ge\delta|y|$ for all $y\in\bar\Omega$.
2. If $x\in\bar\Omega\setminus\Omega$, there is $u\in\bar\Omega\setminus\{0\}$ with $\langle x,u\rangle=0$. If $y\in\Omega$,
   $u\in\bar\Omega$ and $t\ge0$, then $y+tu\in\Omega$.
3. $\Omega$ is the interior of $\bar\Omega$.
4. $G(\Omega)$ is a closed subgroup of $GL(V)$, and $g\in G(\Omega)$ implies $g^T\in G(\Omega)$.

*Proof.* (1) If $0\in\Omega$, then $\Omega$ contains a ball around $0$, hence $\Omega=V$; but $V^*=\emptyset$. For
$x,y\in\Omega$, $\langle x,y\rangle>0$ because $y\in\Omega^*$; by continuity $\langle x,y\rangle\ge0$ on $\bar\Omega$. If the
closed ball $\bar B(x,\delta)$ lies in $\Omega$ and $y\in\bar\Omega\setminus\{0\}$, then
$\langle x-\delta y/|y|,y\rangle\ge0$. (2) $x\notin\Omega=\Omega^*$ gives $u\in\bar\Omega\setminus\{0\}$ with
$\langle u,x\rangle\le0$, and $\langle u,x\rangle\ge0$ by (1). For $x'\in\bar\Omega\setminus\{0\}$,
$\langle x',y+tu\rangle\ge\langle x',y\rangle>0$, so $y+tu\in\Omega^*=\Omega$. (3) If $\bar B(y,\varepsilon)\subseteq\bar\Omega$,
then for $x\in\bar\Omega\setminus\{0\}$, $\langle x,y-\varepsilon x/|x|\rangle\ge0$, so $\langle x,y\rangle\ge\varepsilon|x|>0$ and
$y\in\Omega^*=\Omega$. (4) For $g\in G(\Omega)$, $g$ maps $\bar\Omega\setminus\{0\}$ onto itself, so for $y\in\Omega$ and
$x\in\bar\Omega\setminus\{0\}$, $\langle x,g^Ty\rangle=\langle gx,y\rangle>0$: $g^T\Omega\subseteq\Omega$. The same for $g^{-1}$
gives $g^T\Omega=\Omega$. If $g_k\in G(\Omega)$ converge to $g\in GL(V)$, then $g\Omega\subseteq\bar\Omega$ is open, so
$g\Omega\subseteq\Omega$ by (3); likewise $g^{-1}=\lim g_k^{-1}$ maps $\Omega$ into $\Omega$, so $g\Omega=\Omega$. $\square$

The **characteristic function** of $\Omega$ is $\varphi(x)=\int_{\Omega^*}e^{-\langle x,y\rangle}dy=\int_\Omega e^{-\langle x,y\rangle}dy$.
Put $f=\log\varphi$, let $H(x)$ be the symmetric operator with $\langle H(x)u,v\rangle=D^2f(x)[u,v]$, and
$x^*=-\nabla f(x)$, i.e. $Df(x)[u]=-\langle x^*,u\rangle$.

**Lemma F.9 (characteristic function).** Let $\Omega$ be self-dual.

1. $\varphi$ is finite, positive and $C^\infty$ on $\Omega$, and may be differentiated under the integral sign.
2. $\varphi(tx)=t^{-N}\varphi(x)$ for $t>0$, and $\varphi(gx)=|\det g|^{-1}\varphi(x)$ for $g\in G(\Omega)$.
3. $H(x)>0$ and $x^*\in\Omega$ for $x\in\Omega$.
4. $H(x)x=x^*$ and $D^3f(x)[x,u,v]=-2\langle H(x)u,v\rangle$; for $g\in G(\Omega)$, $(gx)^*=g^{-T}x^*$ and
   $H(gx)=g^{-T}H(x)g^{-1}$.
5. $\varphi(x)\ge c|x|^{-N}$ for a constant $c>0$, and $\varphi(x_k)\to\infty$ if $x_k\in\Omega$ converge to a point of
   $\bar\Omega\setminus\Omega$.

*Proof.* (1) Let $x\in\Omega$ and $\delta$ as in Lemma F.8(1). For $|x'-x|<\delta/2$ and $y\in\Omega$,
$\langle x',y\rangle\ge\frac\delta2|y|$, so the $k$-th derivatives of the integrand in $x'$,
$(-1)^k\prod_{i}\langle u_i,y\rangle e^{-\langle x',y\rangle}$, are bounded by $\prod_i|u_i|\,|y|^ke^{-\delta|y|/2}$, which is
integrable ($|y|\ge N^{-1/2}\sum_i|y_i|$ and Fubini). Dominated convergence allows differentiation under the integral
sign. (2) Substitute $y=y'/t$, which maps $\Omega$ onto $\Omega$; and $y'=g^Ty$, which maps $\Omega$ onto $\Omega$ by Lemma
F.8(4), with $\langle gx,y\rangle=\langle x,y'\rangle$ and $dy=|\det g|^{-1}dy'$.
(3) Let $\mu_x$ be the probability measure $e^{-\langle x,y\rangle}dy/\varphi(x)$ on $\Omega$. By (1),
$Df(x)[u]=-\mathbb E\langle u,y\rangle$ and $D^2f(x)[u,u]=\mathbb E\langle u,y\rangle^2-(\mathbb E\langle u,y\rangle)^2$, the
variance of $\langle u,y\rangle$. It vanishes only if $\langle u,y\rangle$ is constant $\mu_x$-almost everywhere, i.e. on
almost all of the open set $\Omega$, which forces $u=0$. And $x^*=\mathbb E\,y$, so $\langle z,x^*\rangle=\mathbb E\langle z,y\rangle>0$
for $z\in\bar\Omega\setminus\{0\}$: $x^*\in\Omega^*=\Omega$.
(4) By (2), $f(tx)=f(x)-N\log t$, so $Df(x)[x]=-N$ on $\Omega$. Differentiating in $x$ along $u$:
$D^2f(x)[x,u]+Df(x)[u]=0$, i.e. $\langle H(x)x,u\rangle=\langle x^*,u\rangle$; once more along $v$:
$D^3f(x)[x,u,v]+2D^2f(x)[u,v]=0$. Differentiating $f(gx)=f(x)-\log|\det g|$ once and twice gives
$\langle(gx)^*,gu\rangle=\langle x^*,u\rangle$ and $\langle H(gx)gu,gv\rangle=\langle H(x)u,v\rangle$.
(5) $\langle x,y\rangle\le|x||y|\le1$ for $|y|\le1/|x|$, so $\varphi(x)\ge e^{-1}\operatorname{vol}(\Omega\cap B(0,1/|x|))=c|x|^{-N}$
with $c=e^{-1}\operatorname{vol}(\Omega\cap B(0,1))>0$. Let $x_k\to x_0\in\bar\Omega\setminus\Omega$. By Fatou's lemma,
$\liminf\varphi(x_k)\ge\int_\Omega e^{-\langle x_0,y\rangle}dy$. Take $u$ as in Lemma F.8(2), a ball $B(y_1,r)\subseteq\Omega$
and the disjoint balls $B_l=B(y_1+2lru/|u|,r)$, $l\ge0$; they lie in $\Omega$ by Lemma F.8(2), and on $B_l$,
$\langle x_0,y\rangle=\langle x_0,y-2lru/|u|\rangle\le|x_0|(|y_1|+r)$. So the integral is infinite. $\square$

**Lemma F.10 (one-parameter groups).** Let $G$ be a closed subgroup of $GL(V)$ and
$\mathfrak g=\{Z:e^{tZ}\in G\text{ for all }t\in\mathbb R\}$.

1. If $c:[0,s_0)\to G$ satisfies $c(s)=I+sZ+o(s)$ as $s\to0^+$, then $Z\in\mathfrak g$.
2. $\mathfrak g$ is a Lie algebra: a linear subspace closed under $[Z,W]=ZW-WZ$.

*Proof.* (1) Fix $t>0$ and put $A_n=c(t/n)$, $B_n=e^{tZ/n}$. Then $\|A_n-B_n\|=o(1/n)$ and
$\|A_n\|,\|B_n\|\le1+K/n$ for some $K$ and large $n$, so
$$
\|A_n^n-e^{tZ}\|=\Big\|\sum_{r=0}^{n-1}A_n^r(A_n-B_n)B_n^{n-1-r}\Big\|\le n\,e^K\|A_n-B_n\|\to0 .
$$
$A_n^n\in G$ and $G$ is closed in $GL(V)$, so $e^{tZ}\in G$; also $e^{-tZ}=(e^{tZ})^{-1}\in G$. (2) Multiples are clear.
Apply (1) to $c(s)=e^{sZ}e^{sW}=I+s(Z+W)+O(s^2)$, and to $c(s)=e^{aZ}e^{aW}e^{-aZ}e^{-aW}$ with $a=\sqrt s$: expanding each
factor to second order, the terms of order $a$ cancel and those of order $a^2$ add up to $ZW-WZ$, so
$c(s)=I+s[Z,W]+O(s^{3/2})$. $\square$

**Proposition F.11 (the base point).** Let $\Omega$ be a symmetric cone, $G=G(\Omega)$, $\mathfrak g$ as in Lemma F.10.

1. $x^{**}=x$ and $H(x^*)=H(x)^{-1}$ for $x\in\Omega$.
2. There is $e\in\Omega$ with $e^*=e$, and every such $e$ has $H(e)=I$.
3. For $x\in\Omega$ and any $g\in G$ with $ge=x$: $H(x)^{-1}=gg^T\in G$.
4. If $Z\in\mathfrak g$ and $Ze=0$, then $Z^T=-Z$.

*Proof.* (1) The map $x\mapsto x^*$ is smooth from $\Omega$ to $\Omega$ (Lemma F.9(3)), with derivative $-H(x)$. Since
$g^{-T}\in G$ (Lemma F.8(4)), Lemma F.9(2,4) gives
$f(gx)+f((gx)^*)=f(x)-\log|\det g|+f(g^{-T}x^*)=f(x)+f(x^*)$, so $F(x)=f(x)+f(x^*)$ is constant on the orbit
$G x=\Omega$. Hence $0=DF(x)[u]=-\langle x^*,u\rangle+\langle x^{**},H(x)u\rangle$, i.e. $H(x)x^{**}=x^*=H(x)x$
(Lemma F.9(4)), and $x^{**}=x$ as $H(x)$ is invertible. Differentiating $x^{**}=x$ gives $H(x^*)H(x)=I$.
(2) $F_2(x)=f(x)+\frac12|x|^2$ satisfies $F_2(x)\ge\log c-N\log|x|+\frac12|x|^2$ (Lemma F.9(5)), which tends to $\infty$ as
$|x|\to\infty$. So $S=\{x\in\Omega:F_2(x)\le F_2(x_1)\}$ ($x_1\in\Omega$) is bounded; it is closed in $V$, because a limit
point in $\bar\Omega\setminus\Omega$ would make $f$ unbounded on $S$ (Lemma F.9(5)). The continuous $F_2$ attains its minimum
on the compact set $S$ at some $e$ in the open set $\Omega$, where $0=\nabla F_2(e)=-e^*+e$. Then $H(e)=H(e^*)=H(e)^{-1}$
by (1), so $H(e)^2=I$, and $H(e)=I$ since $H(e)>0$.
(3) $H(x)=H(ge)=g^{-T}g^{-1}$ by Lemma F.9(4) and (2), and $gg^T\in G$ by Lemma F.8(4).
(4) $e^{tZ}\in G$ fixes $e$, so $I=H(e)=e^{-tZ^T}H(e)e^{-tZ}=e^{-tZ^T}e^{-tZ}$; differentiate at $t=0$. $\square$

**Theorem F.12 (Koecher–Vinberg).** Let $\Omega$ be a nonempty open convex cone in a Euclidean space $V$. Then $\Omega$ is
a symmetric cone if and only if there is a product $\circ$ on $V$ such that $(V,\circ,\langle\cdot,\cdot\rangle)$ is a
Euclidean Jordan algebra and $\Omega$ is the interior of its cone of squares. For a symmetric cone one may take as unit any
$e\in\Omega$ with $e^*=e$, and the product
$$
\langle u\circ v,w\rangle=-\tfrac12\,D^3\log\varphi(e)[u,v,w].
$$

*Proof.* "If" is Theorems F.5(2) and F.6(2). Conversely let $\Omega$ be symmetric, $e$ as in Proposition F.11(2), and
define $L(u)$ by $\langle L(u)v,w\rangle=-\frac12D^3f(e)[u,v,w]$ and $u\circ v=L(u)v$.

*Step 1: commutativity, symmetry, unit.* $D^3f(e)$ is symmetric in its three arguments, so $u\circ v=v\circ u$ and
$L(u)^T=L(u)$. By Lemma F.9(4) and $H(e)=I$, $\langle e\circ v,w\rangle=\langle H(e)v,w\rangle=\langle v,w\rangle$: $e$ is a unit.

*Step 2: $L(u)\in\mathfrak g$.* $P(x)=H(x)^{-1}$ lies in $G$ (Proposition F.11(3)), is smooth, and $P(e)=I$. Its
derivative is $DP(e)[u]=-H(e)^{-1}DH(e)[u]H(e)^{-1}=-DH(e)[u]=2L(u)$, since $\langle DH(e)[u]v,w\rangle=D^3f(e)[u,v,w]$.
Lemma F.10(1) for $c(s)=P(e+su)$ gives $2L(u)\in\mathfrak g$.

*Step 3: uniqueness.* If $Y\in\mathfrak g$ is symmetric, then $Y=L(Ye)$: the symmetric $Y-L(Ye)\in\mathfrak g$ kills $e$
($L(w)e=w\circ e=w$), so it is also antisymmetric by Proposition F.11(4), hence $0$.

*Step 4: inner derivations.* Let $D=[L(x),L(y)]$. It lies in $\mathfrak g$ (Lemma F.10(2)), $D^T=-D$, and
$De=x\circ y-y\circ x=0$. For $z\in V$, $[D,L(z)]\in\mathfrak g$ is symmetric and $[D,L(z)]e=Dz-L(z)De=Dz$, so
$[D,L(z)]=L(Dz)$ by Step 3; applied to $w$: $D(z\circ w)=Dz\circ w+z\circ Dw$.

*Step 5: the Jordan identity.* For $D=[L(y),L(x)]$, $Dx=y\circ x^2-x\circ(x\circ y)$ and
$D(x^2)=y\circ(x\circ x^2)-x\circ(y\circ x^2)$. Step 4 gives $D(x^2)=2x\circ Dx$, i.e.
$y\circ(x\circ x^2)=3x\circ(x^2\circ y)-2x\circ(x\circ(x\circ y))$ for all $y$:
$$
L(x\circ x^2)=3L(x)L(x^2)-2L(x)^3 .
$$
All $L$'s are symmetric, so transposing gives $L(x\circ x^2)=3L(x^2)L(x)-2L(x)^3$, and $[L(x),L(x^2)]=0$. So $V$ is a
Euclidean Jordan algebra; let $\Omega_J$ be the interior of its cone of squares.

*Step 6: $\Omega=\Omega_J$.* For $b\in V$, $e^{L(b)}\in G$ by Step 2, so $\exp b=e^{L(b)}e\in\Omega$ (Theorem F.6(1)). By
Theorem F.6(2), $\Omega_J\subseteq\Omega$. Both cones are self-dual for $\langle\cdot,\cdot\rangle$ (Theorem F.5(2)), so
$\Omega=\Omega^*\subseteq\Omega_J^*=\Omega_J$. $\square$

The proof follows the structure of J. Faraut, A. Korányi, *Analysis on Symmetric Cones* (Oxford 1994), Chapters I–III,
but replaces their use of Lie group theory by Lemma F.10 and Proposition F.11(3), which produce the needed one-parameter
groups directly. M. Koecher (1957) and E. B. Vinberg (1960) proved the theorem.

**Example F.13 (the cone of semidefinite programming).** Let $V=\operatorname{Herm}_n(\mathbb R)$, $\langle x,y\rangle=\operatorname{Tr}(xy)$,
$N=n(n+1)/2$, and $\Omega$ the positive definite matrices, symmetric by Theorems F.5 and F.6. Then
$\varphi(x)=\varphi(1)\det(x)^{-(n+1)/2}$, $e=\kappa1$ with $\kappa=\sqrt{(n+1)/2}$, and the product of Theorem F.12 is
$u\circ_\kappa v=\frac1{2\kappa}(uv+vu)$, so $x\mapsto x/\kappa$ is an isomorphism onto the usual Jordan algebra.

*Proof.* For $a\in\Omega$, $P(a)x=axa$ lies in $G(\Omega)$ (Theorem F.6(3)). In an orthonormal eigenbasis of $a$ the
matrices $E_{ii}$ and $(E_{ij}+E_{ji})/\sqrt2$ ($i<j$) form an orthonormal basis of eigenvectors of $P(a)$ with
eigenvalues $a_ia_j$, so $\det P(a)=\prod_{i\le j}a_ia_j=\det(a)^{n+1}$ (each $a_i$ occurs $n+1$ times). Lemma F.9(2) gives
$\varphi(a^2)=\varphi(P(a)1)=\det(a)^{-(n+1)}\varphi(1)$. With $D\log\det(x)[u]=\operatorname{Tr}(x^{-1}u)$:
$x^*=\kappa^2x^{-1}$, $e=\kappa1$, $D^2f(x)[v,w]=\kappa^2\operatorname{Tr}(x^{-1}vx^{-1}w)$, and
$D^3f(e)[u,v,w]=-\kappa^2\kappa^{-3}\operatorname{Tr}(uvw+vuw)=-\frac2\kappa\operatorname{Tr}((u\circ v)w)$. $\square$

The code computes $D^3f$ at $e$ by automatic differentiation of $f(x)=-\frac{n+1}2\log\det x$ in orthonormal
coordinates of $\operatorname{Herm}_3(\mathbb R)$ and checks the claims of Proposition F.11 and Theorem F.12:

```python
from torch.autograd.functional import jacobian, hessian

n = 3
pairs = [(i, j) for i in range(n) for j in range(i, n)]
Nd = len(pairs)                                                    # N = 6
Bs = torch.zeros(Nd, n, n, dtype=f64)                              # orthonormal basis for Tr(xy)
for k, (i, j) in enumerate(pairs):
    Bs[k, i, j] = Bs[k, j, i] = 1.0 if i == j else 1 / math.sqrt(2)
mat = lambda v: torch.einsum('k,kij->ij', v, Bs)
vec = lambda A: torch.einsum('kij,ij->k', Bs, A)
f = lambda v: -(n + 1) / 2 * torch.logdet(mat(v))                  # log phi up to a constant (Example F.13)
kappa = math.sqrt((n + 1) / 2)
e = vec(kappa * torch.eye(n, dtype=f64))
assert close(-jacobian(f, e), e)                                   # e* = e
assert close(hessian(f, e), torch.eye(Nd, dtype=f64))             # H(e) = I (Proposition F.11(2))
T = jacobian(lambda v: hessian(f, v, create_graph=True), e)        # T[a, b, c] = D^3 f(e)[B_a, B_b, B_c]
prod = lambda u, v: torch.einsum('a,b,abc->c', u, v, -0.5 * T)     # Theorem F.12
u, v = torch.randn(2, Nd, dtype=f64)
assert close(prod(e, v), v)                                        # unit
assert close(prod(u, v), vec(mat(u) @ mat(v) + mat(v) @ mat(u)) / (2 * kappa))
u2 = prod(u, u)
assert close(prod(u, prod(u2, v)), prod(u2, prod(u, v)))           # Jordan identity
```

## F.4 Order isomorphisms: Kadison's theorem in finite dimensions

Let $\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$ and $\mathcal N=\bigoplus_dM_{k'_d}(\mathbb C)$, as block-diagonal matrices;
$x\ge0$ means $x=x^*$ with nonnegative eigenvalues (equivalently, every block is positive semidefinite), and
$[0,1]=\{a\in\mathcal M_{sa}:0\le a\le1\}$. A **projection** is $p=p^*=p^2$. A linear map $\Phi$ is **positive** if
$x\ge0$ implies $\Phi(x)\ge0$; a **Jordan $*$-isomorphism** is a linear bijection with $\Phi(x^*)=\Phi(x)^*$ and
$\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$ (Section 13.5).

**Lemma F.14.** In $\mathcal M$:

1. The extreme points of the convex set $[0,1]$ are exactly the projections.
2. For projections $p,q$: $pq=0$ iff $p+q\le1$ iff $p\circ q=0$.
3. For projections $p,q$: $p\le q$ iff $pq=p$ iff $p\circ q=p$. The minimal nonzero projections are the rank-one
   projections (in some block).
4. A projection $p$ is central iff $pxp+(1-p)x(1-p)=x$ for all $x\in\mathcal M$. The central projections are the sums of
   the block units $P_c$, and the minimal nonzero ones are the $P_c$.

*Proof.* (1) Let $p$ be a projection and $p=\frac12(a+b)$ with $a,b\in[0,1]$. For $\xi\in\ker p$,
$0=\langle\xi,a\xi\rangle+\langle\xi,b\xi\rangle$ with both terms $\ge0$, so $\|a^{1/2}\xi\|^2=0$ and $a\xi=0$; for $\xi$ in
the range of $p$, the same with $1-a,1-b$ gives $a\xi=\xi$. So $a=p$. Conversely, if $a\in[0,1]$ is not a projection,
it has an eigenvalue $\lambda\in(0,1)$; with the spectral projection $q_\lambda\in\mathcal M$ (a polynomial in $a$) and
$\varepsilon=\min(\lambda,1-\lambda)$, $a\pm\varepsilon q_\lambda\in[0,1]$ and $a$ is their midpoint.
(2) If $pq=0$ then $qp=(pq)^*=0$, $p+q$ is a projection, and $1-(p+q)=(1-p-q)^2\ge0$; also $p\circ q=0$. If
$p+q\le1$, then for every $\xi$, $\|qp\xi\|^2=\langle p\xi,qp\xi\rangle\le\langle p\xi,(1-p)p\xi\rangle=0$, so $qp=0=pq$. If
$p\circ q=0$, then $pq=-qp$, so $pq=p(pq)=-pqp=(qp)p=qp$, and $pq=-pq=0$.
(3) $pq=p$ implies $qp=p$, $p\circ q=p$, and $q-p=(q-p)^2\ge0$. If $p\circ q=p$, multiplying $pq+qp=2p$ by $p$ from
the left and from the right gives $pq+pqp=2p=pqp+qp$, so $pq=qp=p$. If $p\le q$, then $p+(1-q)\le1$ and (2) gives
$p(1-q)=0$. A projection of rank $\ge2$ dominates the projection onto a line in its range (inside one block).
(4) $x=pxp+px(1-p)+(1-p)xp+(1-p)x(1-p)$, and $px(1-p)+(1-p)xp=0$ iff both terms vanish (multiply by $p$ and $1-p$),
iff $px=pxp=xp$. A central element commutes with all matrix units of each block, so it is a scalar $z_c$ on block
$c$, and $z_c\in\{0,1\}$ for a projection. $\square$

**Theorem F.15 (Kadison, finite dimensions).** Let $\Phi:\mathcal M\to\mathcal N$ be a linear bijection with
$\Phi(1)=1$ such that $\Phi$ and $\Phi^{-1}$ are positive. Then $\Phi$ is a Jordan $*$-isomorphism. Conversely, every
Jordan $*$-isomorphism is unital, and it and its inverse are positive.

*Proof.* Every $x\in\mathcal M_{sa}$ is $x_+-x_-$ with $x_\pm\ge0$, so $\Phi(x)$ is self-adjoint; writing $x=h+ik$ with
$h,k$ self-adjoint, $\Phi(x^*)=\Phi(h)-i\Phi(k)=\Phi(x)^*$. The same holds for $\Phi^{-1}$, so $\Phi$ is a real-linear
bijection $\mathcal M_{sa}\to\mathcal N_{sa}$ with $x\ge0\iff\Phi(x)\ge0$. As $\Phi(1)=1$, $\Phi$ maps $[0,1]$ onto $[0,1]$.
A linear bijection between convex sets maps extreme points onto extreme points (if $\Phi(a)=\frac12(b+b')$, then
$a=\frac12(\Phi^{-1}b+\Phi^{-1}b')$), so by Lemma F.14(1) $\Phi$ maps projections onto projections, and by Lemma
F.14(2), $pq=0\iff p+q\le1\iff\Phi(p)+\Phi(q)\le1\iff\Phi(p)\Phi(q)=0$. For $x=\sum_i\lambda_ip_i\in\mathcal M_{sa}$ (spectral
decomposition, $p_ip_j=0$ for $i\ne j$), $\Phi(x)^2=\sum_i\lambda_i^2\Phi(p_i)=\Phi(x^2)$. Polarisation,
$x\circ y=\frac12((x+y)^2-x^2-y^2)$, gives $\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$ on $\mathcal M_{sa}$, and complex
bilinearity on $\mathcal M=\mathcal M_{sa}+i\mathcal M_{sa}$.

Conversely, $\Phi(1)\circ\Phi(y)=\Phi(y)$ for all $y$ and $\Phi$ is onto, so $\Phi(1)$ is a unit for $\circ$, and
$\Phi(1)=\Phi(1)\circ1=1$. Every $x\ge0$ is $h^2$ with $h=\sqrt x$ self-adjoint (Theorem 1.6), and
$\Phi(h^2)=\Phi(h)^2\ge0$ since $\Phi(h)$ is self-adjoint; $\Phi^{-1}$ is a Jordan $*$-isomorphism too. $\square$

Kadison proved the theorem for all unital C\*-algebras (R. V. Kadison, *A generalized Schwarz inequality and algebraic
invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503); that version is not used in the course.

```python
def blocks_rand(ks, psd=False):
    out = []
    for k in ks:
        A = torch.randn(k, k, dtype=dt)
        out.append(A @ A.mH if psd else A + A.mH)
    return out

ks = (2, 2, 3)                                                     # M = M_2 (+) M_2 (+) M_3
Us = [torch.linalg.qr(torch.randn(k, k, dtype=dt))[0] for k in ks]
def Phi(x):                                                        # swap the M_2 blocks, transpose two blocks
    x1, x2, x3 = x
    return [Us[0] @ x2.mT @ Us[0].mH, Us[1] @ x1 @ Us[1].mH, Us[2] @ x3.mT @ Us[2].mH]
def Phi_inv(y):
    y1, y2, y3 = y
    return [Us[1].mH @ y2 @ Us[1], (Us[0].mH @ y1 @ Us[0]).mT, (Us[2].mH @ y3 @ Us[2]).mT]
is_psd = lambda x: all(torch.linalg.eigvalsh(b)[0] > -1e-10 for b in x)
one = [torch.eye(k, dtype=dt) for k in ks]
assert all(close(a, b) for a, b in zip(Phi(one), one))            # unital
for _ in range(50):
    p = blocks_rand(ks, psd=True)
    assert is_psd(Phi(p)) and is_psd(Phi_inv(p))                   # positive in both directions
x, y = blocks_rand(ks), blocks_rand(ks)
assert all(close(a, b) for a, b in zip(Phi_inv(Phi(x)), x))
jordan = lambda a, b: [0.5 * (s @ t + t @ s) for s, t in zip(a, b)]
assert all(close(a, b) for a, b in zip(Phi(jordan(x, y)), jordan(Phi(x), Phi(y))))   # Theorem F.15
assert not close(Phi([s @ t for s, t in zip(x, y)])[0], (Phi(x)[0] @ Phi(y)[0]))      # not multiplicative
```

## F.5 Jordan $*$-isomorphisms of $\bigoplus_cM_{k_c}$

**Lemma F.16.** Let $\Phi:\mathcal M\to\mathcal N$ be a Jordan $*$-isomorphism. Then $\Phi(1)=1$,
$\Phi(xyx)=\Phi(x)\Phi(y)\Phi(x)$, and $\Phi$ maps the projections of $\mathcal M$ bijectively onto those of $\mathcal N$,
preserving and reflecting $pq=0$, $p\le q$, minimality, centrality, and minimality among central projections.

*Proof.* $\Phi(1)=1$ as in Theorem F.15. In every associative algebra $xyx=2x\circ(x\circ y)-(x\circ x)\circ y$, which
gives the second claim. $p=p^*=p\circ p$ is preserved, and so are the characterisations of Lemma F.14(2,3,4), which
use only $\circ$, the triple product $xyx$, $1$, and quantifiers over the whole algebra; $\Phi^{-1}$ is a Jordan
$*$-isomorphism too. $\square$

**Theorem F.17 (Jordan $*$-automorphisms of $M_n$).** Every Jordan $*$-isomorphism $\Phi:M_n(\mathbb C)\to M_n(\mathbb C)$
is $x\mapsto UxU^*$ or $x\mapsto Ux^TU^*$ for a unitary $U$. For $n\ge2$ the first maps are multiplicative and the second
antimultiplicative, so no map is of both forms.

*Proof.* For $n=1$, $\Phi=\mathrm{id}$. Let $n\ge2$ and $E_{ij}$ the matrix units. By Lemma F.16, $\Phi(E_{ii})$ are
pairwise orthogonal rank-one projections, $\Phi(E_{ii})=\xi_i\xi_i^*$ with orthonormal $\xi_i$; with the unitary $W$,
$We_i=\xi_i$, the Jordan $*$-isomorphism $\Phi_1(x)=W^*\Phi(x)W$ satisfies $\Phi_1(E_{ii})=E_{ii}$.

Let $i\ne j$. The space $S_{ij}=\operatorname{span}\{E_{ij},E_{ji}\}$ is $\{x:pxp=x,\ E_{ii}xE_{ii}=E_{jj}xE_{jj}=0\}$ with
$p=E_{ii}+E_{jj}$, so $\Phi_1(S_{ij})\subseteq S_{ij}$ by Lemma F.16. Write $\Phi_1(E_{ij})=\alpha E_{ij}+\beta E_{ji}$. From
$E_{ij}\circ E_{ij}=0$: $0=\Phi_1(E_{ij})^2=\alpha\beta(E_{ii}+E_{jj})$, so exactly one of $\alpha,\beta$ is nonzero
($\Phi_1$ is injective). Let $\varepsilon(i,j)=+$ if $\Phi_1(E_{ij})\in\mathbb CE_{ij}$ and $-$ if $\Phi_1(E_{ij})\in\mathbb CE_{ji}$;
$\varepsilon(j,i)=\varepsilon(i,j)$ since $\Phi_1(E_{ji})=\Phi_1(E_{ij})^*$. For distinct $i,j,k$,
$E_{ij}\circ E_{jk}=\frac12E_{ik}\ne0$; if $\varepsilon(i,j)=+$ and $\varepsilon(j,k)=-$, then
$\Phi_1(E_{ij})\circ\Phi_1(E_{jk})\in\mathbb C\,E_{ij}\circ E_{kj}=\{0\}$, a contradiction; likewise for $(-,+)$. So
$\varepsilon(i,j)=\varepsilon(j,k)$ whenever $i,j,k$ are distinct, and $\varepsilon$ is constant: two pairs either share an
index or are linked through a third pair. If $\varepsilon=-$, replace $\Phi_1$ by $\Phi_1\circ T$ ($T$ the transpose, a
Jordan $*$-automorphism by Proposition 13.9), which has $\varepsilon=+$ and still fixes the $E_{ii}$.

Now $\Phi_1(E_{ij})=\alpha_{ij}E_{ij}$ for all $i,j$, with $\alpha_{ii}=1$ and $\alpha_{ji}=\bar\alpha_{ij}$. From
$E_{ij}\circ E_{ji}=\frac12(E_{ii}+E_{jj})$: $|\alpha_{ij}|^2=1$; from $E_{ij}\circ E_{jk}=\frac12E_{ik}$ ($i,j,k$ distinct):
$\alpha_{ij}\alpha_{jk}=\alpha_{ik}$. With $u_i=\alpha_{i1}$ this gives $\alpha_{ij}=\alpha_{i1}\alpha_{1j}=u_i\bar u_j$ (also for
$i=j$, $i=1$ or $j=1$). So $\Phi_1(x)=DxD^*$ with $D=\operatorname{diag}(u_i)$, and $\Phi(x)=WDxD^*W^*$ or
$\Phi(x)=WDx^TD^*W^*$. If both forms agreed, $T$ would be multiplicative, but $(E_{12}E_{21})^T=E_{11}\ne E_{22}=E_{12}^TE_{21}^T$.
$\square$

This proves the result of N. Jacobson, C. E. Rickart (Trans. Amer. Math. Soc. 69 (1950) 479–502) and I. N.
Herstein (Trans. Amer. Math. Soc. 81 (1956) 331–341) quoted in Section 13.5, in the finite-dimensional case, for all
$n$; Proposition 13.11 is the case $n=2$. A Jordan automorphism of $\operatorname{Herm}_n(\mathbb C)$ extends to a
Jordan $*$-automorphism of $M_n$ by $\Phi(x+iy)=\Phi(x)+i\Phi(y)$, so it is also of one of the two forms.

**Theorem F.18 (direct sums).** Let $\Phi:\bigoplus_{c=1}^CM_{k_c}\to\bigoplus_{d=1}^{C'}M_{k'_d}$ be a Jordan
$*$-isomorphism. Then $C'=C$, and there are a bijection $\sigma$ with $k'_{\sigma(c)}=k_c$, unitaries $U_c\in M_{k_c}$ and
a set $S$ of blocks such that
$$
\Phi(x)_{\sigma(c)}=U_cx_cU_c^*\ (c\notin S),\qquad \Phi(x)_{\sigma(c)}=U_cx_c^TU_c^*\ (c\in S).
$$
Conversely every such map is a Jordan $*$-isomorphism, and $S$ is unique up to blocks with $k_c=1$.

*Proof.* By Lemma F.16 and Lemma F.14(4), $\Phi(P_c)=Q_{\sigma(c)}$ for a bijection $\sigma$ between the minimal central
projections, so $C=C'$. If $x=P_cxP_c$, then $\Phi(x)=Q_{\sigma(c)}\Phi(x)Q_{\sigma(c)}$, so $\Phi$ maps block $c$ into block
$\sigma(c)$, and $\Phi^{-1}$ maps block $\sigma(c)$ into block $c$; hence $k_c^2=k'^2_{\sigma(c)}$, and the restriction is a
Jordan $*$-isomorphism $M_{k_c}\to M_{k_c}$, of one of the two forms by Theorem F.17. Since $x=\sum_cP_cxP_c$, this
determines $\Phi$. The converse holds blockwise (Proposition 13.9), and uniqueness of $S$ follows from Theorem F.17. $\square$

**Corollary F.19 (unital order isomorphisms).** The unital linear bijections $\Phi:\mathcal M\to\mathcal N$ such that
$\Phi$ and $\Phi^{-1}$ are positive are exactly the maps of Theorem F.18. Such a $\Phi$ is $2$-positive iff it is
completely positive iff $S$ contains no block with $k_c\ge2$, iff $\Phi$ is a $*$-isomorphism.

*Proof.* The first sentence is Theorems F.15 and F.18. If $S$ has no block with $k_c\ge2$, every block map is
$x\mapsto U_cxU_c^*$ (for $k_c=1$, $x^T=x$), so $\Phi$ is completely positive (Proposition 6.4) and multiplicative. If
$c\in S$ with $k_c\ge2$ and $\Phi$ is $2$-positive, then so is the block map $\Phi_{\sigma(c)c}(a)=U_ca^TU_c^*$
(Proposition 6.4), and so is $a\mapsto U_c^*\Phi_{\sigma(c)c}(a)U_c=a^T$, since $[A_{ij}]\ge0$ implies
$[U^*A_{ij}U]=(1\otimes U)^*[A_{ij}](1\otimes U)\ge0$; this contradicts Proposition 6.1. Such a $\Phi$ is not multiplicative
either, by Theorem F.17. $\square$

So a unital order isomorphism is determined, up to the unitaries, by which equal blocks it exchanges and which blocks it
transposes; only the transposes distinguish it from a reversible quantum operation. The code builds the Jordan
automorphism $x\mapsto Ux^TU^*$ of $M_3$ and recovers $U$ (up to a phase) by the steps of the proof of Theorem F.17:

```python
k = 3
E = lambda i, j: torch.zeros(k, k, dtype=dt).index_put_((torch.tensor(i), torch.tensor(j)), torch.tensor(1.0, dtype=dt))
U0 = torch.linalg.qr(torch.randn(k, k, dtype=dt))[0]
Psi = lambda A: U0 @ A.mT @ U0.mH                                  # the unknown Jordan *-automorphism

W = torch.stack([torch.linalg.eigh(Psi(E(i, i)))[1][:, -1] for i in range(k)], dim=1)   # Psi(E_ii) = xi_i xi_i^*
Psi1 = lambda A: W.mH @ Psi(A) @ W
assert all(close(Psi1(E(i, i)), E(i, i)) for i in range(k))
transposed = Psi1(E(0, 1))[1, 0].abs().item() > 0.5               # epsilon = - : Psi1(E_01) is a multiple of E_10
Psi2 = (lambda A: Psi1(A.mT)) if transposed else Psi1
alpha = torch.stack([Psi2(E(i, 0))[i, 0] for i in range(k)])       # u_i = alpha_{i1}
for i in range(k):
    for j in range(k):
        assert close(Psi2(E(i, j)), alpha[i] * alpha[j].conj() * E(i, j))
Ufound = W @ torch.diag(alpha)
X = torch.randn(k, k, dtype=dt)
assert transposed and close(Psi(X), Ufound @ X.mT @ Ufound.mH)     # Theorem F.17, second form
phase = (U0.mH @ Ufound)[0, 0]
assert close(U0.mH @ Ufound, phase * torch.eye(k, dtype=dt)) and abs(abs(phase) - 1) < 1e-12
print("U recovered up to the phase", complex(phase))
```

## Exercises

1. Show that $\mathbb R^N$ with the coordinatewise product and the standard inner product is a Euclidean Jordan algebra
   whose cone of squares is $\mathbb R^N_{\ge0}$, and deduce from Corollary F.19 that the unital order isomorphisms of
   $\mathbb C^N=\bigoplus_{c=1}^NM_1$ are the permutations of coordinates.
2. (**Positive is not enough.**) For $0<t<1$ let $\Psi_t(x)=(1-t)x+\frac tn\operatorname{Tr}(x)1$ on $M_n$, $n\ge2$. Show that
   $\Psi_t$ is a unital positive bijection, that $\Psi_t^{-1}$ is not positive, and that $\Psi_t$ does not map
   projections to projections.

<details>
<summary>Solutions</summary>

1. The product is commutative and associative, hence Jordan, with unit $(1,\dots,1)$, and $L(x)$ is the diagonal matrix
   $\operatorname{diag}(x)$, which is symmetric. The squares are the vectors with nonnegative coordinates. By
   Corollary F.19 with all $k_c=1$, $\Phi(x)_{\sigma(c)}=U_cx_cU_c^*=x_c$ (unitaries of $M_1$ are phases), so
   $\Phi$ permutes the coordinates; conversely permutations are unital order isomorphisms.

2. $\Psi_t(1)=(1-t)1+t1=1$; for $x\ge0$ both terms are $\ge0$. Since $\operatorname{Tr}\Psi_t(x)=\operatorname{Tr}x$,
   $\Psi_t^{-1}(y)=\frac1{1-t}\big(y-\frac tn\operatorname{Tr}(y)1\big)$. For a rank-one projection $p$, $\Psi_t^{-1}(p)$
   has the eigenvalue $-\frac{t}{n(1-t)}<0$ on $\ker p$, and $\Psi_t(p)$ has the eigenvalue $\frac tn\notin\{0,1\}$. This
   is consistent with Theorem F.15: positivity of the inverse is what forces projections to projections.

```python
nn, t = 3, 0.25
p = E(0, 0)
Psi_t = lambda A: (1 - t) * A + t / nn * torch.trace(A) * torch.eye(nn, dtype=dt)
Psi_t_inv = lambda A: (A - t / nn * torch.trace(A) * torch.eye(nn, dtype=dt)) / (1 - t)
assert close(Psi_t_inv(Psi_t(X)), X)
assert abs(torch.linalg.eigvalsh(Psi_t_inv(p))[0].item() + t / (nn * (1 - t))) < 1e-12
assert not close(Psi_t(p) @ Psi_t(p), Psi_t(p))                   # not a projection
```

</details>

**Back to:** [Chapter 13](13_jordan_algebras.md), Sections 13.4–13.6, where Theorems 13.8 and 13.10, Proposition 13.11
and Proposition 13.14 use the results of this appendix.
