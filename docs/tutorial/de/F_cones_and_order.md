[English](../F_cones_and_order.md) · [Español](../es/F_cones_and_order.md) · [Français](../fr/F_cones_and_order.md) · **Deutsch** · [中文](../zh/F_cones_and_order.md) · [日本語](../ja/F_cones_and_order.md) · [Русский](../ru/F_cones_and_order.md)

# Anhang F. Symmetrische Kegel, Ordnungsisomorphismen und Jordan-Abbildungen

Kapitel 13 hat die Ordnung der Observablen durch drei Sätze mit ihrem Jordan-Produkt verbunden: den Satz von
Koecher–Vinberg (Satz 13.8), den Satz von Kadison (Satz 13.10) und die Beschreibung der Jordan-Automorphismen von
$M_n(\mathbb C)$ (Proposition 13.11 beweist sie für $n=2$). Dieser Anhang beweist den ersten in voller Allgemeinheit
und die beiden anderen für die Algebren $\bigoplus_cM_{k_c}(\mathbb C)$ des Kurses, für alle Blockgrößen. Nebenbei
beweist er die Tatsachen über euklidische Jordan-Algebren, die Abschnitt 13.6 benötigt: einen Spektralsatz, die
Selbstdualität des Kegels der Quadrate und $P(y)Q\subseteq Q$.

**Sie lernen**, warum jede euklidische Jordan-Algebra potenzassoziativ ist und einen Spektralsatz besitzt; warum ihr
Kegel der Quadrate selbstdual und homogen ist; wie sich ein Jordan-Produkt allein aus einem selbstdualen homogenen
Kegel zurückgewinnen lässt, über die dritte Ableitung des Logarithmus seiner charakteristischen Funktion; warum eine
unitale lineare Bijektion, die die Positivität in beide Richtungen erhält, Quadrate erhält; und warum jede solche
Abbildung von $\bigoplus_cM_{k_c}$ gleiche Blöcke permutiert und auf jedem Block durch $UxU^*$ oder durch $Ux^TU^*$
wirkt.

Neben linearer Algebra, Differential- und Integralrechnung und den mit Nummer zitierten Ergebnissen der Kapitel 1, 6
und 13 verwendet Abschnitt F.3 die Lebesgue-Integration (das Lemma von Fatou, majorisierte Konvergenz und
Differentiation unter dem Integralzeichen, Fubini–Tonelli, die Formel für lineare Variablentransformationen); mehr wird
nicht vorausgesetzt. Der Spektralsatz aus Anhang E mit Jordan-Rahmen verfeinert den unten stehenden Satz F.3; nichts hier
hängt von ihm ab. Ein
**euklidischer Raum** ist ein endlichdimensionaler reeller Vektorraum $V\ne\{0\}$ mit einem Skalarprodukt
$\langle\cdot,\cdot\rangle$; $\|x\|$ ist die Norm (einfache Striche $|x|$ an einem Element von $V$ bezeichnen den Betrag
aus Satz F.3, nie die Norm), $\|A\|$ die Operatornorm einer linearen Abbildung $A$ von $V$, $A^T$ ihre Adjungierte
($\langle Ax,y\rangle=\langle x,A^Ty\rangle$), und $A\ge0$ ($A>0$) bedeutet, dass $A=A^T$ positiv semidefinit (definit)
ist. $e^A=\exp A=\sum_{k\ge0}A^k/k!$.

## F.1 Euklidische Jordan-Algebren

Wie in Abschnitt 13.1 ist eine **Jordan-Algebra** ein reeller Vektorraum mit einem bilinearen, kommutativen Produkt
$\circ$, sodass $x\circ(x^2\circ y)=x^2\circ(x\circ y)$, d. h. $[L(x),L(x^2)]=0$, wobei $L(x)y=x\circ y$ und $x^2=x\circ x$;
ferner $P(x)=2L(x)^2-L(x^2)$.

**Definition.** Eine **euklidische Jordan-Algebra** ist eine endlichdimensionale Jordan-Algebra $V$ mit einer Eins $e$
($e\circ x=x$) und einem Skalarprodukt, für das jedes $L(x)$ symmetrisch ist: $\langle x\circ y,z\rangle=\langle y,x\circ z\rangle$.
Ihr **Kegel der Quadrate** ist $Q=\{x^2:x\in V\}$.

Beispiele: $M_{sa}$ (der selbstadjungierte Teil von $\bigoplus_cM_{k_c}(\mathbb C)$) und $\operatorname{Herm}_n(\mathbb R)$
mit $\langle x,y\rangle=\operatorname{Tr}(xy)$, das reell ist und $\operatorname{Tr}(x^2)=\sum_{ij}|x_{ij}|^2$ erfüllt (Lemma 13.1;
$\operatorname{Tr}((x\circ y)z)=\frac12\operatorname{Tr}(xyz+yxz)=\operatorname{Re}\operatorname{Tr}(xyz)$ ist, da
$\operatorname{Tr}(yxz)=\overline{\operatorname{Tr}((yxz)^*)}=\overline{\operatorname{Tr}(zxy)}=\overline{\operatorname{Tr}(xyz)}$,
invariant unter zyklischen Vertauschungen und aus demselben Grund unter $x\leftrightarrow z$, also symmetrisch in
$x,y,z$); die Spinfaktoren $V_m$ (Jordan-Algebren nach Lemma 13.3) mit dem Skalarprodukt
$\langle x,y\rangle=\operatorname{tr}(x\circ y)=2(ab+u\cdot v)$ für $x=(a,u)$, $y=(b,v)$, da
$\operatorname{tr}((x\circ y)\circ z)=2(abc+c\,u\cdot v+a\,v\cdot w+b\,u\cdot w)$ für $z=(c,w)$ symmetrisch in
$x,y,z$ ist; und $\mathbb R^N$ mit dem koordinatenweisen Produkt (Aufgabe 1).

**Lemma F.1.** In jeder Jordan-Algebra gilt für alle $a,b,d,x,y$:

1. $[L(a),L(b\circ d)]+[L(b),L(d\circ a)]+[L(d),L(a\circ b)]=0$;
2. $L(x^2\circ y)=L(x^2)L(y)-2L(x)L(y)L(x)+2L(x\circ y)L(x)$.

*Beweis.* (1) Für $s,t,r\in\mathbb R$ und $x=sa+tb+rd$ ist $[L(x),L(x^2)]=0$ eine polynomiale Identität in $(s,t,r)$,
also verschwindet jeder Koeffizient; der Koeffizient von $str$ ist das Doppelte der linken Seite von (1). (2) Man setze
$a=b=x$ in (1): $2[L(x),L(x\circ d)]+[L(d),L(x^2)]=0$. Angewandt auf $y$ lautet das
$2x\circ((x\circ d)\circ y)-2(x\circ d)\circ(x\circ y)+d\circ(x^2\circ y)-x^2\circ(d\circ y)=0$, d. h.
$\big(2L(x)L(y)L(x)-2L(x\circ y)L(x)+L(x^2\circ y)-L(x^2)L(y)\big)d=0$ für alle $d$. $\square$

**Proposition F.2 (Potenzassoziativität).** Sei $V$ eine Jordan-Algebra mit Eins $e$, $x\in V$, $x^0=e$ und
$x^{k+1}=x\circ x^k$. Dann gilt $x^p\circ x^q=x^{p+q}$ für alle $p,q\ge0$, und jedes $L(x^p)$ ist ein Polynom in den
kommutierenden Operatoren $L(x)$ und $L(x^2)$. Folglich gilt $p(x)\circ q(x)=(pq)(x)$ für reelle Polynome $p,q$, und die
Operatoren $L(p(x))$ kommutieren miteinander.

*Beweis.* $L(x)$ und $L(x^2)$ kommutieren nach der Jordan-Identität. Wir zeigen durch Induktion über $N\ge2$: (A$_N$)
$x^i\circ x^j=x^{i+j}$ für $i+j\le N$, und (B$_N$) $L(x^i)$ ist ein Polynom in $L(x),L(x^2)$ für $i\le N$. Beides ist
für $N=2$ klar. Es gelte beides für $N$, und sei $i+j=N+1$ mit $i,j\ge1$ (für $i=0$ oder $j=0$ ist nichts zu zeigen).
Da $i\le N$, kommutiert $L(x^i)$ mit $L(x)$, also nach (A$_N$)
$$
x^i\circ x^j=L(x^i)L(x)x^{j-1}=L(x)L(x^i)x^{j-1}=x\circ x^N=x^{N+1}.
$$
Das ist (A$_{N+1}$). Lemma F.1(2) mit $y=x^{N-1}$ und $x^2\circ x^{N-1}=x^{N+1}$ liefert
$L(x^{N+1})=L(x^2)L(x^{N-1})-2L(x)L(x^{N-1})L(x)+2L(x^N)L(x)$, und das ist (B$_{N+1}$). Der letzte Satz folgt durch
Bilinearität. $\square$

**Satz F.3 (Spektralsatz).** Sei $V$ eine euklidische Jordan-Algebra und $x\in V$. Es gibt eindeutig bestimmte
$k\ge1$, reelle Zahlen $\lambda_1>\dots>\lambda_k$ und von null verschiedene $c_1,\dots,c_k\in V$ mit
$$
c_i\circ c_j=\delta_{ij}c_i,\qquad \textstyle\sum_ic_i=e,\qquad x=\sum_i\lambda_ic_i .
$$
Die $c_i$ sind Polynome in $x$, und $p(x)=\sum_ip(\lambda_i)c_i$ für jedes reelle Polynom $p$.

Wir nennen $(c_i)$ ein **vollständiges System orthogonaler Idempotente** und die $\lambda_i$ die **Eigenwerte** von $x$
und setzen $f(x)=\sum_if(\lambda_i)c_i$ für $f:\{\lambda_1,\dots,\lambda_k\}\to\mathbb R$; zum Beispiel $|x|$,
$x_\pm=\max(\pm x,0)$, $\sqrt x$ (falls alle $\lambda_i\ge0$), $\exp x$ und $\log x$ (falls alle $\lambda_i>0$).

*Beweis.* Existenz. $\mathbb R[x]=\{p(x)\}$ ist ein Unterraum, der $e$ enthält und nach Proposition F.2
$p(x)\circ q(x)=(pq)(x)$ erfüllt. Also bildet $M=L(x)|_{\mathbb R[x]}$ den Raum $\mathbb R[x]$ in sich ab, ist symmetrisch
und erfüllt $p(M)q(x)=p(x)\circ q(x)$. Ein symmetrischer Operator ist diagonalisierbar mit reellen Eigenwerten; seien
$\lambda_1>\dots>\lambda_k$ die von $M$ und $m(t)=\prod_i(t-\lambda_i)$ sein Minimalpolynom. Dann ist $m(x)=m(M)e=0$. Mit
den Lagrange-Polynomen $\ell_i(t)=\prod_{j\ne i}\frac{t-\lambda_j}{\lambda_i-\lambda_j}$ verschwinden die Polynome
$\ell_i\ell_j-\delta_{ij}\ell_i$, $\sum_i\ell_i-1$ und $\sum_i\lambda_i\ell_i-t$ in jedem $\lambda_l$, sind also Vielfache
von $m$; Auswerten in $x$ zeigt, dass die $c_i=\ell_i(x)$ die drei Eigenschaften haben. Ist $c_i=0$, so ist
$\ell_i(M)q(x)=\ell_i(x)\circ q(x)=0$ für alle $q$, also $\ell_i(M)=0$, obwohl $\deg\ell_i<\deg m$; daher $c_i\ne0$.
Schließlich $p(x)=\sum_ip(\lambda_i)c_i$: Das gilt für $p=1$ und $p=t$, und $x\circ\sum_i\lambda_i^rc_i=\sum_i\lambda_i^{r+1}c_i$.

Eindeutigkeit. Sei $x=\sum_j\mu_jd_j$ mit verschiedenen $\mu_j$ und einem vollständigen System von null verschiedener
orthogonaler Idempotente $d_j$. Wie oben ist $p(x)=\sum_jp(\mu_j)d_j$. Die $d_j$ sind linear unabhängig (man
multipliziere $\sum_ja_jd_j=0$ mit $d_i$), also ist $p(x)=0$ genau dann, wenn $p(\mu_j)=0$ für alle $j$. Somit liefert
$m(x)=0$ die Inklusion $\{\mu_j\}\subseteq\{\lambda_i\}$, und
$\tilde m(t)=\prod_j(t-\mu_j)$ erfüllt $\tilde m(M)q(x)=\tilde m(x)\circ q(x)=0$, also teilt $m$ das Polynom $\tilde m$, und $\{\lambda_i\}\subseteq\{\mu_j\}$.
Dann ist $d_j=\ell_j(x)=c_j$. $\square$

**Lemma F.4 (Peirce-Zerlegung).** Sei $V$ eine euklidische Jordan-Algebra und $c$ ein Idempotent ($c^2=c$).

1. $L(c)(2L(c)-I)(L(c)-I)=0$. Daher ist $V$ die orthogonale Summe der Eigenräume $V(c,\nu)=\{z:c\circ z=\nu z\}$,
   $\nu\in\{0,\frac12,1\}$, und $L(c)\ge0$.
2. $V(c,1)\circ V(c,0)=\{0\}$; $V(c,\nu)\circ V(c,\nu)\subseteq V(c,\nu)$ und $V(c,\nu)\circ V(c,\frac12)\subseteq V(c,\frac12)$
   für $\nu\in\{0,1\}$; $V(c,\frac12)\circ V(c,\frac12)\subseteq V(c,0)+V(c,1)$.
3. $s=2c-e$ erfüllt $s^2=e$, und $P(s)$ ist ein Jordan-Automorphismus von $V$: Er ist $+1$ auf $V(c,0)+V(c,1)$ und $-1$
   auf $V(c,\frac12)$, und $P(s)(z\circ w)=P(s)z\circ P(s)w$.
4. Sei $x=\sum_{i=1}^k\lambda_ic_i$ wie in Satz F.3. Die $L(c_i)$ kommutieren, und $V$ ist die orthogonale Summe der
   Unterräume $V_{ii}=\{z:c_i\circ z=z\}$ und $V_{ij}=\{z:c_i\circ z=c_j\circ z=\frac12z\}$ ($i<j$). Für
   $a=\sum_la_lc_l$ ($a_l\in\mathbb R$) wirkt $L(a)$ auf $V_{ij}$ ($i\le j$) als der Skalar $\frac12(a_i+a_j)$ und $P(a)$
   als der Skalar $a_ia_j$.

*Beweis.* (1) Lemma F.1(2) mit $x=y=c$ liefert $L(c)=L(c)^2-2L(c)^3+2L(c)^2$, d. h. $2L(c)^3-3L(c)^2+L(c)=0$, und das
ist die angegebene Faktorisierung. Also liegen die Eigenwerte des symmetrischen Operators $L(c)$ in $\{0,\frac12,1\}$.

(2) Lemma F.1(1) mit $a=b=c$, $d=z\in V(c,\nu)$ liefert $2\nu[L(c),L(z)]+[L(z),L(c)]=0$, also kommutiert $L(z)$ mit
$L(c)$, falls $\nu\in\{0,1\}$, und erhält dann jedes $V(c,\mu)$. Das ergibt die zweite und dritte Inklusion sowie
$V(c,1)\circ V(c,0)\subseteq V(c,0)\cap V(c,1)=\{0\}$. Für $z,w\in V(c,\frac12)$ wende man Lemma F.1(1) mit $(a,b,d)=(z,w,c)$
auf $c$ an: Da $L(w\circ c)=\frac12L(w)$ und $L(c\circ z)=\frac12L(z)$, heben sich die ersten beiden Terme weg, und
$[L(c),L(u)]c=L(c)^2u-L(c)u=0$ für $u=z\circ w$. Schreibt man $u=u_1+u_{1/2}+u_0$ gemäß (1), so ist das $-\frac14u_{1/2}=0$.

(3) $s^2=4c-4c+e=e$, $L(s)=2L(c)-I$ ist $1,0,-1$ auf $V(c,1),V(c,\frac12),V(c,0)$, und $P(s)=2L(s)^2-I$. Man schreibe
$z=z'+z''$, $w=w'+w''$ mit $z',w'\in V(c,0)+V(c,1)$ und $z'',w''\in V(c,\frac12)$. Nach (2) liegt $z'\circ w'+z''\circ w''$ in
$V(c,0)+V(c,1)$ und $z'\circ w''+z''\circ w'$ in $V(c,\frac12)$, also
$P(s)(z\circ w)=(z'-z'')\circ(w'-w'')=P(s)z\circ P(s)w$.

(4) Die $c_l$ sind Polynome in $x$, also kommutieren die $L(c_l)$ (Proposition F.2); da sie symmetrisch sind, besitzen sie
eine gemeinsame Orthonormalbasis aus Eigenvektoren. Auf einem gemeinsamen Eigenvektor hat $L(c_l)$ nach (1) einen
Eigenwert $\varepsilon_l\in\{0,\frac12,1\}$, und $\sum_l\varepsilon_l=1$, weil $\sum_lL(c_l)=L(e)=I$. Also ist entweder ein
$\varepsilon_i=1$ und die übrigen verschwinden, oder $\varepsilon_i=\varepsilon_j=\frac12$ für ein Paar $i<j$ und die übrigen
verschwinden. Somit ist $V$ die orthogonale Summe der gemeinsamen Eigenräume $W_{ii}$ (wo $L(c_i)=1$, $L(c_l)=0$ für
$l\ne i$) und $W_{ij}$ (wo $L(c_i)=L(c_j)=\frac12$, sonst $L(c_l)=0$), und $W_{ij}=V_{ij}$: Ein Vektor mit $c_i\circ z=z$
(bzw. $c_i\circ z=c_j\circ z=\frac12z$) hat keine Komponente in den anderen gemeinsamen Eigenräumen, auf denen $L(c_i)$
(bzw. das Paar $(L(c_i),L(c_j))$) einen anderen Wert annimmt. Da
$a^2=\sum_la_l^2c_l$, ist auf $V_{ij}$ ($i<j$) $L(a)=\frac12(a_i+a_j)$, $L(a^2)=\frac12(a_i^2+a_j^2)$ und
$P(a)=\frac12(a_i+a_j)^2-\frac12(a_i^2+a_j^2)=a_ia_j$; auf $V_{ii}$ ist $L(a)=a_i$ und $P(a)=2a_i^2-a_i^2$. $\square$

## F.2 Der Kegel der Quadrate ist symmetrisch

Sei $\Omega$ ein nichtleerer offener **konvexer Kegel** in einem euklidischen Raum $V$ ($tx\in\Omega$ für $x\in\Omega$,
$t>0$). Sein (offener) **dualer Kegel** ist $\Omega^*=\{y:\langle x,y\rangle>0\text{ für alle }x\in\bar\Omega\setminus\{0\}\}$;
$\Omega$ heißt **selbstdual**, wenn $\Omega^*=\Omega$. $G(\Omega)=\{g\in GL(V):g\Omega=\Omega\}$, und $\Omega$ heißt
**homogen**, wenn $G(\Omega)$ transitiv auf $\Omega$ wirkt. Ein **symmetrischer Kegel** ist ein selbstdualer homogener
Kegel. Ist $\Omega_1\subseteq\Omega_2$, so ist $\Omega_2^*\subseteq\Omega_1^*$.

**Satz F.5 (Selbstdualität).** Sei $V$ eine euklidische Jordan-Algebra und $y=\sum_i\lambda_ic_i$ (Satz F.3).

1. Äquivalent sind: (i) $y\in Q$; (ii) alle $\lambda_i\ge0$; (iii) $L(y)\ge0$; (iv) $\langle y,x\rangle\ge0$
   für alle $x\in Q$. Insbesondere ist $Q$ ein abgeschlossener konvexer Kegel.
2. Äquivalent sind: (i) $y\in\Omega:=\operatorname{int}Q$; (ii) alle $\lambda_i>0$; (iii) $L(y)>0$; (iv)
   $\langle y,x\rangle>0$ für alle $x\in Q\setminus\{0\}$. Außerdem ist $e\in\Omega$, $\bar\Omega=Q$ und $\Omega^*=\Omega$:
   $\Omega$ ist ein selbstdualer offener konvexer Kegel.
3. Sind $x,y\in Q$ und $\langle x,y\rangle=0$, so ist $x\circ y=0$.
4. Aus $x^2=0$ folgt $x=0$, und $Q\cap(-Q)=\{0\}$.

*Beweis.* Wir verwenden $\langle y,x^2\rangle=\langle y\circ x,x\rangle=\langle L(y)x,x\rangle$ und
$\langle x^2,e\rangle=\langle x,x\circ e\rangle=\|x\|^2$; Letzteres liefert zusammen mit (1)(ii) die Aussage (4), da die
Eigenwerte von $-y$ wegen der Eindeutigkeit die $-\lambda_i$ sind.

(1) (i)$\Rightarrow$(iii): Ist $y=x^2$ mit $x=\sum_i\mu_id_i$, so ist $y=\sum_i\mu_i^2d_i$ und
$L(y)=\sum_i\mu_i^2L(d_i)\ge0$ nach Lemma F.4(1). (iii)$\Leftrightarrow$(iv) nach der ersten Identität.
(iii)$\Rightarrow$(ii): $0\le\langle L(y)c_i,c_i\rangle=\lambda_i\|c_i\|^2$. (ii)$\Rightarrow$(i):
$y=(\sum_i\sqrt{\lambda_i}c_i)^2$. Nach (iv) ist $Q$ ein Durchschnitt abgeschlossener Halbräume durch $0$.

(2) (ii)$\Rightarrow$(iii): $L(y)=\sum_i\lambda_iL(c_i)\ge(\min_i\lambda_i)\sum_iL(c_i)=(\min_i\lambda_i)I$.
(iii)$\Rightarrow$(iv): Ein von null verschiedenes $x\in Q$ ist $z^2$ mit $z\ne0$, und $\langle y,x\rangle=\langle L(y)z,z\rangle>0$.
(iv)$\Rightarrow$(ii): $\langle y,c_i\rangle=\lambda_i\|c_i\|^2$ und $c_i=c_i^2\in Q\setminus\{0\}$. (iii)$\Rightarrow$(i): Ist
$L(y)\ge mI$ mit $m>0$, so ist $L(y')\ge(m-\|L\|\,\|y-y'\|)I$, wobei $\|L\|$ die Norm der linearen Abbildung
$b\mapsto L(b)$ ist; also besteht eine Umgebung von $y$ aus Elementen mit
$L>0$, die nach (1) in $Q$ liegen. (i)$\Rightarrow$(iii): $L(y)\ge0$ nach (1); ist $\langle L(y)z,z\rangle=0$ mit $z\ne0$,
so ist $x=z^2\ne0$ und $\langle y,x\rangle=0$; für kleines $\varepsilon>0$ ist $y-\varepsilon x\in Q$, aber
$\langle y-\varepsilon x,x\rangle=-\varepsilon\|x\|^2<0$, im Widerspruch zu (1)(iv). Weiter ist $L(e)=I>0$, und für $y\in Q$
strebt $y+\varepsilon e\in\Omega$ gegen $y$, also $\bar\Omega=Q$; dann ist $\Omega^*=\{y:\langle y,x\rangle>0\ \forall x\in Q\setminus\{0\}\}=\Omega$
nach (iv). $\Omega=\{L>0\}$ ist konvex.

(3) Sei $x=\sum_i\lambda_ic_i$ mit $\lambda_i\ge0$. Dann ist $0=\sum_i\lambda_i\langle c_i,y\rangle$, wobei nach (1) alle
Terme $\ge0$ sind. Ist $\lambda_i>0$, so ist $0=\langle y,c_i\rangle=\langle L(y)c_i,c_i\rangle$ mit $L(y)\ge0$, also
$L(y)c_i=0$. Daher ist $x\circ y=\sum_i\lambda_i\,c_i\circ y=0$. $\square$

**Satz F.6 (Homogenität und quadratische Darstellung).** Sei $V$ eine euklidische Jordan-Algebra, $Q$ und
$\Omega$ wie in Satz F.5.

1. Für jedes $b\in V$: $e^{L(b)}Q=Q$, $e^{L(b)}\Omega=\Omega$ und $e^{L(b)}e=\exp b$.
2. Jedes $y\in\Omega$ ist gleich $\exp(\log y)$. Daher ist $\Omega$ homogen und ein symmetrischer Kegel.
3. Für $y\in\Omega$ ist $P(y)=e^{2L(\log y)}$; also $P(y)Q=Q$, $P(y)\Omega=\Omega$ und $P(y)e=y^2$.
4. $P(y)Q\subseteq Q$ für jedes $y\in V$.

*Beweis.* (1) Sei $0\le s<1/\|L(b)\|$ (beliebiges $s\ge0$, falls $L(b)=0$). Der symmetrische Operator $I-sL(b)$ hat
Eigenwerte $\ge1-s\|L(b)\|>0$. Sei $x\in Q$ und $z=(I-sL(b))^{-1}x=\sum_j\nu_jd_j$ (Satz F.3). Mit $z\circ d_j=\nu_jd_j$,
$\langle d_l,d_j\rangle=\langle d_l,d_j\circ d_j\rangle=\langle d_l\circ d_j,d_j\rangle=\delta_{lj}\|d_j\|^2$ und
$\langle b,d_j\rangle=\langle b,d_j^2\rangle=\langle L(b)d_j,d_j\rangle\le\|L(b)\|\,\|d_j\|^2$ gilt
$$
0\le\langle x,d_j\rangle=\langle z-s\,b\circ z,d_j\rangle=\nu_j\|d_j\|^2-s\langle b,z\circ d_j\rangle
=\nu_j\big(\|d_j\|^2-s\langle b,d_j\rangle\big),
$$
wobei die linke Ungleichung Satz F.5(1)(iv) ist. Die Klammer ist $\ge(1-s\|L(b)\|)\|d_j\|^2>0$, also sind alle
$\nu_j\ge0$ und $z\in Q$. Daher gilt für $t\ge0$ und $n>t\|L(b)\|$: $(I-\frac tnL(b))^{-n}Q\subseteq Q$. Für $n\to\infty$
strebt dieser Operator gegen $e^{tL(b)}$ (in einer Orthonormalbasis aus Eigenvektoren von $L(b)$:
$(1-t\beta/n)^{-n}\to e^{t\beta}$), und $Q$ ist abgeschlossen, also $e^{tL(b)}Q\subseteq Q$. Mit $-b$ statt $b$ folgt
$e^{-L(b)}Q\subseteq Q$, also $e^{L(b)}Q=Q$, und der Homöomorphismus $e^{L(b)}$ bildet $\operatorname{int}Q$ auf sich ab.
Schließlich ist $L(b)^re=b^r$ (Jordan-Potenzen), also $e^{L(b)}e=\sum_rb^r/r!=\sum_ie^{\beta_i}d_i=\exp b$ für
$b=\sum_i\beta_id_i$, nach Satz F.3.

(2) Ist $y=\sum_i\lambda_ic_i$ mit $\lambda_i>0$, so ist $\log y=\sum_i(\log\lambda_i)c_i$ und $\exp\log y=y$. Also ist
$y=e^{L(\log y)}e$ mit $e^{L(\log y)}\in G(\Omega)$ nach (1); die Selbstdualität ist Satz F.5(2).

(3) Nach Lemma F.4(4) für die $c_i$ von $y$ wirken sowohl $P(y)$ als auch $e^{2L(\log y)}$ auf $V_{ij}$ als
$\lambda_i\lambda_j=e^{\log\lambda_i+\log\lambda_j}$. Man verwende (1) und $P(y)e=2y\circ y-y^2=y^2$.

(4) Sei $y=\sum_i\lambda_ic_i$, $c=\sum_{\lambda_i\ge0}c_i$ (ein Idempotent), $s=2c-e=\sum_i\varepsilon_ic_i$ mit
$\varepsilon_i=\pm1$, und $|y|=\sum_i|\lambda_i|c_i$. Nach Lemma F.4(4) wirken $P(s)$, $P(|y|)$ und $P(y)$ auf $V_{ij}$ als
$\varepsilon_i\varepsilon_j$, $|\lambda_i\lambda_j|$ und $\lambda_i\lambda_j$, also $P(y)=P(s)P(|y|)$. Für $\eta>0$ ist
$|y|+\eta e\in\Omega$, also $P(|y|+\eta e)Q\subseteq Q$ nach (3); $P(\cdot)$ ist ein Polynom und $Q$ abgeschlossen, also
$P(|y|)Q\subseteq Q$ ($\eta\to0$). $P(s)$ ist ein Jordan-Automorphismus (Lemma F.4(3)), also $P(s)x^2=(P(s)x)^2\in Q$. $\square$

Abschnitt 13.6 hat $\operatorname{tr}|P(y)x|\le\operatorname{tr}P(y)|x|$ in $V_m$ aus drei Eigenschaften bewiesen. Sie
gelten in jeder euklidischen Jordan-Algebra für die Linearform $\tau(x)=\langle x,e\rangle$, die auf $M_{sa}$ gleich
$\operatorname{Tr}$ und auf $V_m$ gleich $\operatorname{tr}$ ist (mit den obigen Skalarprodukten):

**Korollar F.7 (Proposition 13.14 in jeder euklidischen Jordan-Algebra).** Für $x,y$ in einer euklidischen
Jordan-Algebra gilt $\tau(|P(y)x|)\le\tau(P(y)|x|)$, mit Gleichheit, falls $x\in Q$ oder $-x\in Q$.

*Beweis.* $x_\pm\in Q$, $x=x_+-x_-$, $|x|=x_++x_-$ (Satz F.3). $A=P(y)x_+$ und $B=P(y)x_-$ liegen in $Q$ (Satz F.6(4)),
und $z=P(y)x=A-B$. Sei $z=\sum_i\mu_ic_i$, $c=\sum_{\mu_i>0}c_i$; dann ist $|z|=z\circ(2c-e)$, und
$\tau(a\circ b)=\langle a,b\rangle$. Da $\langle Q,Q\rangle\ge0$ (Satz F.5(1)) und $c,e-c\in Q$, gilt
$$
\tau(|z|)=\langle A,c\rangle-\langle A,e-c\rangle-\langle B,c\rangle+\langle B,e-c\rangle\le\tau(A)+\tau(B)=\tau(P(y)|x|).
$$
Ist $x\in Q$, so ist $x_-=0$ und $z=A\in Q$, also $|z|=z=P(y)|x|$; ebenso für $-x\in Q$. $\square$

Der Code prüft Satz F.5 für $\operatorname{Herm}_4(\mathbb C)$ und für den Lorentz-Kegel $Q\subseteq V_4$ (die negative
Paarung ist der Zeuge $c_i$ aus (1), $\langle y,c_i\rangle=\lambda_i\|c_i\|^2<0$), Satz F.6 in $V_4$, wo die Operatoren $L(b)$ und $P(y)$ in den
Koordinaten $(a,u)$ symmetrische $5\times5$-Matrizen sind, und Korollar F.7 in
$\operatorname{Herm}_5(\mathbb R)$:

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
pair = (Z @ c1).diagonal(dim1=-2, dim2=-1).sum(-1).real            # <Z, c1> = lambda_min ||c1||^2
assert close(pair, w[:, 0]) and (pair < 0).float().mean() > 0.9   # witnesses of Z not in Q

V = SpinFactor(4)                                                  # <x, y> = tr(x o y) = V.inner
x, y = V.random(2000, cone=True), V.random(2000, cone=True)
assert (V.inner(x, y) > -1e-12).all()
z = V.random(2000)
lam, idem = V.spectral_decomposition(z)                            # (lambda-, lambda+), (c-, c+); ||c-||^2 = 1
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

## F.3 Der Satz von Koecher–Vinberg

Abschnitt F.2 hat bewiesen, dass das Innere des Kegels der Quadrate einer euklidischen Jordan-Algebra symmetrisch ist.
Die Umkehrung rekonstruiert das Jordan-Produkt aus dem Kegel. Bis Satz F.12 ist $\Omega$ ein selbstdualer offener
konvexer Kegel in einem euklidischen Raum $V$ der Dimension $N$, und $dy$ ist das Lebesgue-Maß in Orthonormalkoordinaten
(orthogonale Koordinatenwechsel haben $|\det|=1$).

**Lemma F.8.** Sei $\Omega$ selbstdual.

1. $0\notin\Omega$; $\langle x,y\rangle\ge0$ für $x,y\in\bar\Omega$; zu jedem $x\in\Omega$ gibt es $\delta>0$ mit
   $\langle x,y\rangle\ge\delta\|y\|$ für alle $y\in\bar\Omega$.
2. Ist $x\in\bar\Omega\setminus\Omega$, so gibt es $u\in\bar\Omega\setminus\{0\}$ mit $\langle x,u\rangle=0$. Ist $y\in\Omega$,
   $u\in\bar\Omega$ und $t\ge0$, so ist $y+tu\in\Omega$.
3. $\Omega$ ist das Innere von $\bar\Omega$.
4. $G(\Omega)$ ist eine abgeschlossene Untergruppe von $GL(V)$, und aus $g\in G(\Omega)$ folgt $g^T\in G(\Omega)$.

*Beweis.* (1) Ist $0\in\Omega$, so enthält $\Omega$ eine Kugel um $0$, also $\Omega=V$; aber $V^*=\emptyset$. Für
$x,y\in\Omega$ ist $\langle x,y\rangle>0$, weil $y\in\Omega^*$; aus Stetigkeitsgründen ist $\langle x,y\rangle\ge0$ auf
$\bar\Omega$. Liegt die abgeschlossene Kugel $\bar B(x,\delta)$ in $\Omega$ und ist $y\in\bar\Omega\setminus\{0\}$, so ist
$\langle x-\delta y/\|y\|,y\rangle\ge0$. (2) $x\notin\Omega=\Omega^*$ liefert $u\in\bar\Omega\setminus\{0\}$ mit
$\langle u,x\rangle\le0$, und $\langle u,x\rangle\ge0$ nach (1). Für $x'\in\bar\Omega\setminus\{0\}$ ist
$\langle x',y+tu\rangle\ge\langle x',y\rangle>0$, also $y+tu\in\Omega^*=\Omega$. (3) Ist $\bar B(y,\varepsilon)\subseteq\bar\Omega$,
so gilt für $x\in\bar\Omega\setminus\{0\}$: $\langle x,y-\varepsilon x/\|x\|\rangle\ge0$, also $\langle x,y\rangle\ge\varepsilon\|x\|>0$
und $y\in\Omega^*=\Omega$. (4) Für $g\in G(\Omega)$ bildet $g$ die Menge $\bar\Omega\setminus\{0\}$ auf sich ab, also gilt für
$y\in\Omega$ und $x\in\bar\Omega\setminus\{0\}$: $\langle x,g^Ty\rangle=\langle gx,y\rangle>0$, d. h. $g^T\Omega\subseteq\Omega$.
Dasselbe für $g^{-1}$ liefert $g^T\Omega=\Omega$. Konvergieren $g_k\in G(\Omega)$ gegen $g\in GL(V)$, so ist
$g\Omega\subseteq\bar\Omega$ offen, also $g\Omega\subseteq\Omega$ nach (3); ebenso bildet $g^{-1}=\lim g_k^{-1}$ die Menge
$\Omega$ in $\Omega$ ab, also $g\Omega=\Omega$. $\square$

Die **charakteristische Funktion** von $\Omega$ ist $\varphi(x)=\int_{\Omega^*}e^{-\langle x,y\rangle}dy=\int_\Omega e^{-\langle x,y\rangle}dy$.
Man setze $f=\log\varphi$, sei $H(x)$ der symmetrische Operator mit $\langle H(x)u,v\rangle=D^2f(x)[u,v]$ und
$x^*=-\nabla f(x)$, d. h. $Df(x)[u]=-\langle x^*,u\rangle$.

**Lemma F.9 (charakteristische Funktion).** Sei $\Omega$ selbstdual.

1. $\varphi$ ist endlich, positiv und $C^\infty$ auf $\Omega$ und darf unter dem Integralzeichen differenziert werden.
2. $\varphi(tx)=t^{-N}\varphi(x)$ für $t>0$, und $\varphi(gx)=|\det g|^{-1}\varphi(x)$ für $g\in G(\Omega)$.
3. $H(x)>0$ und $x^*\in\Omega$ für $x\in\Omega$.
4. $H(x)x=x^*$ und $D^3f(x)[x,u,v]=-2\langle H(x)u,v\rangle$; für $g\in G(\Omega)$ gilt $(gx)^*=g^{-T}x^*$ und
   $H(gx)=g^{-T}H(x)g^{-1}$.
5. $\varphi(x)\ge c\|x\|^{-N}$ mit einer Konstanten $c>0$, und $\varphi(x_k)\to\infty$, wenn $x_k\in\Omega$ gegen einen Punkt
   von $\bar\Omega\setminus\Omega$ konvergieren.

*Beweis.* (1) Sei $x\in\Omega$ und $\delta$ wie in Lemma F.8(1). Für $\|x'-x\|<\delta/2$ und $y\in\Omega$ ist
$\langle x',y\rangle\ge\frac\delta2\|y\|$, also sind die $k$-ten Ableitungen des Integranden nach $x'$,
$(-1)^k\prod_{i}\langle u_i,y\rangle e^{-\langle x',y\rangle}$, beschränkt durch $\prod_i\|u_i\|\,\|y\|^ke^{-\delta\|y\|/2}$, und
das ist integrierbar ($\|y\|\ge N^{-1/2}\sum_i|y_i|$ und Fubini). Majorisierte Konvergenz erlaubt die Differentiation
unter dem Integralzeichen. (2) Man substituiere $y=y'/t$, was $\Omega$ auf $\Omega$ abbildet, und $y'=g^Ty$, was $\Omega$
nach Lemma F.8(4) auf $\Omega$ abbildet, mit $\langle gx,y\rangle=\langle x,y'\rangle$ und $dy=|\det g|^{-1}dy'$.
(3) Sei $\mu_x$ das Wahrscheinlichkeitsmaß $e^{-\langle x,y\rangle}dy/\varphi(x)$ auf $\Omega$. Nach (1) ist
$Df(x)[u]=-\mathbb E\langle u,y\rangle$ und $D^2f(x)[u,u]=\mathbb E\langle u,y\rangle^2-(\mathbb E\langle u,y\rangle)^2$, die
Varianz von $\langle u,y\rangle$. Sie verschwindet nur, wenn $\langle u,y\rangle$ $\mu_x$-fast überall gleich einer
Konstanten $\gamma$ ist, d. h. Lebesgue-fast überall auf $\Omega$ (die Dichte von $\mu_x$ ist positiv); die Menge
$\{y\in\Omega:\langle u,y\rangle\ne\gamma\}$ ist offen und eine Nullmenge, also leer, sodass die lineare Funktion
$\langle u,\cdot\rangle$ auf einer Kugel konstant ist, was $u=0$ erzwingt. Und $x^*=\mathbb E\,y$, also
$\langle z,x^*\rangle=\mathbb E\langle z,y\rangle>0$ für $z\in\bar\Omega\setminus\{0\}$: $x^*\in\Omega^*=\Omega$.
(4) Nach (2) ist $f(tx)=f(x)-N\log t$, also $Df(x)[x]=-N$ auf $\Omega$. Differentiation nach $x$ in Richtung $u$:
$D^2f(x)[x,u]+Df(x)[u]=0$, d. h. $\langle H(x)x,u\rangle=\langle x^*,u\rangle$; noch einmal in Richtung $v$:
$D^3f(x)[x,u,v]+2D^2f(x)[u,v]=0$. Ein- und zweimaliges Differenzieren von $f(gx)=f(x)-\log|\det g|$ liefert
$\langle(gx)^*,gu\rangle=\langle x^*,u\rangle$ und $\langle H(gx)gu,gv\rangle=\langle H(x)u,v\rangle$.
(5) $\langle x,y\rangle\le\|x\|\|y\|\le1$ für $\|y\|\le1/\|x\|$, also $\varphi(x)\ge e^{-1}\operatorname{vol}(\Omega\cap B(0,1/\|x\|))=c\|x\|^{-N}$
mit $c=e^{-1}\operatorname{vol}(\Omega\cap B(0,1))>0$. Sei $x_k\to x_0\in\bar\Omega\setminus\Omega$. Nach dem Lemma von Fatou ist
$\liminf\varphi(x_k)\ge\int_\Omega e^{-\langle x_0,y\rangle}dy$. Man nehme $u$ wie in Lemma F.8(2), eine Kugel
$B(y_1,r)\subseteq\Omega$ und die disjunkten Kugeln $B_l=B(y_1+2lru/\|u\|,r)$, $l\ge0$; sie liegen nach Lemma F.8(2) in
$\Omega$, und auf $B_l$ ist $\langle x_0,y\rangle=\langle x_0,y-2lru/\|u\|\rangle\le\|x_0\|(\|y_1\|+r)$. Also ist das Integral
unendlich. $\square$

**Lemma F.10 (Einparametergruppen).** Sei $G$ eine abgeschlossene Untergruppe von $GL(V)$ und
$\mathfrak g=\{Z:e^{tZ}\in G\text{ für alle }t\in\mathbb R\}$.

1. Erfüllt $c:[0,s_0)\to G$ die Bedingung $c(s)=I+sZ+o(s)$ für $s\to0^+$, so ist $Z\in\mathfrak g$.
2. $\mathfrak g$ ist eine Lie-Algebra: ein linearer Unterraum, der unter $[Z,W]=ZW-WZ$ abgeschlossen ist.

*Beweis.* (1) Man fixiere $t>0$ und setze $A_n=c(t/n)$, $B_n=e^{tZ/n}$. Dann ist $\|A_n-B_n\|=o(1/n)$ und
$\|A_n\|,\|B_n\|\le1+K/n$ für ein $K$ und großes $n$, also
$$
\|A_n^n-e^{tZ}\|=\Big\|\sum_{r=0}^{n-1}A_n^r(A_n-B_n)B_n^{n-1-r}\Big\|\le n\,e^K\|A_n-B_n\|\to0 .
$$
$A_n^n\in G$, und $G$ ist abgeschlossen in $GL(V)$, also $e^{tZ}\in G$; ebenso $e^{-tZ}=(e^{tZ})^{-1}\in G$. (2) Vielfache
sind klar. Man wende (1) auf $c(s)=e^{sZ}e^{sW}=I+s(Z+W)+O(s^2)$ an und auf $c(s)=e^{aZ}e^{aW}e^{-aZ}e^{-aW}$ mit
$a=\sqrt s$: Entwickelt man jeden Faktor bis zur zweiten Ordnung, so heben sich die Terme der Ordnung $a$ weg, und die der
Ordnung $a^2$ summieren sich zu $ZW-WZ$, also $c(s)=I+s[Z,W]+O(s^{3/2})$. $\square$

**Proposition F.11 (der Basispunkt).** Sei $\Omega$ ein symmetrischer Kegel, $G=G(\Omega)$, $\mathfrak g$ wie in Lemma F.10.

1. $x^{**}=x$ und $H(x^*)=H(x)^{-1}$ für $x\in\Omega$.
2. Es gibt $e\in\Omega$ mit $e^*=e$, und jedes solche $e$ erfüllt $H(e)=I$.
3. Für $x\in\Omega$ und jedes $g\in G$ mit $ge=x$ gilt $H(x)^{-1}=gg^T\in G$.
4. Ist $Z\in\mathfrak g$ und $Ze=0$, so ist $Z^T=-Z$.

*Beweis.* (1) Die Abbildung $x\mapsto x^*$ ist glatt von $\Omega$ nach $\Omega$ (Lemma F.9(3)), mit Ableitung $-H(x)$. Da
$g^{-T}\in G$ (Lemma F.8(4)), liefert Lemma F.9(2,4)
$f(gx)+f((gx)^*)=f(x)-\log|\det g|+f(g^{-T}x^*)=f(x)+f(x^*)$, also ist $F(x)=f(x)+f(x^*)$ konstant auf der Bahn
$G x=\Omega$. Daher ist $0=DF(x)[u]=-\langle x^*,u\rangle+\langle x^{**},H(x)u\rangle$, d. h. $H(x)x^{**}=x^*=H(x)x$
(Lemma F.9(4)), und $x^{**}=x$, da $H(x)$ invertierbar ist. Differenzieren von $x^{**}=x$ liefert $H(x^*)H(x)=I$.
(2) $F_2(x)=f(x)+\frac12\|x\|^2$ erfüllt $F_2(x)\ge\log c-N\log\|x\|+\frac12\|x\|^2$ (Lemma F.9(5)), was für
$\|x\|\to\infty$ gegen $\infty$ strebt. Also ist $S=\{x\in\Omega:F_2(x)\le F_2(x_1)\}$ ($x_1\in\Omega$) beschränkt; es ist
abgeschlossen in $V$, weil ein Häufungspunkt in $\bar\Omega\setminus\Omega$ die Funktion $f$ auf $S$ unbeschränkt machen
würde (Lemma F.9(5)). Die stetige Funktion $F_2$ nimmt ihr Minimum auf der kompakten Menge $S$ in einem Punkt $e$ der
offenen Menge $\Omega$ an, wo $0=\nabla F_2(e)=-e^*+e$. Dann ist $H(e)=H(e^*)=H(e)^{-1}$ nach (1), also $H(e)^2=I$, und
$H(e)=I$, da $H(e)>0$.
(3) $H(x)=H(ge)=g^{-T}g^{-1}$ nach Lemma F.9(4) und (2), und $gg^T\in G$ nach Lemma F.8(4).
(4) $e^{tZ}\in G$ fixiert $e$, also $I=H(e)=e^{-tZ^T}H(e)e^{-tZ}=e^{-tZ^T}e^{-tZ}$; man differenziere bei $t=0$. $\square$

**Satz F.12 (Koecher–Vinberg).** Sei $\Omega$ ein nichtleerer offener konvexer Kegel in einem euklidischen Raum $V$.
Dann ist $\Omega$ genau dann ein symmetrischer Kegel, wenn es ein Produkt $\circ$ auf $V$ gibt, sodass
$(V,\circ,\langle\cdot,\cdot\rangle)$ eine euklidische Jordan-Algebra und $\Omega$ das Innere ihres Kegels der Quadrate ist.
Für einen symmetrischen Kegel kann man als Eins jedes $e\in\Omega$ mit $e^*=e$ nehmen, und als Produkt
$$
\langle u\circ v,w\rangle=-\tfrac12\,D^3\log\varphi(e)[u,v,w].
$$

*Beweis.* Die Richtung „wenn“ ergibt sich aus den Sätzen F.5(2) und F.6(2). Umgekehrt sei $\Omega$ symmetrisch, $e$ wie in Proposition F.11(2),
und man definiere $L(u)$ durch $\langle L(u)v,w\rangle=-\frac12D^3f(e)[u,v,w]$ und $u\circ v=L(u)v$.

*Schritt 1: Kommutativität, Symmetrie, Eins.* $D^3f(e)$ ist symmetrisch in seinen drei Argumenten, also $u\circ v=v\circ u$
und $L(u)^T=L(u)$. Nach Lemma F.9(4) und $H(e)=I$ ist $\langle e\circ v,w\rangle=\langle H(e)v,w\rangle=\langle v,w\rangle$:
$e$ ist eine Eins.

*Schritt 2: $L(u)\in\mathfrak g$.* $R(x)=H(x)^{-1}$ liegt in $G$ (Proposition F.11(3)), ist glatt, und $R(e)=I$. Seine
Ableitung ist $DR(e)[u]=-H(e)^{-1}DH(e)[u]H(e)^{-1}=-DH(e)[u]=2L(u)$, da $\langle DH(e)[u]v,w\rangle=D^3f(e)[u,v,w]$.
Lemma F.10(1) für $c(s)=R(e+su)$ ($s\ge0$ klein, sodass $e+su\in\Omega$), mit der abgeschlossenen Untergruppe $G$
(Lemma F.8(4)), liefert $2L(u)\in\mathfrak g$, also $L(u)\in\mathfrak g$.

*Schritt 3: Eindeutigkeit.* Ist $Y\in\mathfrak g$ symmetrisch, so ist $Y=L(Ye)$: Das symmetrische $Y-L(Ye)\in\mathfrak g$
annulliert $e$ ($L(w)e=w\circ e=w$), ist also nach Proposition F.11(4) auch antisymmetrisch und damit $0$.

*Schritt 4: innere Derivationen.* Sei $D=[L(x),L(y)]$. Es liegt in $\mathfrak g$ (Lemma F.10(2)), $D^T=-D$, und
$De=x\circ y-y\circ x=0$. Für $z\in V$ ist $[D,L(z)]\in\mathfrak g$ symmetrisch und $[D,L(z)]e=Dz-L(z)De=Dz$, also
$[D,L(z)]=L(Dz)$ nach Schritt 3; angewandt auf $w$: $D(z\circ w)=Dz\circ w+z\circ Dw$.

*Schritt 5: die Jordan-Identität.* Für $D=[L(y),L(x)]$ ist $Dx=y\circ x^2-x\circ(x\circ y)$ und
$D(x^2)=y\circ(x\circ x^2)-x\circ(y\circ x^2)$. Schritt 4 liefert $D(x^2)=2x\circ Dx$, d. h.
$y\circ(x\circ x^2)=3x\circ(x^2\circ y)-2x\circ(x\circ(x\circ y))$ für alle $y$:
$$
L(x\circ x^2)=3L(x)L(x^2)-2L(x)^3 .
$$
Alle $L$ sind symmetrisch, also liefert Transponieren $L(x\circ x^2)=3L(x^2)L(x)-2L(x)^3$, und $[L(x),L(x^2)]=0$. Somit
ist $V$ eine euklidische Jordan-Algebra; sei $\Omega_J$ das Innere ihres Kegels der Quadrate.

*Schritt 6: $\Omega=\Omega_J$.* Für $b\in V$ ist $e^{L(b)}\in G$ nach Schritt 2, also $\exp b=e^{L(b)}e\in\Omega$
(Satz F.6(1)). Nach Satz F.6(2) ist $\Omega_J\subseteq\Omega$. Beide Kegel sind selbstdual bezüglich
$\langle\cdot,\cdot\rangle$ (Satz F.5(2)), also $\Omega=\Omega^*\subseteq\Omega_J^*=\Omega_J$. $\square$

Der Beweis folgt dem Aufbau von J. Faraut, A. Korányi, *Analysis on Symmetric Cones* (Oxford 1994), Kapitel I–III;
wo dort die Lie-Gruppen-Struktur von $G(\Omega)$ verwendet wird, liefern Lemma F.10 und Proposition F.11(3) die
benötigten Einparametergruppen direkt. M. Koecher (1957) und E. B. Vinberg (1960) haben den Satz bewiesen.

**Beispiel F.13 (der Kegel der semidefiniten Programmierung).** Sei $V=\operatorname{Herm}_n(\mathbb R)$, $\langle x,y\rangle=\operatorname{Tr}(xy)$,
$N=n(n+1)/2$ und $\Omega$ die Menge der positiv definiten Matrizen, symmetrisch nach den Sätzen F.5 und F.6. Dann ist
$\varphi(x)=\varphi(1)\det(x)^{-(n+1)/2}$, $e=\kappa1$ mit $\kappa=\sqrt{(n+1)/2}$, und das Produkt aus Satz F.12 ist
$u\circ_\kappa v=\frac1{2\kappa}(uv+vu)$, also ist $x\mapsto x/\kappa$ ein Isomorphismus auf die übliche Jordan-Algebra.

*Beweis.* Die Spektralzerlegung $x=\sum_i\lambda_ip_i$ einer symmetrischen Matrix in ihre Eigenprojektionen ist die
Zerlegung aus Satz F.3 (Eindeutigkeit; $p_i\circ p_j=p_ip_j$, da die $p_i$ kommutieren), also ist nach Satz F.5(2) das
Innere des Kegels der Quadrate tatsächlich die Menge der positiv definiten Matrizen. Für $a\in\Omega$ liegt
$P(a)x=axa$ (Proposition 13.2(1)) in $G(\Omega)$ (Satz F.6(3)). In einer Orthonormalbasis aus Eigenvektoren von $a$
bilden die Matrizen $E_{ii}$ und $(E_{ij}+E_{ji})/\sqrt2$ ($i<j$) eine Orthonormalbasis aus Eigenvektoren von $P(a)$ mit
den Eigenwerten $a_ia_j$, also $\det P(a)=\prod_{i\le j}a_ia_j=\det(a)^{n+1}$ (jedes $a_i$ kommt $n+1$-mal vor). Lemma
F.9(2) liefert $\varphi(a^2)=\varphi(P(a)1)=\det(a)^{-(n+1)}\varphi(1)$, und jedes $x\in\Omega$ ist $a^2$ mit $a=\sqrt x$. Mit
$D\log\det(x)[u]=\operatorname{Tr}(x^{-1}u)$ (aus $\det(x+su)=\det x\,\det(1+sx^{-1}u)=\det x\,(1+s\operatorname{Tr}(x^{-1}u)+O(s^2))$)
und $D(x^{-1})[u]=-x^{-1}ux^{-1}$ folgt:
$x^*=\kappa^2x^{-1}$; also bedeutet $e^*=e$, dass $e^2=\kappa^2 1$ mit positiv definitem $e$, d. h. $e=\kappa1$; $D^2f(x)[v,w]=\kappa^2\operatorname{Tr}(x^{-1}vx^{-1}w)$, und
$D^3f(e)[u,v,w]=-\kappa^2\kappa^{-3}\operatorname{Tr}(uvw+vuw)=-\frac2\kappa\operatorname{Tr}((u\circ v)w)$. $\square$

Der Code berechnet $D^3f$ in $e$ durch automatisches Differenzieren von $f(x)=-\frac{n+1}2\log\det x$ in
Orthonormalkoordinaten von $\operatorname{Herm}_3(\mathbb R)$ und prüft die Aussagen von Proposition F.11 und Satz F.12:

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

## F.4 Ordnungsisomorphismen: der Satz von Kadison in endlicher Dimension

Seien $\mathcal M=\bigoplus_cM_{k_c}(\mathbb C)$ und $\mathcal N=\bigoplus_dM_{k'_d}(\mathbb C)$ als blockdiagonale Matrizen
gegeben; $x\ge0$ bedeutet $x=x^*$ mit nichtnegativen Eigenwerten (gleichbedeutend: jeder Block ist positiv
semidefinit), und $[0,1]=\{a\in\mathcal M_{sa}:0\le a\le1\}$. Eine **Projektion** ist ein $p=p^*=p^2$. Eine lineare
Abbildung $\Phi$ heißt **positiv**, wenn aus $x\ge0$ folgt $\Phi(x)\ge0$; ein **Jordan-$*$-Isomorphismus** ist eine
lineare Bijektion mit $\Phi(x^*)=\Phi(x)^*$ und $\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$ (Abschnitt 13.5). Ein Punkt $a$ einer
konvexen Menge $K$ heißt **Extremalpunkt**, wenn er nicht der Mittelpunkt $\frac12(b+b')$ zweier verschiedener Punkte
$b,b'\in K$ ist (gleichbedeutend: nicht von der Form $tb+(1-t)b'$ mit $b\ne b'$ und $0<t<1$; ein solches $a$ ist der
Mittelpunkt von $a\pm\varepsilon(b-b')\in K$ für kleines $\varepsilon>0$).

**Lemma F.14.** In $\mathcal M$ gilt:

1. Die Extremalpunkte der konvexen Menge $[0,1]$ sind genau die Projektionen.
2. Für Projektionen $p,q$: $pq=0$ genau dann, wenn $p+q\le1$, genau dann, wenn $p\circ q=0$.
3. Für Projektionen $p,q$: $p\le q$ genau dann, wenn $pq=p$, genau dann, wenn $p\circ q=p$. Die minimalen von null
   verschiedenen Projektionen sind die Rang-eins-Projektionen (in einem Block).
4. Eine Projektion $p$ ist genau dann zentral, wenn $pxp+(1-p)x(1-p)=x$ für alle $x\in\mathcal M$. Die zentralen
   Projektionen sind die Summen der Blockeinheiten $P_c$, und die minimalen von null verschiedenen unter ihnen sind die
   $P_c$.

*Beweis.* (1) Sei $p$ eine Projektion und $p=\frac12(a+b)$ mit $a,b\in[0,1]$. Für $\xi\in\ker p$ ist
$0=\langle\xi,a\xi\rangle+\langle\xi,b\xi\rangle$ mit zwei Termen $\ge0$, also $\|a^{1/2}\xi\|^2=0$ und $a\xi=0$; für $\xi$ im
Bild von $p$ liefert dasselbe mit $1-a,1-b$ die Gleichung $a\xi=\xi$. Also $a=p$. Ist umgekehrt $a\in[0,1]$ keine
Projektion, so hat es einen Eigenwert $\lambda\in(0,1)$; mit der Spektralprojektion $q_\lambda\in\mathcal M$ (einem
Polynom in $a$) und $\varepsilon=\min(\lambda,1-\lambda)$ ist $a\pm\varepsilon q_\lambda\in[0,1]$, und $a$ ist ihr
Mittelpunkt.
(2) Ist $pq=0$, so ist $qp=(pq)^*=0$, $p+q$ ist eine Projektion, und $1-(p+q)=(1-p-q)^2\ge0$; außerdem $p\circ q=0$. Ist
$p+q\le1$, so gilt für jedes $\xi$: $\|qp\xi\|^2=\langle p\xi,qp\xi\rangle\le\langle p\xi,(1-p)p\xi\rangle=0$, also $qp=0=pq$.
Ist $p\circ q=0$, so ist $pq=-qp$, also $pq=p(pq)=-pqp=(qp)p=qp$, und $pq=-pq=0$.
(3) Aus $pq=p$ folgt $qp=p$, $p\circ q=p$ und $q-p=(q-p)^2\ge0$. Ist $p\circ q=p$, so liefert Multiplikation von $pq+qp=2p$
mit $p$ von links und von rechts $pq+pqp=2p=pqp+qp$, also $pq=qp=p$. Ist $p\le q$, so ist $p+(1-q)\le1$, und (2) liefert
$p(1-q)=0$. Eine Projektion vom Rang $\ge2$ dominiert echt die Projektion $\xi\xi^*$ auf eine Gerade in ihrem Bild
innerhalb eines Blocks ($\xi$ ein Einheitsvektor im Bild eines Blocks von $p$). Eine Rang-eins-Projektion $p$ liegt in
einem Block, und ist $0\ne q\le p$, so hat $q=pq$ sein Bild im Bild von $p$, also $q=p$.
(4) $x=pxp+px(1-p)+(1-p)xp+(1-p)x(1-p)$, und $px(1-p)+(1-p)xp=0$ gilt genau dann, wenn beide Terme verschwinden (man
multipliziere mit $p$ und $1-p$), also genau dann, wenn $px=pxp=xp$. Ein zentrales Element kommutiert mit allen
Matrixeinheiten jedes Blocks, ist also auf Block $c$ ein Skalar $z_c$, und $z_c\in\{0,1\}$ für eine Projektion. $\square$

**Satz F.15 (Kadison, endliche Dimension).** Sei $\Phi:\mathcal M\to\mathcal N$ eine lineare Bijektion mit
$\Phi(1)=1$, sodass $\Phi$ und $\Phi^{-1}$ positiv sind. Dann ist $\Phi$ ein Jordan-$*$-Isomorphismus. Umgekehrt ist jeder
Jordan-$*$-Isomorphismus unital, und er und seine Inverse sind positiv.

*Beweis.* Jedes $x\in\mathcal M_{sa}$ ist $x_+-x_-$ mit $x_\pm\ge0$, also ist $\Phi(x)$ selbstadjungiert; schreibt man
$x=h+ik$ mit selbstadjungierten $h,k$, so ist $\Phi(x^*)=\Phi(h)-i\Phi(k)=\Phi(x)^*$. Dasselbe gilt für $\Phi^{-1}$, also
ist $\Phi$ eine reell-lineare Bijektion $\mathcal M_{sa}\to\mathcal N_{sa}$ mit $x\ge0\iff\Phi(x)\ge0$. Da $\Phi(1)=1$,
bildet $\Phi$ das Intervall $[0,1]$ auf $[0,1]$ ab. Eine lineare Bijektion zwischen konvexen Mengen bildet
Extremalpunkte auf Extremalpunkte ab (ist $\Phi(a)=\frac12(b+b')$, so ist $a=\frac12(\Phi^{-1}b+\Phi^{-1}b')$), also bildet
$\Phi$ nach Lemma F.14(1) Projektionen auf Projektionen ab, und nach Lemma F.14(2) gilt
$pq=0\iff p+q\le1\iff\Phi(p)+\Phi(q)\le1\iff\Phi(p)\Phi(q)=0$. Für $x=\sum_i\lambda_ip_i\in\mathcal M_{sa}$ (Spektralzerlegung,
$p_ip_j=0$ für $i\ne j$) ist $\Phi(x)^2=\sum_i\lambda_i^2\Phi(p_i)=\Phi(x^2)$. Polarisation,
$x\circ y=\frac12((x+y)^2-x^2-y^2)$, liefert $\Phi(x\circ y)=\Phi(x)\circ\Phi(y)$ auf $\mathcal M_{sa}$, und komplexe
Bilinearität auf $\mathcal M=\mathcal M_{sa}+i\mathcal M_{sa}$.

Umgekehrt ist $\Phi(1)\circ\Phi(y)=\Phi(y)$ für alle $y$, und $\Phi$ ist surjektiv, also ist $\Phi(1)$ eine Eins für
$\circ$, und $\Phi(1)=\Phi(1)\circ1=1$. Jedes $x\ge0$ ist $h^2$ mit selbstadjungiertem $h=\sqrt x$ (Satz 1.6), und
$\Phi(h^2)=\Phi(h)^2\ge0$, da $\Phi(h)$ selbstadjungiert ist; $\Phi^{-1}$ ist ebenfalls ein Jordan-$*$-Isomorphismus. $\square$

Kadison hat den Satz für alle unitalen C\*-Algebren bewiesen (R. V. Kadison, *A generalized Schwarz inequality and
algebraic invariants for operator algebras*, Ann. of Math. 56 (1952) 494–503); diese Fassung wird im Kurs nicht
verwendet.

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

## F.5 Jordan-$*$-Isomorphismen von $\bigoplus_cM_{k_c}$

**Lemma F.16.** Sei $\Phi:\mathcal M\to\mathcal N$ ein Jordan-$*$-Isomorphismus. Dann ist $\Phi(1)=1$,
$\Phi(xyx)=\Phi(x)\Phi(y)\Phi(x)$, und $\Phi$ bildet die Projektionen von $\mathcal M$ bijektiv auf die von $\mathcal N$ ab,
wobei $pq=0$, $p\le q$, Minimalität, Zentralität und Minimalität unter den zentralen Projektionen erhalten und
reflektiert werden.

*Beweis.* $\Phi(1)=1$ wie in Satz F.15. In jeder assoziativen Algebra ist $xyx=2x\circ(x\circ y)-(x\circ x)\circ y$, was die
zweite Aussage liefert. $p=p^*=p\circ p$ bleibt erhalten, ebenso die Charakterisierungen aus Lemma F.14(2,3,4), die nur
$\circ$, das Tripelprodukt $xyx$, $1$ und Quantoren über die ganze Algebra verwenden; $\Phi^{-1}$ ist ebenfalls ein
Jordan-$*$-Isomorphismus. $\square$

**Satz F.17 (Jordan-$*$-Automorphismen von $M_n$).** Jeder Jordan-$*$-Isomorphismus $\Phi:M_n(\mathbb C)\to M_n(\mathbb C)$
ist $x\mapsto UxU^*$ oder $x\mapsto Ux^TU^*$ mit einem unitären $U$. Für $n\ge2$ sind die ersten Abbildungen
multiplikativ und die zweiten antimultiplikativ, also ist keine Abbildung von beiden Formen.

*Beweis.* Für $n=1$ ist $\Phi=\mathrm{id}$. Sei $n\ge2$ und seien $E_{ij}$ die Matrixeinheiten. Nach Lemma F.16 sind die
$\Phi(E_{ii})$ paarweise orthogonale Rang-eins-Projektionen, $\Phi(E_{ii})=\xi_i\xi_i^*$ mit orthonormalen $\xi_i$; mit dem
Unitären $W$, $We_i=\xi_i$, erfüllt der Jordan-$*$-Isomorphismus $\Phi_1(x)=W^*\Phi(x)W$ die Gleichungen
$\Phi_1(E_{ii})=E_{ii}$.

Sei $i\ne j$. Der Raum $S_{ij}=\operatorname{span}\{E_{ij},E_{ji}\}$ ist $\{x:pxp=x,\ E_{ii}xE_{ii}=E_{jj}xE_{jj}=0\}$ mit
$p=E_{ii}+E_{jj}$, also $\Phi_1(S_{ij})\subseteq S_{ij}$ nach Lemma F.16. Man schreibe $\Phi_1(E_{ij})=\alpha E_{ij}+\beta E_{ji}$.
Aus $E_{ij}\circ E_{ij}=0$ folgt $0=\Phi_1(E_{ij})^2=\alpha\beta(E_{ii}+E_{jj})$, also ist genau eines von $\alpha,\beta$
von null verschieden ($\Phi_1$ ist injektiv). Sei $\varepsilon(i,j)=+$, falls $\Phi_1(E_{ij})\in\mathbb CE_{ij}$, und $-$,
falls $\Phi_1(E_{ij})\in\mathbb CE_{ji}$; $\varepsilon(j,i)=\varepsilon(i,j)$, da $\Phi_1(E_{ji})=\Phi_1(E_{ij})^*$. Für
verschiedene $i,j,k$ ist $E_{ij}\circ E_{jk}=\frac12E_{ik}\ne0$; ist $\varepsilon(i,j)=+$ und $\varepsilon(j,k)=-$, so ist
$\Phi_1(E_{ij})\circ\Phi_1(E_{jk})\in\mathbb C\,E_{ij}\circ E_{kj}=\{0\}$, ein Widerspruch; ebenso für $(-,+)$. Also ist
$\varepsilon(i,j)=\varepsilon(j,k)$, wann immer $i,j,k$ verschieden sind, und $\varepsilon$ ist konstant: Zwei Paare haben
entweder einen Index gemeinsam oder sind über ein drittes Paar verbunden. Ist $\varepsilon=-$, so ersetze man $\Phi_1$
durch $\Phi_1\circ T$ ($T$ die Transposition, ein Jordan-$*$-Automorphismus nach Proposition 13.9), das $\varepsilon=+$
hat und die $E_{ii}$ weiterhin fixiert.

Nun ist $\Phi_1(E_{ij})=\alpha_{ij}E_{ij}$ für alle $i,j$, mit $\alpha_{ii}=1$ und $\alpha_{ji}=\bar\alpha_{ij}$. Aus
$E_{ij}\circ E_{ji}=\frac12(E_{ii}+E_{jj})$ folgt $|\alpha_{ij}|^2=1$; aus $E_{ij}\circ E_{jk}=\frac12E_{ik}$ ($i,j,k$
verschieden) folgt $\alpha_{ij}\alpha_{jk}=\alpha_{ik}$. Mit $u_i=\alpha_{i1}$ ergibt das $\alpha_{ij}=\alpha_{i1}\alpha_{1j}=u_i\bar u_j$
(auch für $i=j$, $i=1$ oder $j=1$). Also ist $\Phi_1(x)=DxD^*$ mit $D=\operatorname{diag}(u_i)$, und $\Phi(x)=WDxD^*W^*$
oder $\Phi(x)=WDx^TD^*W^*$. Stimmten beide Formen überein, so wäre $T$ multiplikativ, aber $(E_{12}E_{21})^T=E_{11}\ne E_{22}=E_{12}^TE_{21}^T$.
$\square$

Das beweist das in Abschnitt 13.5 zitierte Ergebnis von N. Jacobson, C. E. Rickart (Trans. Amer. Math. Soc. 69 (1950)
479–502) und I. N. Herstein (Trans. Amer. Math. Soc. 81 (1956) 331–341) im endlichdimensionalen Fall, für alle $n$;
Proposition 13.11 ist der Fall $n=2$. Ein Jordan-Automorphismus von $\operatorname{Herm}_n(\mathbb C)$ setzt sich durch
$\Phi(x+iy)=\Phi(x)+i\Phi(y)$ zu einem Jordan-$*$-Automorphismus von $M_n$ fort, ist also ebenfalls von einer der beiden
Formen.

**Satz F.18 (direkte Summen).** Sei $\Phi:\bigoplus_{c=1}^CM_{k_c}\to\bigoplus_{d=1}^{C'}M_{k'_d}$ ein
Jordan-$*$-Isomorphismus. Dann ist $C'=C$, und es gibt eine Bijektion $\sigma$ mit $k'_{\sigma(c)}=k_c$, Unitäre
$U_c\in M_{k_c}$ und eine Menge $S$ von Blöcken, sodass
$$
\Phi(x)_{\sigma(c)}=U_cx_cU_c^*\ (c\notin S),\qquad \Phi(x)_{\sigma(c)}=U_cx_c^TU_c^*\ (c\in S).
$$
Umgekehrt ist jede solche Abbildung ein Jordan-$*$-Isomorphismus, und $S$ ist bis auf Blöcke mit $k_c=1$ eindeutig.

*Beweis.* Nach Lemma F.16 und Lemma F.14(4) ist $\Phi(P_c)=P'_{\sigma(c)}$ für eine Bijektion $\sigma$ zwischen den
minimalen zentralen Projektionen, also $C=C'$. Ist $x=P_cxP_c$, so ist $\Phi(x)=P'_{\sigma(c)}\Phi(x)P'_{\sigma(c)}$, also
bildet $\Phi$ Block $c$ in Block $\sigma(c)$ ab, und $\Phi^{-1}$ bildet Block $\sigma(c)$ in Block $c$ ab; daher ist
$k_c^2=k'^2_{\sigma(c)}$, und die Einschränkung ist ein Jordan-$*$-Isomorphismus $M_{k_c}\to M_{k_c}$, nach Satz F.17 von
einer der beiden Formen. Da $x=\sum_cP_cxP_c$, bestimmt das $\Phi$. Die Umkehrung gilt blockweise (Proposition 13.9), und
die Eindeutigkeit von $S$ folgt aus Satz F.17. $\square$

**Korollar F.19 (unitale Ordnungsisomorphismen).** Die unitalen linearen Bijektionen $\Phi:\mathcal M\to\mathcal N$, für
die $\Phi$ und $\Phi^{-1}$ positiv sind, sind genau die Abbildungen aus Satz F.18. Ein solches $\Phi$ ist genau dann
$2$-positiv, wenn es vollständig positiv ist, genau dann, wenn $S$ keinen Block mit $k_c\ge2$ enthält, genau dann,
wenn $\Phi$ ein $*$-Isomorphismus ist.

*Beweis.* Der erste Satz folgt aus den Sätzen F.15 und F.18. Enthält $S$ keinen Block mit $k_c\ge2$, so ist jede
Blockabbildung $x\mapsto U_cxU_c^*$ (für $k_c=1$ ist $x^T=x$), also ist $\Phi$ vollständig positiv (Proposition 6.4) und
multiplikativ. Ist $c\in S$ mit $k_c\ge2$ und $\Phi$ $2$-positiv, so ist es auch die Blockabbildung
$\Phi_{\sigma(c)c}(a)=U_ca^TU_c^*$ (Proposition 6.4), und ebenso $a\mapsto U_c^*\Phi_{\sigma(c)c}(a)U_c=a^T$, da aus
$[A_{ij}]\ge0$ folgt $[U^*A_{ij}U]=(1\otimes U)^*[A_{ij}](1\otimes U)\ge0$; das widerspricht Proposition 6.1. Ein solches
$\Phi$ ist nach Satz F.17 auch nicht multiplikativ. $\square$

Ein unitaler Ordnungsisomorphismus ist also bis auf die Unitären dadurch bestimmt, welche gleichen Blöcke er vertauscht
und welche Blöcke er transponiert, und nach Korollar F.19 sind es genau die transponierten Blöcke, die vollständige
Positivität verhindern. Der Code konstruiert den Jordan-Automorphismus $x\mapsto Ux^TU^*$ von $M_3$ und gewinnt $U$
(bis auf eine Phase) mit den Schritten des Beweises von Satz F.17 zurück:

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

## Aufgaben

1. Man zeige, dass $\mathbb R^N$ mit dem koordinatenweisen Produkt und dem Standardskalarprodukt eine euklidische
   Jordan-Algebra ist, deren Kegel der Quadrate $\mathbb R^N_{\ge0}$ ist, und folgere aus Korollar F.19, dass die unitalen
   Ordnungsisomorphismen von $\mathbb C^N=\bigoplus_{c=1}^NM_1$ die Permutationen der Koordinaten sind.
2. (**Positivität genügt nicht.**) Für $0<t<1$ sei $\Psi_t(x)=(1-t)x+\frac tn\operatorname{Tr}(x)1$ auf $M_n$, $n\ge2$.
   Man zeige, dass $\Psi_t$ eine unitale positive Bijektion ist, dass $\Psi_t^{-1}$ nicht positiv ist und dass $\Psi_t$
   Projektionen nicht auf Projektionen abbildet.

<details>
<summary>Lösungen</summary>

1. Das Produkt ist kommutativ und assoziativ, also ein Jordan-Produkt, mit Eins $(1,\dots,1)$, und $L(x)$ ist die
   Diagonalmatrix $\operatorname{diag}(x)$, die symmetrisch ist. Die Quadrate sind die Vektoren mit nichtnegativen
   Koordinaten. Nach Korollar F.19 mit allen $k_c=1$ ist $\Phi(x)_{\sigma(c)}=U_cx_cU_c^*=x_c$ (Unitäre in $M_1$ sind
   Phasen), also permutiert $\Phi$ die Koordinaten; umgekehrt sind Permutationen unitale Ordnungsisomorphismen.

2. $\Psi_t(1)=(1-t)1+t1=1$; für $x\ge0$ sind beide Terme $\ge0$. Da $\operatorname{Tr}\Psi_t(x)=\operatorname{Tr}x$, ist
   $\Psi_t^{-1}(y)=\frac1{1-t}\big(y-\frac tn\operatorname{Tr}(y)1\big)$. Für eine Rang-eins-Projektion $p$ hat
   $\Psi_t^{-1}(p)$ auf $\ker p$ den Eigenwert $-\frac{t}{n(1-t)}<0$, und $\Psi_t(p)$ hat den Eigenwert
   $\frac tn\notin\{0,1\}$. Das ist mit Satz F.15 verträglich: Die Positivität der Inversen ist es, die Projektionen auf
   Projektionen zwingt.

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

**Zurück zu:** [Kapitel 13](13_jordan_algebras.md), Abschnitte 13.4–13.6: Satz F.12 beweist Satz 13.8, Satz F.15
beweist Satz 13.10 für die Algebren $\bigoplus_cM_{k_c}(\mathbb C)$, Satz F.17 erweitert Proposition 13.11 auf alle
$n$, und Korollar F.7 erweitert Proposition 13.14 auf jede euklidische Jordan-Algebra.
