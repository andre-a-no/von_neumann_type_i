[English](../08_dynamics.md) · [Español](../es/08_dynamics.md) · [Français](../fr/08_dynamics.md) · **Deutsch** · [中文](../zh/08_dynamics.md) · [日本語](../ja/08_dynamics.md) · [Русский](../ru/08_dynamics.md)

# 8. Dynamik

In der Quantenmechanik wird ein abgeschlossenes System durch eine einparametrige Gruppe von Unitären
modelliert, ein offenes – ein System, das an eine nicht beobachtete Umgebung gekoppelt ist – durch
eine Halbgruppe von Kanälen. In einer Algebra mit Sektoren haben beide eine Eigenschaft, die auf einer
einzelnen Matrixalgebra unsichtbar bleibt: Ein Hamiltonoperator in $`\mathcal M`$ kann den Sektor
niemals ändern, die Ladungen sind also erhalten, während ein offenes System zwischen Sektoren
springen kann, etwa indem es ein Teilchen verliert. Dieses Kapitel beweist die grundlegenden
Tatsachen über beide Arten der Zeitentwicklung, stellt die Lindblad-Gleichung mit Sprüngen zwischen
Sektoren auf und zeigt, dass „kein Teilchen ging verloren“ eine Messung der Ladung ist, deren
postselektierter Zustand sich gemäß einem nicht-hermiteschen Hamiltonoperator entwickelt.

**Sie lernen**, warum die Sektorgewichte Erhaltungsgrößen sind, warum Generatoren in
Gorini-Kossakowski-Sudarshan-Lindblad-Form (GKSL-Form) Kanäle erzeugen, wie Sprünge zwischen
Sektoren in die Lindblad-Gleichung eingehen und wann die Sektorgewichte einer klassischen
Ratengleichung genügen, was die sprungfreie Entwicklung und die Postselektion bedeuten und wie man
stationäre Zustände und Relaxationsraten bestimmt. Die numerischen Kontrollen dieses Kapitels stehen im
[Tutorial 8](../../tutorial/de/08_dynamics.md).

Im Folgenden wirkt $`\mathcal M=\bigoplus_{c=1}^CM_{k_c}(\mathbb C)`$ auf
$`\mathcal H=\bigoplus_c\mathbb C^{k_c}`$, es ist $`D=\sum_ck_c`$, die $`P_c`$ sind die minimalen
zentralen Projektionen, Zustände werden durch ihre ungewichteten Dichten $`\rho`$ mit den
Sektorgewichten $`p_c=\mathrm{Tr}\,\rho_c`$ angegeben (Kapitel 3), $`[A,B]=AB-BA`$, $`\{A,B\}=AB+BA`$
und $`\hbar=1`$. Ein **Superoperator** ist eine lineare Abbildung auf $`\mathcal M`$ oder auf
$`M_D(\mathbb C)`$; alle Normen von Superoperatoren sind Operatornormen und daher submultiplikativ.
Für eine lineare Abbildung $`B`$ auf einem endlichdimensionalen Raum ist
$`e^{tB}=\sum_{n\ge0}t^nB^n/n!`$.

<a id="s-8-1"></a>

## 8.1 Abgeschlossene Systeme

Ein **Hamiltonoperator** ist ein selbstadjungiertes $`H\in\mathcal M`$. Ein Vektor entwickelt sich
gemäß der **Schrödinger-Gleichung** $`i\psi'(t)=H\psi(t)`$, eine Dichte gemäß der
**von-Neumann-Gleichung** $`\rho'(t)=-i[H,\rho(t)]`$. Beide sind lineare Gleichungen mit konstanten
Koeffizienten auf einem endlichdimensionalen Raum, und das folgende Lemma ist alles, was wir aus der
Theorie der Differentialgleichungen brauchen.

<a id="r-8-1"></a>

**Lemma 8.1.** Sei $`B`$ eine lineare Abbildung auf einem endlichdimensionalen normierten Raum $`X`$.
Dann gilt $`\frac d{dt}e^{tB}=Be^{tB}=e^{tB}B`$, $`e^{(s+t)B}=e^{sB}e^{tB}`$,
$`\|e^{tB}\|\le e^{|t|\,\|B\|}`$, und für jedes $`x_0\in X`$ ist $`x(t)=e^{tB}x_0`$ die einzige
differenzierbare Funktion $`x`$ auf einem Intervall $`I\ni0`$ mit $`x'=Bx`$ und $`x(0)=x_0`$.

*Beweis.* Die Reihe konvergiert absolut, mit der angegebenen Normschranke, und als Potenzreihe in
$`t`$ darf sie gliedweise differenziert werden; das liefert die erste Formel. Dann ist
$`\frac d{dt}(e^{tB}e^{-tB})=Be^{tB}e^{-tB}-e^{tB}Be^{-tB}=0`$, also $`e^{tB}e^{-tB}=1`$. Gilt $`x'=Bx`$,
so ist $`\frac d{dt}(e^{-tB}x(t))=e^{-tB}(x'-Bx)=0`$, also $`e^{-tB}x(t)=x_0`$ und $`x(t)=e^{tB}x_0`$.
Schließlich lösen $`t\mapsto e^{(s+t)B}x_0`$ und $`t\mapsto e^{tB}e^{sB}x_0`$ die Gleichung $`y'=By`$ mit
demselben Anfangswert. $`\square`$

<a id="r-8-2"></a>

**Proposition 8.2.** Sei $`H=H^*\in\mathcal M`$ und $`U(t)=e^{-itH}`$, definiert durch den
Funktionalkalkül aus Abschnitt 1.3 mit $`f_t(\lambda)=e^{-it\lambda}`$.

1. $`U(t)=\bigoplus_ce^{-itH_c}`$ ist ein unitäres Element von $`\mathcal M`$, $`U(s+t)=U(s)U(t)`$ und
   $`U'(t)=-iHU(t)`$; $`U(t)`$ stimmt mit der Exponentialreihe von $`-itH`$ überein.
2. $`\psi(t)=U(t)\psi_0`$ ist die einzige Lösung von $`i\psi'=H\psi`$ mit $`\psi(0)=\psi_0`$, und
   $`\rho(t)=U(t)\rho_0U(t)^*`$ ist die einzige Lösung in $`\mathcal M`$ von $`\rho'=-i[H,\rho]`$ mit
   $`\rho(0)=\rho_0`$. Ist $`\rho_0`$ ein Zustand, so auch jedes $`\rho(t)`$.
3. (*Heisenberg-Bild*) $`\mathrm{Tr}(\rho(t)A)=\mathrm{Tr}(\rho_0A(t))`$ mit
   $`A(t)=U(t)^*AU(t)`$, und $`A'=i[H,A]`$.

*Beweis.* (1) Man diagonalisiere jeden Block, $`H_c=V_c\Lambda_cV_c^*`$; dann diagonalisiert
$`V=\bigoplus_cV_c`$ die Matrix $`H`$, und $`U(t)=Ve^{-it\Lambda}V^*`$ ist blockdiagonal (nach Satz 1.5(a)
hängt das Ergebnis nicht von der diagonalisierenden unitären Matrix ab). Wegen $`|f_t|^2=1`$ und
$`f_sf_t=f_{s+t}`$ liefert Satz 1.5(c) die Unitarität und das Gruppengesetz. Eintragsweises
Differenzieren von $`e^{-it\Lambda}`$ ergibt $`U'=V(-i\Lambda)e^{-it\Lambda}V^*=-iHU`$, also ist $`U(t)`$
nach Lemma 8.1 die Reihe $`e^{-itH}`$ ($`X=M_D`$, $`B`$ die Linksmultiplikation mit $`-iH`$). (2)
$`(U\psi_0)'=-iHU\psi_0`$, und $`(U\rho_0U^*)'=-iHU\rho_0U^*+U\rho_0U^*(iH)`$, da $`(U')^*=iU^*H`$. Die
Eindeutigkeit ist Lemma 8.1 mit $`B=-iH`$ auf $`\mathcal H`$ bzw. $`B=-i[H,\cdot\,]`$ auf $`\mathcal M`$.
Positivität: $`\langle x,U\rho_0U^*x\rangle=\langle U^*x,\rho_0U^*x\rangle\ge0`$; die Spur wegen der
Zyklizität. (3) Zyklizität der Spur, und $`(U^*AU)'=iU^*HAU-iU^*AHU=i[H,A(t)]`$. $`\square`$

*In der Bibliothek:* $`U(t)=e^{-itH}`$, Block für Block, ist `dynamics.propagator`; die Lösungen aus Proposition 8.2(2) sind `dynamics.schrodinger` und `dynamics.von_neumann` — [Tutorial, Abschnitt 8.1](../../tutorial/de/08_dynamics.md#s-8-1).

<a id="r-8-3"></a>

**Satz 8.3 (Ladungserhaltung).** Sei $`t\mapsto H(t)`$ eine Abbildung auf einem Intervall $`I`$ mit
Werten in den selbstadjungierten Elementen von $`\mathcal M`$.

1. Ist $`\psi:I\to\mathcal H`$ differenzierbar mit $`i\psi'=H(t)\psi`$, so ist jedes $`\|P_c\psi(t)\|`$
   konstant.
2. Ist $`\rho:I\to\mathcal M`$ differenzierbar mit $`\rho'=-i[H(t),\rho]`$, so ist
   $`\mathrm{Tr}(Z\rho(t))`$ für jedes zentrale $`Z`$ konstant; insbesondere ist jedes
   Sektorgewicht $`p_c(t)=\mathrm{Tr}(P_c\rho(t))`$ konstant.
3. Ist $`H`$ konstant und kommutiert $`A\in\mathcal M`$ mit $`H`$, so ist $`\mathrm{Tr}(A\rho(t))`$
   konstant.

*Beweis.* (1) Das Skalarprodukt ist im ersten Argument antilinear, also gilt
$`\frac d{dt}\langle\psi,P_c\psi\rangle=\langle-iH\psi,P_c\psi\rangle+\langle\psi,-iP_cH\psi\rangle =i\langle\psi,HP_c\psi\rangle-i\langle\psi,P_cH\psi\rangle=0`$,
da $`P_c`$ zentral ist. (2)
$`\frac d{dt}\mathrm{Tr}(Z\rho)=-i\mathrm{Tr}(Z[H,\rho])=-i\mathrm{Tr}([Z,H]\rho)=0`$
wegen der Zyklizität. (3) Dieselbe Rechnung mit $`A`$. $`\square`$

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Ein Hamiltonoperator
> in $`\mathcal M`$ kommutiert mit den Ladungen (Teilchenzahl, Magnetisierung, eine
> Symmetriequantenzahl), daher ist die Wahrscheinlichkeit jeder Ladung eine Erhaltungsgröße, selbst
> für getriebene, zeitabhängige Hamiltonoperatoren.

Auf $`\mathcal M`$ definiert ein Einheitsvektor $`\psi=\bigoplus_c\psi_c`$ den Zustand der Mischung
$`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$: Jedes $`A\in\mathcal M`$ ist blockdiagonal, also ist
$`\langle\psi,A\psi\rangle=\sum_c\langle\psi_c,A_c\psi_c\rangle=\mathrm{Tr}\big(A\bigoplus_c|\psi_c\rangle\langle\psi_c|\big)`$
(wie in Kapitel 3, Aufgabe 3). Löst $`\psi(t)`$ die Gleichung $`i\psi'=H\psi`$, so ist
$`\bigoplus_c|\psi_c(t)\rangle\langle\psi_c(t)|=U(t)\big(\bigoplus_c|\psi_c(0)\rangle\langle\psi_c(0)|\big)U(t)^*`$, da $`U(t)`$
blockdiagonal ist.

*Rechnung:* Satz 8.3 und die Übereinstimmung der Entwicklungen von $`\psi`$ und von $`\bigoplus_c|\psi_c\rangle\langle\psi_c|`$, mit den exakten Lösern und mit dem Runge-Kutta-Löser `dynamics.schrodinger_rk4` für einen getriebenen Hamiltonoperator — [Tutorial, Abschnitt 8.1](../../tutorial/de/08_dynamics.md#s-8-1).

<a id="s-8-2"></a>

## 8.2 Quantendynamische Halbgruppen

CP- und TP-Abbildungen sowie Kanäle sind in Abschnitt 6.1 definiert ([Kapitel
6](06_quantum_channels.md)). Neben Satz 6.2 und Proposition 6.4 (eine Abbildung
$`A\mapsto\sum_lK_lAK_l^*`$, auch zwischen Sektoren, ist CP) benötigen wir drei
Abgeschlossenheitseigenschaften.

<a id="r-8-4"></a>

**Lemma 8.4.** (a) Summen, nichtnegative Vielfache und Kompositionen von CP-Abbildungen sind CP. Sei
$`\Phi`$ CP auf $`\mathcal A=M_D`$ oder auf $`\mathcal A=\mathcal M`$, und sei $`\mathcal N\subseteq\mathcal A`$
gleich $`\mathcal M`$ (falls $`\mathcal A=M_D`$) oder eine der Ecken $`P_c\mathcal MP_c\cong M_{k_c}`$ (mit
Einselement $`P_c`$). Bildet $`\Phi`$ die Algebra $`\mathcal N`$ in sich ab, so ist $`\Phi|_{\mathcal N}`$ CP.
(b) Konvergieren CP-Abbildungen $`\Phi_m`$ gegen $`\Phi`$ (d. h. $`\Phi_m(A)\to\Phi(A)`$ für alle
$`A`$), so ist $`\Phi`$ CP. (c) Ist $`s\mapsto\Psi(s)`$ eine stetige Abbildung von $`[a,b]`$ in die
CP-Abbildungen, so ist $`\int_a^b\Psi(s)\,ds`$ CP.

*Beweis.* (a) $`\Phi\mapsto\Phi_n`$ ist linear und $`(\Phi\Psi)_n=\Phi_n\Psi_n`$; Summen positiver
Matrizen sind positiv. Zur Einschränkung: $`M_n(\mathcal N)\subseteq M_n(\mathcal A)`$, und die
Positivität ist in beiden dieselbe Bedingung. Denn nach Abschnitt 6.1 ist ein Element $`[A_{ij}]`$ von
$`M_n(\mathcal M)`$ genau dann positiv, wenn jede Blockmatrix $`[(A_{ij})_c]`$ es ist; für
$`\mathcal A=M_D`$ ist das die Positivität von $`[A_{ij}]\in M_{nD}`$, weil $`[A_{ij}]`$ nach einer
Permutation der Basis die direkte Summe dieser Blöcke ist (alle $`A_{ij}`$ sind blockdiagonal); und ein
Element von $`M_n(P_c\mathcal MP_c)`$ hat den einzigen von null verschiedenen Block $`[(A_{ij})_c]`$,
dessen Positivität die Positivität in $`M_n(M_{k_c})`$ ist. (b) Die positiven Matrizen bilden eine
abgeschlossene Menge ($`\langle x,Xx\rangle\ge0`$ für jedes $`x`$ ist eine abgeschlossene Bedingung), und
$`(\Phi_m)_n(X)\to\Phi_n(X)`$. (c) Die Riemann-Summen sind nach (a) CP und konvergieren gegen das
Integral (Lemma A.4(a)); man wende (b) an. $`\square`$

**Definition.** Eine **quantendynamische Halbgruppe** auf $`\mathcal M`$ ist eine Familie
$`(T_t)_{t\ge0}`$ von Kanälen $`\mathcal M\to\mathcal M`$ mit $`T_0=\mathrm{id}`$, $`T_{s+t}=T_sT_t`$ und
stetigem $`t\mapsto T_t`$.

<a id="r-8-5"></a>

**Lemma 8.5.** Jede stetige Halbgruppe $`(T_t)_{t\ge0}`$ von Superoperatoren ($`T_0=\mathrm{id}`$,
$`T_{s+t}=T_sT_t`$) ist von der Form $`T_t=e^{t\mathcal L}`$ für genau einen Superoperator $`\mathcal L`$,
ihren **Generator**, und es gilt $`\mathcal L=\lim_{t\to0}(T_t-\mathrm{id})/t`$. Dasselbe gilt, mit
demselben Beweis, für stetige Halbgruppen linearer Abbildungen auf einem beliebigen
endlichdimensionalen Raum (zum Beispiel von Matrizen).

*Beweis.* (Integrale stetiger superoperatorwertiger Funktionen sind die aus Lemma A.4.) Es gilt
$`\|V_h-\mathrm{id}\|\le\max_{0\le s\le h}\|T_s-\mathrm{id}\|\to0`$ für $`V_h=\frac1h\int_0^hT_s\,ds`$
(Lemma A.4(c)), also ist $`\det V_h\neq0`$ für ein festes kleines $`h`$, wegen der Stetigkeit der
Determinante. Nun ist $`T_tV_h=\frac1h\int_0^hT_{t+s}\,ds=\frac1h\int_t^{t+h}T_s\,ds`$ (Lemma
A.4(b),(e)) in $`t`$ differenzierbar
mit der Ableitung $`\frac1h(T_{t+h}-T_t)=T_t\frac{T_h-\mathrm{id}}h`$. Folglich ist $`T_t`$
differenzierbar und $`T_t'=T_t\mathcal L`$ mit $`\mathcal L=\frac{T_h-\mathrm{id}}hV_h^{-1}`$, und nach
Lemma 8.1 (mit $`B`$ der Komposition von rechts mit $`\mathcal L`$) ist $`T_t=e^{t\mathcal L}`$.
Eindeutigkeit: $`\mathcal L=T'(0)`$. $`\square`$

<a id="r-8-6"></a>

**Satz 8.6 (GKSL-Generatoren).** Seien $`H=H^*`$ und $`L_1,\dots,L_J`$ in $`M_D`$, und sei
```math
\mathcal L(\rho)=-i[H,\rho]+\sum_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big).
```
Dann ist $`(e^{t\mathcal L})_{t\ge0}`$ eine quantendynamische Halbgruppe auf $`M_D`$. Genauer: Man setze
$`H_{\mathrm{eff}}=H-\frac i2\sum_jL_j^*L_j`$,
$`\mathcal L_0(\rho)=-i(H_{\mathrm{eff}}\rho-\rho H_{\mathrm{eff}}^*)`$ und
$`\Phi(\rho)=\sum_jL_j\rho L_j^*`$, sodass $`\mathcal L=\mathcal L_0+\Phi`$, und definiere
```math
S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},\qquad
S_n(t)=\int_0^tS_0(t-s)\,\Phi\,S_{n-1}(s)\,ds\quad(n\ge1).
```
Dann ist jedes $`S_n(t)`$ CP, und $`e^{t\mathcal L}=\sum_{n\ge0}S_n(t)`$ (die **Dyson-Entwicklung**),
mit $`\|S_n(t)\|\le e^{t\|\mathcal L_0\|}(t\|\Phi\|)^n/n!`$.

*Beweis.* *Schritt 1.* $`X(t)=e^{-itH_{\mathrm{eff}}}\rho e^{itH_{\mathrm{eff}}^*}`$ erfüllt
$`X'=\mathcal L_0(X)`$ und $`X(0)=\rho`$, also ist $`S_0(t)=e^{t\mathcal L_0}`$ nach Lemma 8.1; es ist CP
mit dem einzigen Kraus-Operator $`e^{-itH_{\mathrm{eff}}}`$. *Schritt 2.* Man schreibe
$`a=\|\mathcal L_0\|`$, $`b=\|\Phi\|`$. Durch Induktion sind die $`S_n`$ stetig, und es gilt
$`\|S_n(t)\|\le\int_0^te^{(t-s)a}\,b\,e^{sa}\frac{(bs)^{n-1}}{(n-1)!}ds=e^{ta}\frac{(bt)^n}{n!}`$; die
Reihe konvergiert also gleichmäßig auf jedem $`[0,T]`$. *Schritt 3.* Da
$`S_n(t)=S_0(t)\int_0^te^{-s\mathcal L_0}\Phi S_{n-1}(s)\,ds`$, liefern die Produktregel und der
Hauptsatz der Differential- und Integralrechnung (Lemma A.4(e)) $`S_n'=\mathcal L_0S_n+\Phi S_{n-1}`$ ($`n\ge1`$) und
$`S_0'=\mathcal L_0S_0`$. Wegen der Schranken konvergiert die Reihe der Ableitungen gleichmäßig auf
$`[0,T]`$, also ist $`Y(t)=\sum_nS_n(t)`$ differenzierbar mit $`Y'=\mathcal L_0Y+\Phi Y=\mathcal LY`$,
$`Y(0)=\mathrm{id}`$, und nach Lemma 8.1 ist $`Y(t)=e^{t\mathcal L}`$. *Schritt 4.* $`\Phi`$ ist CP (Satz
6.2), also ist durch Induktion der Integrand $`s\mapsto S_0(t-s)\Phi S_{n-1}(s)`$ eine stetige Familie
von CP-Abbildungen und $`S_n(t)`$ CP (Lemma 8.4); ebenso die Partialsummen und ihr Grenzwert
$`e^{t\mathcal L}`$. *Schritt 5.* Wegen der Zyklizität ist
$`\mathrm{Tr}\,\mathcal L(\rho)=\sum_j\big(\mathrm{Tr}(L_j\rho L_j^*)-\mathrm{Tr}(L_j^*L_j\rho)\big)=0`$,
folglich
$`\frac d{dt}\mathrm{Tr}\,e^{t\mathcal L}\rho=\mathrm{Tr}\,\mathcal L(e^{t\mathcal L}\rho)=0`$:
Jedes $`e^{t\mathcal L}`$ ist TP. Halbgruppengesetz und Stetigkeit folgen aus Lemma 8.1. $`\square`$

*Rechnung:* Satz 8.6 für beliebige Sprünge auf $`M_3`$, einschließlich der vollständigen Positivität des Sprunganteils $`\sum_{n\ge1}S_n(t)`$ — [Tutorial, Abschnitt 8.2](../../tutorial/de/08_dynamics.md#s-8-2). *In der Bibliothek:* $`e^{t\mathcal L}`$ ist `dynamics.lindblad_channel`, seine Matrix $`\mathcal L`$ ist `dynamics.lindblad_superoperator` — [Tutorial, Abschnitt 8.2](../../tutorial/de/08_dynamics.md#s-8-2).

<a id="r-8-7"></a>

**Satz 8.7 (Gorini-Kossakowski-Sudarshan; Lindblad).** (1) Umgekehrt hat der Generator jeder
quantendynamischen Halbgruppe auf $`M_D`$ die Form aus Satz 8.6, mit $`J\le D^2-1`$ Sprüngen. Man kann
$`\mathrm{Tr}\,H=0`$, $`\mathrm{Tr}\,L_j=0`$ und $`\mathrm{Tr}(L_j^*L_k)=0`$ für $`j\neq k`$
wählen; mit diesen Normierungen sind $`H`$ und $`\Phi=\sum_jL_j\cdot L_j^*`$ durch $`\mathcal L`$ bestimmt
und die $`L_j`$ bis auf $`L'_j=\sum_ku_{jk}L_k`$ mit einer Unitären $`u`$ (nachdem man die kürzere Familie
mit Nullen aufgefüllt hat).

(2) Der Generator jeder quantendynamischen Halbgruppe auf $`\mathcal M=\bigoplus_cM_{k_c}`$ ist
```math
\mathcal L(\rho)=\Phi(\rho)-\tfrac12\{\Phi^*(1),\rho\}-i[H,\rho]
```
mit einer CP-Abbildung $`\Phi:\mathcal M\to\mathcal M`$ (Kraus-Operatoren zwischen den Sektoren,
Proposition 6.4) und $`H=H^*\in\mathcal M`$: die Form aus Abschnitt 8.3, mit einem einzigen Sprung
$`\Phi_1=\Phi`$ zwischen Sektoren und der Rate $`1`$.

*Beweis.* (1) ist Satz D.4 zusammen mit den Eindeutigkeitsaussagen von Satz D.5(b),(c), und (2) ist
Satz D.7; beide werden in [Anhang D](D_gksl.md) bewiesen. $`\square`$

Der entscheidende Schritt (Satz D.3) ist, dass die Choi-Matrix des Generators einer Halbgruppe von
CP-Abbildungen auf dem orthogonalen Komplement des maximal verschränkten Vektors positiv ist; die
Halbgruppeneigenschaft allein, mit positiven Abbildungen $`T_t`$, würde nicht genügen (das Beispiel
nach Satz D.3). Für (2) sind die Abbildungen $`T_t\circ E`$ mit der bedingten Erwartung
$`E(X)=\sum_cP_cXP_c`$ CP auf $`M_D`$; die Choi-Matrix der Ableitung $`\mathcal L\circ E`$ bei $`t=0`$ ist
positiv auf dem orthogonalen Komplement der Vektoren $`(1\otimes P_c)\Omega`$, was eine Zerlegung von
$`\mathcal L\circ E`$ liefert (Lemma D.2), und diese wird durch $`E`$ wieder auf $`\mathcal M`$
komprimiert. Im Heisenberg-Bild (Korollar D.8) ist der Generator einer stetigen Halbgruppe unitaler
CP-Abbildungen auf $`\mathcal M`$ die Abbildung $`A\mapsto\Psi(A)+K^*A+AK`$ mit $`\Psi`$ CP und
$`K\in\mathcal M`$. Historisch geht (1) auf V. Gorini, A. Kossakowski, E. C. G. Sudarshan [*J. Math.
Phys.* 17, 821 (1976)] und G. Lindblad [*Commun. Math. Phys.* 48, 119 (1976)] zurück, und die Form im
Heisenberg-Bild, für normstetige Halbgruppen auf allgemeinen von-Neumann-Algebren, auf E. Christensen
und D. E. Evans [*J. London Math. Soc.* 20, 358 (1979)].

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Von einem System, das
> schwach mit einer großen Umgebung wechselwirkt, deren Korrelationen schnell abklingen (die
> Born-Markov-Näherung), erwartet man, dass es das Gedächtnis an seine Vergangenheit verliert: Sein
> Zustand zur Zeit $`s+t`$ geht dann aus dem Zustand zur Zeit $`s`$ durch einen Kanal hervor, der nur von
> $`t`$ abhängt. *Wenn* das so ist, besagt Satz 8.7, dass die Dynamik aus einem Hamiltonschen Anteil und
> „Sprüngen“ $`L_j`$ (Photonenemission, Teilchenverlust, Dephasierung) besteht, mit Raten, die durch die
> Normierung der $`L_j`$ festgelegt sind.

<a id="s-8-3"></a>

## 8.3 Lindblad-Gleichungen in einer Algebra mit Sektoren

Auf $`\mathcal M`$ gibt es zwei Arten von Sprüngen. Ein Sprung **innerhalb der Sektoren** ist ein
Operator $`L_j\in\mathcal M`$. Ein Sprung **zwischen Sektoren** ist eine CP-Abbildung
$`\Phi:\mathcal M\to\mathcal M`$ mit Kraus-Operatoren $`K^{dc}_l:\mathbb C^{k_c}\to\mathbb C^{k_d}`$
(Proposition 6.4, ein `InterSectorChannel`): $`\Phi(\rho)_d=\sum_{c,l}K^{dc}_l\rho_cK^{dc*}_l`$. Nach
Proposition 6.8 erfüllt ihre ungewichtete duale Abbildung $`\Phi^*(1)_c=\sum_{d,l}K^{dc*}_lK^{dc}_l`$,
ein blockdiagonales Element von $`\mathcal M`$. Mit Raten $`g_j,g_m>0`$ lautet der Generator von
`lindblad_rhs`
```math
\mathcal L(\rho)=-i[H,\rho]+\sum_jg_j\Big(L_j\rho L_j^*-\tfrac12\{L_j^*L_j,\rho\}\Big)
+\sum_mg_m\Big(\Phi_m(\rho)-\tfrac12\{\Phi_m^*(1),\rho\}\Big).
```
Nach der folgenden Proposition 8.8(1) erzeugt er eine quantendynamische Halbgruppe auf $`\mathcal M`$;
umgekehrt hat nach Satz 8.7(2) (Satz D.7) jede quantendynamische Halbgruppe auf $`\mathcal M`$ einen
Generator dieser Form.

*In der Bibliothek:* Ein Sprung zwischen Sektoren ist ein `InterSectorChannel`, und $`\mathcal L(\rho)`$ ist `dynamics.lindblad_rhs` — [Tutorial, Abschnitt 8.3](../../tutorial/de/08_dynamics.md#s-8-3).

<a id="r-8-8"></a>

**Proposition 8.8.** Sei $`J^{dc}_l\in M_D`$ die Matrix mit dem Block $`K^{dc}_l`$ an der Position
$`(d,c)`$ und Nullen sonst.

1. $`\mathcal L`$ bildet $`\mathcal M`$ in $`\mathcal M`$ ab und ist die Einschränkung auf $`\mathcal M`$
   des GKSL-Generators $`\hat{\mathcal L}`$ auf $`M_D`$ mit den Sprüngen $`\sqrt{g_j}L_j`$ und
   $`\sqrt{g_m}J^{dc}_{m,l}`$ (alle $`m,l,d,c`$). Folglich ist $`(e^{t\mathcal L})_{t\ge0}`$ eine
   quantendynamische Halbgruppe auf $`\mathcal M`$.
2. Sei stattdessen $`\tilde{\mathcal L}`$ der GKSL-Generator auf $`M_D`$, in dem die Blöcke eines
   Kraus-Index zu einem einzigen Sprung summiert werden, $`L_{m,l}=\sum_{d,c}J^{dc}_{m,l}`$ (mit
   demselben Hamiltonoperator, denselben Sprüngen $`\sqrt{g_j}L_j`$ und den Sprüngen
   $`\sqrt{g_m}L_{m,l}`$). Dann gilt
   $`\mathcal L(\rho)=E\big(\tilde{\mathcal L}(\rho)\big)`$ für $`\rho\in\mathcal M`$, wobei
   $`E(X)=\sum_dP_dXP_d`$ die bedingte Erwartung von $`M_D`$ auf $`\mathcal M`$ ist; im Allgemeinen ist
   $`\tilde{\mathcal L}(\rho)\notin\mathcal M`$.
3. Hat jedes $`L_{m,l}`$ in jeder Blockzeile und in jeder Blockspalte höchstens einen von null
   verschiedenen Block (etwa wenn es nur den Sektor $`c`$ in den Sektor $`c-1`$ abbildet, wie die
   Vernichtungsoperatoren im folgenden Beispiel), so gilt $`\tilde{\mathcal L}(\rho)=\mathcal L(\rho)`$
   für $`\rho\in\mathcal M`$.

*Beweis.* Da $`\rho=\sum_cP_c\rho P_c`$ und $`J^{dc}P_{c'}=\delta_{cc'}J^{dc}`$, ist
$`J^{dc}\rho J^{dc*}`$ der Block $`K^{dc}\rho_cK^{dc*}`$ an der Position $`(d,d)`$, und $`J^{dc*}J^{dc}`$
ist $`K^{dc*}K^{dc}`$ an der Position $`(c,c)`$. Summation über $`d,l`$ liefert $`\Phi(\rho)`$ und
$`\Phi^*(1)`$; das beweist die Identität in (1), und alle Terme liegen in $`\mathcal M`$. Also bilden
$`\hat{\mathcal L}`$ und alle seine Potenzen $`\mathcal M`$ in $`\mathcal M`$ ab, es gilt
$`e^{t\hat{\mathcal L}}|_{\mathcal M}=e^{t\mathcal L}`$, und dies ist nach Satz 8.6 und Lemma 8.4(a)
ein Kanal von $`\mathcal M`$. (2) Es ist $`L\rho L^*=\sum_{d,d',c}J^{dc}\rho J^{d'c*}`$ für
$`L=\sum_{d,c}J^{dc}`$, und $`E`$ behält die Terme mit $`d=d'`$. Für $`B\in\mathcal M`$ gilt
$`P_dXBP_d=P_dXP_dB`$, also $`E(XB)=E(X)B`$ und $`E(BX)=BE(X)`$; somit ist
$`E\{L^*L,\rho\}=\{E(L^*L),\rho\}`$, und $`E(L^*L)=\sum_cP_cL^*LP_c=\sum_{c,d}J^{dc*}J^{dc}`$, weil
$`J^{dc*}J^{d'c}=0`$ für $`d\ne d'`$. Schließlich ist $`E(-i[H,\rho])=-i[H,\rho]`$, und die Terme der
Sprünge $`\sqrt{g_j}L_j\in\mathcal M`$ liegen in $`\mathcal M`$. (3) Sei $`L=\sum_{(d,c)\in S}J^{dc}`$, wobei
die Menge $`S`$ von Paaren höchstens ein Paar mit gegebenem $`d`$ und höchstens eines mit gegebenem $`c`$
enthält. Für $`\rho\in\mathcal M`$ verschwindet $`J^{dc}\rho J^{d'c'*}=J^{dc}P_c\rho P_{c'}J^{d'c'*}`$,
außer wenn $`c=c'`$, und dann ist $`d=d'`$; und $`J^{dc*}J^{d'c'}`$ verschwindet, außer wenn $`d=d'`$, und
dann ist $`c=c'`$. Also sind $`L\rho L^*`$ und $`L^*L`$ blockdiagonal, alle Terme von
$`\tilde{\mathcal L}(\rho)`$ liegen in $`\mathcal M`$, und nach (2) ist
$`\tilde{\mathcal L}(\rho)=E(\tilde{\mathcal L}(\rho))=\mathcal L(\rho)`$. $`\square`$

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Teil (2) besagt,
> dass die Lindblad-Gleichung auf $`\mathcal M`$ das ist, was von einer Gleichung auf $`M_D`$ übrig
> bleibt, wenn die von den Sprüngen erzeugten Kohärenzen zwischen Sektoren sofort gelöscht werden –
> der Gehalt einer Superauswahlregel. Beim Teilchenverlust wird nichts gelöscht (Teil (3)): Die
> übliche Lindblad-Gleichung mit den Vernichtungsoperatoren als Sprüngen ist bereits eine Gleichung
> auf $`\mathcal M`$.

*Rechnung:* Proposition 8.8(2) für eine zufällige Abbildung zwischen den Sektoren und Proposition 8.8(3) für den Teilchenverlust, im Vergleich mit GKSL-Generatoren auf $`M_6`$ — [Tutorial, Abschnitt 8.3](../../tutorial/de/08_dynamics.md#s-8-3).

<a id="r-8-9"></a>

**Proposition 8.9 (Sektorgewichte).** Sei $`\rho(t)=e^{t\mathcal L}\rho_0`$, und seien $`p(t)`$ seine
Sektorgewichte.

1. Gibt es keine Sprünge zwischen Sektoren ($`K^{dc}_{m,l}=0`$ für alle $`m,l`$ und $`d\ne c`$), so ist
   jedes $`p_c(t)`$ konstant.
2. Gilt für jedes $`m`$ und jedes Paar $`d\ne c`$ die Beziehung
   $`\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}=r_{m,dc}1_{k_c}`$ mit Zahlen $`r_{m,dc}`$, so ist $`p'=Qp`$ mit
   $`Q_{dc}=\sum_mg_mr_{m,dc}`$ für $`d\ne c`$ und $`Q_{cc}=-\sum_{d\ne c}Q_{dc}`$. Folglich hängt
   $`p(t)=e^{tQ}p(0)`$ weder von $`H`$ noch von den $`L_j`$ noch von den Zuständen innerhalb der Sektoren
   ab.

*Beweis.* $`p_d'=\mathrm{Tr}(P_d\mathcal L(\rho))`$ nach Lemma 8.1. Der Hamiltonsche Term verschwindet wie in
Satz 8.3. Für $`L_j\in\mathcal M`$ ist
$`\mathrm{Tr}(P_dL_j\rho L_j^*)-\mathrm{Tr}(P_dL_j^*L_j\rho) =\mathrm{Tr}(L_j^*P_dL_j\rho)-\mathrm{Tr}(L_j^*L_jP_d\rho)=0`$,
da $`P_d`$ zentral ist. Für einen Sprung $`\Phi`$ zwischen Sektoren ist
$`\mathrm{Tr}(P_d\Phi(\rho))=\sum_c\mathrm{Tr}\big(\sum_lK^{dc*}_lK^{dc}_l\rho_c\big)`$
und
$`\mathrm{Tr}(P_d\tfrac12\{\Phi^*(1),\rho\})=\sum_{d'}\mathrm{Tr}\big(\sum_lK^{d'd*}_lK^{d'd}_l\rho_d\big)`$;
die Terme mit $`c=d`$ und $`d'=d`$ heben sich weg. Das beweist (1), und unter der Voraussetzung von (2)
sind die verbleibenden Terme $`\sum_{c\ne d}r_{dc}p_c-\sum_{d'\ne d}r_{d'd}p_d`$. Lemma 8.1 liefert
$`p(t)=e^{tQ}p(0)`$. $`\square`$

Unter der Voraussetzung von (2) ist für $`d\ne c`$ die Zahl $`r_{m,dc}`$ für jeden Zustand $`\rho`$ mit
$`p_c>0`$ der Eintrag $`T_{dc}(\rho)`$ der Übergangsmatrix von $`\Phi_m`$ (definiert vor Proposition 6.12),
da
$`\mathrm{Tr}\,\Phi_{m,dc}(\rho_c)=\mathrm{Tr}\big(\sum_lK^{dc*}_{m,l}K^{dc}_{m,l}\rho_c\big)=r_{m,dc}p_c`$.
Für einen
einzigen Sprung mit der Rate $`1`$ ist also $`Q=R-\mathrm{diag}(\sum_dR_{dc})_c`$ mit $`R`$ =
`transition_matrix()` (die Diagonale von $`R`$ hebt sich weg).

*In der Bibliothek:* $`R`$ ist `InterSectorChannel.transition_matrix` — [Tutorial, Abschnitt 8.3](../../tutorial/de/08_dynamics.md#s-8-3).

**Beispiel: Teilchenverlust.** Zwei bosonische Moden mit höchstens zwei Teilchen; die Sektoren sind
$`N=0,1,2`$ mit den Basen $`|00\rangle`$; $`|10\rangle,|01\rangle`$; $`|20\rangle,|11\rangle,|02\rangle`$,
also $`\mathcal M=M_1\oplus M_2\oplus M_3`$. Der Bose-Hubbard-Dimer
$`H=J(a_1^*a_2+a_2^*a_1)+\frac U2\sum_jn_j(n_j-1)`$ kommutiert mit $`N`$, und jede Mode verliert
Teilchen mit der Rate $`\gamma`$: ein Sprung zwischen Sektoren mit den Kraus-Operatoren $`a_1,a_2`$,
eingeschränkt auf $`N\to N-1`$. Da $`a_1^*a_1+a_2^*a_2=N`$, ist Proposition 8.9(2) mit
$`Q_{N-1,N}=\gamma N`$ anwendbar, und ausgehend von $`N=2`$ gilt
```math
p_2(t)=q^2,\qquad p_1(t)=2q(1-q),\qquad p_0(t)=(1-q)^2,\qquad q=e^{-\gamma t},
```
die Binomialverteilung mit den Parametern $`2`$ und $`q`$ (die Verteilung der Anzahl der Überlebenden,
wenn jedes Teilchen unabhängig mit der Wahrscheinlichkeit $`q`$ überlebt), was auch immer das Tunneln
(Hopping) und die Wechselwirkung bewirken: Beide Seiten lösen $`p'=Qp`$ (differenzieren) mit
$`p(0)=(0,0,1)`$, stimmen also nach Lemma 8.1 überein.

*Rechnung:* Das Beispiel, mit $`\Phi^*(1)=0\oplus1\oplus2`$ und $`p(t)`$ geprüft gegen beide Formeln — [Tutorial, Abschnitt 8.3](../../tutorial/de/08_dynamics.md#s-8-3).

<a id="s-8-4"></a>

## 8.4 Postselektion und die sprungfreie Entwicklung

Man fasse die Sprünge des Generators aus Abschnitt 8.3 zu
$`\Phi_{\mathrm{tot}}(\rho)=\sum_jg_jL_j\rho L_j^*+\sum_mg_m\Phi_m(\rho)`$ zusammen, mit
$`A=\Phi_{\mathrm{tot}}^*(1)=\sum_jg_jL_j^*L_j+\sum_mg_m\Phi_m^*(1)\ge0`$ in $`\mathcal M`$ und
$`H_{\mathrm{eff}}=H-\frac i2A`$. Man wende Satz 8.6 auf den Generator $`\hat{\mathcal L}`$ auf $`M_D`$ aus
Proposition 8.8(1) an: Nach dem Beweis von Proposition 8.8 bildet sein Sprunganteil $`\mathcal M`$ in
$`\mathcal M`$ ab und stimmt dort mit $`\Phi_{\mathrm{tot}}`$ überein, und seine Summe
$`\sum\hat L^*\hat L`$ über die Sprünge ist $`A`$, also ist sein $`H_{\mathrm{eff}}`$ gleich
$`H-\frac i2A\in\mathcal M`$. Folglich bilden seine $`S_n(t)`$ den Raum $`\mathcal M`$ in $`\mathcal M`$ ab
(durch Induktion über $`n`$), und für $`\rho\in\mathcal M`$ gilt
```math
\rho(t)=\sum_{n\ge0}S_n(t)\rho,\qquad S_0(t)\rho=e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*},
\qquad S_n(t)=\int_0^tS_0(t-s)\,\Phi_{\mathrm{tot}}\,S_{n-1}(s)\,ds .
```

<a id="r-8-10"></a>

**Proposition 8.10.** Für einen Zustand $`\rho`$ sei $`w_n(t)=\mathrm{Tr}\,S_n(t)\rho`$. Dann ist
$`w_n(t)\ge0`$, $`\sum_nw_n(t)=1`$, und $`w_0`$ ist monoton nicht wachsend mit
$`w_0'(t)=-\mathrm{Tr}\big(A\,S_0(t)\rho\big)`$. Für einen reinen Zustand
$`\rho=|\psi\rangle\langle\psi|`$ ist $`w_0(t)=\|e^{-itH_{\mathrm{eff}}}\psi\|^2`$.

*Beweis.* $`S_n(t)`$ ist CP, also positiv, und $`\sum_nw_n=\mathrm{Tr}\,\rho(t)=1`$. Mit
$`X=S_0(t)\rho\ge0`$ ist
$`w_0'=\mathrm{Tr}(-iH_{\mathrm{eff}}X+iXH_{\mathrm{eff}}^*)=-\mathrm{Tr}(AX)`$, und
$`\mathrm{Tr}(AX)=\mathrm{Tr}(X^{1/2}AX^{1/2})\ge0`$. Die letzte Behauptung ist
$`\mathrm{Tr}|\phi\rangle\langle\phi|=\|\phi\|^2`$. $`\square`$

> **Physikalische Bedeutung** (physikalische Deutung, kein Satz dieses Kurses). Ein Detektor registriere jeden Sprung – jedes emittierte Photon,
> jedes verlorene Atom. Dann ist $`w_n(t)`$ die Wahrscheinlichkeit für genau $`n`$ Klicks bis zur Zeit
> $`t`$, $`S_n(t)\rho/w_n(t)`$ der Zustand unter der Bedingung von $`n`$ Klicks und $`S_0(t)\rho/w_0(t)`$
> der Zustand unter der Bedingung, dass es **keinen Klick** gab: Auch ein Nullergebnis ändert den
> Zustand, und zwar über das nicht-hermitesche $`H_{\mathrm{eff}}`$. Diese Deutung als Zählprozess [E.
> B. Davies, *Quantum Theory of Open Systems*, 1976; M. D. Srinivas und E. B. Davies, *Optica Acta*
> 28, 981 (1981)] ist die Grundlage der Methode der Quantentrajektorien
> (Monte-Carlo-Wellenfunktionen) [J. Dalibard, Y. Castin, K. Mølmer, *Phys. Rev. Lett.* 68, 580
> (1992)]; hier verwenden wir nur ihren mathematischen Teil.

Bei Verlusten sind die Dyson-Terme die Sektorkomponenten von $`\rho(t)`$ (in der Sprache des Kastens:
Das Zählen der Klicks ist dasselbe wie das Messen der Ladung). Die **Lüders-Aktualisierung** eines
Zustands $`\rho`$ mit einer Projektion $`P`$, für die $`\mathrm{Tr}(P\rho)>0`$ ist, ist der Zustand
$`P\rho P/\mathrm{Tr}(P\rho P)`$ (man beachte $`\mathrm{Tr}(P\rho P)=\mathrm{Tr}(P\rho)`$).

*In der Bibliothek:* Die Lüders-Aktualisierung ist `states.lueders_update` — [Tutorial, Abschnitt 8.4](../../tutorial/de/08_dynamics.md#s-8-4).

<a id="r-8-11"></a>

**Proposition 8.11.** Die Sektoren seien mit $`N=0,\dots,C-1`$ bezeichnet, es gebe keine Sprünge
innerhalb der Sektoren, jeder Sprung senke die Ladung um eins ($`K^{dc}_{m,l}=0`$ außer für $`d=c-1`$),
und es sei $`\rho=P_{N_0}\rho P_{N_0}`$. Dann gilt $`S_n(t)\rho=P_{N_0-n}\,\rho(t)\,P_{N_0-n}`$ für
$`0\le n\le N_0`$ und $`S_n(t)\rho=0`$ für $`n>N_0`$. Insbesondere ist $`p_{N_0}(t)=w_0(t)>0`$, und die
Lüders-Aktualisierung von
$`\rho(t)`$ mit $`P_{N_0}`$ („kein Teilchen verloren“) der sprungfreie Zustand:
```math
\frac{P_{N_0}\rho(t)P_{N_0}}{p_{N_0}(t)}=\frac{e^{-itH_{\mathrm{eff}}}\rho\,e^{itH_{\mathrm{eff}}^*}}{w_0(t)},
\qquad p_{N_0}(t)=w_0(t).
```

*Beweis.* $`S_0(s)`$ bildet Operatoren mit Träger im Sektor $`N`$ auf Operatoren mit Träger im Sektor
$`N`$ ab, da $`H_{\mathrm{eff}}\in\mathcal M`$; $`\Phi_{\mathrm{tot}}`$ bildet den Sektor $`N\ge1`$ auf den
Sektor $`N-1`$ ab und annulliert den Sektor $`0`$. Durch Induktion hat $`S_n(s)\rho`$ für alle $`s`$ seinen
Träger im Sektor $`N_0-n`$ und ist $`0`$ für $`n>N_0`$. Die Kompression von $`\rho(t)=\sum_nS_n(t)\rho`$ mit
$`P_{N_0-n}`$ behält genau den $`n`$-ten Term. Schließlich ist $`w_0(t)>0`$, weil $`e^{-itH_{\mathrm{eff}}}`$
invertierbar ist (Lemma 8.1), sodass $`S_0(t)\rho`$ eine von null verschiedene positive Matrix ist.
$`\square`$

*Rechnung:* Proposition 8.11 für $`n=0`$ und $`n=1`$, für zwei Moden, die mit verschiedenen Raten Teilchen verlieren — [Tutorial, Abschnitt 8.4](../../tutorial/de/08_dynamics.md#s-8-4).

<a id="s-8-5"></a>

## 8.5 Stationäre Zustände und Relaxation

Ein Zustand $`\sigma`$ heißt **stationär**, wenn $`e^{t\mathcal L}\sigma=\sigma`$ für alle $`t\ge0`$ gilt,
äquivalent $`\mathcal L(\sigma)=0`$ (man differenziere bei $`t=0`$; umgekehrt verwende man die Reihe).

<a id="r-8-12"></a>

**Proposition 8.12.** Jede quantendynamische Halbgruppe auf $`\mathcal M`$ besitzt einen stationären
Zustand.

*Beweis.* Man fixiere einen Zustand $`\rho`$ und setze
$`\sigma_\tau=\frac1\tau\int_0^\tau T_s\rho\,ds`$. Dies ist nach Lemma A.4(d) ein Zustand (die Zustände
bilden eine abgeschlossene konvexe Menge). Nach Lemma A.4(b),(e) ist
$`T_t\sigma_\tau-\sigma_\tau=\frac1\tau\big(\int_\tau^{\tau+t}-\int_0^t\big)T_s\rho\,ds`$, und dessen
Spurnorm ist nach Lemma A.4(c) höchstens $`2t/\tau`$, da $`\|T_s\rho\|_1=1`$. Die Zustände bilden eine beschränkte
abgeschlossene Menge in einem endlichdimensionalen Raum, also konvergiert eine Folge
$`\sigma_{\tau_m}`$, $`\tau_m\to\infty`$, gegen einen Zustand $`\sigma`$ (Bolzano-Weierstraß), und
$`T_t\sigma-\sigma=\lim_m(T_t\sigma_{\tau_m}-\sigma_{\tau_m})=0`$. $`\square`$

<a id="r-8-13"></a>

**Korollar 8.13.** Gibt es keine Sprünge zwischen Sektoren (wie in Proposition 8.9(1)), so spannen
die stationären Zustände einen Raum der Dimension mindestens $`C`$ auf: Jeder Sektor trägt seinen
eigenen stationären Zustand.

*Beweis.* $`\mathcal L`$ bildet $`P_c\mathcal MP_c\cong M_{k_c}`$ in sich ab ($`H`$ und die $`L_j`$ liegen in
$`\mathcal M`$, und $`\Phi_m(X)=\sum_lK^{cc}_{m,l}X_cK^{cc*}_{m,l}`$ und $`\Phi_m^*(1)`$ sind
blockdiagonal), also auch jedes $`e^{t\mathcal L}`$. Seine Einschränkung ist CP (Lemma 8.4(a)) und TP,
also eine quantendynamische Halbgruppe auf $`M_{k_c}`$, mit einem stationären Zustand $`\sigma_c`$ nach
Proposition 8.12; als Element von $`\mathcal M`$ ist $`\sigma_c`$ ein stationärer Zustand von
$`e^{t\mathcal L}`$. Die $`\sigma_c`$ liegen in verschiedenen Blöcken und sind daher linear unabhängig.
$`\square`$

Satz D.6 (der endlichdimensionale Fall eines Kriteriums von A. Frigerio) gibt eine
hinreichende Bedingung für die Eindeutigkeit: Hat ein GKSL-Generator auf $`M_k`$ mit den Daten $`(H;L_1,\dots,L_r)`$ einen
treuen stationären Zustand $`\sigma>0`$ und sind die einzigen Matrizen, die mit $`H`$ und allen $`L_j`$,
$`L_j^*`$ kommutieren, die Vielfachen von $`1`$, so ist $`\sigma`$ der einzige stationäre Zustand. Ohne
Sprünge zwischen Sektoren ist die Einschränkung von $`\mathcal L`$ auf
$`P_c\mathcal MP_c\cong M_{k_c}`$ der GKSL-Generator mit den Daten
$`(H_c;\sqrt{g_1}L_{1,c},\sqrt{g_2}L_{2,c},\dots)`$, also ist der Satz sektorweise anwendbar. Aufgabe 4 zeigt,
dass die Eindeutigkeit selbst in einem Faktor verletzt
sein kann: Dort ist die zweite Voraussetzung verletzt.

*Rechnung:* Korollar 8.13 und die Voraussetzungen von Satz D.6 für zufällige Sprünge innerhalb der Sektoren von $`M_2\oplus M_3`$ — [Tutorial, Abschnitt 8.5](../../tutorial/de/08_dynamics.md#s-8-5).

<a id="r-8-14"></a>

**Proposition 8.14 (spontane Emission).** Auf $`M_2`$ mit der Basis $`|g\rangle,|e\rangle`$ seien
$`H=\omega|e\rangle\langle e|`$, ein Sprung $`L=|g\rangle\langle e|`$ und die Rate $`\gamma>0`$ gegeben.
Dann gilt
```math
\rho_{ee}(t)=e^{-\gamma t}\rho_{ee}(0),\qquad\rho_{ge}(t)=e^{(i\omega-\gamma/2)t}\rho_{ge}(0),
\qquad e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}},
```
wobei $`\mathcal U_t(\rho)=e^{-itH}\rho e^{itH}`$ und $`\mathcal A_p`$ der Amplitudendämpfungskanal mit
Parameter $`p`$ ist (Abschnitt 6.4). Das Spektrum von $`\mathcal L`$ ist
$`\{0,-\gamma,-\gamma/2\pm i\omega\}`$, $`|g\rangle\langle g|`$ ist der einzige stationäre Zustand, und
es gilt
$`\tfrac12\|\rho(t)-|g\rangle\langle g|\|_1=\big(\rho_{ee}(0)^2e^{-2\gamma t}+|\rho_{ge}(0)|^2e^{-\gamma t}\big)^{1/2}`$.

*Beweis.* Hier ist $`L^*L=|e\rangle\langle e|`$ und $`L\rho L^*=\rho_{ee}|g\rangle\langle g|`$, also
eintragsweise $`\rho_{ee}'=-\gamma\rho_{ee}`$, $`\rho_{gg}'=\gamma\rho_{ee}`$ und
$`\rho_{ge}'=(i\omega-\gamma/2)\rho_{ge}`$, woraus die Lösung folgt. In der Basis
$`(\rho_{gg},\rho_{ge},\rho_{eg},\rho_{ee})`$ ist die Matrix von $`\mathcal L`$ eine obere
Dreiecksmatrix mit der Diagonale $`(0,-\gamma/2+i\omega,-\gamma/2-i\omega,-\gamma)`$.
$`\mathcal L(\sigma)=0`$ erzwingt $`\sigma_{ee}=\sigma_{ge}=0`$. Der Amplitudendämpfungskanal
multipliziert $`\rho_{ee}`$ mit $`1-p`$ und $`\rho_{ge}`$ mit $`\sqrt{1-p}`$ (seine Kraus-Operatoren sind
$`|g\rangle\langle g|+\sqrt{1-p}\,|e\rangle\langle e|`$ und $`\sqrt p\,|g\rangle\langle e|`$), und
$`\mathcal U_t`$ multipliziert $`\rho_{ge}`$ mit $`e^{i\omega t}`$; beide sind TP, was $`\rho_{gg}`$
festlegt. Schließlich ist $`\rho(t)-|g\rangle\langle g|`$ spurfrei und selbstadjungiert mit der
Diagonale $`(-\rho_{ee}(t),\rho_{ee}(t))`$, seine Eigenwerte sind also
$`\pm(\rho_{ee}(t)^2+|\rho_{ge}(t)|^2)^{1/2}`$. $`\square`$

Die Kohärenz zerfällt mit der halben Rate der Besetzung: Mit $`T_1=1/\gamma`$ und $`T_2=2/\gamma`$ ist
das die Beziehung $`\frac1{T_2}=\frac1{2T_1}+\frac1{T_\varphi}`$ aus Abschnitt 6.4 ohne den
Dephasierungsterm. Die langsamste von null verschiedene Zerfallsrate, die **Spektrallücke**
$`\min\{-\mathrm{Re}\,\mu:\mu\ne0\text{ ein Eigenwert von }\mathcal L\}=\gamma/2`$, ist hier die
Relaxationsrate: Nach der letzten Formel ist der Abstand zu $`|g\rangle\langle g|`$ höchstens
$`e^{-\gamma t/2}`$ (da $`|\rho_{ge}|^2\le\rho_{gg}\rho_{ee}`$ wegen der Positivität,
$`\rho_{ee}^2+|\rho_{ge}|^2\le\rho_{ee}\le1`$), und $`e^{\gamma t/2}`$ mal dieser Abstand strebt für
$`t\to\infty`$ gegen $`|\rho_{ge}(0)|`$.

*Rechnung:* Das Spektrum von $`\mathcal L`$, die Faktorisierung $`e^{t\mathcal L}=\mathcal U_t\circ\mathcal A_{1-e^{-\gamma t}}`$ und der Spurabstand aus Proposition 8.14, mit `dynamics.lindblad_superoperator` und `dynamics.lindblad_channel` — [Tutorial, Abschnitt 8.5](../../tutorial/de/08_dynamics.md#s-8-5).

<a id="s-8-6"></a>

## 8.6 Womit man vergleichen sollte

*Rechnung:* Der Vergleich mit anderer Software (QuTiP, dynamiqs, QuantumOptics.jl, SciPy) — [Tutorial, Abschnitt 8.6](../../tutorial/de/08_dynamics.md#s-8-6).

<a id="x-8"></a>

## Aufgaben

1. (Heisenberg-Bild) Man zeige
   $`\mathrm{Tr}(\mathcal L(\rho)A)=\mathrm{Tr}(\rho\,\mathcal L^*(A))`$ für den Generator
   aus Abschnitt 8.3, mit
   $`\mathcal L^*(A)=i[H,A]+\sum_jg_j(L_j^*AL_j-\frac12\{L_j^*L_j,A\})+\sum_mg_m(\Phi_m^*(A)-\frac12\{\Phi_m^*(1),A\})`$,
   ferner $`\mathcal L^*(1)=0`$ sowie $`\mathcal L^*(P_c)=0`$, falls es keine Sprünge zwischen Sektoren
   gibt.
2. (Exzeptioneller Punkt) Im Sektor $`N=1`$ aus Abschnitt 8.4 ist
   $`H_{\mathrm{eff}}=\begin{pmatrix}-i\gamma_1/2&J\\J&-i\gamma_2/2\end{pmatrix}`$. Man bestimme seine
   Eigenwerte und zeige, dass er für $`J=|\gamma_1-\gamma_2|/4>0`$ nicht diagonalisierbar ist.
3. $`N_0`$ Bosonen in $`m`$ Moden mögen jeweils mit der Rate $`\gamma`$ verloren gehen (Sprünge
   $`a_1,\dots,a_m`$, eingeschränkt auf $`N\to N-1`$). Man zeige $`p_N(t)=\binom{N_0}Nq^N(1-q)^{N_0-N}`$
   mit $`q=e^{-\gamma t}`$ für jeden Hamiltonoperator in $`\mathcal M`$.
4. (Reine Dephasierung) Auf $`M_2`$ seien $`H=\omega|e\rangle\langle e|`$ und
   $`L=\sigma_z=\mathrm{diag}(1,-1)`$ mit der Rate $`\gamma`$. Man zeige, dass die Besetzungen
   konstant sind, dass $`\rho_{ge}(t)=e^{(i\omega-2\gamma)t}\rho_{ge}(0)`$ gilt und dass die
   stationären Zustände den zweidimensionalen Raum der Diagonalmatrizen aufspannen, obwohl $`C=1`$
   ist.

<details>
<summary>Lösungen</summary>

1. Wegen der Zyklizität gilt $`\mathrm{Tr}(L\rho L^*A)=\mathrm{Tr}(\rho L^*AL)`$,
   $`\mathrm{Tr}(\{L^*L,\rho\}A)=\mathrm{Tr}(\rho\{L^*L,A\})`$ und
   $`\mathrm{Tr}(-i[H,\rho]A)=\mathrm{Tr}(\rho\,i[H,A])`$; für Sprünge zwischen Sektoren
   ist $`\mathrm{Tr}(\Phi(\rho)A)=\mathrm{Tr}(\rho\Phi^*(A))`$ die Definition von $`\Phi^*`$
   (Abschnitt 6.3). Dann ist
   $`\mathcal L^*(1)=\sum_jg_j(L_j^*L_j-L_j^*L_j)+\sum_mg_m(\Phi_m^*(1)-\Phi_m^*(1))=0`$, was noch
   einmal die Spurerhaltung ist. Ohne Sprünge zwischen Sektoren kommutiert $`P_c`$ mit $`H`$ und allen
   $`L_j`$, also ist $`L_j^*P_cL_j=L_j^*L_jP_c=\frac12\{L_j^*L_j,P_c\}`$ und $`\mathcal L^*(P_c)=0`$; dann
   ist $`p_c'=\mathrm{Tr}(\rho\mathcal L^*(P_c))=0`$, ein zweiter Beweis von Proposition 8.9(1).

2. $`H_{\mathrm{eff}}=-\frac i4(\gamma_1+\gamma_2)1+\begin{pmatrix}-i\delta&J\\J&i\delta\end{pmatrix}`$
   mit $`\delta=(\gamma_1-\gamma_2)/4`$; die zweite Matrix hat die Spur $`0`$ und die Determinante
   $`\delta^2-J^2`$, die Eigenwerte sind also $`-\frac i4(\gamma_1+\gamma_2)\pm\sqrt{J^2-\delta^2}`$:
   zwei Zerfallsraten für $`J<|\delta|`$, eine Rate und zwei Frequenzen für $`J>|\delta|`$. Für
   $`J=|\delta|>0`$ erfüllt die zweite Matrix $`Y`$ die Gleichung $`Y^2=(J^2-\delta^2)1=0`$, aber $`Y\ne0`$;
   daher hat $`H_{\mathrm{eff}}`$ den einzigen Eigenwert $`-\frac i4(\gamma_1+\gamma_2)`$ und ist kein
   Vielfaches von $`1`$: Er ist nicht diagonalisierbar.

3. Auf dem Sektor $`N`$ gilt $`\sum_ja_j^*a_j=N1`$, und das einzige Sektorpaar ist $`(N-1,N)`$, also
   liefert Proposition 8.9(2) $`p_N'=-\gamma Np_N+\gamma(N+1)p_{N+1}`$. Mit $`q'=-\gamma q`$ ist die
   Ableitung des Binomialterms gleich
   $`-\gamma Np_N+\gamma(N_0-N)\binom{N_0}Nq^{N+1}(1-q)^{N_0-N-1}`$, und
   $`(N+1)\binom{N_0}{N+1}=(N_0-N)\binom{N_0}N`$ zeigt, dass der zweite Term gleich
   $`\gamma(N+1)p_{N+1}`$ ist. Bei $`t=0`$ ist dies die Punktmasse in $`N_0`$; die Eindeutigkeit ist Lemma
   8.1.

4. $`L\rho L^*-\frac12\{L^*L,\rho\}=\sigma_z\rho\sigma_z-\rho`$ multipliziert die
   Nichtdiagonaleinträge mit $`-2`$ und annulliert die Diagonale, und
   $`-i[H,\rho]_{ge}=i\omega\rho_{ge}`$. Also gilt $`\mathcal L(\sigma)=0`$ genau dann, wenn
   $`\sigma_{ge}=0`$ (für $`\gamma>0`$), und das Spektrum ist $`\{0,0,-2\gamma\pm i\omega\}`$. Mit
   $`U_t=e^{-itH}`$ gilt
   $`e^{t\mathcal L}(\rho)=\frac{1+e^{-2\gamma t}}2U_t\rho U_t^*+\frac{1-e^{-2\gamma t}}2U_t\sigma_z\rho\sigma_zU_t^*`$:
   Da $`\sigma_z`$ mit $`U_t`$ kommutiert und $`\sigma_z^2=1`$ ist, zeigt Differenzieren, dass die rechte
   Seite $`R(t)`$ die Gleichungen $`R'=-i[H,R]+\gamma(\sigma_zR\sigma_z-R)=\mathcal L(R)`$ und $`R(0)=\rho`$
   erfüllt, also ist Lemma 8.1 anwendbar. Für $`t>0`$ ist dies ein unitaler Kanal, dessen Fixpunkte die
   mit seinen Kraus-Operatoren kommutierenden Matrizen sind (Satz 6.10), also die Diagonalmatrizen.

</details>

*Rechnung:* Aufgaben 2 und 4 numerisch — [Tutorial, Aufgaben](../../tutorial/de/08_dynamics.md#x-8).

**Weiter:** [Kapitel 9](09_spin_chains.md) behandelt Spinketten: Sektoren von Magnetisierung und
Impuls, Lanczos- und Krylov-Zeitentwicklung für große Sektoren sowie Unordnung.
