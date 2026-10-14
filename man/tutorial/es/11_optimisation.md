[English](../11_optimisation.md) · **Español** · [Français](../fr/11_optimisation.md) · [Deutsch](../de/11_optimisation.md) · [中文](../zh/11_optimisation.md) · [日本語](../ja/11_optimisation.md) · [Русский](../ru/11_optimisation.md)

# Tutorial 11. Optimización y desigualdades de trazas

Este tutorial acompaña al [capítulo 11](../../lecture_notes/es/11_optimisation.md) de los apuntes de clase, que demuestra los extremos de tres
funcionales de traza y explica por qué el muestreo de unitarias no los encuentra. Aquí las parametrizaciones de
conjuntos de restricciones de la biblioteca (`optimize.PositiveParam`, `UnitaryParam`, `SelfAdjointParam`) y su
optimizador multiarranque (`optimize.extremize`) se contrastan con esos teoremas. Con $`\Delta`$ el contraste de Michelson
y $`\ell(\delta)=(1-\delta)/(1+\delta)`$, los tres funcionales son:

| | funcional | variables |
|---|---|---|
| Exp. 1 | $`z=\lvert\mathrm{Tr}(XUY)\rvert-\mathrm{Tr}(XY)`$ | $`X,Y\ge0`$, $`U\in\mathcal M`$ unitaria |
| Exp. 2 | $`z=\mathrm{Tr}\lvert XY\rvert-\mathrm{Tr}(\lvert X\rvert\,\lvert Y\rvert)`$ | $`X=UX_0`$, $`Y=VY_0`$, $`X_0,Y_0\ge0`$, $`U,V`$ unitarias |
| Exp. 3 | $`z=\mathrm{Tr}(Y\lvert X\rvert Y)-\mathrm{Tr}\lvert YXY\rvert`$ | $`X=X^*`$, $`Y\ge0`$ |

con $`|X|`$, $`|Y|`$ normalizados a $`\lambda_{\max}=1`$ y sus contrastes como parámetros. Las **envolventes** son
$`\overline z(\delta)=\sup\{z:\Delta(|X|)=\delta\}`$ y $`\underline z(\delta)=\inf\{z:\Delta(|X|)=\delta\}`$, con todas las
demás variables libres.

<a id="s-11-1"></a>

## 11.1 Conjuntos de restricciones y sus parametrizaciones

*Apuntes:* [sección 11.1](../../lecture_notes/es/11_optimisation.md#s-11-1); [lema 11.1](../../lecture_notes/es/11_optimisation.md#r-11-1), [proposición 11.2](../../lecture_notes/es/11_optimisation.md#r-11-2).

`optimize.PositiveParam` con contraste prescrito $`\delta`$ implementa la parametrización de la
[sección 11.1](../../lecture_notes/es/11_optimisation.md#s-11-1) de los apuntes de clase: en cada arranque sortea dos casillas distintas $`s_\pm`$ (`spectrum.pin_max`,
`spectrum.pin_min`; todo par ordenado tiene probabilidad $`1/(D(D-1))`$) y una unitaria aleatoria de Haar $`U_0`$
(`unitary.base`, sorteada por `random_unitary_operator`, [lema 11.10](../../lecture_notes/es/11_optimisation.md#r-11-10)), y envía los parámetros $`\theta`$
(`spectrum.theta`) y $`G`$ (`unitary.G`) a
```math
X(\theta,G)=U\mathrm{diag}(\lambda)U^*,\qquad U=U_0\,e^{G-G^*},\qquad
\lambda_{s_+}=1,\quad\lambda_{s_-}=\ell(\delta),\quad\lambda_s=\ell+(1-\ell)\tfrac{1+\sin\theta_s}2\ \text{otherwise}.
```
`UnitaryParam` es la aplicación $`G\mapsto U_0e^{G-G^*}`$ sola, y `SelfAdjointParam` añade signos,
$`X=U\mathrm{diag}(\sigma\lambda)U^*`$. Por la [proposición 11.2](../../lecture_notes/es/11_optimisation.md#r-11-2), la imagen es el conjunto de los operadores de
la capa de contraste $`\mathcal P_\delta`$ con $`1`$ y $`\ell(\delta)`$ en las casillas sorteadas. Un contraste libre en $`[a,b]`$ es
un parámetro más; el valor por defecto de la biblioteca es $`[a,b]=[0,0.999]`$, y el código de abajo lo mantiene.

El código comprueba la restricción, y la proposición 11.2 de forma constructiva: los parámetros de la demostración
reproducen un objetivo.

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

<a id="s-11-2"></a>

## 11.2 La desigualdad de trazas de von Neumann

*Apuntes:* [sección 11.2](../../lecture_notes/es/11_optimisation.md#s-11-2); [teorema 11.4](../../lecture_notes/es/11_optimisation.md#r-11-4), [corolario 11.5](../../lecture_notes/es/11_optimisation.md#r-11-5), [proposición 11.6](../../lecture_notes/es/11_optimisation.md#r-11-6), [corolario 11.7](../../lecture_notes/es/11_optimisation.md#r-11-7).

Según [apuntes de clase, teorema 11.4](../../lecture_notes/es/11_optimisation.md#r-11-4) (desigualdad de trazas de von Neumann),
$`|\mathrm{Tr}(AB)|\le\sum_is_i(A)s_i(B)`$, y por el [corolario 11.5](../../lecture_notes/es/11_optimisation.md#r-11-5), $`\sum_ip^\downarrow_iq^\uparrow_i\le\mathrm{Tr}(PQ)\le\sum_ip^\downarrow_iq^\downarrow_i`$ para $`P,Q`$ autoadjuntas con valores
propios $`p,q`$, y $`\mathrm{Tr}(PQ)\le\|PQ\|_1\le\sum_ip^\downarrow_iq^\downarrow_i`$ para $`P,Q\ge0`$. Para Exp. 2 con espectros fijos, la
[proposición 11.6](../../lecture_notes/es/11_optimisation.md#r-11-6) da $`\sup z=-\inf z=\sum_c\sum_ix^\downarrow_{c,i}(y^\downarrow_{c,i}-y^\uparrow_{c,i})`$, que se alcanza para $`X_0=\mathrm{diag}(x^\downarrow)`$, $`V=J`$ la permutación de inversión e $`Y_0=\mathrm{diag}(y^\uparrow)`$ (resp.
$`\mathrm{diag}(y^\downarrow)`$), y el [corolario 11.7](../../lecture_notes/es/11_optimisation.md#r-11-7) da la envolvente $`\sup z=-\inf z=\sum_c\lfloor k_c/2\rfloor(1-\ell_X)(1-\ell_Y)`$.

Aquí $`1-\ell_X=2\delta/(1+\delta)`$ para $`\Delta(|X|)=\delta`$; para $`k_c`$ impar, el valor propio central es irrelevante.
El código comprueba el teorema 11.4, el corolario 11.5 y las configuraciones extremales de la proposición 11.6 para
$`k=3`$, $`C=2`$.

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

<a id="s-11-3"></a>

## 11.3 Una desigualdad de trazas con condición de igualdad

*Apuntes:* [sección 11.3](../../lecture_notes/es/11_optimisation.md#s-11-3); [lema 11.8](../../lecture_notes/es/11_optimisation.md#r-11-8), [proposición 11.9](../../lecture_notes/es/11_optimisation.md#r-11-9).

Según [apuntes de clase, proposición 11.9](../../lecture_notes/es/11_optimisation.md#r-11-9), $`z=\mathrm{Tr}(Y|X|Y)-\mathrm{Tr}|YXY|\ge0`$ en Exp. 3, con igualdad si y solo si
$`YX_+Y^2X_-Y=0`$, en particular si $`XY=YX`$. Para $`X=\mathrm{diag}(1,-1,0)`$ e $`Y`$ con filas $`(1,0,a)`$, $`(0,1,0)`$, $`(a,0,1)`$, $`-1<a<1`$, la igualdad se cumple
aunque $`XY\ne YX`$ para $`a\ne0`$. El código comprueba ambas cosas (con $`a=0.5`$); `X.abs()` es $`|X|`$ y `trace_norm()` es
$`\|\cdot\|_1`$.

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

<a id="s-11-4"></a>

## 11.4 Muestreo frente a optimización

*Apuntes:* [sección 11.4](../../lecture_notes/es/11_optimisation.md#s-11-4); [lema 11.10](../../lecture_notes/es/11_optimisation.md#r-11-10), [proposición 11.11](../../lecture_notes/es/11_optimisation.md#r-11-11), [corolario 11.12](../../lecture_notes/es/11_optimisation.md#r-11-12), [lema 11.13](../../lecture_notes/es/11_optimisation.md#r-11-13).

¿Qué tamaño tiene típicamente $`|\mathrm{Tr}(AU)|`$, comparado con su máximo $`\|A\|_1`$? La biblioteca muestrea $`U`$
como el factor $`Q`$ de la descomposición QR $`Z=QR`$, $`R_{ii}>0`$, de una **matriz de Ginibre** $`Z`$ (entradas
independientes con densidad $`\pi^{-1}e^{-|z|^2}`$; la receta de Mezzadri), de forma independiente en cada bloque.

Según [apuntes de clase, proposición 11.11](../../lecture_notes/es/11_optimisation.md#r-11-11), $`\mathbb E|\mathrm{Tr}(AU)|^2=\|A\|_2^2/k`$ para $`U\in U(k)`$ de Haar, y por el
[corolario 11.12](../../lecture_notes/es/11_optimisation.md#r-11-12), para $`X,Y\ge0`$ con espectros en $`[a,b]`$, $`a>0`$, la media cuadrática de $`|\mathrm{Tr}(XUY)|`$ es a lo
sumo $`\min(1/\sqrt k,b^2/(a^2k))`$ veces su máximo sobre $`U`$, y $`\mathbb P\big(|\mathrm{Tr}(XUY)|\ge\mathrm{Tr}(XY)\big)\le b^4/(a^4k^2)`$.

Así, $`N`$ muestras de Haar encuentran una violación de $`|\mathrm{Tr}(XUY)|\le\mathrm{Tr}(XY)`$ con probabilidad
a lo sumo $`Nb^4/(a^4k^2)`$ (una unión de $`N`$ sucesos), aunque existe una para todo par que no conmuta
(proposición 3.9). La cota solo es informativa para $`k>\sqrt N\,b^2/a^2`$ (para $`b/a=6`$ y $`N=2000`$: $`k>1609`$);
numéricamente, las muestras pasan por alto la violación ya para $`k`$ mucho menores (véase abajo). El código comprueba
la proposición 11.11 por Monte Carlo (la media muestral está a menos de cinco errores estándar del valor exacto, en
esta ejecución con semilla fija) e imprime, para espectros uniformes en $`[0.2,1.2]`$ ($`b/a=6`$), el cociente del
corolario 11.12, $`\sup z`$ y el mejor $`z`$ de $`2000`$ muestras, ambos relativos a $`\mathrm{Tr}(XY)`$.

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

*Observación numérica* (la tabla impresa; afirmada solo para esta ejecución con semilla fija): el cociente se
mantiene por debajo de la cota demostrada $`\min(1/\sqrt k,b^2/(a^2k))`$ y se reduce aproximadamente a la mitad cada
vez que $`k`$ se duplica, es decir, decrece como $`1/k`$, y a partir de $`k=4`$ la mejor de $`2000`$ muestras ni siquiera
alcanza $`z=0`$. Una búsqueda por gradiente, en cambio, se acerca a ambos extremos en la ejecución de abajo (de nuevo,
una observación numérica). Para el par `Xf, Yf` se conocen ambos extremos sobre $`U`$:

Según [apuntes de clase, lema 11.13](../../lecture_notes/es/11_optimisation.md#r-11-13), $`\min_U|\mathrm{Tr}(AU)|=0`$ para toda $`A\in M_k(\mathbb C)`$, $`k\ge2`$.

Así, $`\inf_Uz=-\mathrm{Tr}(XY)`$ y $`\sup_Uz=\|YX\|_1-\mathrm{Tr}(XY)`$ (proposición 3.9). Una sola llamada a
`extremize` aproxima ambos, con $`12`$ arranques que maximizan $`z`$ y $`12`$ que maximizan $`-z`$ (con una precisión relativa
de $`10^{-3}`$ en esta ejecución, afirmada abajo). Por la proposición 11.11, la media cuadrática de
$`|\mathrm{Tr}(XUY)|`$ es $`\|YX\|_2/\sqrt k`$; con las cotas de la demostración del corolario 11.12 (espectros en
$`[a,b]`$), $`\mathbb E|\mathrm{Tr}(XUY)|^2\le b^4`$ y $`\mathrm{Tr}(XY)\ge ka^2`$, de modo que la desigualdad de
Markov da $`\mathbb P\big(z\ge-(1-\varepsilon)\mathrm{Tr}(XY)\big)\le b^4/(\varepsilon^2a^4k^2)`$ para
$`\varepsilon>0`$: para $`k`$ grande, $`z`$ está cerca de su ínfimo con probabilidad alta.
*Observación numérica* (afirmada para esta ejecución con semilla fija): el muestreo llega a menos de un $`1\,\%`$ del
ínfimo, pero ni siquiera alcanza $`z=0`$, mientras que el supremo es positivo:

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

<a id="s-11-5"></a>

## 11.5 Búsqueda multiarranque de envolventes

*Apuntes:* [sección 11.5](../../lecture_notes/es/11_optimisation.md#s-11-5); [proposición 11.14](../../lecture_notes/es/11_optimisation.md#r-11-14), [corolario 11.7](../../lecture_notes/es/11_optimisation.md#r-11-7), [proposición 11.9](../../lecture_notes/es/11_optimisation.md#r-11-9).

`optimize.extremize(objective, params, maximize)` ejecuta Adam en `rounds` rondas de `steps` pasos (tasa de
aprendizaje `lr * decay**round`) sobre la pérdida $`\mp\sum_b\text{objective}_b`$, recentra las cartas unitarias tras
cada ronda y devuelve, **para cada arranque**, el mayor (el menor, para `maximize=False`) valor del objetivo entre los
calculados en los iterados del gradiente y una vez en el punto final. Cada elemento del lote es un arranque
independiente con sus propias casillas, signos, $`U_0`$ y contraste: Adam actúa coordenada a coordenada, así que, si la
$`b`$-ésima entrada del objetivo depende solo de los $`b`$-ésimos parámetros (como en todos los objetivos de abajo), los
arranques no interactúan. Los signos de `SelfAdjointParam` son discretos: entre rondas, `discrete_step` prueba ocho
cambios de signo individuales por arranque (cada uno en una casilla uniformemente aleatoria), uno tras otro, y
conserva cada cambio que mejora el objetivo.

En aritmética exacta todo iterado es admisible (la proposición 11.2 y las observaciones que la siguen), así que todo
valor devuelto es el valor de $`z`$ en un punto admisible; en coma flotante es ese valor salvo errores de redondeo (las
aserciones de abajo admiten $`10^{-12}`$). Por tanto, una maximización da una **cota inferior** del supremo, una
minimización una cota superior del ínfimo, y nada más: una búsqueda local puede detenerse en un extremo local. La
otra dirección requiere una demostración, como el corolario 11.7 y la proposición 11.9, o una relajación convexa
(sección 11.6). Una envolvente calculada se contrasta (i) con formas cerradas; (ii) con configuraciones construidas a
mano, como sumas directas de soluciones más pequeñas, por debajo de las cuales no debe quedar; (iii) con muestras
estructuradas (espectros de dos niveles, bases de permutación), cuyos valores son otros valores alcanzados.

Según [apuntes de clase, proposición 11.14](../../lecture_notes/es/11_optimisation.md#r-11-14), si todos los $`k_c\ge2`$, la envolvente inferior de Exp. 1 es $`\underline z(\delta)=-\big(D-2\delta/(1+\delta)\big)`$.

El código comprueba la proposición 11.14 en $`M_2\oplus M_2`$, el corolario 11.7 en $`M_3\oplus M_3`$ (ambos contrastes
prescritos, $`k=3`$ impar) y la proposición 11.9; supremos e ínfimos comparten una llamada, como antes. Las aserciones
unilaterales (ningún valor devuelto más allá de un extremo demostrado) son consecuencias de estos teoremas; las
aserciones de que la búsqueda llega a la tolerancia indicada del extremo son observaciones numéricas de esta
ejecución con semilla fija. La palabra clave `is_self_adjoint` selecciona el algoritmo para la norma traza (valores
singulares, o valores propios para $`YXY`$) sin una comprobación en cada paso. La última comprobación es del tipo (ii):
el $`z`$ de Exp. 3 es aditivo respecto de las sumas directas, y el contraste de una suma directa es el mayor contraste
de sus sumandos (ejercicio 3), luego, para las envolventes de Exp. 3 con $`\Delta(Y)`$ libre en $`[0,0.999]`$ como en el
código,
```math
\overline z_{M_2\oplus M_2}(\delta)\ \ge\ \overline z_{M_2}(\delta)+\sup_{\delta'\le\delta}\overline z_{M_2}(\delta');
```
en particular, la suma de dos valores alcanzados en $`M_2`$ con contrastes $`\delta`$ y $`\delta'\le\delta`$ es un valor
alcanzado en $`M_2\oplus M_2`$ con contraste $`\delta`$.

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

*Observaciones numéricas* (los valores impresos de esta ejecución). La búsqueda en $`M_2\oplus M_2`$ se queda
alrededor de un $`0.2\,\%`$ por debajo de la suma directa: no ha convergido del todo, y la mejor cota inferior de
$`\overline z_{M_2\oplus M_2}(0.5)`$ es el número mayor; detectar esto es el propósito de la comprobación (ii). Los
valores para $`M_2`$ coinciden hasta $`10^{-4}`$ con $`(1-\delta)\big(2d/(1+d)\big)^2=(1-\delta)(1-\ell(d))^2`$, donde
$`d=0.999`$ es el mayor contraste admitido de $`Y`$. Esta forma cerrada **no** es un teorema de este curso: el artículo
que acompaña a la biblioteca la observó numéricamente (en álgebras reales) y no la demostró, y el código la imprime
sin afirmar nada sobre ella. Para bloques mayores no conocemos ninguna forma cerrada de las envolventes superiores
de Exp. 1 y Exp. 3; las curvas calculadas son entonces cotas inferiores (valores alcanzados), nada más.

<a id="s-11-6"></a>

## 11.6 Con qué comparar

**Manopt**, **pymanopt** y **geoopt** proporcionan variedades matriciales (variedades de Stiefel, grupos
ortogonales y unitarios) con retracciones como la QR, la polar o la de Cayley, y resolvedores riemannianos
(gradientes conjugados, regiones de confianza, Adam riemanniano), una alternativa a Adam en una carta recentrada (no
hacemos aquí ninguna comparación de rendimiento). **cvxpy** pasa problemas convexos en variables operador a
resolvedores SDP, que devuelven un óptimo junto con un certificado dual, ambos salvo la tolerancia del resolvedor;
una relajación —un problema convexo sobre un conjunto factible mayor— acota el supremo por arriba, la dirección que
una búsqueda local no puede dar. **NCTSSOS**, **ncpol2sdpa** y
**NCSOStools** acotan valores propios y trazas de polinomios no conmutativos en todas las dimensiones a la vez
(sección 0.3). La biblioteca añade parametrizaciones exactas de las restricciones espectrales y de contraste de
$`\bigoplus_cM_{k_c}`$ (los extremos pueden estar en cualquier sector), unitarias que respetan los bloques, y un único
cálculo por lotes para muchos arranques y valores de la restricción, en CPU o en GPU. Produce candidatos a
extremizadores y cotas inferiores de supremos; el resto deben aportarlo demostraciones o relajaciones.

<a id="x-11"></a>

## Ejercicios

[Ejercicios 1 y 2](../../lecture_notes/es/11_optimisation.md#x-11) de los apuntes de clase, numéricamente: la transformada de Cayley $`C(S)=(1-S)(1+S)^{-1}`$ de una
$`S`$ antiautoadjunta es unitaria, no tiene valor propio $`-1`$, y $`S=(1-C)(1+C)^{-1}`$ (ejercicio 1); para $`O`$ distribuida según la
medida de Haar en $`O(k)`$, muestreada con `random_unitary_operator` de un álgebra real (`complex_valued=False`),
$`\mathbb E\mathrm{Tr}(AO)^2=\|A\|_2^2/k`$ (ejercicio 2).

<details>
<summary>Soluciones</summary>

```python
Sk = torch.randn(4, 4, dtype=torch.complex128)
Sk = Sk - Sk.mH                                                   # skew-adjoint
I4 = torch.eye(4, dtype=torch.complex128)
Cay = (I4 - Sk) @ torch.linalg.inv(I4 + Sk)
assert torch.allclose(Cay.mH @ Cay, I4, atol=1e-12)
assert torch.allclose((I4 - Cay) @ torch.linalg.inv(I4 + Cay), Sk, atol=1e-10)
assert torch.linalg.eigvals(Cay).add(1).abs().min() > 1e-6        # -1 is never an eigenvalue
```

```python
Mr = TypeIAlgebra([5], [5], complex_valued=False, precision='double')
Ar = Mr.operator(torch.randn(1, 1, 5, 5, dtype=torch.float64))
tr2 = (Ar @ Mr.random_unitary_operator(20_000)).trace ** 2
assert abs(tr2.mean().item() - Ar.frobenius_norm().item() ** 2 / 5) < 5 * tr2.std().item() / math.sqrt(20_000)
```

</details>

**Siguiente:** [Tutorial 12](12_equivariant_learning.md): las capas equivariantes como elementos de un conmutante, en la biblioteca.
