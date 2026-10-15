**English** · [Español](es/10_quantum_computing.md) · [Français](fr/10_quantum_computing.md) · [Deutsch](de/10_quantum_computing.md) · [中文](zh/10_quantum_computing.md) · [日本語](ja/10_quantum_computing.md) · [Русский](ru/10_quantum_computing.md)

# Tutorial 10. Quantum computing

This tutorial accompanies [Chapter 10](../lecture_notes/10_quantum_computing.md) of the lecture notes, which describes noise on a quantum
computer by the $`*`$-algebra of its Kraus operators. We compute fixed-point algebras and noiseless subsystems
([Theorem 10.2](../lecture_notes/10_quantum_computing.md#r-10-2)), check the Knill–Laflamme theorem ([Theorem 10.5](../lecture_notes/10_quantum_computing.md#r-10-5)) and the correction of a
sector label ([Proposition 10.6](../lecture_notes/10_quantum_computing.md#r-10-6)) on the repetition code, and simulate an XY-mixer ansatz with amplitude
damping and post-selection in the algebra $`\bigoplus_wM_{\binom nw}`$ of Hamming-weight sectors
([Propositions 10.8](../lecture_notes/10_quantum_computing.md#r-10-8)–[10.10](../lecture_notes/10_quantum_computing.md#r-10-10)).

<a id="s-10-1"></a>

## 10.1 Noise algebras, noiseless subsystems and fixed points

*Lecture notes:* [Section 10.1](../lecture_notes/10_quantum_computing.md#s-10-1); [Theorem 10.2](../lecture_notes/10_quantum_computing.md#r-10-2), [Corollary 10.3](../lecture_notes/10_quantum_computing.md#r-10-3), [Proposition 10.4](../lecture_notes/10_quantum_computing.md#r-10-4).

By [lecture notes, Corollary 10.3](../lecture_notes/10_quantum_computing.md#r-10-3), collective dephasing $`\Phi(\rho)=\tfrac12(U\rho U^*+U^*\rho U)`$,
$`U=e^{-i\theta Z_{\mathrm{tot}}}`$, on $`n`$ qubits has the interaction algebra $`\mathcal A=\mathrm{span}\{P_0,\dots,P_n\}`$
($`P_w`$ the projection onto the span $`H_w`$ of the basis vectors of Hamming weight $`w`$), and
$`\mathrm{Fix}(\Phi)=\mathcal A'=\bigoplus_wB(H_w)`$ ([Theorem 10.2](../lecture_notes/10_quantum_computing.md#r-10-2)).

Commutants are null spaces of $`T\mapsto(AT-TA)_A`$ (Chapter 2), fixed points that of
$`\sum_iE_i\otimes\bar E_i-1`$ (the matrix of $`\Phi-\mathrm{id}`$ on row-major vectorised matrices). For $`n=4`$
the three spaces of the corollary have dimensions $`5`$, $`70=1+16+36+16+1`$ and $`70`$:

```python
import itertools
import math
import torch

torch.manual_seed(0)
dt = torch.complex128
I2, X = torch.eye(2, dtype=dt), torch.tensor([[0, 1], [1, 0]], dtype=dt)
Y, Z = torch.tensor([[0, -1j], [1j, 0]], dtype=dt), torch.tensor([[1, 0], [0, -1]], dtype=dt)

def on(op, j, n):                             # op on qubit j of n (j = 0 is the leftmost factor)
    out = torch.ones(1, 1, dtype=dt)
    for i in range(n):
        out = torch.kron(out, op if i == j else I2)
    return out

def null_space(L, tol=1e-9):                  # orthonormal columns spanning {x : L x = 0}
    w, V = torch.linalg.eigh(L.mH @ L)
    return V[:, w < tol]

def commutant_basis(ops):
    d, ops = ops[0].shape[0], [A.resolve_conj().contiguous() for A in ops]
    eye = torch.eye(d, dtype=dt)
    L = torch.cat([torch.kron(A, eye) - torch.kron(eye, A.T.contiguous()) for A in ops])
    return [v.reshape(d, d) for v in null_space(L).T]

def fixed_points(kraus):                      # vec(E X E^*) = (E (x) conj E) vec(X)
    d = kraus[0].shape[0]
    S = sum(torch.kron(E.resolve_conj(), E.conj().resolve_conj()) for E in kraus)
    return [v.reshape(d, d) for v in null_space(S - torch.eye(d * d, dtype=dt)).T]

n, theta = 4, 0.3
U = torch.linalg.matrix_exp(-1j * theta * sum(on(Z, j, n) for j in range(n)))
kraus = [U / math.sqrt(2), U.mH / math.sqrt(2)]                  # unital: U U^* = U^* U = 1
comm = commutant_basis(kraus + [E.mH for E in kraus])             # A'
A_int, fix = commutant_basis(comm), fixed_points(kraus)           # A = A'' (Theorem 2.2), Fix(Phi)
print(len(A_int), len(comm), len(fix))
assert (len(A_int), len(comm), len(fix)) == (5, 70, 70)
span = torch.linalg.svdvals(torch.stack([F.reshape(-1) for F in comm + fix]))
assert int((span > 1e-9).sum()) == 70                             # Fix(Phi) = A' (Theorem 10.2)
```

$`\bigoplus_wM_{\binom nw}`$ is the algebra of a `SpinChain` (Chapter 9), whose sectors are the numbers
of up spins, i.e. Hamming weights. There the Kraus operators are central,
$`U=\sum_we^{-i\theta(n-2w)}P_w`$, and the channel is the identity: the library stores the noiseless part.

```python
from torch_vn_algebra import SpinChain, TypeIAlgebra, DensityMatrix, channels, states, dynamics

M = SpinChain(n).algebra                                          # (+)_w M_binom(4, w)
print(M.k_factors)
assert sum(k * k for k in M.k_factors) == len(fix)
phases = [torch.exp(torch.tensor(-1j * theta * (n - 2 * w), dtype=dt)) for w in range(n + 1)]
Phi_M = channels.kraus_channel([M.central([z / math.sqrt(2) for z in phases]),
                                M.central([z.conj() / math.sqrt(2) for z in phases])])
rho = states.random_density_matrix(M, batch_size=3)
assert (Phi_M(rho).matrix - rho.matrix).abs().max() < 1e-12      # no state of M is disturbed
```

Collective rotations $`e^{-i\theta\,\mathbf n\cdot\mathbf S}`$, $`S_\alpha=\tfrac12\sum_j\sigma_{\alpha,j}`$, lie in the $`*`$-algebra
$`\mathcal B`$ generated by $`1,S_x,S_y,S_z`$. For $`n=4`$, [lecture notes, Section 10.1](../lecture_notes/10_quantum_computing.md#s-10-1) shows $`\dim\mathcal B=35`$,
$`\dim\mathcal B'=14`$ and $`\dim K=2`$ for $`K=\ker S_x\cap\ker S_y\cap\ker S_z`$; the code finds the same numbers and
checks [Proposition 10.4](../lecture_notes/10_quantum_computing.md#r-10-4): a channel with random collective rotations as Kraus operators fixes a
state supported in $`K`$, and moves a generic one.

```python
S = [sum(on(P, j, n) for j in range(n)) / 2 for P in (X, Y, Z)]
comm_S = commutant_basis(S)
print(len(comm_S), len(commutant_basis(comm_S)))
assert len(comm_S) == 14 and len(commutant_basis(comm_S)) == 35
K = null_space(torch.cat(S))                                      # columns: a basis of K
assert K.shape[1] == 2
axes = torch.nn.functional.normalize(torch.randn(3, 3, dtype=torch.float64), dim=1)
rotations = [torch.linalg.matrix_exp(-1j * t * sum(c * Sa for c, Sa in zip(ax, S))) / math.sqrt(3)
             for t, ax in zip(torch.rand(3, dtype=torch.float64) * math.pi, axes)]
apply = lambda kr, r: sum(E @ r @ E.mH for E in kr)
for G in (K @ torch.randn(2, 2, dtype=dt), torch.randn(16, 16, dtype=dt)):  # supported in K / generic
    r = G @ G.mH / (G @ G.mH).trace()
    assert ((apply(rotations, r) - r).abs().max() < 1e-12) == (G.shape[1] == 2)
```

<a id="s-10-2"></a>

## 10.2 Error correction: the Knill–Laflamme theorem

*Lecture notes:* [Section 10.2](../lecture_notes/10_quantum_computing.md#s-10-2); [Theorem 10.5](../lecture_notes/10_quantum_computing.md#r-10-5).

The bit-flip example of [lecture notes, Section 10.2](../lecture_notes/10_quantum_computing.md#s-10-2): for $`\mathcal C=\mathrm{span}\{|000\rangle,|111\rangle\}`$ and
"at most one bit flip" the Knill–Laflamme condition of [Theorem 10.5](../lecture_notes/10_quantum_computing.md#r-10-5) holds with
$`\Lambda=\mathrm{diag}(p_0,\dots,p_3)`$, and the recovery $`\mathcal R(\sigma)=\sum_{k=0}^3PX_k\sigma X_kP`$ restores every
state of the code. If every qubit flips independently with probability $`p`$, the condition fails, and the
fidelity of $`|0_L\rangle`$ after recovery is $`1-3p^2+2p^3=0.972`$ for $`p=0.1`$. The helper `kl_matrix` returns
$`\Lambda`$, or `None` if the condition fails. In the library $`M_8`$ is an algebra with one block, and `R @ Phi`
is the composition $`\mathcal R\circ\Phi`$:

```python
M8 = TypeIAlgebra([8], [8], precision='double')
op8 = lambda T: M8.operator(T.reshape(1, 1, 8, 8))
state8 = lambda v: DensityMatrix(M8, matrix=torch.outer(v, v.conj()).reshape(1, 1, 8, 8))
ket = lambda bits: torch.eye(8, dtype=dt)[int(bits, 2)]
def kl_matrix(kraus, P, tol=1e-12):
    """(lambda_ij) with P E_i^* E_j P = lambda_ij P, or None if the Knill-Laflamme condition fails."""
    lam = torch.tensor([[(P @ a.mH @ b @ P).trace() / P.trace() for b in kraus] for a in kraus])
    ok = all((P @ a.mH @ b @ P - lam[i, j] * P).abs().max() < tol
             for i, a in enumerate(kraus) for j, b in enumerate(kraus))
    return lam if ok else None
def independent(op, prob):                    # op on each qubit independently with probability prob
    return [math.sqrt(prob ** sum(s) * (1 - prob) ** (3 - sum(s)))
            * torch.linalg.multi_dot([on(op, j, 3) if s[j] else torch.eye(8, dtype=dt) for j in range(3)])
            for s in itertools.product([0, 1], repeat=3)]

P = torch.outer(ket('000'), ket('000')) + torch.outer(ket('111'), ket('111'))
flips = [torch.eye(8, dtype=dt)] + [on(X, j, 3) for j in range(3)]
E_one = [math.sqrt(pk) * Xk for pk, Xk in zip((0.7, 0.1, 0.1, 0.1), flips)]
assert torch.allclose(kl_matrix(E_one, P), torch.diag(torch.tensor([0.7, 0.1, 0.1, 0.1], dtype=dt)))
Phi_one = channels.kraus_channel([op8(E) for E in E_one])
R = channels.kraus_channel([op8(P @ Xk) for Xk in flips])        # syndrome measurement + correction
assert R.is_trace_preserving()
a, b = torch.nn.functional.normalize(torch.randn(2, dtype=dt), dim=0)
rho_L = state8(a * ket('000') + b * ket('111'))
assert abs((R @ Phi_one)(rho_L).fidelity(rho_L).item() - 1) < 1e-12       # Theorem 10.5
p = 0.1                                                           # independent flips: KL fails
E_ind = independent(X, p)
assert kl_matrix(E_ind, P) is None
zero_L = state8(ket('000'))
F = (R @ channels.kraus_channel([op8(E) for E in E_ind]))(zero_L).fidelity(zero_L).item()
print(f"fidelity after recovery {F:.4f}, exact {1 - 3 * p**2 + 2 * p**3:.4f}")
assert abs(F - (1 - 3 * p**2 + 2 * p**3)) < 1e-12
```

<a id="s-10-3"></a>

## 10.3 Operator algebra error correction

*Lecture notes:* [Section 10.3](../lecture_notes/10_quantum_computing.md#s-10-3); [Proposition 10.6](../lecture_notes/10_quantum_computing.md#r-10-6).

The example of [lecture notes, Section 10.3](../lecture_notes/10_quantum_computing.md#s-10-3): the repetition code under independent phase flips
(probability $`q`$ per qubit) followed by at most one bit flip. The qubit is not protected (the Knill–Laflamme
condition fails), but the classical bit is: condition 2 of [Proposition 10.6](../lecture_notes/10_quantum_computing.md#r-10-6) holds for
$`P_0=|000\rangle\langle000|`$, $`P_1=|111\rangle\langle111|`$, and the majority vote $`N_0=\sum_{|b|\le1}|b\rangle\langle b|`$,
$`N_1=1-N_0`$ reads the bit. After the recovery of Section 10.2, the fidelity of $`|{+_L}\rangle`$ is
$`(1+(1-2q)^3)/2=0.8645`$ for $`q=0.05`$.

```python
q = 0.05
Phi2 = Phi_one @ channels.kraus_channel([op8(E) for E in independent(Z, q)])   # dephasing, then a flip
E2 = [Phi2.kraus[0, i, 0] for i in range(Phi2.kraus_rank)]        # its 32 Kraus operators
P0, P1 = torch.outer(ket('000'), ket('000')), torch.outer(ket('111'), ket('111'))
assert kl_matrix(E2, P) is None                                   # the qubit is not protected
assert max((P0 @ e.mH @ f @ P1).abs().max().item() for e in E2 for f in E2) == 0      # condition 2
N0 = torch.diag(torch.tensor([float(bin(i).count('1') <= 1) for i in range(8)], dtype=dt))
for N, Pd in ((N0, P0), (torch.eye(8, dtype=dt) - N0, P1)):       # P Phi^*(N_d) P = P_d
    assert (P @ sum(e.mH @ N @ e for e in E2) @ P - Pd).abs().max() < 1e-12
votes = states.born_probabilities(Phi2(rho_L), [op8(N0), op8(torch.eye(8, dtype=dt) - N0)])[0]
assert torch.allclose(votes, torch.stack([abs(a) ** 2, abs(b) ** 2]))      # the classical bit survives
plus_L = state8((ket('000') + ket('111')) / math.sqrt(2))
F_plus = (R @ Phi2)(plus_L).fidelity(plus_L).item()
print(f"fidelity of |+_L> {F_plus:.4f}, exact {(1 + (1 - 2 * q)**3) / 2:.4f}")
assert abs(F_plus - (1 + (1 - 2 * q) ** 3) / 2) < 1e-12
```

<a id="s-10-4"></a>

## 10.4 Constrained mixers

*Lecture notes:* [Section 10.4](../lecture_notes/10_quantum_computing.md#s-10-4); [Proposition 10.8](../lecture_notes/10_quantum_computing.md#r-10-8), [Proposition 10.9](../lecture_notes/10_quantum_computing.md#r-10-9), [Proposition 10.10](../lecture_notes/10_quantum_computing.md#r-10-10).

The quantum alternating operator ansatz alternates $`e^{-i\gamma H_C}`$, $`H_C`$ a diagonal cost Hamiltonian,
and $`e^{-i\beta H_M}`$, $`H_M=\tfrac14\sum_{(i,j)}(X_iX_j+Y_iY_j)`$ the XY mixer on a ring. By
[lecture notes, Proposition 10.8](../lecture_notes/10_quantum_computing.md#r-10-8) all of them lie in $`\mathcal M=\{P_0,\dots,P_n\}'\cong\bigoplus_wM_{\binom nw}`$,
so the ansatz never leaves the Hamming-weight sector of its initial state.

In the library $`\mathcal M`$ is `SpinChain(n).algebra` (bit $`1`$ = spin up, so $`S^z_i=-\tfrac12Z_i`$), and $`H_M`$
is its XXZ Hamiltonian $`J\sum_{(i,j)}\big(\tfrac12(S^+_iS^-_j+S^-_iS^+_j)+\Delta S^z_iS^z_j\big)`$ with $`J=1`$,
$`\Delta=0`$, as $`S^+_iS^-_j+S^-_iS^+_j=\tfrac12(X_iX_j+Y_iY_j)`$; we check it against the dense
$`64\times64`$ matrix. For $`H_C`$ we take graph bisection on the triangles $`\{0,1,2\}`$, $`\{3,4,5\}`$ joined
by the edges $`(2,3)`$, $`(5,0)`$: $`H_C=\sum_{(i,j)\in E}\tfrac12(1-Z_iZ_j)`$, $`Z_iZ_j=4S^z_iS^z_j`$, counts the cut
edges. The best bisection cuts $`2`$ edges (the two triangles; any other set of three vertices splits both
triangles, and a split triangle has two cut edges); we start from $`|101010\rangle`$, which cuts $`6`$. Finally the whole
ansatz is recomputed with dense $`64\times64`$ matrices, independently of the library.

```python
n, w = 6, 3
chain = SpinChain(n, 'periodic')
M = chain.algebra                                                 # (+)_w M_binom(6, w), padded to 20
H_M = chain.xxz(J=1.0, Delta=0.0)
H_dense = sum(on(X, i, n) @ on(X, j, n) + on(Y, i, n) @ on(Y, j, n) for i, j in chain.bonds) / 4
spec = torch.sort(torch.cat([e[0] for e in H_M.eigenvalues()]))[0]
assert torch.allclose(spec, torch.linalg.eigvalsh(H_dense))       # same spectrum, sector by sector
edges = chain.bonds + [(0, 2), (3, 5)]                            # two triangles joined by (2,3), (5,0)
H_C = M.zero()
for i, j in edges:
    H_C = H_C + (M.identity() - chain.sz(i) @ chain.sz(j) * 4.0) * 0.5
layers = [dynamics.propagator(H_M, beta) @ dynamics.propagator(H_C, gamma)
          for gamma, beta in ((0.6, 1.0), (1.4, 1.4))]
rho0 = ideal = chain.basis_state('101010')
for U_l in layers:
    ideal = channels.unitary_channel(U_l)(ideal)
assert abs(ideal.sector_probabilities()[0, w].item() - 1) < 1e-12           # Proposition 10.8

HC_dense = sum(torch.eye(64, dtype=dt) - on(Z, i, n) @ on(Z, j, n) for i, j in edges) / 2
psi = torch.eye(64, dtype=dt)[int('101010', 2)]                   # the same circuit on 6 full qubits
for gamma, beta in ((0.6, 1.0), (1.4, 1.4)):
    psi = torch.linalg.matrix_exp(-1j * beta * H_dense) @ torch.linalg.matrix_exp(-1j * gamma * HC_dense) @ psi
weight = torch.tensor([bin(b).count('1') for b in range(64)])
assert abs(psi[weight == w].norm().item() - 1) < 1e-12           # Proposition 10.8, dense
assert abs((psi.conj() @ HC_dense @ psi).real.item() - states.expectation(ideal, H_C).real.item()) < 1e-10
```

**Amplitude damping** of qubit $`i`$ ($`|1\rangle\to|0\rangle`$ with probability $`g`$) has the Kraus operators
$`K_0^{(i)}=|0\rangle\langle0|_i+\sqrt{1-g}\,|1\rangle\langle1|_i`$ and $`K_1^{(i)}=\sqrt g\,|0\rangle\langle1|_i`$, which lowers
the weight: a map between sectors, `chain.site_amplitude_damping(i, g)` (an `InterSectorChannel`).

By [lecture notes, Proposition 10.9](../lecture_notes/10_quantum_computing.md#r-10-9), if every layer of the circuit is followed by the damping
of every qubit, the weight $`w`$ of the initial state is still found with probability $`(1-g)^{wp}`$ after $`p`$ layers,
and the post-selected state is the noiseless one.

An `InterSectorChannel` keeps only the diagonal blocks of its output; here there are no others (the proof
of Proposition 10.9 applies verbatim to the damping of one qubit, with $`s\in\{0,1\}`$), so nothing is lost. With $`g=0.05`$ and $`p=2`$ the acceptance is $`0.95^6=0.7351`$:

```python
g = 0.05
damping = [chain.site_amplitude_damping(i, g) for i in range(n)]
noisy = rho0
for U_l in layers:
    noisy = channels.unitary_channel(U_l)(noisy)
    for D in damping:                                             # one qubit after the other
        noisy = D(noisy)
print(noisy.sector_probabilities()[0])                            # weight leaks to w = 2, 1, 0
P_w = M.central([1.0 if v == w else 0.0 for v in range(n + 1)])
post, accept = noisy.condition_on(P_w)                            # measure the weight, keep w = 3
cut = lambda r: states.expectation(r, H_C).real.item()
print(f"P(weight 3) {accept.item():.4f}, fidelity {noisy.fidelity(ideal).item():.4f}, "
      f"post-selected {post.fidelity(ideal).item():.4f}; "
      f"cut: ideal {cut(ideal):.4f}, noisy {cut(noisy):.4f}, post-selected {cut(post):.4f}")
assert abs(accept.item() - (1 - g) ** (w * len(layers))) < 1e-12              # Proposition 10.9
assert abs(post.fidelity(ideal).item() - 1) < 1e-10 and abs(cut(post) - cut(ideal)) < 1e-10
```

Numerical observation (the printed values): the infeasible part lowers the expected cut from $`4.3789`$ to
$`4.3123`$, a spurious improvement that would mislead an optimiser of the angles; after post-selection the
value is exact (Proposition 10.9).

By [lecture notes, Proposition 10.10](../lecture_notes/10_quantum_computing.md#r-10-10), under the Lindblad equation with amplitude damping
(rate $`\gamma_1`$) and dephasing $`L_k=S^z_k`$ (rate $`\gamma_2`$), a state starting in $`H_w`$ keeps the weight $`w`$ with
probability $`e^{-\gamma_1wt}`$, and the post-selected state does not depend on $`\gamma_1`$; for $`\gamma_2=0`$
post-selection is perfect. The companion paper
(`examples/paper/ex9_qaoa_xy_mixer.py`) evolves $`|111000\rangle`$ under $`H_M`$ for $`T=2`$; we reproduce its
row $`\gamma_1=\gamma_2=0.1`$: $`P(w=3)=e^{-0.6}=0.5488`$ (Proposition 10.10), and, as a numerical observation of the
code below, fidelity $`0.4805`$ and post-selected fidelity $`0.8755`$; by Proposition 10.10 the latter is the same as
for $`\gamma_1=0`$.

```python
T = 2.0
jumps = [chain.lowering(i) for i in range(n)] + [chain.sz(i) for i in range(n)]
rho0 = chain.basis_state('111000')                                # the initial state of the paper
ideal_T = dynamics.von_neumann(rho0, H_M, [T])[0]
post_T = {}
for g1, g2 in ((0.1, 0.0), (0.0, 0.1), (0.1, 0.1)):          # the last one is printed
    rho_T = dynamics.lindblad_evolve(rho0, H_M, jumps, [0.0, T], rates=[g1] * n + [g2] * n, substeps=80)[-1]
    post_T[g1, g2], accept_T = rho_T.condition_on(P_w)
    assert abs(accept_T.item() - math.exp(-w * g1 * T)) < 1e-8                 # Proposition 10.10
print(f"P(weight 3) {accept_T.item():.4f}, fidelity {rho_T.fidelity(ideal_T).item():.4f}, "
      f"post-selected {post_T[0.1, 0.1].fidelity(ideal_T).item():.4f}")
assert abs(post_T[0.1, 0.0].fidelity(ideal_T).item() - 1) < 1e-6               # gamma_2 = 0: perfect
assert post_T[0.1, 0.1].trace_distance(post_T[0.0, 0.1]).item() < 1e-6       # independent of gamma_1
```

<a id="s-10-5"></a>

## 10.5 What to compare with

Qiskit (`qiskit.quantum_info`, Aer noise models, `QAOAAnsatz` with a custom mixer), Cirq, PennyLane
(`qml.qaoa.xy_mixer`, `default.mixed`) and QuTiP/`qutip-qip` simulate all of this on $`n`$ full qubits
and are the right tools for circuits; stim samples Clifford circuits with Pauli-type noise on thousands
of qubits (amplitude damping is not a Pauli channel: Pauli channels are unital, amplitude damping is not). We know of no function in these packages that
computes the interaction algebra of a channel, its commutant or its whole fixed-point algebra (QuTiP's
`steadystate` returns one stationary state of a Lindbladian); for numerical block decompositions of matrix $`*`$-algebras see
Murota et al., *Japan J. Indust. Appl. Math.* 27 (2010). The library adds storage in
$`\bigoplus_wM_{\binom nw}`$, maps between sectors and post-selection as a central projection, batched.
One sector needs $`\binom nw^2`$ entries instead of $`4^n`$ ($`1.7\cdot10^8`$ versus $`4.3\cdot10^9`$ for $`n=16`$,
$`w=8`$), but the full algebra is padded to $`(n+1)\binom n{\lfloor n/2\rfloor}^2`$ entries ($`2.8\cdot10^9`$), a
modest saving. It is not a circuit simulator and has no stabilizer formalism.

<a id="x-10"></a>

## Exercises

[Exercises 1–3](../lecture_notes/10_quantum_computing.md#x-10) of the lecture notes, numerically: for amplitude damping on one qubit,
$`\mathcal A'`$, $`\mathrm{Fix}(\Phi)`$ and $`\mathrm{Fix}(\Phi^*)`$ are one-dimensional (Exercise 1); for the channel of
[Corollary 10.3](../lecture_notes/10_quantum_computing.md#r-10-3) with $`\theta=\pi/2`$ and $`n=4`$, $`\dim\mathrm{Fix}(\Phi)=128`$ (Exercise 2); the recovery
$`\mathcal R`$ of the bit-flip code corrects the coherent error $`e^{-i\phi X_1}`$ (Exercise 3).

<details>
<summary>Solutions</summary>

```python
E_ad = [torch.tensor([[1, 0], [0, math.sqrt(0.7)]], dtype=dt), torch.tensor([[0, math.sqrt(0.3)], [0, 0]], dtype=dt)]
fix_ad = fixed_points(E_ad)
assert len(fix_ad) == 1 and len(commutant_basis(E_ad + [E.mH for E in E_ad])) == 1
assert torch.allclose(fix_ad[0] / fix_ad[0][0, 0], torch.tensor([[1, 0], [0, 0]], dtype=dt))
assert len(fixed_points([E.mH for E in E_ad])) == 1                     # Fix(Phi^*) = C 1
```

```python
U2 = torch.linalg.matrix_exp(-1j * (math.pi / 2) * sum(on(Z, j, 4) for j in range(4)))
assert len(fixed_points([U2 / math.sqrt(2), U2.mH / math.sqrt(2)])) == 128
assert sum(math.comb(4, v) ** 2 for v in range(5)) == math.comb(8, 4) == 70
```

```python
for phi in (0.2, 1.0, 2.5):
    G1 = torch.linalg.matrix_exp(-1j * phi * on(X, 0, 3))
    out = (R @ channels.kraus_channel([op8(G1)]))(rho_L)
    assert (out.matrix - rho_L.matrix).abs().max() < 1e-12       # corrected exactly
```

</details>

**Next:** [Tutorial 11](11_optimisation.md): optimisation over the unitary groups of an algebra, and trace inequalities.
