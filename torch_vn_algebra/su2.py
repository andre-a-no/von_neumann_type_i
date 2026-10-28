"""
SU(2) symmetry of spin-1/2 chains: spin multiplets, highest-weight vectors and the Wigner-Eckart theorem.

Mathematics
-----------
The total spin S = sum_i S_i of L spins 1/2 (hbar = 1) acts on H = (C^2)^{(x) L}. By complete reducibility and
Schur's lemma (in the form of Schur-Weyl duality),

    H  =  (+)_S  C^{m_S} (x) V_S,       dim V_S = 2S + 1,   S = L/2, L/2 - 1, ..., >= 0,

with V_S the spin-S irreducible representation and m_S its multiplicity. An operator commutes with SU(2) if and
only if it acts as A_S (x) 1_{2S+1} on every summand, so the commutant of SU(2) is the Type I algebra

    SU(2)'  =  (+)_S  M_{m_S}(C)      (block S of size m_S, one block per total spin),

and an SU(2)-invariant Hamiltonian, e.g. the Heisenberg model H = sum_<ij> J_ij S_i . S_j, is H = (+)_S H_S (x) 1:
the spectrum of H is the union of the spectra of the H_S, each level 2S+1 times degenerate. Compared with the
U(1) sectors N (number of up spins, dimension binom(L, N)) of `chains.SpinChain`, every block shrinks from
binom(L, L/2 + S) to m_S, and the 2S+1 copies of each multiplet are computed once.

Highest-weight vectors. The multiplicity space C^{m_S} is realised inside the U(1) sector N = L/2 + S
(S^z_total = N - L/2 = S) as the space of highest-weight vectors

    H_S^hw = ker S^+  cap  {S^z = S},      dim H_S^hw = m_S,

because a spin-S' multiplet has a state with S^z = S iff S' >= S, and that state is annihilated by S^+ iff
S' = S. With an orthonormal basis B_S (a binom(L, N) x m_S isometry) the block of an invariant operator A is
A_S = B_S^* A_N B_S, where A_N is the U(1) block of A in sector N: A_N maps H_S^hw to itself because A commutes
with S^+ and S^z. The map A -> (A_S)_S is a *-isomorphism of SU(2)' onto (+)_S M_{m_S}.

Multiplicities. Every multiplet with S' >= |M| contributes exactly one state to S^z = M, so the U(1) sector
dimensions are d_{L/2+M} = sum_{S' >= |M|} m_S'. Taking differences,

    m_S = d_{L/2+S} - d_{L/2+S+1} = binom(L, L/2 - S) - binom(L, L/2 - S - 1),      sum_S (2S+1) m_S = 2^L.

(Equivalently: S^+ maps sector N onto sector N + 1 for N >= L/2, and H_S^hw is its kernel.) m_S is also the
number of standard Young tableaux of the two-row shape (L/2 + S, L/2 - S).

Construction of B_S (Yamanouchi basis). Adding one site at a time, a spin-S multiplet of n sites arises from a
spin-(S -+ 1/2) multiplet |S', M'> of the first n - 1 sites coupled with the new spin 1/2 by Clebsch-Gordan
coefficients (Condon-Shortley phases, old (x) new):

    S' = S - 1/2:  |S, S> = |S', S'> (x) |up>
    S' = S + 1/2:  |S, S> = sqrt(2S'/(2S'+1)) |S', S'> (x) |down>  -  sqrt(1/(2S'+1)) |S', S'-1> (x) |up>,
                   |S', S'-1> = S^- |S', S'> / sqrt(2S').

Each basis vector is labelled by its path of intermediate spins (S_1, ..., S_L) = (1/2, ..., S), S_n >= 0,
|S_n - S_{n-1}| = 1/2 (a Yamanouchi symbol / standard Young tableau). Different paths give orthogonal vectors and
each step is an isometry, so B_S is orthonormal by construction; no null space, SVD or Gram-Schmidt is needed and
the cost is O(L * sum_S m_S binom(L, L/2+S)) for all S together (sparse lowering operators only). The blocks
B_S^* H B_S then cost O(m_S^2 binom(L, L/2+S)) dense flops. The basis of the n-site sector N is the
concatenation [(sector N-1 of n-1 sites) (x) up; (sector N of n-1 sites) (x) down], which is exactly the order of
the iterated fused product used by SpinChain (labels in decreasing order), so B_S is expressed directly in the
U(1) basis of `SpinChain(L)`.

Wigner-Eckart theorem. The site spins form a rank-1 tensor operator with spherical components
T_{+1} = -S^+_i / sqrt 2, T_0 = S^z_i, T_{-1} = S^-_i / sqrt 2. For multiplets |S, M, a> (a = Yamanouchi label),

    <S', M', a | T_q | S, M, b>  =  <S M; 1 q | S' M'>  <S', a || S_i || S, b>  /  sqrt(2S' + 1)

(Edmonds' / Racah's normalisation of the reduced matrix element, for which
sum_{M, M', q} |<S' M' a|T_q|S M b>|^2 = |<S' a||S_i||S b>|^2; Sakurai divides by sqrt(2S+1) instead). Only
S' - S in {-1, 0, 1} occur. The reduced elements are computed from highest-weight vectors only, choosing M = S,
q = S' - S, M' = S'. Scalar products of two rank-1 operators then follow without any M sum on the U(1) level:

    <S a| S_i . S_j |S b> = sum_{S''} c(S, S'') sum_c <S a||S_i||S'' c> <S'' c||S_j||S b>,
    c(S, S'') = (-1)^{S - S''} / (2S + 1),

with S_i . S_j = sum_q (-1)^q T^{(i)}_q T^{(j)}_{-q}; `spin_correlation` evaluates c(S, S'') from Clebsch-Gordan
coefficients.

Conventions
-----------
* Basis and sectors as in `chains.SpinChain`: bit i of a label is site i, 1 = up, S^z_i = bit - 1/2, sector N =
  number of up spins, S^z_total = N - L/2; the highest-weight vectors of spin S live in sector N = L/2 + S.
* Sectors of the algebra are ordered by increasing S (S_min = 0 or 1/2, ..., L/2) and carry the charges 2S
  (integers, exact for half-integer spins); `spins` gives the S values. Functions accept spins as int, float
  (multiples of 1/2) or fractions.Fraction.
* H = sum_bonds J_ij S_i . S_j + h sum_i S^z_i, J > 0 antiferromagnetic (the convention of SpinChain.xxz with
  Delta = 1). A uniform field commutes with S^2 and S^z but splits the multiplets: block S of `heisenberg` is H
  restricted to the highest-weight vectors (M = S), i.e. H_S^(h=0) + h S 1, and the multiplet levels are
  E + h (M - S), M = -S..S (`full_spectrum` adds them).

References: E. Lieb and D. Mattis, J. Math. Phys. 3, 749 (1962) (ordering of energy levels); A. R. Edmonds,
Angular Momentum in Quantum Mechanics (Princeton, 1957), ch. 3-5; J. J. Sakurai, Modern Quantum Mechanics, sec.
3.8-3.11; G. Racah, Phys. Rev. 62, 438 (1942) (closed formula for the Clebsch-Gordan coefficients);
H. Weyl / I. Schur (Schur-Weyl duality); T. Yamanouchi, Proc. Phys.-Math. Soc. Japan 19, 436 (1937).

    chain = SU2Chain(L=12, boundary='periodic')
    chain.multiplicities                     # [132, 297, 275, 154, 54, 11, 1] for S = 0..6
    H = chain.heisenberg(J=1.0, J2=0.2)      # Operator of chain.algebra = (+)_S M_{m_S}
    E = chain.full_spectrum(H)               # all 2^L levels, sorted
    gs = chain.ground_state(H)               # energy, spin S, S^z and highest-weight coefficients
    psi = chain.embed(gs.spin[0], gs.coefficients[0, :chain.multiplicity(gs.spin[0])])   # SpinChain vector

Limits: the bases of all spins together hold sum_S m_S binom(L, L/2+S) numbers in float64 (L = 14: 7.2e6,
55 MB; L = 16: 1.0e8, 0.75 GB; L = 18: 1.4e9, 10.6 GB), and the largest block has m_S = 1001 (L = 14), 3640
(L = 16), 13260 (L = 18) rows, at S = 1 or 2, against binom(L, L/2) = 3432, 12870, 48620 for U(1). The
padded algebra of all spins holds C m_max^2 numbers per operator (L = 16: 1.2e8, 0.9 GB real, 1.8 GB complex);
`heisenberg(..., spin=S)` builds one block in `spin_algebra(S)` without padding, and the bases of a single spin
are built by a pruned recursion (L = 18, S = 0: 4862 x 48620, 1.9 GB). Dense methods are practical up to
L = 16 on a workstation (L = 14: about 1 s for the full spectrum on one CPU core, against 10 s for the dense
U(1) sectors); `cost.check_memory` guards the large allocations.
"""
import math
from fractions import Fraction
from numbers import Integral, Rational
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple, Union

import torch

from .algebra import TypeIAlgebra
from . import cost
from .chains import SpinChain, _labels_tensor

Operator = TypeIAlgebra.Operator

Spin = Union[int, float, Fraction]

MAX_L = 18                               # the bases of all spins need 0.75 GB at L = 16, 10.6 GB at L = 18


# ======================================================================================================
# angular-momentum arithmetic (exact, with doubled integers)
# ======================================================================================================
def _twice(j: Spin, name: str = 'j') -> int:
    """2j as an int for j an integer or half-integer (int, float, Fraction); ValueError otherwise."""
    if isinstance(j, bool):
        raise TypeError(f"{name} must be a number, got a bool")
    if isinstance(j, Integral):
        return 2 * int(j)
    if isinstance(j, Rational):
        t = 2 * Fraction(j)
        if t.denominator != 1:
            raise ValueError(f"{name} = {j} is not a multiple of 1/2")
        return int(t)
    if isinstance(j, torch.Tensor):
        j = j.item()                                   # API boundary: an explicit value is requested
    t = 2 * float(j)
    r = round(t)
    if abs(t - r) > 1e-9:
        raise ValueError(f"{name} = {j} is not a multiple of 1/2")
    return int(r)


def _cg_selection(a: int, b: int, c: int, d: int, e: int, f: int) -> bool:
    """Selection rules for <j1 m1 j2 m2 | J M> in doubled integers (a = 2 j1, b = 2 m1, ...)."""
    if min(a, c, e) < 0:
        raise ValueError("angular momenta must be non-negative")
    if b + d != f:
        return False
    for j, m in ((a, b), (c, d), (e, f)):
        if abs(m) > j or (j - m) % 2:
            return False
    return abs(a - c) <= e <= a + c and (a + c + e) % 2 == 0


def clebsch_gordan_exact(j1: Spin, m1: Spin, j2: Spin, m2: Spin, J: Spin, M: Spin) -> Tuple[int, Fraction]:
    """
    Clebsch-Gordan coefficient <j1 m1; j2 m2 | J M> in exact arithmetic, returned as (sign, square) with sign in
    {-1, 0, 1} and square a Fraction, so that the coefficient is sign * sqrt(square). Racah's formula with the
    Condon-Shortley phase convention (<j1 j1; j2 (J - j1) | J J> > 0):

        <j1 m1 j2 m2|J M> = delta_{M, m1+m2} sqrt[(2J+1) (J+j1-j2)! (J-j1+j2)! (j1+j2-J)! / (j1+j2+J+1)!]
                            * sqrt[(J+M)! (J-M)! (j1-m1)! (j1+m1)! (j2-m2)! (j2+m2)!]
                            * sum_k (-1)^k / [k! (j1+j2-J-k)! (j1-m1-k)! (j2+m2-k)! (J-j2+m1+k)! (J-j1-m2+k)!].

    Arguments are integers or half-integers (int, float or Fraction); coefficients violating the selection rules
    (triangle, |m| <= j, j - m integer, M = m1 + m2) are 0.
    """
    a, b, c, d, e, f = (_twice(x, n) for x, n in zip((j1, m1, j2, m2, J, M), ('j1', 'm1', 'j2', 'm2', 'J', 'M')))
    if not _cg_selection(a, b, c, d, e, f):
        return 0, Fraction(0)
    fact = math.factorial
    # all of the following are integers by the selection rules
    j1pj2mJ, Jpj1mj2, Jmj1pj2 = (a + c - e) // 2, (e + a - c) // 2, (e - a + c) // 2
    total = (a + c + e) // 2 + 1
    JpM, JmM = (e + f) // 2, (e - f) // 2
    j1mm1, j1pm1, j2mm2, j2pm2 = (a - b) // 2, (a + b) // 2, (c - d) // 2, (c + d) // 2
    Jmj2pm1, Jmj1mm2 = (e - c + b) // 2, (e - a - d) // 2
    pref = Fraction((e + 1) * fact(Jpj1mj2) * fact(Jmj1pj2) * fact(j1pj2mJ), fact(total))
    pref *= fact(JpM) * fact(JmM) * fact(j1mm1) * fact(j1pm1) * fact(j2mm2) * fact(j2pm2)
    kmin = max(0, -Jmj2pm1, -Jmj1mm2)
    kmax = min(j1pj2mJ, j1mm1, j2pm2)
    s = Fraction(0)
    for k in range(kmin, kmax + 1):
        den = (fact(k) * fact(j1pj2mJ - k) * fact(j1mm1 - k) * fact(j2pm2 - k) * fact(Jmj2pm1 + k)
               * fact(Jmj1mm2 + k))
        s += Fraction(-1 if k % 2 else 1, den)
    if s == 0:
        return 0, Fraction(0)
    return (1 if s > 0 else -1), pref * s * s


def clebsch_gordan(j1: Spin, m1: Spin, j2: Spin, m2: Spin, J: Spin, M: Spin) -> float:
    """<j1 m1; j2 m2 | J M> as a float64 (exact rational arithmetic, one rounding at the end); see
    `clebsch_gordan_exact`."""
    sign, square = clebsch_gordan_exact(j1, m1, j2, m2, J, M)
    return sign * math.sqrt(float(square))         # Fraction -> float is correctly rounded


def wigner_eckart(reduced, j_out: Spin, m_out: Spin, k: Spin, q: Spin, j_in: Spin, m_in: Spin):
    """
    Matrix element <j_out m_out | T^k_q | j_in m_in> = <j_in m_in; k q | j_out m_out> <j_out||T^k||j_in> /
    sqrt(2 j_out + 1) from a reduced matrix element (number or tensor) in Edmonds' normalisation.
    """
    cg = clebsch_gordan(j_in, m_in, k, q, j_out, m_out)
    return reduced * (cg / math.sqrt(_twice(j_out) + 1))


def multiplicity(L: int, S: Spin) -> int:
    """m_S = binom(L, L/2 - S) - binom(L, L/2 - S - 1): the number of spin-S multiplets of L spins 1/2
    (0 if S is not one of L/2, L/2 - 1, ..., >= 0)."""
    t = _twice(S, 'S')
    if t < 0 or t > L or (L - t) % 2:
        return 0
    n = (L - t) // 2
    return math.comb(L, n) - (math.comb(L, n - 1) if n >= 1 else 0)


def two_spins(L: int) -> List[int]:
    """The values 2S of the total spins of L spins 1/2, increasing: L mod 2, L mod 2 + 2, ..., L."""
    return list(range(L % 2, L + 1, 2))


class SU2GroundState(NamedTuple):
    energy: torch.Tensor          # (batch,) lowest energy
    spin: torch.Tensor            # (batch,) total spin S of the ground multiplet (float; smallest S on ties)
    sz: torch.Tensor              # (batch,) S^z of the ground state (-S for h > 0, +S for h <= 0)
    coefficients: torch.Tensor    # (batch, k_max): coordinates in the highest-weight basis B_S, zero padded


# ======================================================================================================
class SU2Chain:
    """
    Spin-1/2 chain of L sites with its SU(2) symmetry: the algebra (+)_S M_{m_S} of SU(2)-invariant operators,
    orthonormal highest-weight bases B_S in the U(1) sectors of `SpinChain`, invariant Heisenberg Hamiltonians
    (nearest and next-nearest neighbours or arbitrary bonds), the full spectrum with its degeneracies, reduced
    matrix elements of the site spins (Wigner-Eckart) and spin correlations.

    The U(1) chain with the same basis conventions is `self.u1` (a SpinChain); `embed`, `multiplet_basis`,
    `restrict` and `lift` translate between the two descriptions.
    """

    def __init__(self, L: int, boundary: str = 'open', complex_valued: bool = True, precision: str = 'double',
                 device=None):
        if L < 1:
            raise ValueError("L must be at least 1")
        if L > MAX_L:
            gb = sum(multiplicity(L, t / 2) * math.comb(L, (L + t) // 2) for t in two_spins(L)) * 8 / 2 ** 30
            raise ValueError(f"L > {MAX_L} is not supported: the highest-weight bases of all spins need about "
                             f"{gb:.0f} GB (dense SU(2) methods reach L = 16 for all spins, L = 18 for single "
                             "spins; use the sparse Lanczos solvers of krylov for ground states beyond)")
        self.u1 = SpinChain(L, boundary, complex_valued, precision, device)
        self.L, self.boundary = L, boundary
        self.complex_valued, self.precision = complex_valued, precision
        self.device, self.dtype, self.real_dtype = self.u1.device, self.u1.dtype, self.u1.real_dtype
        self.two_spins = two_spins(L)                                  # 2S, increasing
        self.spins = [t / 2 for t in self.two_spins]                   # S as floats (exact: multiples of 1/2)
        self.multiplicities = [multiplicity(L, t / 2) for t in self.two_spins]
        self._algebra = None
        self._spin_algebras = {}
        self._bases = {}                                               # {2S: (d_N, m_S) float64}
        self._paths = {}                                               # {2S: list of tuples of 2 S_n}
        self._cast = {}
        self._reduced = None

    # ------------------------------------------------------------------ structure
    @property
    def algebra(self) -> TypeIAlgebra:
        """(+)_S M_{m_S}, sectors ordered by increasing S, charges 2S."""
        if self._algebra is None:
            m = self.multiplicities
            self._algebra = TypeIAlgebra(m, m, complex_valued=self.complex_valued, precision=self.precision,
                                         device=self.device, charges=list(self.two_spins))
        return self._algebra

    def multiplicity(self, S: Spin) -> int:
        return multiplicity(self.L, S)

    def sector_index(self, S: Spin) -> int:
        """Position of spin S among the sectors of `algebra`."""
        t = _twice(S, 'S')
        if t not in self.two_spins:
            raise ValueError(f"S = {S} is not a total spin of {self.L} spins 1/2")
        return self.two_spins.index(t)

    def u1_sector(self, S: Spin, M: Optional[Spin] = None) -> int:
        """U(1) sector N = L/2 + M (default M = S) holding the states |S, M, a>."""
        t = _twice(S, 'S')
        tm = t if M is None else _twice(M, 'M')
        if abs(tm) > t or (t - tm) % 2:
            raise ValueError(f"M = {M} is not a magnetic quantum number of spin S = {S}")
        self.sector_index(S)
        return (self.L + tm) // 2

    @property
    def bonds(self) -> List[Tuple[int, int]]:
        """Nearest-neighbour bonds, as SpinChain.bonds (periodic L = 2: one bond)."""
        return self.u1.bonds

    @property
    def next_nearest_bonds(self) -> List[Tuple[int, int]]:
        """Distinct pairs {i, i+2}: open chains i = 0..L-3; periodic chains i + 2 mod L, every unordered pair once
        (L = 4: the two diagonals once each; L = 3: the pairs coincide with the nearest-neighbour bonds)."""
        L = self.L
        if self.boundary == 'open' or L < 3:
            return [(i, i + 2) for i in range(L - 2)]
        seen, out = set(), []
        for i in range(L):
            p = (min(i, (i + 2) % L), max(i, (i + 2) % L))
            if p not in seen:
                seen.add(p)
                out.append((i, (i + 2) % L))
        return out

    # ------------------------------------------------------------------ highest-weight bases
    def _build_bases(self, targets: Sequence[int]):
        """
        Yamanouchi recursion (module docstring) for n = 1, ..., L sites, float64, for the spins 2S in `targets`
        (stored in the cache). Paths change 2S_n by +-1 per site, so at n sites only the spins within L - n of a
        target are needed: the bases of a single spin cost much less memory than those of all spins.
        """
        L, dev = self.L, self.device
        targets = sorted(set(targets))
        need = sum(math.comb(L, (L + t) // 2) * multiplicity(L, Fraction(t, 2)) for t in targets)
        cost.check_memory(2.5 * need * 8, dev, f"SU(2) highest-weight bases (L={L}, spins {[t / 2 for t in targets]}, "
                                               f"{need:.2e} entries)")
        f64 = torch.float64
        prev = {1: torch.ones(1, 1, dtype=f64, device=dev)}            # one site: |up> spans spin 1/2
        paths = {1: [(1,)]}
        for n in range(2, L + 1):
            labels = _labels_tensor(n - 1, dev)
            index = torch.empty(2 ** (n - 1), dtype=torch.int64, device=dev)
            for lab in labels:
                index[lab] = torch.arange(len(lab), device=dev)
            new, new_paths = {}, {}
            for t in range(n % 2, n + 1, 2):
                if all(abs(t - T) > L - n for T in targets):            # cannot reach a target spin
                    continue
                N = (n + t) // 2
                d_dn = math.comb(n - 1, N)
                cols, pth = [], []
                if t - 1 in prev:                                       # S' = S - 1/2: |S', S'> (x) |up>
                    v = prev[t - 1]
                    cols.append(torch.cat([v, v.new_zeros(d_dn, v.shape[1])]))
                    pth += [p + (t,) for p in paths[t - 1]]
                if t + 1 in prev:                                       # S' = S + 1/2, coupled antiparallel
                    v = prev[t + 1]
                    tp = t + 1                                          # 2 S'
                    a = math.sqrt(tp / (tp + 1))                        # <S' S'; 1/2 -1/2 | S S>
                    b = -math.sqrt(1 / (tp + 1))                        # <S' S'-1; 1/2 1/2 | S S>
                    lv = self._lower_sector(v, N, n - 1, labels, index) / math.sqrt(tp)    # |S', S'-1>
                    cols.append(torch.cat([b * lv, a * v]))
                    pth += [p + (t,) for p in paths[t + 1]]
                new[t] = torch.cat(cols, dim=1)
                new_paths[t] = pth
            prev, paths = new, new_paths
        for t in targets:
            self._bases[t], self._paths[t] = prev[t], paths[t]

    def _ensure_bases(self, targets: Sequence[int]):
        missing = [t for t in targets if t not in self._bases]
        if missing:
            self._build_bases(missing)

    def build_bases(self, spins: Optional[Sequence[Spin]] = None) -> 'SU2Chain':
        """Build (and cache) the highest-weight bases of the given spins (default: all) in one recursion; otherwise
        they are built on first use. Returns self."""
        ts = self.two_spins if spins is None else [self.two_spins[self.sector_index(S)] for S in spins]
        self._ensure_bases(ts)
        return self

    @staticmethod
    def _lower_sector(X: torch.Tensor, N: int, n: int, labels, index) -> torch.Tensor:
        """S^-_total of an n-site chain applied to the columns of X (sector N -> N-1), by gathers."""
        lab_t = labels[N - 1]
        out = X.new_zeros(len(lab_t), X.shape[1])
        for i in range(n):
            sel = torch.nonzero(((lab_t >> i) & 1) == 0).squeeze(-1)
            out.index_add_(0, sel, X[index[lab_t[sel] | (1 << i)]])
        return out

    @staticmethod
    def _raise_sector(X: torch.Tensor, N: int, n: int, labels, index) -> torch.Tensor:
        """S^+_total of an n-site chain (sector N -> N+1)."""
        lab_t = labels[N + 1]
        out = X.new_zeros(len(lab_t), X.shape[1])
        for i in range(n):
            sel = torch.nonzero(((lab_t >> i) & 1) == 1).squeeze(-1)
            out.index_add_(0, sel, X[index[lab_t[sel] ^ (1 << i)]])
        return out

    def lower(self, X: torch.Tensor, N: int) -> torch.Tensor:
        """S^-_total on vectors of the U(1) sector N of `u1`: X of shape (binom(L, N),) or (binom(L, N), r)."""
        v = X if X.dim() == 2 else X[:, None]
        out = self._lower_sector(v, N, self.L, self.u1.labels, self.u1.index_of)
        return out if X.dim() == 2 else out[:, 0]

    def raise_(self, X: torch.Tensor, N: int) -> torch.Tensor:
        """S^+_total on vectors of the U(1) sector N of `u1` (sector N -> N + 1)."""
        v = X if X.dim() == 2 else X[:, None]
        out = self._raise_sector(v, N, self.L, self.u1.labels, self.u1.index_of)
        return out if X.dim() == 2 else out[:, 0]

    def hw_basis(self, S: Spin, dtype: Optional[torch.dtype] = None) -> torch.Tensor:
        """
        Orthonormal basis of the highest-weight vectors of spin S, as the columns of a real (binom(L, N), m_S)
        matrix in the basis of the U(1) sector N = L/2 + S of `u1` (Yamanouchi order, see `paths`). Computed once
        on first use (float64) and cached; dtype: cast (cached), default float64. The cached tensor is returned:
        clone it before modifying it in place.
        """
        t = _twice(S, 'S')
        self.sector_index(S)
        self._ensure_bases([t])
        B = self._bases[t]
        if dtype is None or dtype == B.dtype:
            return B
        key = (t, dtype)
        if key not in self._cast:
            self._cast[key] = B.to(dtype)
        return self._cast[key]

    def paths(self, S: Spin) -> List[Tuple[int, ...]]:
        """Yamanouchi labels of the columns of hw_basis(S): (2 S_1, ..., 2 S_L) with S_n the spin of sites 0..n-1."""
        t = _twice(S, 'S')
        self.sector_index(S)
        self._ensure_bases([t])
        return self._paths[t]

    def multiplet_basis(self, S: Spin, M: Spin) -> torch.Tensor:
        """The states |S, M, a> (a = 1..m_S) as columns in the U(1) sector N = L/2 + M of `u1`:
        (S^-)^{S-M} B_S, normalised by prod sqrt((S + M')(S - M' + 1)) (Condon-Shortley); float64 (for M = S the
        cached hw_basis(S) itself)."""
        t, tm = _twice(S, 'S'), _twice(M, 'M')
        N = self.u1_sector(S, M)
        X = self.hw_basis(S)
        Nc, tc = (self.L + t) // 2, t                                     # current sector and 2M
        cost.check_memory(2 * X.shape[1] * math.comb(self.L, N) * 8, self.device, f"multiplet basis S={S}, M={M}")
        while tc > tm:
            X = self.lower(X, Nc) / math.sqrt((t + tc) * (t - tc + 2) / 4)
            Nc, tc = Nc - 1, tc - 2
        return X

    def embed(self, S: Spin, coefficients: torch.Tensor, M: Optional[Spin] = None) -> torch.Tensor:
        """
        The vector sum_a c_a |S, M, a> (default M = S, the highest weight) in the U(1) sector N = L/2 + M of `u1`,
        for coefficients of shape (m_S,) or (batch, m_S) (e.g. eigenvectors of block S); dtype of `u1`.
        """
        c = torch.as_tensor(coefficients, device=self.device)
        X = self.multiplet_basis(S, S if M is None else M).to(torch.promote_types(c.dtype, torch.float64))
        out = c.to(X.dtype) @ X.T if c.dim() == 2 else X @ c.to(X.dtype)
        return out.to(self.dtype)

    def multiplet(self, S: Spin, coefficients: torch.Tensor) -> List[Tuple[int, torch.Tensor]]:
        """All 2S+1 states of the multiplet with highest-weight coordinates `coefficients`: [(N, vector)] for
        M = S, S-1, ..., -S (vectors in the U(1) sector N of `u1`)."""
        t = _twice(S, 'S')
        return [((self.L + tm) // 2, self.embed(S, coefficients, Fraction(tm, 2))) for tm in range(t, -t - 1, -2)]

    # ------------------------------------------------------------------ invariant operators
    def _apply_bonds(self, N: int, X: torch.Tensor, bonds: Sequence[Tuple[int, int]], weights=None) -> torch.Tensor:
        """sum_b w_b S_i . S_j applied to the columns of X (U(1) sector N of the L-site chain), sparse."""
        u1 = self.u1
        lab = u1.labels[N]
        if not bonds:
            return torch.zeros_like(X)
        I = torch.tensor([i for i, _ in bonds], device=self.device)
        Jx = torch.tensor([j for _, j in bonds], device=self.device)
        bits = ((lab[:, None] >> torch.arange(self.L, device=self.device)) & 1).to(X.dtype) - 0.5    # S^z_i
        w = (torch.ones(len(bonds), dtype=X.dtype, device=self.device) if weights is None
             else torch.as_tensor(weights, dtype=X.dtype, device=self.device))
        out = ((bits[:, I] * bits[:, Jx]) @ w)[:, None] * X                       # S^z_i S^z_j
        for b, (i, j) in enumerate(bonds):                                         # (S^+_i S^-_j + h.c.) / 2
            sel = torch.nonzero(((lab >> i) & 1) != ((lab >> j) & 1)).squeeze(-1)
            src = u1.index_of[lab[sel] ^ ((1 << i) | (1 << j))]
            out.index_add_(0, sel, X[src] * (w[b] / 2))
        return out

    def _project(self, t: int, bonds, weights=None) -> torch.Tensor:
        """B_S^T (sum_b w_b S_i . S_j) B_S, (m_S, m_S), real dtype of the chain."""
        B = self.hw_basis(Fraction(t, 2), self.real_dtype)
        N = (self.L + t) // 2
        d, m = B.shape
        # columns in chunks: the temporaries of the sparse product stay below about 2^25 entries
        chunk = max(1, min(m, 2 ** 25 // max(d, 1)))
        cost.check_memory(3 * d * chunk * B.element_size(), self.device, f"projection on spin {t / 2} (L={self.L})")
        P = torch.empty(m, m, dtype=B.dtype, device=self.device)
        for k in range(0, m, chunk):
            P[:, k:k + chunk] = B.T @ self._apply_bonds(N, B[:, k:k + chunk], bonds, weights)
        return (P + P.T) / 2                                                         # exact symmetry

    def _batch_param(self, x, name) -> torch.Tensor:
        # numbers are converted directly in the chain's precision (no detour through float32)
        x = torch.as_tensor(x, dtype=self.real_dtype, device=self.device) if not isinstance(x, torch.Tensor) \
            else x.to(device=self.device, dtype=self.real_dtype)
        if x.dim() > 1:
            raise ValueError(f"{name} must be a number or a tensor of shape (batch,)")
        return x.reshape(-1)

    def spin_algebra(self, S: Spin) -> TypeIAlgebra:
        """The single block M_{m_S} (no padding, charge 2S), cached; for large L, where the padded algebra of all
        spins is mostly padding."""
        t = _twice(S, 'S')
        m = self.multiplicities[self.sector_index(S)]
        if t not in self._spin_algebras:
            self._spin_algebras[t] = TypeIAlgebra([m], [m], complex_valued=self.complex_valued,
                                                  precision=self.precision, device=self.device, charges=[t])
        return self._spin_algebras[t]

    def heisenberg(self, J=1.0, h=0.0, J2=0.0, bonds: Optional[Sequence[Tuple[int, int]]] = None,
                   spin: Optional[Spin] = None) -> Operator:
        """
        SU(2)-invariant Heisenberg Hamiltonian in the algebra (+)_S M_{m_S}:

            H = J sum_<ij> S_i . S_j + J2 sum_<<ij>> S_i . S_j + h sum_i S^z_i

        (nearest-neighbour bonds `bonds`, next-nearest `next_nearest_bonds`; J2 != 0 gives the frustrated J1-J2
        chain). With an explicit `bonds` list [(i, j), ...] the exchange part is sum_b J_b S_i . S_j with J a number,
        a tensor (n_bonds,) or (batch, n_bonds) (one Hamiltonian per row, e.g. random couplings), and J2 must be 0.
        J, J2, h may be numbers or tensors (batch,). Block S is H restricted to the highest-weight vectors,
        H_S^(h=0) + h S 1 (the field is recorded as `su2_field` on the result for `full_spectrum`).
        spin=None: operator of `algebra` (all spins); spin=S: the block of spin S only, an operator of
        `spin_algebra(S)`. Cost: building the bases once, then O(m_S^2 binom(L, L/2+S)) per term and spin.
        """
        terms, coefs = [], []
        if bonds is None:
            terms.append((self.bonds, None))
            coefs.append(self._batch_param(J, 'J'))
            if isinstance(J2, torch.Tensor) or J2 != 0:
                terms.append((self.next_nearest_bonds, None))
                coefs.append(self._batch_param(J2, 'J2'))
        else:
            bonds = [(int(i), int(j)) for i, j in bonds]
            for i, j in bonds:
                if not (0 <= i < self.L and 0 <= j < self.L) or i == j:
                    raise ValueError(f"invalid bond ({i}, {j}) for L = {self.L}")
            if not (isinstance(J2, (int, float)) and J2 == 0):
                raise ValueError("with explicit bonds, give all couplings in J (J2 must be 0)")
            Jt = torch.as_tensor(J, dtype=self.real_dtype, device=self.device) if not isinstance(J, torch.Tensor) \
                else J.to(device=self.device, dtype=self.real_dtype)
            if Jt.dim() <= 1 and (Jt.dim() == 0 or Jt.shape[0] == len(bonds)):
                terms.append((bonds, None if Jt.dim() == 0 else Jt))
                coefs.append(torch.ones(1, dtype=self.real_dtype, device=self.device) if Jt.dim() else Jt.reshape(1))
            elif Jt.dim() == 2 and Jt.shape[1] == len(bonds):
                for b in range(len(bonds)):
                    terms.append(([bonds[b]], None))
                    coefs.append(Jt[:, b])
            else:
                raise ValueError(f"J must be a number, (n_bonds,) or (batch, n_bonds) with n_bonds = {len(bonds)}")
        hb = self._batch_param(h, 'h')
        batch = max([hb.shape[0]] + [c.shape[0] for c in coefs])
        blocks = []
        sel = range(len(self.two_spins)) if spin is None else [self.sector_index(spin)]
        self._ensure_bases([self.two_spins[k] for k in sel])                        # one recursion for all
        for c_idx in sel:
            t, m = self.two_spins[c_idx], self.multiplicities[c_idx]
            P = torch.stack([self._project(t, bd, w) for bd, w in terms])            # (n_terms, m, m)
            C = torch.stack([c.expand(batch) for c in coefs], dim=1)                 # (batch, n_terms)
            Hs = torch.einsum('bt,tij->bij', C, P)
            Hs = Hs + (hb.expand(batch) * (t / 2))[:, None, None] * torch.eye(m, dtype=Hs.dtype, device=self.device)
            blocks.append(Hs)
        op = (self.algebra if spin is None else self.spin_algebra(spin)).from_blocks(blocks)
        op._is_self_adjoint = True
        op.su2_field = hb.expand(batch)
        return op

    def restrict(self, A: Operator) -> Operator:
        """
        The SU(2)-invariant operator A of the full U(1) algebra `u1.algebra` (e.g. SpinChain.xxz with Delta = 1)
        as an element of (+)_S M_{m_S}: block S = B_S^* A_N B_S with N = L/2 + S (A restricted to the highest-weight
        vectors). For an A that is not SU(2)-invariant this is a compression, not a homomorphism; for an A that
        commutes with S^2 and S^z (e.g. with a uniform field) it is the restriction to S^z = S.
        """
        if A.algebra.k_factors != self.u1.algebra.k_factors:
            raise ValueError("restrict expects an operator of SpinChain(L).algebra (all U(1) sectors)")
        dt = A._dtype()
        self._ensure_bases(self.two_spins)
        blocks = []
        for t in self.two_spins:
            N = (self.L + t) // 2
            B = self.hw_basis(Fraction(t, 2), A._real_dtype()).to(dt)
            blocks.append(B.mT @ A.block(N) @ B)
        op = self.algebra.from_blocks(blocks)
        op._is_self_adjoint = A._is_self_adjoint
        return op

    def lift(self, A: Operator) -> Operator:
        """
        The inverse of `restrict` on SU(2)' : the element (+)_S A_S of (+)_S M_{m_S} as the operator
        (+)_S A_S (x) 1_{2S+1} of the U(1) algebra `u1.algebra`, block N = sum_{S >= |M|} X_{S,M} A_S X_{S,M}^*
        with X_{S,M} = multiplet_basis(S, M). Memory: one dense block per U(1) sector (sum_N binom(L, N)^2).
        """
        alg = self.u1.algebra
        B_, dt = A.batch_size, A._dtype()
        cost.check_memory(sum(cost.tensor_bytes((B_, k, k), dt) for k in alg.k_factors), self.device,
                          f"lift to the U(1) algebra (L={self.L})")
        self._ensure_bases(self.two_spins)
        blocks = [torch.zeros(B_, k, k, dtype=dt, device=self.device) for k in alg.k_factors]
        for c, (t, m) in enumerate(zip(self.two_spins, self.multiplicities)):
            As = A.block(c)
            for tm in range(-t, t + 1, 2):
                X = self.multiplet_basis(Fraction(t, 2), Fraction(tm, 2)).to(dt)
                blocks[(self.L + tm) // 2] += X @ As @ X.mT
        op = alg.from_blocks(blocks)
        op._is_self_adjoint = A._is_self_adjoint
        return op

    def casimir(self) -> Operator:
        """S^2 = (+)_S S(S+1) 1 (central)."""
        return self.algebra.central([s * (s + 1) for s in self.spins])

    # ------------------------------------------------------------------ spectra
    def _require_full(self, H: Operator, what: str):
        if H.algebra.k_factors != self.multiplicities:
            raise ValueError(f"{what} needs an operator of SU2Chain.algebra (all spins); for a single block "
                             "use its eigenvalues directly")

    def _field(self, H: Operator, h) -> torch.Tensor:
        if h is None:
            h = getattr(H, 'su2_field', None)
        B = H.batch_size
        if h is None:
            return torch.zeros(B, dtype=H._real_dtype(), device=self.device)
        return self._batch_param(h, 'h').to(H._real_dtype()).expand(B)

    def full_spectrum(self, H: Operator, h=None, return_quantum_numbers: bool = False):
        """
        All 2^L energies of H (an operator of `algebra` with blocks H restricted to highest-weight vectors, e.g.
        from `heisenberg`), sorted, shape (batch, 2^L): every eigenvalue E of block S gives the 2S+1 levels
        E + h (M - S), M = -S..S. h: the uniform field (number or (batch,)); default the field recorded by
        `heisenberg`, otherwise 0 (pass it explicitly for operators derived by arithmetic).
        return_quantum_numbers: also return S and M of every level (tensors of the same shape).
        """
        self._require_full(H, 'full_spectrum')
        hb = self._field(H, h)
        evs = H.eigenvalues()
        E, Ss, Ms = [], [], []
        for t, e in zip(self.two_spins, evs):
            M = torch.arange(-t, t + 1, 2, device=self.device, dtype=e.dtype) / 2
            lev = e[:, :, None] + hb[:, None, None] * (M - t / 2)[None, None, :]
            E.append(lev.reshape(e.shape[0], -1))
            if return_quantum_numbers:
                Ss.append(torch.full_like(E[-1], t / 2))
                Ms.append(M.expand(e.shape[1], -1).reshape(-1).expand(e.shape[0], -1))
        E = torch.cat(E, dim=1)
        E, order = torch.sort(E, dim=1)
        if not return_quantum_numbers:
            return E
        return E, torch.gather(torch.cat(Ss, 1), 1, order), torch.gather(torch.cat(Ms, 1), 1, order)

    def lowest_energies(self, H: Operator, h=None) -> torch.Tensor:
        """(batch, C): the lowest energy among the states of total spin S, for every S (with a field h the
        multiplet member M = -S for h > 0)."""
        self._require_full(H, 'lowest_energies')
        hb = self._field(H, h)
        evs = H.eigenvalues()
        return torch.stack([e[:, 0] - t * torch.clamp(hb, min=0) for t, e in zip(self.two_spins, evs)], dim=1)

    def ground_state(self, H: Operator, h=None) -> SU2GroundState:
        """Lowest level of H (an operator of `algebra`): energy, total spin S (the smallest S on exact ties), S^z,
        and the eigenvector of block S in the highest-weight basis, zero padded to k_max (pass
        `coefficients[b, :m_S]` to `embed(S, ..., M=sz)` for the vector in the U(1) sector)."""
        self._require_full(H, 'ground_state')
        hb = self._field(H, h)
        alg = self.algebra
        eig = H.eigh()
        E = torch.stack([w[:, 0] - t * torch.clamp(hb, min=0) for t, (w, _) in zip(self.two_spins, eig)], dim=1)
        vecs = torch.stack([torch.nn.functional.pad(V[:, :, 0], (0, alg.k_max - V.shape[1])) for _, V in eig], dim=1)
        E0, c = E.min(dim=1)
        spins = torch.tensor(self.spins, dtype=E.dtype, device=self.device)
        S = spins[c]
        sz = torch.where(hb > 0, -S, S)
        coef = torch.gather(vecs, 1, c[:, None, None].expand(-1, 1, alg.k_max))[:, 0]
        return SU2GroundState(E0, S, sz, coef)

    def expectation(self, A: Operator, S: Spin, coefficients: torch.Tensor) -> torch.Tensor:
        """<psi|A|psi> for psi = sum_a c_a |S, M, a> (any M: A is SU(2)-invariant), c of shape (m_S,) or
        (batch, m_S); A an operator of `algebra` or of `spin_algebra(S)`. Returns (batch,)."""
        c_idx = self.sector_index(S)
        m = self.multiplicities[c_idx]
        if A.algebra.C == 1 and A.algebra.k_factors == [m] and A.algebra.charges == [self.two_spins[c_idx]]:
            block = A.block(0)
        else:
            self._require_full(A, 'expectation')
            block = A.block(c_idx)
        c = torch.as_tensor(coefficients, device=self.device).to(block.dtype)
        c = c[None] if c.dim() == 1 else c
        val = torch.einsum('bi,bij,bj->b', c.conj(), block, c)
        return val.real if torch.is_complex(val) else val

    # ------------------------------------------------------------------ Wigner-Eckart
    def reduced_matrix_elements(self, sites: Optional[Sequence[int]] = None) -> Dict[Tuple[int, int], torch.Tensor]:
        """
        Reduced matrix elements <S', a || S_i || S, b> (Edmonds' normalisation, module docstring) of the site
        spins, from highest-weight vectors only: <S' S' a | T_q | S S b> with q = S' - S, divided by
        <S S; 1 q | S' S'> / sqrt(2S'+1). Returns {(2S', 2S): tensor (n_sites, m_S', m_S)} for all S' - S in
        {-1, 0, 1} (zero blocks for S = S' = 0), real dtype of the chain. sites: default all, in order.
        Hermiticity: R_{S,S'} = (-1)^{S'-S} R_{S',S}^T.
        """
        sites = list(range(self.L)) if sites is None else [int(i) for i in sites]
        u1, dev, rd = self.u1, self.device, self.real_dtype
        self._ensure_bases(self.two_spins)
        out = {}
        for t in self.two_spins:
            S = Fraction(t, 2)
            N = (self.L + t) // 2
            B = self.hw_basis(S, rd)
            lab = u1.labels[N]
            # q = 0, S' = S
            if t > 0:
                cg = clebsch_gordan(S, S, 1, 0, S, S)
                R = torch.stack([B.T @ ((((lab >> i) & 1).to(rd) - 0.5)[:, None] * B) for i in sites])
                out[(t, t)] = R * (math.sqrt(t + 1) / cg)
            else:
                out[(t, t)] = torch.zeros(len(sites), B.shape[1], B.shape[1], dtype=rd, device=dev)
            # q = +1, S' = S + 1: T_{+1} = -S^+_i / sqrt 2
            if t + 2 <= self.L:
                Bp = self.hw_basis(S + 1, rd)
                lt = u1.labels[N + 1]
                vals = []
                for i in sites:
                    sel = torch.nonzero(((lt >> i) & 1) == 1).squeeze(-1)
                    vals.append(Bp[sel].T @ B[u1.index_of[lt[sel] ^ (1 << i)]])
                cg = clebsch_gordan(S, S, 1, 1, S + 1, S + 1)
                out[(t + 2, t)] = torch.stack(vals) * (-math.sqrt(t + 3) / (math.sqrt(2) * cg))
            # q = -1, S' = S - 1: T_{-1} = S^-_i / sqrt 2
            if t - 2 >= 0:
                Bm = self.hw_basis(S - 1, rd)
                lt = u1.labels[N - 1]
                vals = []
                for i in sites:
                    sel = torch.nonzero(((lt >> i) & 1) == 0).squeeze(-1)
                    vals.append(Bm[sel].T @ B[u1.index_of[lt[sel] | (1 << i)]])
                cg = clebsch_gordan(S, S, 1, -1, S - 1, S - 1)
                out[(t - 2, t)] = torch.stack(vals) * (math.sqrt(t - 1) / (math.sqrt(2) * cg))
        return out

    @staticmethod
    def _scalar_coefficient(t: int, tpp: int) -> float:
        """c(S, S'') = sum_q (-1)^q <S'' S-q; 1 q|S S> <S S; 1 -q|S'' S-q> / sqrt((2S+1)(2S''+1))."""
        S, Spp = Fraction(t, 2), Fraction(tpp, 2)
        c = 0.0
        for q in (-1, 0, 1):
            Mpp = S - q
            if abs(Mpp) > Spp:
                continue
            c += (-1) ** (q % 2) * clebsch_gordan(Spp, Mpp, 1, q, S, S) * clebsch_gordan(S, S, 1, -q, Spp, Mpp)
        return c / math.sqrt((t + 1) * (tpp + 1))

    def spin_correlation(self, pairs: Sequence[Tuple[int, int]], method: str = 'wigner_eckart') -> Operator:
        """
        The SU(2)-invariant operators S_i . S_j for pairs [(i, j), ...] as a batch of operators of `algebra`
        (batch index = pair). method='wigner_eckart': from the reduced matrix elements of the site spins,
        block S = sum_{S''} c(S, S'') R^i_{S,S''} R^j_{S'',S}; method='direct': B_S^T (S_i . S_j) B_S with the
        sparse bond operator. Expectation values in a multiplet: `expectation(C, S, coefficients)`.
        """
        pairs = [(int(i), int(j)) for i, j in pairs]
        blocks = []
        if method == 'direct':
            self._ensure_bases(self.two_spins)
            for t, m in zip(self.two_spins, self.multiplicities):
                bl = []
                for i, j in pairs:
                    if i == j:
                        bl.append(0.75 * torch.eye(m, dtype=self.real_dtype, device=self.device))
                    else:
                        bl.append(self._project(t, [(i, j)]))
                blocks.append(torch.stack(bl))
        elif method == 'wigner_eckart':
            for i, j in pairs:
                if not (0 <= i < self.L and 0 <= j < self.L):
                    raise ValueError(f"invalid pair ({i}, {j}) for L = {self.L}")
            if self._reduced is None:                       # all sites, cached
                self._reduced = self.reduced_matrix_elements()
            R = self._reduced
            I = torch.tensor([i for i, _ in pairs], device=self.device)
            Jx = torch.tensor([j for _, j in pairs], device=self.device)
            for t, m in zip(self.two_spins, self.multiplicities):
                acc = torch.zeros(len(pairs), m, m, dtype=self.real_dtype, device=self.device)
                for tpp in (t - 2, t, t + 2):
                    if (tpp, t) not in R or (t, tpp) not in R:
                        continue
                    c = self._scalar_coefficient(t, tpp)
                    if c == 0:
                        continue
                    acc = acc + c * (R[(t, tpp)][I] @ R[(tpp, t)][Jx])
                blocks.append((acc + acc.mT) / 2)
        else:
            raise ValueError("method must be 'wigner_eckart' or 'direct'")
        op = self.algebra.from_blocks(blocks)
        op._is_self_adjoint = True
        return op

    def __repr__(self) -> str:
        return (f"SU2Chain(L={self.L}, boundary={self.boundary!r}, spins={self.spins}, "
                f"multiplicities={self.multiplicities})")
