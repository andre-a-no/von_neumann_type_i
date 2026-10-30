"""
Unital inclusions of finite-dimensional von Neumann algebras, conditional expectations and index.

Inclusion matrices
------------------
Let N = (+)_{c=1}^{C} M_{k_c} and M = (+)_{d=1}^{D} M_{m_d}. Every unital injective *-homomorphism
iota: N -> M is, up to an inner automorphism of M, of the standard form

    iota_0(x)_d = (+)_c  x_c (x) 1_{lambda_dc}        (summands in the order c = 1, ..., C),

given by the inclusion matrix Lambda = (lambda_dc) in Z_{>=0}^{D x C}: lambda_dc is the multiplicity of
the irreducible representation of the summand M_{k_c} in that of M_{m_d}. Unitality is the condition
m_d = sum_c lambda_dc k_c, injectivity says that every column of Lambda is non-zero. The general
inclusion is iota = Ad(u) o iota_0 with a unitary u of M ([GHJ, Ch. 2]); `Inclusion` stores Lambda
and the optional u (one unitary per sample: a batch). In block d the basis vector
off_dc + a lambda_dc + i (a < k_c, i < lambda_dc, off_dc = sum_{c' < c} lambda_dc' k_c') is e_a (x) e_i of
the i-th copy of C^{k_c}, the convention of composite.kron; the isometries

    V_{dc,i}: C^{k_c} -> C^{m_d},  e_a -> e_{off_dc + a lambda_dc + i},   iota_0(x)_d = sum_{c,i} V_{dc,i} x_c V_{dc,i}^*,

are the Kraus operators of iota as an InterSectorChannel N -> M (`Inclusion.as_channel`). Inclusion
matrices compose by multiplication: N c M c P has matrix Lambda_{PM} Lambda_{MN} (`Inclusion.compose`);
a tower of inclusions is drawn as its Bratteli diagram (`bratteli_diagram`).

Traces and conditional expectations
-----------------------------------
A faithful trace on M is tau = sum_d w_d Tr_d with w_d > 0 (w_d is the trace of a minimal projection of
the summand d); the library's traces 'blunt', 'norm', 'vn' (tau_vN) have w_d = 1, 1/m_d, 1/(D m_d). Its
restriction to N is the trace with weights nu = Lambda^T w (`Inclusion.induced_trace`):
tau(iota(x)) = sum_d w_d sum_c lambda_dc Tr x_c. The map E: M -> N,

    E(y)_c = nu_c^{-1} sum_{d,i} w_d V_{dc,i}^* u_d^* y_d u_d V_{dc,i},                       (E)

is the unique tau-preserving conditional expectation onto N. Proof (u = 1; conjugate by u otherwise).
iota(n)_d V_{dc,i} = V_{dc,i} n_c gives the bimodule property E(iota(n1) y iota(n2)) = n1 E(y) n2; with
V^* V = 1 it gives E(iota(n)) = nu_c^{-1} sum_d w_d lambda_dc n_c = n, so E o iota = id and E(1) = 1.
Summing Tr(n_c V^* y_d V) over (d, c, i) with weights w_d gives tau|_N(n E(y)) = tau(iota(n) y): E is
the adjoint of iota with respect to tau on M and the restricted trace tau|_N (weights nu, not the trace
of the same name on N) on N; n = 1 gives tau o E = tau. If F is another tau-preserving conditional
expectation, tau|_N(n F(y)) = tau(F(iota(n) y)) = tau(iota(n) y) = tau|_N(n E(y)) for all n, so F = E by
faithfulness of tau|_N. E is completely positive and unital: (E) is a Kraus form with non-negative
weights (`ConditionalExpectation.as_channel`).

Index
-----
* Markov trace and index ([Jones 1983], [GHJ, Ch. 2-3]). For a connected inclusion (the
  bipartite graph of Lambda is connected) let t > 0 be the Perron-Frobenius vector of Lambda Lambda^T,
  beta = ||Lambda||^2 its eigenvalue (largest singular value of Lambda squared), s = Lambda^T t, normalised
  by sum_d m_d t_d = 1. Then t is the weight vector of the unique Markov trace of M (s is its
  restriction to N) and Lambda s = beta t. The index [M : N] = beta: for factors N = M_k c M = M_m it is
  lambda^2 = (m/k)^2, for C c M_n it is n^2, for the diagonal masa C^n c M_n (Lambda = (1, ..., 1)) it is n.
  A disconnected inclusion is the direct sum of its connected components; `Inclusion.index` then
  returns max_j ||Lambda_j||^2 = ||Lambda||^2 with a warning, and `markov_trace` raises.
* Watatani index ([Watatani 1990]). A quasi-basis of E is a family {u_j} in M with
  y = sum_j u_j iota(E(u_j^* y)) for all y in M; Ind E = sum_j u_j u_j^* does not depend on it and is
  central. With matrix units of M a quasi-basis is u_{d,r,c,i} = (nu_c / w_d)^{1/2} u e^d_{r, off_dc + i} u^*
  (r < m_d, lambda_dc > 0, i < lambda_dc), which gives

      (Ind E)_d = sum_c lambda_dc nu_c / w_d,     i.e.  Ind E = (Lambda Lambda^T w)_d / w_d,

  equal to [M : N] 1 exactly for the Markov trace.
* Pimsner-Popa constant ([Pimsner-Popa 1986]). lambda(E) = max{l : E(y) >= l y for all y >= 0}. The
  minimum over the extreme rays y = xi xi^* (xi in C^{m_d}) gives, by the Schur-complement criterion
  A >= l xi xi^* <=> l xi^* A^+ xi <= 1 (xi in ran A), with A = iota(E(xi xi^*))_d,

      lambda(E)^{-1} = max_d sum_c (nu_c / w_d) min(k_c, lambda_dc),

  since xi^* A^+ xi = sum_c (nu_c / w_d) rank(sum_i xi_ci xi_ci^*) and the rank is at most min(k_c, lambda_dc),
  attained for generic xi. For C c M_n: lambda(E) = 1/n (not 1/n^2: for type I factors lambda(E)^{-1}
  differs from the index, unlike for II_1 subfactors, where lambda(E) = [M : N]^{-1}). The constant of
  the complete order, max{l : E - l id is completely positive}, is 1 / ||Ind E|| ([Frank-Kirchberg 1998]):
  1/n^2 for C c M_n and 1/m^2 for M_n (x) 1 c M_n (x) M_m. Both are in
  `ConditionalExpectation.pimsner_popa_constant`.

Basic construction and relative commutant
-----------------------------------------
The basic construction M_1 = <M, e_N> = (J N J)' on L^2(M, tau) has the sectors of N, sizes
(Lambda^T m)_c, and M c M_1 has inclusion matrix Lambda^T ([Jones 1983], [GHJ, Ch. 2]); iterating
gives the Jones tower N c M c M_1 c M_2 c ... with matrices Lambda, Lambda^T, Lambda, ... and constant
index. In the standard form of M c M_1 used here the Jones projection is
(e_N)_c = nu_c^{-1} sum_a eta_{c,a} eta_{c,a}^* with
eta_{c,a} = sum_{d,i} w_d^{1/2} e_{(d, off_dc + a lambda_dc + i, i)} (basis of summand d: (row of M_d, copy i))
(conjugated by the image of u), satisfying e_N y e_N = E(y) e_N; for the Markov trace the trace of M_1
with weights s / beta extends tau and tau_1(y e_N) = beta^{-1} tau(y) (Markov property).
The relative commutant is N' cap M = (+)_{(d,c): lambda_dc > 0} M_{lambda_dc}, block d acting as
1_{k_c} (x) y_dc; dim = sum lambda_dc^2.

Standard inclusions
-------------------
from_tensor_product (N (x) 1 c N (x) P, Lambda = diag blocks p_e), diagonal (C^n c M_n), scalars (C c M),
centre (Z(M) c M), and chain_tower: the U(1) spin-chain algebras A_L = (+)_{N=0}^{L} M_{binom(L, N)} of
operators on (C^2)^{(x) L} commuting with total S^z, with A_L c A_{L+1}, x -> x (x) 1_2 (the new site L is
the highest bit). Sector N of A_{L+1} receives N - 1 (new spin up) and N (new spin down) of A_L: Lambda
is the (L + 2) x (L + 1) bidiagonal 0/1 matrix, the incidence matrix of the path graph with 2L + 3
vertices, so [A_{L+1} : A_L] = ||Lambda||^2 = 4 cos^2(pi / (2L + 4)) (the norm of the adjacency matrix of
the path A_n is 2 cos(pi / (n + 1)), [GHJ, Ch. 1]). The basis order coincides with that of chains.SpinChain.

All maps act on batches (batch, C, k_max, k_max) by gathers with index maps built once per inclusion
(no loop over sectors, no host synchronisation), and are differentiable.

References
----------
[Jones 1983] V. F. R. Jones, Index for subfactors, Invent. Math. 72 (1983) 1-25.
[GHJ] F. M. Goodman, P. de la Harpe, V. F. R. Jones, Coxeter graphs and towers of algebras, MSRI Publ. 14,
    Springer (1989).
[Pimsner-Popa 1986] M. Pimsner, S. Popa, Entropy and index for subfactors, Ann. Sci. ENS 19 (1986) 57-106.
[Watatani 1990] Y. Watatani, Index for C*-subalgebras, Mem. Amer. Math. Soc. 83 (1990), no. 424.
[Frank-Kirchberg 1998] M. Frank, E. Kirchberg, On conditional expectations of finite index,
    J. Operator Theory 40 (1998) 87-111.
"""
import math
import warnings
from collections import namedtuple
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np
import torch

from .algebra import TypeIAlgebra, _same_device
from . import cost
from .channels import InterSectorChannel, _check_input, _same_structure
from .states import _trace_weights

Operator = TypeIAlgebra.Operator

__all__ = ['Inclusion', 'ConditionalExpectation', 'MarkovTrace', 'BratteliDiagram', 'bratteli_diagram',
           'from_tensor_product', 'diagonal', 'scalars', 'centre', 'center', 'chain_algebra', 'chain_tower']

MarkovTrace = namedtuple('MarkovTrace', ['t', 's', 'beta'])
MarkovTrace.__doc__ = """Markov trace of a connected inclusion: weights t (D,) on M, s = Lambda^T t (C,) on N
(the trace of a minimal projection of each summand; sum_d m_d t_d = 1), and beta = ||Lambda||^2 = [M : N]
with Lambda Lambda^T t = beta t and Lambda s = beta t."""

_TRACE_NAMES = {'blunt': 'blunt', 'norm': 'norm', 'vn': 'tau_vN', 'tau_vN': 'tau_vN'}

TraceSpec = Union[str, torch.Tensor, Sequence[float], np.ndarray]


def _new_algebra(like: TypeIAlgebra, sizes: Sequence[int], charges=None) -> TypeIAlgebra:
    """A TypeIAlgebra with blocks `sizes` (zero allowed) and the field, precision and device of `like`."""
    sizes = [int(s) for s in sizes]
    return TypeIAlgebra([max(s, 1) for s in sizes], sizes, complex_valued=like.complex_valued,
                        precision=like.precision, device=like.hilbert.device, charges=charges)


def _int_matrix(matrix, shape: Tuple[int, int]) -> np.ndarray:
    a = matrix.detach().cpu().numpy() if isinstance(matrix, torch.Tensor) else np.asarray(matrix)
    if a.ndim != 2 or a.shape != shape:
        raise ValueError(f"the inclusion matrix must have shape (D, C) = {shape} (rows: sectors of M, columns: "
                         f"sectors of N), got {a.shape}")
    if a.size and (not np.all(np.isfinite(a.astype(float))) or np.any(np.mod(a, 1) != 0) or np.any(a < 0)):
        raise ValueError("the inclusion matrix must have non-negative integer entries")
    return a.astype(np.int64)


def _components(lam: np.ndarray, rows: np.ndarray, cols: np.ndarray) -> List[Tuple[List[int], List[int]]]:
    """Connected components of the bipartite graph of lam on the given (non-empty) rows and columns."""
    D = lam.shape[0]
    parent = list(range(D + lam.shape[1]))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for d in rows:
        for c in cols:
            if lam[d, c]:
                parent[find(d)] = find(D + c)
    comps = {}
    for d in rows:
        comps.setdefault(find(d), ([], []))[0].append(int(d))
    for c in cols:
        comps.setdefault(find(D + c), ([], []))[1].append(int(c))
    return sorted(comps.values(), key=lambda rc: (rc[0] or [D], rc[1]))


def _random_operator(alg: TypeIAlgebra, batch: int) -> Operator:
    return alg.operator(torch.randn(batch, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype, device=alg.hilbert.device))


class Inclusion:
    """
    Unital inclusion iota = Ad(u) o iota_0: N -> M of finite-dimensional von Neumann algebras with inclusion
    matrix Lambda (see the module docstring).

    N, M: TypeIAlgebra with the same field, precision and device; the block sizes used are the active
        sizes k_factors (empty sectors, k_c = 0, are not part of the algebra and need zero columns / rows).
    matrix: Lambda, non-negative integers of shape (M.C, N.C) with m_d = sum_c lambda_dc k_c (unital) and a
        non-zero column for every non-empty sector of N (injective).
    unitary: optional Operator of M, batch 1 or B (one inclusion per sample): iota(x) = u iota_0(x) u^*.
    validate: check that `unitary` is unitary (one host synchronisation, at construction).
    """

    def __init__(self, N: TypeIAlgebra, M: TypeIAlgebra, matrix, unitary: Optional[Operator] = None,
                 validate: bool = True):
        if not isinstance(N, TypeIAlgebra) or not isinstance(M, TypeIAlgebra):
            raise TypeError("N and M must be TypeIAlgebra instances")
        if N.complex_valued != M.complex_valued:
            raise ValueError("N and M must both be real or both be complex")
        if N.dtype != M.dtype:
            raise ValueError(f"N and M must have the same precision, got {N.dtype} and {M.dtype}")
        if not _same_device(N.hilbert.device, M.hilbert.device):
            raise ValueError("N and M must be on the same device")
        lam = _int_matrix(matrix, (M.C, N.C))
        k = np.asarray(N.k_factors, dtype=np.int64)
        m = np.asarray(M.k_factors, dtype=np.int64)
        if np.any(lam[:, k == 0]):
            raise ValueError("empty sectors of N (k_c = 0) must have zero columns in the inclusion matrix")
        if not np.array_equal(lam @ k, m):
            raise ValueError(f"not a unital inclusion: the block sizes of M are {m.tolist()}, but Lambda k = "
                             f"{(lam @ k).tolist()} (need m_d = sum_c lambda_dc k_c)")
        dead = [c for c in range(N.C) if k[c] and not lam[:, c].any()]
        if dead:
            raise ValueError(f"not injective: the sectors {dead} of N are mapped to zero (zero columns of Lambda)")
        if unitary is not None:
            if not isinstance(unitary, Operator):
                raise TypeError("unitary must be an Operator of M")
            if not _same_structure(unitary.algebra, M):
                raise ValueError("unitary must be an Operator of M (same sectors, field, precision, device)")
            if validate:
                U = unitary.matrix
                err = (U @ U.mH - M._const('eye', U.dtype)).abs().amax()
                tol = 1e-4 if U.dtype in (torch.float32, torch.complex64) else 1e-10
                if float(err) > tol:
                    raise ValueError(f"unitary is not unitary (|u u^* - 1| = {float(err):.2e})")
        self.N, self.M = N, M
        self._lam = lam
        self.unitary = unitary
        self._cache = {}
        self._basic = None

    # ------------------------------------------------------------------
    @property
    def matrix(self) -> np.ndarray:
        """The inclusion matrix Lambda, an int64 array of shape (D, C) (a copy)."""
        return self._lam.copy()

    @property
    def batch_size(self) -> int:
        """Number of inclusions (batch of the unitary; 1 without unitary)."""
        return 1 if self.unitary is None else self.unitary.shape[0]

    def __repr__(self) -> str:
        u = '' if self.unitary is None else f", unitary batch {self.batch_size}"
        return (f"Inclusion(N={self.N.k_factors} -> M={self.M.k_factors}, Lambda={self._lam.tolist()}{u})")

    def _k(self) -> np.ndarray:
        return np.asarray(self.N.k_factors, dtype=np.int64)

    def _m(self) -> np.ndarray:
        return np.asarray(self.M.k_factors, dtype=np.int64)

    def _offsets(self) -> np.ndarray:
        """off[d, c] = sum_{c' < c} lambda_dc' k_c': first row of the summand c in block d."""
        lam, k = self._lam, self._k()
        return np.cumsum(lam * k[None, :], axis=1) - lam * k[None, :]

    # ------------------------------------------------------------------
    # index maps (built once)
    # ------------------------------------------------------------------
    def _embed_maps(self):
        """For every row p of block d of M: the sector c and row a of N it comes from and the id of its copy
        (-1 for padding); the mask (D, m_max, m_max) of entries inside one copy."""
        maps = self._cache.get('embed')
        if maps is not None:
            return maps
        lam, k = self._lam, self._k()
        D, C = lam.shape
        mx = self.M.k_max
        rc = np.zeros((D, mx), dtype=np.int64)
        ra = np.zeros((D, mx), dtype=np.int64)
        cp = np.full((D, mx), -1, dtype=np.int64)
        off = self._offsets()
        cid = 0
        for d, c in zip(*np.nonzero(lam)):
            n = lam[d, c] * k[c]
            r = np.arange(n)
            rc[d, off[d, c] + r] = c
            ra[d, off[d, c] + r] = r // lam[d, c]
            cp[d, off[d, c] + r] = cid + r % lam[d, c]
            cid += lam[d, c]
        dev = self.M.hilbert.device
        cp_t = torch.as_tensor(cp, device=dev)
        mask = (cp_t[:, :, None] == cp_t[:, None, :]) & (cp_t >= 0)[:, :, None]
        maps = (torch.as_tensor(rc, device=dev), torch.as_tensor(ra, device=dev), mask)
        self._cache['embed'] = maps
        return maps

    def _expect_maps(self):
        """For every sector c of N its copies (d, i), padded to R = max_c sum_d lambda_dc: the sector d (C, R),
        a validity mask (C, R) and the rows off_dc + a lambda_dc + i of the copy in block d (C, R, k_max)."""
        maps = self._cache.get('expect')
        if maps is not None:
            return maps
        lam, k = self._lam, self._k()
        C = lam.shape[1]
        kx = self.N.k_max
        R = max(1, int(lam.sum(axis=0).max()))
        ed = np.zeros((C, R), dtype=np.int64)
        valid = np.zeros((C, R), dtype=bool)
        pos = np.zeros((C, R, kx), dtype=np.int64)
        off = self._offsets()
        for c in range(C):
            j = 0
            for d in np.nonzero(lam[:, c])[0]:
                for i in range(lam[d, c]):
                    ed[c, j] = d
                    valid[c, j] = True
                    pos[c, j, :k[c]] = off[d, c] + np.arange(k[c]) * lam[d, c] + i
                    j += 1
        dev = self.M.hilbert.device
        maps = (torch.as_tensor(ed, device=dev), torch.as_tensor(valid, device=dev), torch.as_tensor(pos, device=dev))
        self._cache['expect'] = maps
        return maps

    def _isometries0(self) -> torch.Tensor:
        """The isometries V_{dc,i} of iota_0 as a Kraus tensor (1, r, D, C, m_max, k_max), r = max lambda_dc."""
        V = self._cache.get('isometries')
        if V is not None:
            return V
        lam, k = self._lam, self._k()
        D, C = lam.shape
        r = max(1, int(lam.max()))
        V = torch.zeros(1, r, D, C, self.M.k_max, self.N.k_max, dtype=self.M.dtype, device=self.M.hilbert.device)
        off = self._offsets()
        ii, dd, cc, pp, aa = [], [], [], [], []
        for d, c in zip(*np.nonzero(lam)):
            a, i = np.meshgrid(np.arange(k[c]), np.arange(lam[d, c]), indexing='ij')
            ii.append(i.ravel()); aa.append(a.ravel())
            dd.append(np.full(a.size, d)); cc.append(np.full(a.size, c))
            pp.append((off[d, c] + a * lam[d, c] + i).ravel())
        idx = [torch.as_tensor(np.concatenate(t), device=V.device) for t in (ii, dd, cc, pp, aa)]
        V[(torch.zeros_like(idx[0]),) + tuple(idx)] = 1
        self._cache['isometries'] = V
        return V

    # ------------------------------------------------------------------
    # the embedding
    # ------------------------------------------------------------------
    def _embed0(self, X: torch.Tensor) -> torch.Tensor:
        """iota_0 on matrices (B, C, k_max, k_max) -> (B, D, m_max, m_max): one gather and one mask."""
        rc, ra, mask = self._embed_maps()
        cost.check_memory(2 * cost.tensor_bytes((X.shape[0], self.M.C, self.M.k_max, self.M.k_max), X.dtype),
                          X.device, f"Inclusion.embed (batch {X.shape[0]}, {self.M.C} sectors up to {self.M.k_max})")
        return X[:, rc[:, :, None], ra[:, :, None], ra[:, None, :]] * mask

    def embed(self, x: Operator) -> Operator:
        """
        iota(x) for an Operator x of N with batch B: an Operator of M with batch max(B, batch of u) (lazy,
        differentiable; iota_0(x)_d = (+)_c x_c (x) 1_{lambda_dc}, then u . u^*). iota is a unital injective
        *-homomorphism, so self-adjointness, normality, positivity, invertibility and projections carry over.
        """
        _check_input(x, self.N, "Inclusion.embed")
        U = self.unitary

        def generator():
            Y = self._embed0(x.matrix)
            if U is not None:
                u = U.matrix
                Y = u @ Y @ u.mH
            return Y
        out = Operator(self.M, generator=generator)
        for flag in ('self_adjoint', 'normal', 'positive', 'invertible', 'projection'):
            setattr(out, f'_is_{flag}', getattr(x, f'_is_{flag}'))
        if x._batch is not None:
            out._batch = max(x._batch, self.batch_size)
        return out

    __call__ = embed

    def isometries(self) -> torch.Tensor:
        """Kraus tensor of iota, shape (batch, r, D, C, m_max, k_max): u_d V_{dc,i} (zero beyond lambda_dc)."""
        V = self._isometries0()
        if self.unitary is None:
            return V
        u = self.unitary.matrix
        return u[:, None, :, None] @ V

    def as_channel(self) -> InterSectorChannel:
        """
        iota as an InterSectorChannel N -> M (Kraus operators u_d V_{dc,i}), for the channel machinery:
        `adjoint('blunt')` of it is the restriction of states (densities with respect to the blunt traces),
        `compose` gives composite maps. iota is unital (as_channel().is_unital()), not trace preserving.
        """
        return InterSectorChannel(self.N, self.M, self.isometries())

    # ------------------------------------------------------------------
    # traces and conditional expectations
    # ------------------------------------------------------------------
    def trace_weights(self, trace: TraceSpec = 'blunt') -> torch.Tensor:
        """
        Weights w (batch, D) of a faithful trace tau = sum_d w_d Tr_d on M: 'blunt' (1), 'norm' (1/m_d),
        'vn' or 'tau_vN' (1/(D m_d)), 'markov' (the Markov trace, connected inclusions only), or explicit
        weights of shape (D,) or (batch, D), positive on the non-empty sectors (checked: one synchronisation).
        """
        M = self.M
        rdt = M.hilbert.real_dtype
        if isinstance(trace, str):
            if trace == 'markov':
                return self.markov_trace().t[None]
            if trace not in _TRACE_NAMES:
                raise ValueError(f"trace must be one of {tuple(_TRACE_NAMES) + ('markov',)} or a weight vector")
            return _trace_weights(M, _TRACE_NAMES[trace])[None]
        w = torch.as_tensor(trace)
        if torch.is_complex(w):
            raise ValueError("trace weights must be real")
        w = w.to(dtype=rdt, device=M.hilbert.device)
        if w.dim() == 1:
            w = w[None]
        if w.dim() != 2 or w.shape[1] != M.C:
            raise ValueError(f"trace weights must have shape ({M.C},) or (batch, {M.C}), got {tuple(w.shape)}")
        act = [d for d, m_d in enumerate(M.k_factors) if m_d]
        if not bool((w.detach()[:, act] > 0).all()):
            raise ValueError("trace weights must be positive on the non-empty sectors of M (a faithful trace)")
        return w

    def induced_trace(self, trace: TraceSpec = 'blunt') -> torch.Tensor:
        """
        Weights nu = Lambda^T w of the restriction tau|_N of tau = sum_d w_d Tr_d: tau(iota(x)) = sum_c nu_c Tr x_c.
        Shape (C,) for a named trace or weights (D,), (batch, C) for weights (batch, D).
        """
        w = self.trace_weights(trace)
        nu = w @ torch.as_tensor(self._lam, dtype=w.dtype, device=w.device)
        unbatched = isinstance(trace, str) or torch.as_tensor(trace).dim() == 1
        return nu[0] if unbatched else nu

    def conditional_expectation(self, trace: TraceSpec = 'blunt') -> 'ConditionalExpectation':
        """
        The unique tau-preserving conditional expectation E: M -> N for the trace `trace` of M (see
        trace_weights; a batch of weights gives a batch of expectations). See ConditionalExpectation.
        """
        return ConditionalExpectation(self, self.trace_weights(trace))

    # ------------------------------------------------------------------
    # index
    # ------------------------------------------------------------------
    def _active(self) -> Tuple[np.ndarray, np.ndarray]:
        return np.nonzero(self._m())[0], np.nonzero(self._k())[0]

    def connected_components(self) -> List[Tuple[List[int], List[int]]]:
        """Connected components of the bipartite graph of Lambda (non-empty sectors only), as pairs
        (sectors of M, sectors of N). The inclusion is the direct sum of the inclusions of its components."""
        rows, cols = self._active()
        return _components(self._lam, rows, cols)

    @property
    def is_connected(self) -> bool:
        """True if the Bratteli diagram of N c M is connected (equivalently Z(M) cap Z(N) = C 1)."""
        return len(self.connected_components()) == 1

    def component_indices(self) -> List[float]:
        """||Lambda_j||^2 of every connected component (float64)."""
        return [float(np.linalg.norm(self._lam[np.ix_(r, c)].astype(np.float64), 2) ** 2)
                for r, c in self.connected_components()]

    def index(self) -> float:
        """
        [M : N] = ||Lambda||^2 (largest singular value of Lambda squared, float64): the Jones index for factors
        ((m/k)^2) and the index of the Markov trace for connected inclusions (Lambda s = beta t). For a
        disconnected inclusion the Markov trace does not exist or is not unique; the value returned is
        max_j ||Lambda_j||^2 over the components, with a warning (component_indices gives all of them).
        """
        if not self.is_connected:
            warnings.warn(f"the inclusion is not connected (components {self.connected_components()}, indices "
                          f"{self.component_indices()}): it is a direct sum of connected inclusions and has no "
                          f"unique Markov trace; returning the largest component index ||Lambda||^2",
                          UserWarning, stacklevel=2)
        return float(np.linalg.norm(self._lam.astype(np.float64), 2) ** 2)

    def markov_trace(self) -> MarkovTrace:
        """
        The Markov trace of a connected inclusion (MarkovTrace(t, s, beta)): t (D,) is the Perron-Frobenius
        vector of Lambda Lambda^T (computed in float64 by a symmetric eigendecomposition), normalised to
        sum_d m_d t_d = 1, s = Lambda^T t (C,), beta = ||Lambda||^2; empty sectors get weight 0. Raises
        ValueError for a disconnected inclusion.
        """
        mt = self._cache.get('markov')
        if mt is not None:
            return mt
        comps = self.connected_components()
        if len(comps) != 1:
            raise ValueError(f"no unique Markov trace: the inclusion is not connected (components {comps} with "
                             f"indices {self.component_indices()}); it is the direct sum of these inclusions, "
                             f"each with its own Markov trace; a Markov trace of the sum exists only if all "
                             f"indices agree, and then is not unique")
        rows, cols = self._active()
        A = self._lam[np.ix_(rows, cols)].astype(np.float64)
        evals, evecs = np.linalg.eigh(A @ A.T)
        beta = float(evals[-1])
        t = np.abs(evecs[:, -1])
        t = t / (self._m()[rows].astype(np.float64) @ t)
        if not np.all(t > 0):                              # Perron-Frobenius for an irreducible matrix
            raise RuntimeError("Perron-Frobenius vector not positive (numerical failure)")
        t_full = np.zeros(self.M.C)
        t_full[rows] = t
        s_full = self._lam.T.astype(np.float64) @ t_full
        rdt, dev = self.M.hilbert.real_dtype, self.M.hilbert.device
        mt = MarkovTrace(torch.as_tensor(t_full, dtype=rdt, device=dev), torch.as_tensor(s_full, dtype=rdt, device=dev),
                         beta)
        self._cache['markov'] = mt
        return mt

    def watatani_index(self, E: Union['ConditionalExpectation', TraceSpec] = 'markov',
                       method: str = 'formula') -> Operator:
        """Ind E for a ConditionalExpectation or a trace of M; see ConditionalExpectation.index."""
        if not isinstance(E, ConditionalExpectation):
            E = self.conditional_expectation(E)
        return E.index(method)

    def pimsner_popa_constant(self, E: Union['ConditionalExpectation', TraceSpec] = 'markov',
                              complete: bool = False) -> torch.Tensor:
        """lambda(E) for a ConditionalExpectation or a trace of M; see ConditionalExpectation.pimsner_popa_constant."""
        if not isinstance(E, ConditionalExpectation):
            E = self.conditional_expectation(E)
        return E.pimsner_popa_constant(complete)

    # ------------------------------------------------------------------
    # relative commutant, composition, basic construction
    # ------------------------------------------------------------------
    def _permutation(self, alg: TypeIAlgebra, sec: np.ndarray, actual: np.ndarray, canon: np.ndarray):
        """The permutation unitary P of `alg` with P e_canon = e_actual in sector sec, or None if P = 1."""
        if np.array_equal(actual, canon):
            return None
        P = torch.zeros(1, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype, device=alg.hilbert.device)
        idx = [torch.as_tensor(a, device=P.device) for a in (sec, actual, canon)]
        P[(torch.zeros_like(idx[0]),) + tuple(idx)] = 1
        op = alg.operator(P)
        op._is_normal = op._is_invertible = True
        return op

    def relative_commutant(self) -> 'Inclusion':
        """
        N' cap M as the inclusion R c M with R = (+)_{(d,c): lambda_dc > 0} M_{lambda_dc} (sectors ordered by d,
        then c; charges (d, c)): sector (d, c) acts in block d as 1_{k_c} (x) y on the copies of C^{k_c},
        conjugated by u. Its inclusion matrix has entry k_c in row d and column (d, c); dim R = sum lambda_dc^2.
        `commutator_error` checks [iota(N), R] = 0 numerically.
        """
        lam, k = self._lam, self._k()
        pairs = [(int(d), int(c)) for d, c in zip(*np.nonzero(lam))]       # row-major: ordered by d, then c
        R = _new_algebra(self.M, [lam[d, c] for d, c in pairs], charges=pairs)
        lam_R = np.zeros((self.M.C, len(pairs)), dtype=np.int64)
        off = self._offsets()
        sec, actual, canon = [], [], []
        for j, (d, c) in enumerate(pairs):
            lam_R[d, j] = k[c]
            a, i = np.meshgrid(np.arange(k[c]), np.arange(lam[d, c]), indexing='ij')
            sec.append(np.full(a.size, d))
            actual.append((off[d, c] + a * lam[d, c] + i).ravel())
            canon.append((off[d, c] + i * k[c] + a).ravel())
        P = self._permutation(self.M, *(np.concatenate(t) for t in (sec, actual, canon)))
        U = P if self.unitary is None else (self.unitary if P is None else self.unitary @ P)
        return Inclusion(R, self.M, lam_R, unitary=U, validate=False)

    def commutator_error(self, other: 'Inclusion', batch_size: int = 4) -> torch.Tensor:
        """
        ||[iota(x), iota'(y)]||_F / (||iota(x)||_F ||iota'(y)||_F) for random (Gaussian) x in N, y in other.N, for two
        inclusions into the same M (e.g. other = self.relative_commutant()). Shape (batch,): batch_size samples,
        or one per inclusion when the unitaries are batched.
        """
        if not _same_structure(self.M, other.M):
            raise ValueError("commutator_error: both inclusions must have the same algebra M")
        B = max(self.batch_size, other.batch_size)
        B = batch_size if B == 1 else B
        X = self.embed(_random_operator(self.N, B)).matrix
        Y = other.embed(_random_operator(other.N, B)).matrix
        comm = X @ Y - Y @ X
        nrm = lambda Z: torch.linalg.vector_norm(Z, dim=(-3, -2, -1))
        return nrm(comm) / (nrm(X) * nrm(Y))

    def compose(self, other: 'Inclusion') -> 'Inclusion':
        """
        The inclusion N c P of a tower N c M c P (self: N c M, other: M c P), iota_P o iota_M, with matrix
        Lambda_PM Lambda_MN and unitary iota_P(u_M) u_P Pi, where the permutation Pi orders the copies of
        C^{k_c} in a block e of P by (d, i, j) (sector d of M, copy i in M_d, copy j of M_d in P_e).
        """
        if not _same_structure(self.M, other.N):
            raise ValueError("compose: self.M must be other.N (self: N c M, other: M c P)")
        lam, mu, k = self._lam, other._lam, self._k()
        lam2 = mu @ lam
        offM, offP = self._offsets(), other._offsets()
        off2 = np.cumsum(lam2 * k[None, :], axis=1) - lam2 * k[None, :]
        sec, actual, canon = [], [], []
        for e in range(mu.shape[0]):
            for c in range(lam.shape[1]):
                base = 0                                      # copies of c in P_e before sector d
                for d in range(lam.shape[0]):
                    if mu[e, d] and lam[d, c] and k[c]:
                        a, i, j = np.meshgrid(np.arange(k[c]), np.arange(lam[d, c]), np.arange(mu[e, d]),
                                              indexing='ij')
                        sec.append(np.full(a.size, e))
                        actual.append((offP[e, d] + (offM[d, c] + a * lam[d, c] + i) * mu[e, d] + j).ravel())
                        canon.append((off2[e, c] + a * lam2[e, c] + base + i * mu[e, d] + j).ravel())
                    base += lam[d, c] * mu[e, d]
        cat = (lambda ts: np.concatenate(ts)) if sec else (lambda ts: np.zeros(0, dtype=np.int64))
        P = self._permutation(other.M, cat(sec), cat(actual), cat(canon))
        U = None
        for f in ((other.embed(self.unitary) if self.unitary is not None else None), other.unitary, P):
            if f is not None:
                U = f if U is None else U @ f
        return Inclusion(self.N, other.M, lam2, unitary=U, validate=False)

    def basic_construction(self) -> 'Inclusion':
        """
        The inclusion M c M_1 of the basic construction M_1 = <M, e_N>: M_1 has one sector per sector c of N
        (charges of N) of size (Lambda^T m)_c, and the inclusion matrix is Lambda^T (in standard form; the
        unitary u of N c M enters only the Jones projection). Cached; iterate for the Jones tower.
        """
        if self._basic is None:
            lamT = np.ascontiguousarray(self._lam.T)
            M1 = _new_algebra(self.M, lamT @ self._m(), charges=list(self.N.charges))
            self._basic = Inclusion(self.M, M1, lamT, validate=False)
        return self._basic

    def jones_tower(self, n: int) -> List['Inclusion']:
        """[N c M, M c M_1, ..., M_{n-2} c M_{n-1}] (n inclusions), matrices Lambda, Lambda^T, Lambda, ..."""
        tower = [self]
        while len(tower) < n:
            tower.append(tower[-1].basic_construction())
        return tower

    def jones_projection(self, trace: TraceSpec = 'markov') -> Operator:
        """
        The Jones projection e_N in M_1 = basic_construction().M for the trace `trace` of M: the projection
        of L^2(M, tau) onto L^2(N, tau). (e_N)_c = nu_c^{-1} sum_a eta_{c,a} eta_{c,a}^* in the standard form (see
        the module docstring), conjugated by j(u). It satisfies e_N j(y) e_N = j(iota(E(y))) e_N for y in M
        (j: M -> M_1), commutes with j(iota(N)), and M_1 = span j(M) e_N j(M). Shape (batch, C, n1_max, n1_max).
        """
        basic = self.basic_construction()
        M1 = basic.M
        w = self.trace_weights(trace)
        nu = w @ torch.as_tensor(self._lam, dtype=w.dtype, device=w.device)
        key = 'jones'
        idx = self._cache.get(key)
        if idx is None:
            lam, k = self._lam, self._k()
            off, off1 = self._offsets(), basic._offsets()          # off1[c, d]: summand d in sector c of M_1
            cc, rr, ss, dx, dy = [], [], [], [], []
            for c in range(lam.shape[1]):
                pos, ds, As = [], [], []
                for d in np.nonzero(lam[:, c])[0]:
                    a, i = np.meshgrid(np.arange(k[c]), np.arange(lam[d, c]), indexing='ij')
                    pos.append((off1[c, d] + (off[d, c] + a * lam[d, c] + i) * lam[d, c] + i).ravel())
                    ds.append(np.full(a.size, d)); As.append(a.ravel())
                if not pos:
                    continue
                pos, ds, As = np.concatenate(pos), np.concatenate(ds), np.concatenate(As)
                x, y = np.nonzero(As[:, None] == As[None, :])          # pairs of entries with the same a
                cc.append(np.full(x.size, c)); rr.append(pos[x]); ss.append(pos[y]); dx.append(ds[x]); dy.append(ds[y])
            idx = tuple(torch.as_tensor(np.concatenate(t), device=M1.hilbert.device) for t in (cc, rr, ss, dx, dy))
            self._cache[key] = idx
        cc, rr, ss, dx, dy = idx
        nu_safe = torch.where(nu > 0, nu, torch.ones_like(nu))
        vals = torch.sqrt(w[:, dx] * w[:, dy]) / nu_safe[:, cc]
        e = torch.zeros(w.shape[0], M1.C, M1.k_max, M1.k_max, dtype=M1.dtype, device=M1.hilbert.device)
        e[:, cc, rr, ss] = vals.to(e.dtype)
        if self.unitary is not None:
            ju = basic.embed(self.unitary).matrix
            e = ju @ e @ ju.mH
        op = M1.operator(e)
        op._is_self_adjoint = op._is_positive = op._is_normal = op._is_projection = True
        return op


class ConditionalExpectation:
    """
    The tau-preserving conditional expectation E: M -> N of an inclusion iota: N c M,

        E(y)_c = nu_c^{-1} sum_{d,i} w_d V_{dc,i}^* u_d^* y_d u_d V_{dc,i},     nu = Lambda^T w,

    for the faithful trace tau = sum_d w_d Tr_d of M (weights (batch, D)). E is unital, completely positive,
    N-bimodular, E o iota = id, tau o E = tau, and tau(iota(x) y) = tau|_N(x E(y)) with tau|_N = sum_c nu_c Tr_c
    (the restricted trace, not a trace of the same name on N). Create with Inclusion.conditional_expectation.
    """

    def __init__(self, inclusion: Inclusion, weights: torch.Tensor):
        self.inclusion = inclusion
        self.weights = weights                                          # (Bw, D)
        lam = torch.as_tensor(inclusion._lam, dtype=weights.dtype, device=weights.device)
        self.nu = weights @ lam                                         # (Bw, C): the trace tau|_N on N

    def __repr__(self) -> str:
        return f"ConditionalExpectation({self.inclusion!r}, weights batch {self.weights.shape[0]})"

    def _nu_safe(self) -> torch.Tensor:
        return torch.where(self.nu > 0, self.nu, torch.ones_like(self.nu))

    def _scale(self) -> torch.Tensor:
        """w_d / nu_c for every copy (d, i) of every sector c, zero for padding: (Bw, C, R)."""
        ed, valid, _ = self.inclusion._expect_maps()
        s = self.weights[:, ed] / self._nu_safe()[:, :, None]
        return s * valid

    def _apply0(self, Y: torch.Tensor) -> torch.Tensor:
        """E for u = 1 on matrices (B, D, m_max, m_max) -> (B, C, k_max, k_max)."""
        inc = self.inclusion
        ed, _, pos = inc._expect_maps()
        C, R, kx = pos.shape
        B = max(Y.shape[0], self.weights.shape[0])
        cost.check_memory(2 * cost.tensor_bytes((B, C, R, kx, kx), Y.dtype), Y.device,
                          f"ConditionalExpectation (batch {B}, {C} sectors, {R} copies of size up to {kx})")
        G = Y[:, ed[:, :, None, None], pos[:, :, :, None], pos[:, :, None, :]]        # (B, C, R, k, k)
        out = (G * self._scale()[:, :, :, None, None].to(G.real.dtype)).sum(2)
        return out * inc.N._const('block_mask', torch.bool)

    def __call__(self, y: Operator) -> Operator:
        """E(y) for an Operator y of M (batch B): an Operator of N with batch max(B, batch of u, batch of w);
        lazy and differentiable. Positive (self-adjoint) y give positive (self-adjoint) E(y)."""
        inc = self.inclusion
        _check_input(y, inc.M, "ConditionalExpectation")
        U = inc.unitary

        def generator():
            Y = y.matrix
            if U is not None:
                u = U.matrix
                Y = u.mH @ Y @ u
            return self._apply0(Y)
        out = Operator(inc.N, generator=generator)
        out._is_self_adjoint = True if y._is_self_adjoint else None
        out._is_positive = True if y._is_positive else None
        if y._batch is not None:
            out._batch = max(y._batch, inc.batch_size, self.weights.shape[0])
        return out

    apply = __call__

    def as_channel(self) -> InterSectorChannel:
        """E as an InterSectorChannel M -> N with Kraus operators (w_d / nu_c)^{1/2} V_{dc,i}^* u_d^*
        (shape (batch, r, C, D, k_max, m_max)); unital, not trace preserving in general."""
        K = self.inclusion.isometries()                                        # (Bu, r, D, C, m, k)
        Kd = K.conj().transpose(-2, -1).transpose(2, 3)                        # (Bu, r, C, D, k, m)
        scale = (self.weights[:, None, :] / self._nu_safe()[:, :, None]).sqrt()   # (Bw, C, D)
        Kd = Kd * scale[:, None, :, :, None, None].to(Kd.dtype)
        return InterSectorChannel(self.inclusion.M, self.inclusion.N, Kd)

    # ------------------------------------------------------------------
    def index(self, method: str = 'formula') -> Operator:
        """
        Watatani index Ind E = sum_j u_j u_j^* (central, positive, invertible), an Operator of M with batch
        max(batch of w, batch of u). method='formula': (Ind E)_d = sum_c lambda_dc nu_c / w_d (exact);
        method='quasi_basis': the sum over the quasi-basis of `quasi_basis` (O(J D m_max^2) memory).
        For the Markov trace Ind E = [M : N] 1.
        """
        inc = self.inclusion
        if method == 'quasi_basis':
            Q = self.quasi_basis()
            mat = (Q @ Q.mH).sum(1)
            op = inc.M.operator(mat)
        elif method == 'formula':
            lam = torch.as_tensor(inc._lam, dtype=self.weights.dtype, device=self.weights.device)
            w = self.weights
            act = inc.M._const('k', w.dtype) > 0
            val = torch.where(act, (self.nu @ lam.T) / torch.where(act, w, torch.ones_like(w)), torch.zeros_like(w))
            B = max(w.shape[0], inc.batch_size)
            eye = inc.M._const('eye')
            mat = (val[:, :, None, None].to(eye.dtype) * eye).expand(B, -1, -1, -1)
            op = inc.M.operator(mat)
        else:
            raise ValueError("method must be 'formula' or 'quasi_basis'")
        op._is_self_adjoint = op._is_normal = op._is_positive = True
        op._is_invertible = True
        return op

    def quasi_basis(self) -> torch.Tensor:
        """
        A quasi-basis (Pimsner-Popa basis) {u_j} of E, y = sum_j u_j iota(E(u_j^* y)) for all y in M:
        u_{d,r,c,i} = (nu_c / w_d)^{1/2} u e^d_{r, off_dc + i} u^* (r < m_d, lambda_dc > 0, i < lambda_dc), as a tensor
        (batch, J, D, m_max, m_max) with J = sum_d m_d sum_c lambda_dc.
        """
        inc = self.inclusion
        lam, m = inc._lam, inc._m()
        off = inc._offsets()
        dd, rr, cols, cc = [], [], [], []
        for d, c in zip(*np.nonzero(lam)):
            r, i = np.meshgrid(np.arange(m[d]), np.arange(lam[d, c]), indexing='ij')
            dd.append(np.full(r.size, d)); cc.append(np.full(r.size, c))
            rr.append(r.ravel()); cols.append((off[d, c] + i).ravel())
        dev = inc.M.hilbert.device
        dd, rr, cols, cc = (torch.as_tensor(np.concatenate(t), device=dev) for t in (dd, rr, cols, cc))
        J = dd.numel()
        w = self.weights
        vals = (self.nu[:, cc] / w[:, dd]).sqrt()                              # (Bw, J)
        Q = torch.zeros(w.shape[0], J, inc.M.C, inc.M.k_max, inc.M.k_max, dtype=inc.M.dtype, device=dev)
        Q[:, torch.arange(J, device=dev), dd, rr, cols] = vals.to(Q.dtype)
        if inc.unitary is not None:
            u = inc.unitary.matrix[:, None]
            Q = u @ Q @ u.mH
        return Q

    def pimsner_popa_constant(self, complete: bool = False) -> torch.Tensor:
        """
        Pimsner-Popa constant lambda(E) = max{l : E(y) >= l y for all y >= 0} = 1 / max_d sum_c (nu_c / w_d)
        min(k_c, lambda_dc) (see the module docstring); with complete=True the best constant of the complete
        order, max{l : E - l id completely positive} = 1 / ||Ind E|| = 1 / max_d sum_c lambda_dc nu_c / w_d.
        Shape (batch,) (batch of the weights). Neither depends on the unitary u.
        """
        inc = self.inclusion
        lam, k = inc._lam, inc._k()
        coef = lam if complete else np.minimum(lam, k[None, :])
        coef = torch.as_tensor(coef, dtype=self.weights.dtype, device=self.weights.device)
        w = self.weights
        act = inc.M._const('k', w.dtype) > 0
        per_block = torch.where(act, (self.nu @ coef.T) / torch.where(act, w, torch.ones_like(w)), torch.zeros_like(w))
        return 1.0 / per_block.amax(dim=-1)


# ----------------------------------------------------------------------
# Bratteli diagrams
# ----------------------------------------------------------------------
class BratteliDiagram:
    """
    Bratteli diagram of a tower A_0 c A_1 c ... c A_n: `levels[l]` are the block sizes of A_l, `charges[l]`
    their labels, `matrices[l]` the inclusion matrix of A_l c A_{l+1} (shape (len(levels[l+1]), len(levels[l]))).
    """

    def __init__(self, levels: List[List[int]], matrices: List[np.ndarray], charges: Optional[List[list]] = None):
        self.levels = [list(map(int, lv)) for lv in levels]
        self.matrices = [np.asarray(a, dtype=np.int64) for a in matrices]
        self.charges = charges if charges is not None else [list(range(len(lv))) for lv in levels]
        if len(self.matrices) != len(self.levels) - 1:
            raise ValueError("need one inclusion matrix between consecutive levels")
        for l, a in enumerate(self.matrices):
            if a.shape != (len(self.levels[l + 1]), len(self.levels[l])):
                raise ValueError(f"inclusion matrix {l} has shape {a.shape}")

    def __repr__(self) -> str:
        return f"BratteliDiagram(levels={self.levels})"

    def edges(self, level: int) -> List[Tuple[int, int, int]]:
        """Edges between level and level + 1 as (c, d, multiplicity), c in A_level, d in A_{level+1}."""
        a = self.matrices[level]
        return [(int(c), int(d), int(a[d, c])) for d, c in zip(*np.nonzero(a))]

    def draw(self, ax=None, show_sizes: bool = True):
        """Draw the diagram with matplotlib (imported on use): level l on the line y = -l, multiplicities
        larger than 1 as edge labels. Returns the axes."""
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(1 + 0.8 * max(len(lv) for lv in self.levels), 1 + 0.9 * len(self.levels)))
        xs = [np.arange(len(lv)) - (len(lv) - 1) / 2 for lv in self.levels]
        for l in range(len(self.matrices)):
            for c, d, mult in self.edges(l):
                x0, x1 = xs[l][c], xs[l + 1][d]
                ax.plot([x0, x1], [-l, -l - 1], color='0.4', lw=0.8 + 0.6 * (mult - 1), zorder=1)
                if mult > 1:
                    ax.text((x0 + x1) / 2, -l - 0.5, str(mult), fontsize=8, ha='center', va='center',
                            bbox=dict(fc='white', ec='none', pad=0.5))
        for l, lv in enumerate(self.levels):
            ax.scatter(xs[l], [-l] * len(lv), s=160, color='white', edgecolors='0.1', zorder=2)
            if show_sizes:
                for x, n in zip(xs[l], lv):
                    ax.text(x, -l, str(n), fontsize=7, ha='center', va='center', zorder=3)
        ax.set_yticks([-l for l in range(len(self.levels))])
        ax.set_yticklabels([f"A_{l}" for l in range(len(self.levels))])
        ax.set_xticks([])
        for s in ('top', 'right', 'bottom'):
            ax.spines[s].set_visible(False)
        return ax


def bratteli_diagram(tower: Union[Inclusion, Sequence[Inclusion]]) -> BratteliDiagram:
    """Bratteli diagram of a tower of inclusions [A_0 c A_1, A_1 c A_2, ...] (consecutive algebras must agree)."""
    tower = [tower] if isinstance(tower, Inclusion) else list(tower)
    if not tower:
        raise ValueError("empty tower")
    for a, b in zip(tower, tower[1:]):
        if not _same_structure(a.M, b.N):
            raise ValueError("not a tower: the larger algebra of an inclusion must be the smaller of the next")
    algs = [tower[0].N] + [inc.M for inc in tower]
    return BratteliDiagram([alg.k_factors for alg in algs], [inc.matrix for inc in tower],
                           [list(alg.charges) for alg in algs])


# ----------------------------------------------------------------------
# standard inclusions
# ----------------------------------------------------------------------
def from_tensor_product(N: TypeIAlgebra, second: Union[int, TypeIAlgebra]) -> Inclusion:
    """
    N (x) 1 c N (x) P with P = M_m (second = m) or any TypeIAlgebra P: M = composite.tensor_product(N, P)
    (sector (c, e) at c * P.C + e), inclusion matrix lambda_{(c,e), c'} = delta_cc' p_e; iota(x) = kron(x, 1),
    and composite.partial_trace / composite.partial_trace_channel apply to M.
    """
    from .composite import tensor_product
    P = second if isinstance(second, TypeIAlgebra) else _new_algebra(N, [int(second)])
    if P.dtype != N.dtype or P.complex_valued != N.complex_valued:
        raise ValueError("both factors must have the same field and precision")
    M = tensor_product(N, P)
    lam = np.zeros((M.C, N.C), dtype=np.int64)
    for c, k_c in enumerate(N.k_factors):
        for e, p_e in enumerate(P.k_factors):
            if k_c:
                lam[c * P.C + e, c] = p_e
    return Inclusion(N, M, lam)


def diagonal(n: int, complex_valued: bool = True, precision: str = 'double', device=None) -> Inclusion:
    """The diagonal masa C^n c M_n (Lambda = (1, ..., 1), shape (1, n)); index n."""
    N = TypeIAlgebra([1] * n, [1] * n, complex_valued=complex_valued, precision=precision, device=device)
    M = _new_algebra(N, [n])
    return Inclusion(N, M, np.ones((1, n), dtype=np.int64))


def scalars(M: TypeIAlgebra) -> Inclusion:
    """C 1 c M (Lambda = (m_d)_d, one column); index sum_d m_d^2 for the Markov trace (connected)."""
    N = _new_algebra(M, [1])
    return Inclusion(N, M, np.asarray(M.k_factors, dtype=np.int64)[:, None])


def centre(M: TypeIAlgebra) -> Inclusion:
    """The centre Z(M) = C^D c M (one sector of size 1 per non-empty sector of M, Lambda = diag(m_d)).
    Disconnected for more than one sector: the components are C c M_{m_d}."""
    act = [d for d, m_d in enumerate(M.k_factors) if m_d]
    N = _new_algebra(M, [1] * len(act), charges=[M.charges[d] for d in act])
    lam = np.zeros((M.C, len(act)), dtype=np.int64)
    for j, d in enumerate(act):
        lam[d, j] = M.k_factors[d]
    return Inclusion(N, M, lam)


center = centre


def chain_algebra(L: int, complex_valued: bool = True, precision: str = 'double', device=None) -> TypeIAlgebra:
    """A_L = (+)_{N=0}^{L} M_{binom(L, N)}: operators on L spins 1/2 commuting with total S^z (charges N); the
    structure of chains.SpinChain(L).algebra (A_0 = C)."""
    dims = [math.comb(L, n) for n in range(L + 1)]
    return TypeIAlgebra(dims, dims, complex_valued=complex_valued, precision=precision, device=device,
                        charges=list(range(L + 1)))


def chain_tower(L_max: int, L_min: int = 0, complex_valued: bool = True, precision: str = 'double',
                device=None) -> List[Inclusion]:
    """
    The U(1) spin-chain tower [A_{L_min} c A_{L_min+1}, ..., A_{L_max-1} c A_{L_max}], A_L c A_{L+1} by
    x -> x (x) 1_2 (a new site L): Lambda[N, N] = 1 (new spin down, N <= L), Lambda[N, N-1] = 1 (new spin up),
    with the basis order of chains.SpinChain (states with the new spin up first). Consecutive inclusions share
    their algebras; [A_{L+1} : A_L] = 4 cos^2(pi / (2L + 4)).
    """
    if not 0 <= L_min < L_max:
        raise ValueError("need 0 <= L_min < L_max")
    algs = [chain_algebra(L, complex_valued, precision, device) for L in range(L_min, L_max + 1)]
    tower = []
    for j, L in enumerate(range(L_min, L_max)):
        lam = np.zeros((L + 2, L + 1), dtype=np.int64)
        lam[np.arange(L + 1), np.arange(L + 1)] = 1
        lam[np.arange(1, L + 2), np.arange(L + 1)] = 1
        tower.append(Inclusion(algs[j], algs[j + 1], lam))
    return tower
