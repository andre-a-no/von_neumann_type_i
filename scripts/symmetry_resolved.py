#!/usr/bin/env python3
"""
Section "Use case: symmetry-resolved entanglement in the XXZ chain" of the paper.

The ground state |psi> of the XXZ chain H = sum_<ij> [S^x_i S^x_j + S^y_i S^y_j + Delta S^z_i S^z_j] at zero
magnetisation (N = L/2 up spins) is a vector in one sector of (+)_N M_binom(L,N). For a block A of ell sites the
reduced state rho_A commutes with the charge Q_A of A (number of up spins in A), i.e. it is an element of

    M_A = (+)_{q=0}^{ell} M_binom(ell, q),

the algebra of SpinChain(ell). The fused tensor product M_A (x)_fused M_B (composite.fused_tensor_product) has the
sector N = L/2 made of pair blocks (q, N - q) of size binom(ell, q) x binom(L - ell, N - q); in this layout the vector
|psi> is a list of matrices Psi_q and the partial trace is rho_A = (+)_q Psi_q Psi_q^*: the sector form of rho_A
without ever forming |psi><psi| (dimension binom(24, 12)^2 = 7.3e12 at L = 24). From the spectrum of each block
(the library's sector-wise eigenvalues) we obtain, for every q,

    p(q) = Tr Pi_q rho_A,   S(q) = S(rho_A|_q / p(q)),   Z_n(q) = Tr Pi_q rho_A^n,   p_n(q) = Z_n(q) / sum_q Z_n(q),

the decomposition S = sum_q p(q) S(q) + H(p) (configurational + number entanglement), the variances b_n of the
charge distributions p_n (b_1 is the full-counting-statistics variance of Q_A), and the widths B_n of the
integer-sampled Gaussian through the two central sectors (q0 = ceil(ell/2), x = q - ell/2, x0 = q0 - ell/2):
B_n = (x1^2 - x0^2) / (2 log(p_n(q0)/p_n(q0+1))), with G_n = (n B_n - B_1)/(n - 1), G_1 = d/dn (n B_n) at n = 1.
Exact identities (no approximation):

    S_n(q) - S_n = [log p_n(q) - n log p(q)] / (1 - n),     S(q) - S = log p(q) - d/dn log p_n(q) |_{n=1}.

If p_n(q) is proportional to exp(-x^2 / (2 B_n)) on the integers, they give S_n(q) - S_n(q0) = -G_n (x^2 - x0^2) /
(2 B_1 B_n) (exact at q0 +- 1 by the definition of B_n, tested at q0 +- 2). Bosonisation: B_n = K/(pi^2 n) log ell_c +
gamma_n + o(1) (Goldstein-Sela), so n B_n has an n-independent slope, G_n = O(1) and equipartition is broken at order
(q - <q>)^2 / log^2 ell. K = pi / (2 (pi - arccos Delta)) is the Luttinger parameter (K = 1 free fermions, K = 1/2 at
the isotropic point), ell_c = (L/pi) sin(pi ell/L) (periodic); for an open chain with A at the boundary the slope is
halved and ell_c = (2L/pi) sin(pi ell/L).

Parts:
  main        L = L_main, periodic and open, a batch of Delta (one batched Lanczos run per chunk: Delta only enters
              the diagonal, so the Hamiltonians share the sparse hopping part), all blocks ell <= L/2; K from b_1
              (periodic) and from n B_n, n = 1, 2, 3, against the Bethe ansatz; equipartition and its breaking.
  validation  (a) XX chain (Delta = 0): Z_n(alpha) = prod_j [nu_j^n e^{i alpha} + (1 - nu_j)^n] over the eigenvalues
              nu_j of the restricted correlation matrix (flux insertion), Z_n(q) by an exact discrete Fourier
              transform on ell + 1 points, S(q) from d/dn Z_n(q); against the library at a small L (dense sector)
              and at L_main (the Lanczos ground state of the main part), periodic and open;
              (b) full Hilbert space: dense 2^L Hamiltonian in NumPy, |psi> reshaped, Pi_q rho_A Pi_q of the
              2^ell x 2^ell matrix, for a batch of Delta, against the library;
              (c) three library routes to rho_A at small L: the pure-state blocks, SpinChain.reduced_state of
              |psi><psi| and composite.partial_trace in the fused product;
              (d) QuSpin (external): ground state and reduced density matrix, projected on the charges.
  scan        half chains ell = L/2 of periodic chains for a range of L.
  ff_inf      infinite XX chain (ell up to a few thousand): b_1 - (log 2 ell + 1 + gamma_E)/pi^2, G_1 (observed to tend
              to 1/(2 pi^2)) and S(q0) - S(q0+1) against (pi^2/4)/(log 2 ell + 1 + gamma_E)^2.
  ff_ring     periodic XX ring, ell = L/2, L up to a few thousand: B_1, G_1 and the Gaussian test at q0 +- 2.

Writes generated/symmetry_resolved.{tex,json} and figures/symmetry_resolved.pdf.

    python scripts/symmetry_resolved.py --mode full --device cpu        # paper values, about 15-30 min on a 4-core CPU

On a shared machine use OMP_NUM_THREADS=2 (oversubscribed threads slow the sparse products down drastically).
"""
import math
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT), str(ROOT / 'scripts')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from torch_vn_algebra import SpinChain, DensityMatrix, composite, krylov, cost   # noqa: E402
from torch_vn_algebra.krylov import SparseSectorHamiltonian                         # noqa: E402

F64 = torch.float64
TOL = 1e-12                      # Lanczos residual tolerance (relative to max(1, |E0|))
P_MIN = 1e-10                    # S(q), S_n(q) are reported for sectors with p(q) > P_MIN


def luttinger_K(delta):
    """Bethe ansatz: K = pi / (2 (pi - arccos Delta)) for -1 < Delta <= 1 (K = 1 at Delta = 0, 1/2 at Delta = 1)."""
    return math.pi / (2 * (math.pi - math.acos(delta)))


def log_chord(L, ell, boundary):
    """log ell_c: ell_c = (L/pi) sin(pi ell/L) (periodic), (2L/pi) sin(pi ell/L) (open, block at the boundary)."""
    f = 1.0 if boundary == 'periodic' else 2.0
    return math.log(f * L / math.pi * math.sin(math.pi * ell / L))


# ----------------------------------------------------------------------------------------------------------------
# ground states: one batched Lanczos run for a batch of Delta
# ----------------------------------------------------------------------------------------------------------------
def xxz_batch(chain, N, deltas):
    """Sector-N XXZ Hamiltonians (J = 1) for a batch of Delta as one krylov.SparseSectorHamiltonian: Delta only
    multiplies the diagonal S^z S^z part, the sparse hopping part is shared."""
    H1 = chain.xxz_sparse(J=1.0, Delta=1.0, sector=N)
    D = torch.as_tensor(deltas, dtype=chain.real_dtype, device=chain.device)
    return SparseSectorHamiltonian(D[:, None] * H1.diag[0][None], H1.offdiag, H1.dtype)


def ground_states(chain, deltas, m, chunk, seed, log=print):
    """Ground states of sector N = L/2 for every Delta: E0 (B,), psi (B, binom(L, L/2)), residuals (B,)."""
    N = chain.L // 2
    Es, psis, ress = [], [], []
    for i in range(0, len(deltas), chunk):
        part = deltas[i:i + chunk]
        t0 = time.time()
        H = xxz_batch(chain, N, part)
        E, psi, res = krylov.ground_state(H, m=min(m, H.dim), tol=TOL, max_restarts=60,
                                          generator=torch.Generator().manual_seed(seed + i))
        Es.append(E)
        psis.append(psi)
        ress.append(res)
        log(f"  Lanczos L={chain.L} {chain.boundary} Delta={part}: dim {H.dim}, {time.time() - t0:.1f} s, "
            f"max residual {res.max().item():.1e}")
        del H
    return torch.cat(Es), torch.cat(psis), torch.cat(ress)


# ----------------------------------------------------------------------------------------------------------------
# the reduced state in sector form, through the fused tensor product
# ----------------------------------------------------------------------------------------------------------------
class Bipartition:
    """
    A = sites 0..ell-1, B = the rest, for vectors of sector N of `chain`. The fused product
    M_A (x)_fused M_B = composite.fused_tensor_product(SpinChain(ell).algebra, SpinChain(L - ell).algebra) fixes the
    layout: its sector N is the list of pair blocks (q, N - q, offset); `perm` maps the configuration basis of the
    chain's sector N to that layout, and `gather[q]` (binom(ell, q) x binom(L - ell, N - q)) holds the positions of
    the amplitudes Psi_q in the chain's basis.
    """

    def __init__(self, chain, N, ell, right_cache=None):
        L, dev = chain.L, chain.device
        self.chain, self.N, self.ell = chain, N, ell
        kw = dict(complex_valued=chain.complex_valued, precision=chain.precision, device=dev)
        self.left = SpinChain(ell, 'open', **kw)
        nB = L - ell
        if right_cache is not None and nB in right_cache:
            right = right_cache[nB]
        else:
            right = SpinChain(nB, 'open', **kw)
            if right_cache is not None:
                right_cache[nB] = right
        with warnings.catch_warnings():               # the product algebra is never stored densely
            warnings.simplefilter('ignore', cost.CostWarning)
            self.alg12 = composite.fused_tensor_product(self.left.algebra, right.algebra)
        self.sector = self.alg12.charges.index(N)
        lab = chain.labels[N]
        a, b = lab & ((1 << ell) - 1), lab >> ell
        qa = torch.zeros_like(a)
        for i in range(ell):
            qa += (a >> i) & 1
        pos = torch.empty_like(lab)
        self.blocks = []
        for q, d, off in self.alg12.product_blocks[self.sector]:
            kq, md = self.left.algebra.k_factors[q], right.algebra.k_factors[d]
            sel = qa == q
            pos[sel] = off + self.left.index_of[a[sel]] * md + right.index_of[b[sel]]
            self.blocks.append((q, off, kq, md))
        self.perm = pos                               # chain index -> fused index
        inv = torch.empty_like(pos)
        inv[pos] = torch.arange(len(pos), device=dev)
        self.gather = {q: inv[off:off + kq * md].reshape(kq, md) for q, off, kq, md in self.blocks}

    def fused_vector(self, psi):
        """psi (B, binom(L, N)) in the basis of sector N of the fused product algebra."""
        v = torch.empty_like(psi)
        v[:, self.perm] = psi
        return v

    def reduced_state(self, psi):
        """rho_A = (+)_q Psi_q Psi_q^* as a DensityMatrix of SpinChain(ell).algebra (batched over psi)."""
        alg = self.left.algebra
        mat = torch.zeros(psi.shape[0], alg.C, alg.k_max, alg.k_max, dtype=psi.dtype, device=psi.device)
        for q, G in self.gather.items():
            P = psi[:, G]                                                  # (B, k_q, m_d)
            mat[:, q, :P.shape[1], :P.shape[1]] = P @ P.mH
        return DensityMatrix(alg, matrix=mat, validate=False)


# ----------------------------------------------------------------------------------------------------------------
# sector-resolved quantities from the block spectra
# ----------------------------------------------------------------------------------------------------------------
def _xlogx(x):
    return torch.special.xlogy(x, x)


def central_width(logZ, ell):
    """
    Width of the integer-sampled Gaussian through the two central sectors: if p_n(q) is proportional to
    exp(-(q - ell/2)^2 / (2 B)) on the integers, then B = (x1^2 - x0^2) / (2 log(p_n(q0) / p_n(q1))) with q0 = ceil(ell/2),
    q1 = q0 + 1, x_i = q_i - ell/2. logZ: (..., ell + 1) logarithms of Z_n(q) (any normalisation). NaN for ell < 2.
    """
    if ell < 2:
        return torch.full(logZ.shape[:-1], math.nan, dtype=logZ.dtype, device=logZ.device)
    q0 = (ell + 1) // 2
    x0 = q0 - ell / 2
    return ((x0 + 1) ** 2 - x0 ** 2) / (2 * (logZ[..., q0] - logZ[..., q0 + 1]))


def sector_statistics(rho, ns=(2, 3)):
    """
    From the eigenvalues of every block of rho_A (rho.eigenvalues(), one tensor per charge q), batched over the
    leading dimension:
      p(q), S(q), S_n(q), S, S_n, H(p), Z_n(q) = Tr Pi_q rho_A^n, the fourth cumulant k4 of p;
      b_n = Var_{p_n}(Q_A) (b_1: full counting statistics);
      B_n = central_width of p_n, G_n = (n B_n - B_1)/(n - 1) and G_1 = d/dn (n B_n) at n = 1, which by the identity
      d/dn log Z_n(q) |_{n=1} = log p(q) - S(q) equals B_1 (S(q0) - S(q1)) / log(p(q0)/p(q1)).
    S(q), S_n(q) are NaN where p(q) <= P_MIN. Assumes <Q_A> = ell/2 (half filling) for B_n.
    """
    ev = [w.clamp(min=0) for w in rho.eigenvalues()]
    dev = ev[0].device
    ell = rho.algebra.C - 1
    q = torch.tensor(rho.algebra.charges, dtype=F64, device=dev)
    p = torch.stack([w.sum(-1) for w in ev], -1)                         # (B, C)
    T = torch.stack([_xlogx(w).sum(-1) for w in ev], -1)                 # sum over the block of lambda log lambda
    S = -T.sum(-1)
    H = -_xlogx(p).sum(-1)
    ok = p > P_MIN
    ps = torch.where(ok, p, torch.ones_like(p))
    Sq = torch.where(ok, -T / ps + torch.log(ps), torch.full_like(p, math.nan))
    mean = (p * q).sum(-1)
    dq2 = (q[None] - mean[:, None]) ** 2
    b1 = (p * dq2).sum(-1)
    k4 = (p * dq2 ** 2).sum(-1) - 3 * b1 ** 2
    logp = torch.log(ps)
    B1 = central_width(logp, ell)
    if ell >= 2:
        q0 = (ell + 1) // 2
        G1 = B1 * (Sq[:, q0] - Sq[:, q0 + 1]) / (logp[:, q0] - logp[:, q0 + 1])
    else:
        G1 = torch.full_like(B1, math.nan)
    out = dict(q=q, p=p, T=T, S=S, H=H, Sq=Sq, mean=mean, b={1: b1}, B={1: B1}, G={1: G1}, k4=k4,
               Sn={}, Snq={}, Z={})
    for n in ns:
        Zq = torch.stack([(w ** n).sum(-1) for w in ev], -1)
        Zn = Zq.sum(-1)
        pn = Zq / Zn[:, None]
        mn = (pn * q).sum(-1)
        out['b'][n] = (pn * (q[None] - mn[:, None]) ** 2).sum(-1)
        Zs = torch.where(ok, Zq, torch.ones_like(Zq))
        out['B'][n] = Bn = central_width(torch.log(Zs), ell)
        out['G'][n] = (n * Bn - B1) / (n - 1)
        out['Sn'][n] = torch.log(Zn) / (1 - n)
        out['Snq'][n] = torch.where(ok, (torch.log(Zs) - n * logp) / (1 - n), torch.full_like(p, math.nan))
        out['Z'][n] = Zq
    return out


def gaussian_breaking(st, n=1):
    """
    S_n(q) - S_n(q0) for p_n(q) proportional to exp(-(q - ell/2)^2 / (2 B_n)) on the integers (all n):
    -G_n ((q - ell/2)^2 - x0^2) / (2 B_1 B_n); exact at q0 +- 1 by the definition of B_n, a prediction elsewhere.
    """
    ell = len(st['q']) - 1
    x0 = (ell + 1) // 2 - ell / 2
    x2 = (st['q'] - ell / 2) ** 2 - x0 ** 2
    return -st['G'][n][:, None] * x2[None] / (2 * st['B'][1][:, None] * st['B'][n][:, None])


# ----------------------------------------------------------------------------------------------------------------
# free fermions (Delta = 0): correlation matrix and flux insertion
# ----------------------------------------------------------------------------------------------------------------
def ff_correlation(L, boundary, N):
    """
    Jordan-Wigner: S^+_i S^-_j -> c_i^* c_j on bonds (i, i+1); on the periodic bond (L-1, 0) the string gives the factor
    (-1)^(N-1) in the sector with N particles. Returns C_ij = <c_i^* c_j> of the N-particle ground state and E0.
    """
    h = np.zeros((L, L))
    for i in range(L - 1):
        h[i, i + 1] = h[i + 1, i] = 0.5
    if boundary == 'periodic' and L > 2:
        h[L - 1, 0] = h[0, L - 1] = 0.5 * (-1) ** (N - 1)
    e, U = np.linalg.eigh(h)
    if N < L and e[N] - e[N - 1] < 1e-9:
        raise ValueError("degenerate free-fermion ground state")
    occ = U[:, :N]
    return occ @ occ.T, float(e[:N].sum())


def ff_flux(C, ell, ns=(2, 3), theta=0.25):
    """
    Charged moments Z_n(alpha) = Tr rho_A^n e^{i alpha Q_A} = prod_j [nu_j^n e^{i alpha} + (1 - nu_j)^n] (nu_j the
    eigenvalues of C restricted to A) on the ell + 1 points alpha_k = 2 pi (k + theta)/(ell + 1); Z_n(q) is the exact
    discrete Fourier transform (Z_n(alpha) is a polynomial of degree ell in e^{i alpha}). For S(q):
    T(q) = Tr Pi_q rho_A log rho_A = FT of dZ_n/dn at n = 1 = Z_1(alpha) sum_j [nu_j e^{i alpha} log nu_j +
    (1 - nu_j) log(1 - nu_j)] / [nu_j e^{i alpha} + 1 - nu_j]. Returns p(q), T(q), S(q), Z_n(q), S_n(q), nu.
    """
    nu = np.clip(np.linalg.eigvalsh(C[:ell, :ell]), 0.0, 1.0)
    k = np.arange(ell + 1)
    alpha = 2 * np.pi * (k + theta) / (ell + 1)
    z = np.exp(1j * alpha)[:, None]
    q = np.arange(ell + 1)

    def ft(Za):                                 # Z(q) = (1/(ell+1)) sum_k e^{-i alpha_k q} Z(alpha_k)
        return (np.exp(-1j * np.outer(q, alpha)) @ Za / (ell + 1)).real

    with np.errstate(divide='ignore', invalid='ignore'):
        lnu = np.where(nu > 0, np.log(np.where(nu > 0, nu, 1)), 0.0)
        l1nu = np.where(nu < 1, np.log(np.where(nu < 1, 1 - nu, 1)), 0.0)
    fac = nu[None] * z + (1 - nu)[None]
    Z1a = np.prod(fac, axis=1)
    dZa = Z1a * np.sum((nu * lnu)[None] * z / fac + ((1 - nu) * l1nu)[None] / fac, axis=1)
    p, T = ft(Z1a), ft(dZa)
    ok = p > P_MIN
    ps = np.where(ok, p, 1.0)
    out = dict(p=p, T=T, Sq=np.where(ok, -T / ps + np.log(ps), np.nan), Z={}, Snq={}, nu=nu)
    for n in ns:
        Zq = ft(np.prod(nu[None] ** n * z + (1 - nu)[None] ** n, axis=1))
        out['Z'][n] = Zq
        with np.errstate(divide='ignore'):                 # Z_n(q) may round to 0 where p(q) is tiny
            out['Snq'][n] = np.where(ok, (np.log(np.where(ok, np.abs(Zq), 1.0)) - n * np.log(ps)) / (1 - n), np.nan)
    return out


def ff_variances(nu, n):
    """Free fermions: Q_A under rho_A^n / Tr is a sum of independent Bernoulli modes with
    sigma_n = nu^n / (nu^n + (1 - nu)^n); b_n = sum sigma_n (1 - sigma_n), and G_1 = d/dn (n b_n)|_1 in closed form."""
    nu = np.clip(nu, 1e-300, 1 - 1e-16)
    eps = np.log((1 - nu) / nu)
    f = 0.25 / np.cosh(np.clip(n * eps, -700, 700) / 2) ** 2
    bn = f.sum()
    f1 = 0.25 / np.cosh(np.clip(eps, -700, 700) / 2) ** 2
    G1 = (f1 - eps * f1 * np.tanh(np.clip(eps, -700, 700) / 2)).sum()
    return float(bn), float(G1)


def ff_central(fl, ell):
    """B_1 and G_1 (as in sector_statistics) from the flux-insertion p(q) and S(q)."""
    q0 = (ell + 1) // 2
    x0 = q0 - ell / 2
    R = math.log(fl['p'][q0] / fl['p'][q0 + 1])
    B1 = ((x0 + 1) ** 2 - x0 ** 2) / (2 * R)
    return B1, B1 * (fl['Sq'][q0] - fl['Sq'][q0 + 1]) / R


# ----------------------------------------------------------------------------------------------------------------
# independent computations for the validation
# ----------------------------------------------------------------------------------------------------------------
def dense_full_space(L, boundary, delta):
    """Ground state of the XXZ chain in the full 2^L space (NumPy, bit labels as in the library: bit i = site i)."""
    dim = 2 ** L
    x = np.arange(dim)
    bonds = [(i, i + 1) for i in range(L - 1)] + ([(L - 1, 0)] if boundary == 'periodic' and L > 2 else [])
    H = np.zeros((dim, dim))
    for i, j in bonds:
        si, sj = ((x >> i) & 1) - 0.5, ((x >> j) & 1) - 0.5
        H[x, x] += delta * si * sj
        flip = ((x >> i) & 1) != ((x >> j) & 1)
        H[x[flip] ^ ((1 << i) | (1 << j)), x[flip]] += 0.5
    w, V = np.linalg.eigh(H)
    if w[1] - w[0] < 1e-8:
        raise ValueError("degenerate ground state in the full space")
    return V[:, 0], w[0]


def full_space_sectors(psi_full, L, ell):
    """p(q), S(q) from the 2^ell x 2^ell matrix rho_A and the projections Pi_q (no sector structure is used); also the
    norm of the part of rho_A between different charges (zero for U(1)-symmetric states)."""
    M = psi_full.reshape(2 ** (L - ell), 2 ** ell)                     # rows: B configurations, columns: A
    rho = M.T @ M.conj()
    pc = np.array([bin(a).count('1') for a in range(2 ** ell)])
    off = rho * (pc[:, None] != pc[None, :])
    p, Sq = np.zeros(ell + 1), np.full(ell + 1, np.nan)
    for q in range(ell + 1):
        P = np.diag((pc == q).astype(float))
        w = np.linalg.eigvalsh(P @ rho @ P)
        p[q] = w.sum()
        if p[q] > P_MIN:
            lam = np.clip(w / p[q], 0, None)
            lam = lam[lam > 0]
            Sq[q] = -(lam * np.log(lam)).sum()
    return p, Sq, float(np.abs(off).max())


def _max_dev(a, b, p=None, pmin=0.0):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    sel = np.isfinite(a) & np.isfinite(b)
    if p is not None:
        sel &= np.asarray(p) > pmin
    return float(np.abs(a[sel] - b[sel]).max()) if sel.any() else 0.0


def compare_free_fermions(ch, psi, E0, how, p_cmp, log=print):
    """Library sector statistics of the XX ground state psi (1, dim) against the flux-insertion formulas, all blocks
    ell <= L/2; maximal deviations."""
    t0 = time.time()
    L, N = ch.L, ch.L // 2
    C, Eff = ff_correlation(L, ch.boundary, N)
    devs = dict(E0=abs(float(E0) - Eff), p=0.0, pS=0.0, Z=0.0, S=0.0, S2=0.0, S3=0.0, b=0.0, G=0.0)
    cache = {}
    for ell in range(1, L // 2 + 1):
        st = sector_statistics(Bipartition(ch, N, ell, cache).reduced_state(psi))
        fl = ff_flux(C, ell)
        p_lib = st['p'][0].cpu().numpy()
        devs['p'] = max(devs['p'], _max_dev(p_lib, fl['p']))
        devs['pS'] = max(devs['pS'], _max_dev(-st['T'][0].cpu().numpy(), -fl['T']))
        devs['S'] = max(devs['S'], _max_dev(st['Sq'][0].cpu().numpy(), fl['Sq'], p_lib, p_cmp))
        for n in (2, 3):
            Zl = st['Z'][n][0].cpu().numpy()
            devs['Z'] = max(devs['Z'], _max_dev(Zl, fl['Z'][n]))
            pn = Zl / Zl.sum()                    # S_n(q) where p(q) and p_n(q) exceed p_cmp
            devs[f'S{n}'] = max(devs[f'S{n}'], _max_dev(st['Snq'][n][0].cpu().numpy(), fl['Snq'][n],
                                                        np.minimum(p_lib, pn), p_cmp))
        b1 = ff_variances(fl['nu'], 1)[0]
        devs['b'] = max(devs['b'], abs(st['b'][1][0].item() - b1),
                        abs(st['b'][2][0].item() - ff_variances(fl['nu'], 2)[0]))
        if ell >= 2:
            B1, G1 = ff_central(fl, ell)
            devs['G'] = max(devs['G'], abs(st['B'][1][0].item() - B1), abs(st['G'][1][0].item() - G1))
    log(f"  XX vs free fermions L={L} {ch.boundary} ({how}): " + ', '.join(f"{k} {v:.1e}" for k, v in devs.items()))
    return dict(L=L, boundary=ch.boundary, how=how, dev=devs, seconds=time.time() - t0)


def part_validation(cfg, dev, seed, log=print):
    """XX chain (dense sector) against free fermions, the full Hilbert space, three library routes, QuSpin."""
    res = {}
    # (a) free fermions, dense sector at small L (the largest L is compared in run(), from the main ground states)
    ff = []
    L = cfg['L_ff_dense']
    for bc in ('periodic', 'open'):
        ch = SpinChain(L, bc, complex_valued=False, precision='double', device=dev)
        w, V = ch.xxz(J=1.0, Delta=0.0, sector=L // 2).eigh()[0]
        ff.append(compare_free_fermions(ch, V[:, :, 0], w[0, 0].item(), 'dense', cfg['p_cmp'], log))
    res['ff'] = ff

    # (b) full Hilbert space (NumPy, no sector structure) for a batch of Delta; (c) three library routes
    L = cfg['L_full']
    deltas = cfg['deltas_val']
    full = []
    for bc in ('periodic', 'open'):
        ch = SpinChain(L, bc, complex_valued=False, precision='double', device=dev)
        N = L // 2
        E0, psi, resid = ground_states(ch, deltas, cfg['m'], len(deltas), seed, log)
        dp, dS, offmax, dE, droutes = 0.0, 0.0, 0.0, 0.0, 0.0
        cache = {}
        for d_i, delta in enumerate(deltas):
            vfull, Efull = dense_full_space(L, bc, delta)
            dE = max(dE, abs(Efull - E0[d_i].item()))
            for ell in range(1, L):
                st = sector_statistics(Bipartition(ch, N, ell, cache).reduced_state(psi[d_i:d_i + 1]))
                p_f, S_f, off = full_space_sectors(vfull, L, ell)
                dp = max(dp, _max_dev(st['p'][0].cpu().numpy(), p_f))
                dS = max(dS, _max_dev(st['Sq'][0].cpu().numpy(), S_f, p_f, cfg['p_cmp']))
                offmax = max(offmax, off)
        # (c) routes to rho_A in the library, from |psi><psi| in sector N (dense, small L only)
        Hd = ch.xxz(J=1.0, Delta=deltas[-1], sector=N)
        w, V = Hd.eigh()[0]
        v = V[:, :, 0]
        rhoN = DensityMatrix(ch.sector_algebra(N), matrix=(v[:, :, None] * v[:, None, :].conj())[:, None], validate=False)
        for ell in range(1, L):
            bp = Bipartition(ch, N, ell, cache)
            r_pure = bp.reduced_state(v).matrix
            r_chain = ch.reduced_state(rhoN, ell, sector=N).matrix
            alg12 = bp.alg12
            mat = torch.zeros(1, alg12.C, alg12.k_max, alg12.k_max, dtype=v.dtype, device=v.device)
            fv = bp.fused_vector(v)
            k = fv.shape[1]
            mat[:, bp.sector, :k, :k] = fv[:, :, None] * fv[:, None, :].conj()
            r_comp = composite.partial_trace(DensityMatrix(alg12, matrix=mat, validate=False), keep=1).matrix
            droutes = max(droutes, (r_pure - r_chain).abs().max().item(), (r_pure - r_comp).abs().max().item())
        full.append(dict(L=L, boundary=bc, deltas=deltas, dE=dE, dp=dp, dS=dS, offdiag=offmax, routes=droutes))
        log(f"  full space L={L} {bc}, Delta={deltas}: E0 {dE:.1e}, p(q) {dp:.1e}, S(q) {dS:.1e}, "
            f"|rho_A between charges| {offmax:.1e}; library routes (pure / reduced_state / partial_trace) {droutes:.1e}")
    res['full'] = full

    # (d) QuSpin
    res['quspin'] = part_quspin(cfg, dev, seed, log)
    return res


def part_quspin(cfg, dev, seed, log=print):
    """QuSpin: ground state and reduced density matrix of A = sites 0..L/2-1, projected on the charges of A."""
    try:
        import quspin                                                          # noqa: F401
        from quspin.basis import spin_basis_1d
        from quspin.operators import hamiltonian
    except ImportError:
        log("  QuSpin not installed: skipped")
        return None
    L, delta = cfg['L_quspin'], cfg['delta_quspin']
    out = []
    for bc in ('periodic', 'open'):
        t0 = time.time()
        basis = spin_basis_1d(L, Nup=L // 2, pauli=0)
        bonds = [[1.0, i, i + 1] for i in range(L - 1)] + ([[1.0, L - 1, 0]] if bc == 'periodic' else [])
        hb = [[0.5, i, j] for _, i, j in bonds]
        zb = [[delta, i, j] for _, i, j in bonds]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            Hq = hamiltonian([["+-", hb], ["-+", hb], ["zz", zb]], [], basis=basis, dtype=np.float64,
                             check_symm=False, check_herm=False, check_pcon=False)
            Eq, Vq = Hq.eigsh(k=1, which='SA', tol=1e-14)
        ch = SpinChain(L, bc, complex_valued=False, precision='double', device=dev)
        E0, psi, _ = ground_states(ch, [delta], cfg['m'], 1, seed, log)
        ell = L // 2
        rdm = basis.ent_entropy(Vq[:, 0], sub_sys_A=list(range(ell)), return_rdm='A', density=False)['rdm_A']
        pc = np.array([bin(a).count('1') for a in range(2 ** ell)])
        p_q, S_q = np.zeros(ell + 1), np.full(ell + 1, np.nan)
        for q in range(ell + 1):
            blk = rdm[np.ix_(pc == q, pc == q)]
            lam = np.linalg.eigvalsh(blk)
            p_q[q] = lam.sum()
            lam = np.clip(lam / p_q[q], 0, None) if p_q[q] > P_MIN else np.zeros(0)
            lam = lam[lam > 0]
            S_q[q] = -(lam * np.log(lam)).sum() if p_q[q] > P_MIN else np.nan
        st = sector_statistics(Bipartition(ch, L // 2, ell).reduced_state(psi))
        # QuSpin may number the sites of A in the opposite order; p(q) and S(q) do not depend on it
        r = dict(L=L, boundary=bc, delta=delta, dE=abs(Eq[0] - E0[0].item()),
                 dp=_max_dev(st['p'][0].cpu().numpy(), p_q), dS=_max_dev(st['Sq'][0].cpu().numpy(), S_q, p_q, cfg['p_cmp']),
                 seconds=time.time() - t0)
        out.append(r)
        log(f"  QuSpin L={L} {bc} Delta={delta}: E0 {r['dE']:.1e}, p(q) {r['dp']:.1e}, S(q) {r['dS']:.1e}")
    return out


# ----------------------------------------------------------------------------------------------------------------
# main computation
# ----------------------------------------------------------------------------------------------------------------
def _np(t):
    return t.detach().cpu().numpy().tolist()


def chain_data(L, bc, deltas, cfg, dev, seed, log=print, ells=None):
    """Ground states for all Delta (batched) and the sector statistics of every block ell (batched over Delta)."""
    t0 = time.time()
    ch = SpinChain(L, bc, complex_valued=False, precision='double', device=dev)
    E0, psi, resid = ground_states(ch, deltas, cfg['m'], cfg['chunk'], seed, log)
    t_gs = time.time() - t0
    N = L // 2
    ells = list(range(1, L // 2 + 1)) if ells is None else ells
    cache = {}
    rows = []
    t1 = time.time()
    for ell in ells:
        st = sector_statistics(Bipartition(ch, N, ell, cache).reduced_state(psi))
        row = dict(ell=ell, logc=log_chord(L, ell, bc), S=_np(st['S']), H=_np(st['H']), mean=_np(st['mean']),
                   k4=_np(st['k4']), p=_np(st['p']), Sq=_np(st['Sq']),
                   b={n: _np(v) for n, v in st['b'].items()}, B={n: _np(v) for n, v in st['B'].items()},
                   G={n: _np(v) for n, v in st['G'].items()},
                   Sn={n: _np(v) for n, v in st['Sn'].items()}, Snq={n: _np(v) for n, v in st['Snq'].items()},
                   gauss1=_np(gaussian_breaking(st, 1)), gauss2=_np(gaussian_breaking(st, 2)))
        rows.append(row)
    t_ent = time.time() - t1
    keep = None
    if 0.0 in list(deltas):
        i0 = list(deltas).index(0.0)
        keep = (ch, psi[i0:i0 + 1].clone(), E0[i0].item())
    del psi
    log(f"  L={L} {bc}: ground states {t_gs:.0f} s, sector statistics for {len(ells)} blocks {t_ent:.1f} s")
    return dict(L=L, boundary=bc, deltas=list(deltas), K=[luttinger_K(d) for d in deltas], E0=_np(E0),
                residual=_np(resid), rows=rows, seconds_gs=t_gs, seconds_entanglement=t_ent, _xx=keep)




def fit_K_fcs(data, ell_min):
    """
    Periodic chains: least-squares fit of the full-counting-statistics variance
    b_1(ell) = (K/pi^2) log ell_c + gamma + C (-1)^ell ell_c^{-2K} over ell_min <= ell <= L/2 (oscillating term as in
    Song et al.), with the exponent fixed self-consistently by the fitted K. Returns K per Delta.
    """
    rows = [r for r in data['rows'] if r['ell'] >= ell_min]
    X = np.array([r['logc'] for r in rows])
    sgn = np.array([(-1) ** r['ell'] for r in rows], dtype=float)
    out = []
    for d in range(len(data['deltas'])):
        y = np.array([r['b'][1][d] for r in rows])
        K = np.polyfit(X, y, 1)[0] * math.pi ** 2
        for _ in range(500):                                            # fixed point for the exponent
            M = np.stack([X, np.ones_like(X), sgn * np.exp(-2 * K * X)], 1)
            K_new = np.linalg.lstsq(M, y, rcond=None)[0][0] * math.pi ** 2
            if abs(K_new - K) < 1e-12:
                break
            K = 0.5 * (K + K_new)
        out.append(float(K))
    return out


def fit_K_width(data, n, ell_min):
    """
    Least-squares fit n B_n(ell) = (K/(f pi^2)) log ell_c + gamma + C (-1)^ell over ell_min <= ell <= L/2 (f = 1
    periodic, f = 2 open: one entangling point); the constant C absorbs the even/odd effect. Returns K per Delta.
    """
    f = 1 if data['boundary'] == 'periodic' else 2
    rows = [r for r in data['rows'] if r['ell'] >= ell_min]
    X = np.array([r['logc'] for r in rows])
    sgn = np.array([(-1) ** r['ell'] for r in rows], dtype=float)
    M = np.stack([X, np.ones_like(X), sgn], 1)
    out = []
    for d in range(len(data['deltas'])):
        y = np.array([n * r['B'][n][d] for r in rows])
        out.append(float(np.linalg.lstsq(M, y, rcond=None)[0][0] * f * math.pi ** 2))
    return out


def part_scan(cfg, dev, seed, log=print):
    """Half chains ell = L/2 of periodic chains for a range of L."""
    out = []
    for L in cfg['L_scan']:
        d = chain_data(L, 'periodic', cfg['deltas'], cfg, dev, seed, log, ells=[L // 2])
        out.append(dict(L=L, row=d['rows'][0], E0=d['E0']))
    return out


def _inf_correlation(ell):
    i = np.arange(ell)
    d = i[:, None] - i[None, :]
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(d == 0, 0.5, np.sin(np.pi * d / 2) / (np.pi * np.where(d == 0, 1, d)))


def part_ff_inf(cfg, log=print):
    """
    Infinite XX chain, C_ij = sin(pi (i-j)/2)/(pi (i-j)): b_1 - (log 2 ell + 1 + gamma_E)/pi^2, G_1 = d/dn (n b_n)|_1
    (closed form over the entanglement spectrum), and S(q0) - S(q0+1) (flux insertion) against the Gaussian value
    G_1 / (2 b_1^2) -> (pi^2/4) / (log 2 ell + 1 + gamma_E)^2 for G_1 -> 1/(2 pi^2).
    """
    out = []
    for ell in cfg['ff_inf_ells']:
        C = _inf_correlation(ell)
        nu = np.linalg.eigvalsh(C)
        b1, G1 = ff_variances(nu, 1)
        b2, _ = ff_variances(nu, 2)
        X = math.log(2 * ell) + 1 + np.euler_gamma
        ref = X / math.pi ** 2
        r = dict(ell=ell, b1=b1, b1_dev=b1 - ref, G1=G1, G2=2 * b2 - b1, dS_asym=math.pi ** 2 / (4 * X ** 2))
        if ell % 2 == 0 and ell <= cfg['ff_inf_flux_max']:
            fl = ff_flux(C, ell, ns=())
            q0 = ell // 2
            r['dS1'] = float(fl['Sq'][q0] - fl['Sq'][q0 + 1])
            r['ratio'] = float((fl['Sq'][q0] - fl['Sq'][q0 + 2]) / r['dS1'])
        out.append(r)
        log(f"  infinite XX chain ell={ell}: b_1 - (log 2ell + 1 + gamma_E)/pi^2 = {b1 - ref:.2e}, G_1 = {G1:.9f}, "
            f"G_2 = {2 * b2 - b1:.6f}" + (f", S(q0)-S(q0+1) = {r['dS1']:.5f} vs (pi^2/4)/X^2 = {r['dS_asym']:.5f}, "
                                          f"ratio {r['ratio']:.3f}" if 'dS1' in r else ''))
    return out


def part_ff_ring(cfg, log=print):
    """Periodic XX ring, half chain ell = L/2, L up to a few thousand (flux insertion): B_1, G_1, S(q0) - S(q0+1) and
    the ratio [S(q0) - S(q0+2)] / [S(q0) - S(q0+1)] (4 for a Gaussian p_n)."""
    out = []
    for L in cfg['ff_ring_L']:
        C, _ = ff_correlation(L, 'periodic', L // 2)
        ell = L // 2
        fl = ff_flux(C, ell, ns=())
        B1, G1 = ff_central(fl, ell)
        q0 = (ell + 1) // 2
        d1, d2 = fl['Sq'][q0] - fl['Sq'][q0 + 1], fl['Sq'][q0] - fl['Sq'][q0 + 2]
        out.append(dict(L=L, B1=B1, G1=G1, dS1=float(d1), ratio=float(d2 / d1)))
        log(f"  XX ring L={L}, ell=L/2: B_1 = {B1:.5f}, G_1 = {G1:.5f}, S(q0)-S(q0+1) = {d1:.5f}, ratio {d2 / d1:.3f}")
    return out


CONFIGS = {
    'quick': dict(L_main=16, L_scan=[8, 10, 12, 14], deltas=[-0.5, 0.0, 0.5, 1.0], chunk=4, m=60,
                  L_ff_dense=10, L_full=8, deltas_val=[-0.5, 0.0, 0.5, 1.0], L_quspin=12, delta_quspin=0.5,
                  ell_min_b=3, ell_min_B=4, ff_inf_ells=[100, 200, 400], ff_inf_flux_max=400, ff_ring_L=[24, 96], p_cmp=1e-6),
    'check': dict(L_main=20, L_scan=[8, 12, 16], deltas=[-0.5, 0.0, 0.5, 1.0], chunk=4, m=50,
                  L_ff_dense=12, L_full=10, deltas_val=[-0.5, 0.0, 0.5, 1.0], L_quspin=16, delta_quspin=0.5,
                  ell_min_b=3, ell_min_B=4, ff_inf_ells=[100, 400, 1600], ff_inf_flux_max=1600, ff_ring_L=[24, 96, 400], p_cmp=1e-6),
    'full': dict(L_main=24, L_scan=[8, 10, 12, 14, 16, 18, 20, 22],
                 deltas=[-0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0], chunk=4, m=40,
                 L_ff_dense=14, L_full=12, deltas_val=[-0.5, 0.0, 0.5, 1.0], L_quspin=20, delta_quspin=0.5,
                 ell_min_b=3, ell_min_B=4, ff_inf_ells=[100, 200, 400, 800, 1600, 3200], ff_inf_flux_max=3200,
                 ff_ring_L=[24, 96, 400, 1600, 3200], p_cmp=1e-6),
}


def run(mode, dev='cpu', seed=0, log=print):
    cfg = CONFIGS[mode]
    warnings.simplefilter('ignore', cost.CostWarning)
    t0 = time.time()
    res = dict(config=cfg)
    log("main")
    res['main'] = {}
    xx = []
    for bc in ('periodic', 'open'):
        d = chain_data(cfg['L_main'], bc, cfg['deltas'], cfg, dev, seed, log)
        xx.append(d.pop('_xx'))
        if bc == 'periodic':
            d['K_fcs'] = fit_K_fcs(d, cfg['ell_min_b'])
        for n in (1, 2, 3):
            d[f'K_width_{n}'] = fit_K_width(d, n, cfg['ell_min_B'])
        for i, delta in enumerate(d['deltas']):
            log(f"  {bc} Delta={delta:+.2f}: K = {d['K'][i]:.4f}; "
                + (f"from b_1: {d['K_fcs'][i]:.4f}; " if bc == 'periodic' else '')
                + "from n B_n, n=1,2,3: " + ', '.join(f"{d[f'K_width_{n}'][i]:.4f}" for n in (1, 2, 3)))
        res['main'][bc] = d
    log("validation")
    res['validation'] = part_validation(cfg, dev, seed, log)
    for item in xx:                                       # the XX ground states of the main run (Lanczos)
        if item is not None:
            ch, psi, E0 = item
            res['validation']['ff'].append(compare_free_fermions(ch, psi, E0, 'lanczos', cfg['p_cmp'], log))
    del xx
    log("scan")
    res['scan'] = part_scan(cfg, dev, seed, log)
    pm = res['main']['periodic']
    res['scan'].append(dict(L=cfg['L_main'], row=pm['rows'][-1], E0=pm['E0']))
    log("free fermions, large sizes")
    res['ff_inf'] = part_ff_inf(cfg, log)
    res['ff_ring'] = part_ff_ring(cfg, log)
    res['cost'] = cost_estimate(cfg['L_main'])
    res['seconds'] = time.time() - t0
    log(f"total {res['seconds']:.0f} s")
    return res


def cost_estimate(L):
    """Memory and flops: sector-resolved computation versus dense full-space ones (double precision, real)."""
    N, ell = L // 2, L // 2
    kN = math.comb(L, N)
    blocks = [math.comb(ell, q) for q in range(ell + 1)]
    return dict(L=L, sector_dim=kN, full_dim=2 ** L,
                sector_vector_GB=kN * 8 / 1e9, full_vector_GB=2 ** L * 8 / 1e9,
                dense_sector_H_TB=kN ** 2 * 8 / 1e12, dense_full_H_TB=4 ** L * 8 / 1e12,
                rhoA_dense_flops=float((2 ** ell) ** 3), rhoA_blocks_flops=float(sum(k ** 3 for k in blocks)),
                rhoA_dense_entries=4 ** ell, rhoA_block_entries=sum(k * k for k in blocks))


# ----------------------------------------------------------------------------------------------------------------
# output
# ----------------------------------------------------------------------------------------------------------------
def _sampled_gaussian_entropy(B, ell):
    x = np.arange(ell + 1) - ell / 2
    w = np.exp(-x ** 2 / (2 * B))
    w /= w.sum()
    return float(-(w * np.log(w)).sum())


def summary(res):
    """Numbers quoted in the text, from the main data at ell = L/2 (and ell_ref = L/4 for the stability of G_1)."""
    out = {}
    for bc in ('periodic', 'open'):
        d = res['main'][bc]
        L = d['L']
        r = d['rows'][-1]
        ell = r['ell']
        q0 = (ell + 1) // 2
        ell_ref = 2 * (L // 8)
        rref = next(x for x in d['rows'] if x['ell'] == ell_ref)
        per = []
        for i, delta in enumerate(d['deltas']):
            Sq, S2q = r['Sq'][i], r['Snq'][2][i]
            B1 = r['B'][1][i]
            dS1 = Sq[q0] - Sq[q0 + 1]
            dS2 = Sq[q0] - Sq[q0 + 2]
            nB = [n * r['B'][n][i] for n in (1, 2, 3)]
            per.append(dict(delta=delta, K=d['K'][i], S=r['S'][i], H=r['H'][i], frac_H=r['H'][i] / r['S'][i],
                            H_sampled=_sampled_gaussian_entropy(B1, ell), b1=r['b'][1][i], B1=B1, nB=nB,
                            nB_spread=(max(nB) - min(nB)) / B1, G1=r['G'][1][i], G1_ref=rref['G'][1][i],
                            B1_ref=rref['B'][1][i],
                            G2=r['G'][2][i], k4=r['k4'][i], dS1=dS1, dS2=dS2, ratio=dS2 / dS1,
                            dS1_renyi2=S2q[q0] - S2q[q0 + 1],
                            K_fcs=d['K_fcs'][i] if 'K_fcs' in d else None,
                            K_width=[d[f'K_width_{n}'][i] for n in (1, 2, 3)]))
        # log-slopes of G_1 and B_1 over the even blocks ell >= 4 (bosonisation: G_1 = O(1), B_1 ~ K log ell_c / pi^2)
        ev = [x for x in d['rows'] if x['ell'] >= 4 and x['ell'] % 2 == 0]
        X = np.array([x['logc'] for x in ev])
        for i, row in enumerate(per):
            sG = np.polyfit(X, [_get(x['G'], 1)[i] for x in ev], 1)[0]
            sB = np.polyfit(X, [_get(x['B'], 1)[i] for x in ev], 1)[0]
            row.update(slope_G1=float(sG), slope_B1=float(sB))
        out[bc] = dict(L=L, ell=ell, ell_ref=ell_ref, rows=per)
    return out


def write_outputs(args, res):
    from common import save_json, write_tex, env_macro, sci
    cfg = res['config']
    res['summary'] = summary(res)
    save_json(args, 'symmetry_resolved', res)
    sm = res['summary']
    val = res['validation']

    def mac(name, value):
        return f"\\newcommand{{\\SR{name}}}{{{value}}}\n"

    tex = env_macro(args, 'SR')
    tex += mac('L', cfg['L_main']) + mac('Half', cfg['L_main'] // 2) + mac('Quarter', sm['periodic']['ell_ref'])
    tex += mac('Dim', f"{res['cost']['sector_dim']:,}".replace(',', '{,}'))
    tex += mac('NDelta', len(cfg['deltas'])) + mac('Chunk', cfg['chunk']) + mac('LanczosM', cfg['m'])
    tex += mac('Deltas', ', '.join(f'{d:g}' for d in cfg['deltas']))
    tex += mac('ScanL', f"{cfg['L_scan'][0]}, {cfg['L_scan'][1]}, \\dots, {cfg['L_main']}")
    tex += mac('ellMinb', cfg['ell_min_b']) + mac('ellMinB', cfg['ell_min_B'])
    tex += mac('Tol', sci(TOL, 1))
    res_max = max(max(res['main'][bc]['residual']) for bc in ('periodic', 'open'))
    tex += mac('Residual', sci(res_max, 1))
    gs = sum(res['main'][bc]['seconds_gs'] for bc in ('periodic', 'open'))
    ent = sum(res['main'][bc]['seconds_entanglement'] for bc in ('periodic', 'open'))
    tex += mac('SecondsGS', f"{gs:.0f}") + mac('SecondsEnt', f"{ent:.0f}") + mac('Seconds', f"{res['seconds']:.0f}")
    # validation
    keys = ['p', 'pS', 'Z', 'b', 'G']
    ffd = [r for r in val['ff'] if r['how'] == 'dense']
    ffl = [r for r in val['ff'] if r['how'] == 'lanczos']

    def worst(rows, ks):
        return max(max(r['dev'][k] for k in ks) for r in rows)
    tex += mac('FFLDense', ffd[0]['L'])
    tex += mac('FFDenseDev', sci(max(worst(ffd, keys), 1e-17), 1))
    tex += mac('FFDenseDevS', sci(max(worst(ffd, ['S', 'S2', 'S3']), 1e-17), 1))
    tex += mac('FFLanczosDev', sci(max(worst(ffl, keys), 1e-17), 1))
    tex += mac('FFLanczosDevS', sci(max(worst(ffl, ['S', 'S2', 'S3']), 1e-17), 1))
    tex += mac('FFEDev', sci(max(max(r['dev']['E0'] for r in val['ff']), 1e-17), 1))
    tex += mac('PCmp', sci(cfg['p_cmp'], 1))
    fu = val['full']
    tex += mac('FullL', fu[0]['L']) + mac('FullDeltas', ', '.join(f'{d:g}' for d in fu[0]['deltas']))
    tex += mac('FullDev', sci(max(max(r['dp'], r['dS']) for r in fu), 1))
    tex += mac('FullOff', sci(max(max(r['offdiag'] for r in fu), 1e-18), 1))
    tex += mac('RoutesDev', sci(max(max(r['routes'] for r in fu), 1e-18), 1))
    if val['quspin']:
        q = val['quspin']
        tex += mac('QuSpinL', q[0]['L']) + mac('QuSpinDelta', f"{q[0]['delta']:g}")
        tex += mac('QuSpinDev', sci(max(max(r['dp'], r['dS']) for r in q), 1))
        tex += mac('QuSpinE', sci(max(max(r['dE'] for r in q), 1e-17), 1))
    else:
        tex += mac('QuSpinL', '--') + mac('QuSpinDelta', '--') + mac('QuSpinDev', '--') + mac('QuSpinE', '--')
    fi = res['ff_inf'][-1]
    tex += mac('InfEll', fi['ell']) + mac('InfGOne', f"{fi['G1']:.7f}") + mac('InfGRef', f"{1 / (2 * math.pi ** 2):.7f}")
    tex += mac('InfBDev', sci(abs(fi['b1_dev']), 1))
    fd = [r for r in res['ff_inf'] if 'dS1' in r]
    tex += mac('InfdSEll', fd[-1]['ell']) + mac('InfdSRatio', f"{fd[-1]['dS1'] / fd[-1]['dS_asym']:.3f}")
    tex += mac('InfdSRatioFirst', f"{fd[0]['dS1'] / fd[0]['dS_asym']:.3f}") + mac('InfdSEllFirst', fd[0]['ell'])
    tex += mac('InfRatio', f"{fd[-1]['ratio']:.2f}")
    rings = res['ff_ring']
    tex += mac('RingLmax', rings[-1]['L']) + mac('RingGOne', f"{rings[-1]['G1']:.4f}")
    tex += mac('RingRatio', f"{rings[-1]['ratio']:.2f}") + mac('RingRatioFirst', f"{rings[0]['ratio']:.2f}")
    tex += mac('RingLfirst', rings[0]['L'])
    # tables
    per, opn = sm['periodic']['rows'], sm['open']['rows']
    tex += "\\newcommand{\\SRKRows}{%\n"
    for rp, ro in zip(per, opn):
        tex += (f"${rp['delta']:+.2f}$ & {rp['K']:.3f} & {rp['K_fcs']:.3f} & "
                + ' & '.join(f"{k:.3f}" for k in rp['K_width']) + f" & {ro['K_width'][0]:.3f} \\\\\n")
    tex += "}\n"
    tex += "\\newcommand{\\SREqRows}{%\n"
    for r in per:
        tex += (f"${r['delta']:+.2f}$ & {r['S']:.3f} & {r['H']:.3f} & {r['H_sampled']:.3f} & {r['B1']:.3f} & "
                f"{100 * r['nB_spread']:.1f} & {r['G1_ref']:.4f} & {r['G1']:.4f} & {r['dS1']:.4f} & "
                f"{r['ratio']:.2f} & {r['dS1_renyi2']:.4f} \\\\\n")
    tex += "}\n"
    tex += mac('FracHMin', f"{min(r['frac_H'] for r in per):.2f}") + mac('FracHMax', f"{max(r['frac_H'] for r in per):.2f}")
    tex += mac('HSampledDevMax', f"{max(abs(r['H'] - r['H_sampled']) for r in per):.3f}")
    tex += mac('KFourMax', f"{max(abs(r['k4']) for r in per):.3f}")
    tex += mac('SpreadMax', f"{100 * max(r['nB_spread'] for r in per):.1f}")
    mid = [r for r in per if -0.5 <= r['delta'] <= 0.75]
    tex += mac('KFcsDevPct', f"{100 * max(abs(r['K_fcs'] / r['K'] - 1) for r in mid):.1f}")
    tex += mac('KFcsMidRange', f"{min(r['delta'] for r in mid):g}")
    tex += mac('KFcsMidRangeHi', f"{max(r['delta'] for r in mid):g}")
    tex += mac('KWidthDevPct', f"{100 * max(abs(r['K_width'][0] / r['K'] - 1) for r in per if r['delta'] <= 0.75):.0f}")
    half = next(r for r in per if r['delta'] == 0.5)
    tex += mac('KWidthHalfTwo', f"{100 * (half['K_width'][1] / half['K'] - 1):.0f}")
    tex += mac('KWidthHalfThree', f"{100 * (half['K_width'][2] / half['K'] - 1):.0f}")
    tex += mac('SpreadMin', f"{100 * min(r['nB_spread'] for r in per):.1f}")
    tex += mac('BGrowMin', f"{100 * min(r['B1'] / r['B1_ref'] - 1 for r in per):.0f}")
    heis = next(r for r in per if r['delta'] == 1.0)
    tex += mac('HeisKWidthThree', f"{heis['K_width'][2]:.3f}")
    tex += mac('HeisKFcs', f"{heis['K_fcs']:.3f}")
    tex += mac('RatioMin', f"{min(r['ratio'] for r in per):.2f}") + mac('RatioMax', f"{max(r['ratio'] for r in per):.2f}")
    tex += mac('dSMin', f"{min(r['dS1'] for r in per):.3f}") + mac('dSMax', f"{max(r['dS1'] for r in per):.3f}")
    tex += mac('GSlopeRatioMax', f"{100 * max(r['slope_G1'] / r['slope_B1'] for r in per):.0f}")
    tex += mac('GRelChangeMax', f"{100 * max(abs(r['G1'] / r['G1_ref'] - 1) for r in per):.0f}")
    xx = next(r for r in per if r['delta'] == 0.0)
    tex += mac('XXGOne', f"{xx['G1']:.4f}") + mac('XXdS', f"{xx['dS1']:.4f}")
    c = res['cost']
    tex += mac('CostSectorH', f"{c['dense_sector_H_TB']:.0f}") + mac('CostFullH', f"{c['dense_full_H_TB'] / 1e3:.1f}")
    tex += mac('CostVecRatio', f"{c['full_dim'] / c['sector_dim']:.1f}")
    tex += mac('CostRhoRatio', f"{c['rhoA_dense_flops'] / c['rhoA_blocks_flops']:.0f}")
    tex += mac('CostVecGB', f"{c['sector_vector_GB']:.3f}")
    write_tex(args, 'symmetry_resolved', tex)
    make_figure(args, res)


def _get(dct, n):
    """dct[n] for dictionaries keyed by n (in memory) or by str(n) (read back from the JSON file)."""
    return dct[n] if n in dct else dct[str(n)]


def make_figure(args, res):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.5))
    d = res['main']['periodic']
    L = d['L']
    cfg = res['config']
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(d['deltas'])))
    ax = axes[0]
    rows = [r for r in d['rows'] if r['ell'] >= 2]
    X = np.array([r['logc'] for r in rows])
    ells = np.array([r['ell'] for r in rows])
    show = [i for i, x in enumerate(d['deltas']) if x in (-0.5, 0.0, 0.5, 1.0)]
    for i in show:
        for n, mk in ((1, 'o'), (2, 's'), (3, '^')):
            y = np.array([n * _get(r['B'], n)[i] for r in rows])
            ax.plot(X, y, mk, ms=2.6 if n == 1 else 2.2, mfc='none' if n > 1 else cols[i], mew=0.6, color=cols[i],
                    label=(f"$\\Delta={d['deltas'][i]:g}$" if n == 1 else None))
        K = d['K'][i]
        y1 = np.array([_get(r['B'], 1)[i] for r in rows])
        sel = ells >= cfg['ell_min_B']
        c0 = np.mean(y1[sel] - K / math.pi ** 2 * X[sel])
        ax.plot(X, K / math.pi ** 2 * X + c0, '-', lw=0.7, color=cols[i])
    ax.set_xlabel('$\\log \\ell_c$')
    ax.set_ylabel('$nB_n$')
    ax.set_title(f'widths, $L={L}$ periodic', fontsize=9)
    ax.legend(frameon=False, fontsize=5.5, loc='upper left')
    ax = axes[1]
    D = np.linspace(-0.8, 1.0, 200)
    ax.plot(D, [luttinger_K(x) for x in D], color='0.3', lw=1.0, label='Bethe ansatz')
    ax.plot(d['deltas'], d['K_fcs'], 'D', ms=3, mfc='none', color='C3', label='$b_1$ (FCS), periodic')
    for n, mk in ((1, 'o'), (2, 's'), (3, '^')):
        ax.plot(d['deltas'], d[f'K_width_{n}'], mk, ms=3, mfc='none', color=f'C{n - 1}', label=f'$nB_{n}$, periodic')
    do = res['main']['open']
    ax.plot(do['deltas'], do['K_width_1'], 'x', ms=3, color='C4', label='$B_1$, open')
    ax.set_xlabel('$\\Delta$')
    ax.set_ylabel('$K$')
    ax.set_title('Luttinger parameter', fontsize=9)
    ax.legend(frameon=False, fontsize=5.5)
    ax = axes[2]
    r = d['rows'][-1]
    ell = r['ell']
    q = np.arange(ell + 1)
    q0 = (ell + 1) // 2
    for i, delta in enumerate(d['deltas']):
        Sq = np.array(r['Sq'][i], dtype=float)
        sel = np.array(r['p'][i]) > 1e-4
        ax.plot(q[sel] - ell / 2, Sq[sel] - Sq[q0], 'o-', ms=2.5, lw=0.6, color=cols[i])
        g = np.array(r['gauss1'][i])
        ax.plot(q[sel] - ell / 2, g[sel], ':', lw=0.8, color=cols[i])
    lo = min(min(np.array(r['Sq'][i], dtype=float)[np.array(r['p'][i]) > 1e-4] - r['Sq'][i][q0])
             for i in range(len(d['deltas'])))
    ax.set_ylim(1.15 * lo, 0.05 * abs(lo))
    ax.axhline(0, color='0.8', lw=0.6)
    ax.set_xlabel('$q-\\ell/2$')
    ax.set_ylabel('$S(q)-S(\\ell/2)$')
    ax.set_title(f'equipartition breaking, $\\ell={ell}$', fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{args.figdir}/symmetry_resolved.pdf", metadata={"CreationDate": None})


def main():
    from common import parse_args
    args = parse_args(__doc__)
    res = run(args.mode, args.device, args.seed, log=lambda s: print(s, flush=True))
    write_outputs(args, res)
    print("written:", args.out, args.figdir)


if __name__ == '__main__':
    main()
