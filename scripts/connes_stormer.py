#!/usr/bin/env python3
"""
Section "Use case: the Connes-Stormer entropy of two MASAs" of the paper.

For subalgebras P, Q of M = M_n with the trace tau = Tr/n, eta(t) = -t log t and the tau-preserving conditional
expectations E_P, E_Q, the Connes-Stormer relative entropy is

    H(P|Q) = sup  sum_j [ tau eta(E_Q x_j) - tau eta(E_P x_j) ]   over finite partitions of unity x_j >= 0, sum x_j = 1.

We compute H(Delta_n | U Delta_n U^*) for the diagonal MASA Delta_n and a rotated copy with torch_vn_algebra:
inclusions.Inclusion (the MASA C^n c M_n, with the unitary U for the rotated copy) and its conditional expectation for
the trace 'norm' (= Tr/n), the spectra of E(x_j) for tau eta, and optimize.extremize (batched multistart ascent, one
start per batch element, batched over the unitaries U as well) over

  * general partitions with m elements, x_j = S^{-1/2} Z_j Z_j^* S^{-1/2}, S = sum_j Z_j Z_j^* (Z_j free complex n x n;
    the power S^{-1/2} by the library's functional calculus), and
  * orthonormal bases, x_j = w_j w_j^* with W = (w_1 ... w_n) an optimize.UnitaryParam;
    the value is H(b(U^* W)) - H(b(W)), b(V) = (|V_ij|^2) the unistochastic matrix, H(b) = (1/n) sum_ij eta(b_ij).

Parts (paper, Section "Use case"):
  n = 2   the closed form of the Proposition (one-dimensional maximisation over the angle alpha) on a grid of
          c = |u_11|^2, against the library optimum over general partitions for complex U with random phases;
  n = 3,4 Haar-random U: optimum over general partitions (m = n^2) versus the optimum over bases, Choda's value
          h(U) = H(b(U)) (the partition into the minimal projections of Delta_n), and the naive upper bound
          G* = max_rho [S(E_Q rho) - S(E_P rho)] of the Lemma (rho over all states);
  n = 5   the same on fewer U; for the U where general partitions beat the best basis found, the basis search is
          repeated with many more starts (also at bases built from the directions of the best partition), and the
          values at both explicit points are recomputed in mpmath with 50 digits (the partition is re-normalised in
          mpmath, so the certified value is a lower bound for H at an exactly admissible point; the U used is the
          polar unitarisation of the double-precision U, stored in the JSON file). No upper bound over bases is
          certified: "general partitions beat bases" remains an observation;
  checks  H(P|P) = 0; complementary MASAs (Fourier U) give log n; two non-abelian examples, M_2 (+) C c M_3 and
          M_2 (x) 1 c M_4 against Ad(u) of themselves for seeded Haar u (m = n, 2n, n^2 elements).

Writes generated/connes_stormer.tex, generated/connes_stormer.json and figures/connes_stormer.pdf.

    python scripts/connes_stormer.py --mode full --device cpu        # paper values, about 15-25 min on a CPU
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

from torch_vn_algebra import TypeIAlgebra, cost, inclusions as inc, optimize as opt   # noqa: E402

DT = torch.complex128
TINY = 1e-300                       # eta(t) is evaluated at max(t, TINY): eta'(TINY) is finite


# ----------------------------------------------------------------------------------------------------------------
# the functional through the library
# ----------------------------------------------------------------------------------------------------------------
def eta(t):
    t = t.clamp_min(TINY)
    return -t * torch.log(t)


class Pair:
    """
    P = iota_0(N) and Q = Ad(U) P inside M = M_n, with the tau-preserving conditional expectations E_P, E_Q for
    tau = Tr/n (the trace 'norm' of the library). N = (+)_c M_{k_c} with inclusion matrix `lam` (one row: M is a
    factor); by default N = C^n, lam = (1, ..., 1), i.e. the diagonal MASA Delta_n. U: tensor (B, n, n), one unitary
    per batch element of the partitions it is applied to (or (1, n, n)).
    """

    def __init__(self, U, sizes=None, lam=None, device='cpu'):
        n = U.shape[-1]
        sizes = [1] * n if sizes is None else list(sizes)
        lam = [[1] * n] if lam is None else lam
        self.n = n
        self.N = TypeIAlgebra(sizes, sizes, precision='double', device=device)
        self.M = TypeIAlgebra([n], [n], precision='double', device=device)
        u = self.M.operator(U.reshape(-1, 1, n, n).to(device=device, dtype=DT))
        self.P = inc.Inclusion(self.N, self.M, lam)
        self.Q = inc.Inclusion(self.N, self.M, lam, unitary=u, validate=False)
        self.EP = self.P.conditional_expectation('norm')
        self.EQ = self.Q.conditional_expectation('norm')

    @staticmethod
    def tau_eta(E, X):
        """tau(eta(E(x))) = sum_c nu_c Tr eta(E(x)_c), nu = Lambda^T w the trace tau restricted to N (batch of x)."""
        Ex = E(X)
        nu = E.nu
        return sum(nu[:, c] * eta(w).sum(-1) for c, w in enumerate(Ex.eigenvalues()))

    def terms(self, x):
        """F(x) = tau eta(E_Q x) - tau eta(E_P x) for positive matrices x (batch, n, n)."""
        X = self.M.operator(x.reshape(-1, 1, self.n, self.n), is_self_adjoint=True, is_positive=True)
        return self.tau_eta(self.EQ, X) - self.tau_eta(self.EP, X)


class PartitionParam(opt.Param):
    """Partitions of unity of M_n with m elements: x_j = S^{-1/2} Z_j Z_j^* S^{-1/2}, S = sum_j Z_j Z_j^*."""

    def __init__(self, M, batch, m, init=None):
        n = M.k_max
        self.M, self.m, self.n = M, m, n
        Z = torch.randn(batch, m, n, n, dtype=DT, device=M.hilbert.device) if init is None else init.clone()
        self.Z = Z.requires_grad_(True)

    def parameters(self):
        return [self.Z]

    def elements(self):
        G = self.Z @ self.Z.mH                                                  # (B, m, n, n)
        S = self.M.operator(G.sum(1)[:, None], is_self_adjoint=True, is_positive=True)
        Si = S.power(-0.5).matrix[:, 0]
        return Si[:, None] @ G @ Si[:, None]


class BasisParam(opt.UnitaryParam):
    """Projective partitions x_j = w_j w_j^* onto the columns of a unitary W (optimize.UnitaryParam)."""

    def elements(self):
        W = self.matrix()[:, 0]
        return torch.einsum('bij,bkj->bjik', W, W.conj())                      # x_j = w_j w_j^*, (B, n, n, n)


class StateParam(opt.Param):
    """States rho = Z Z^* / Tr(Z Z^*) of M_n (for the bound G* of the Lemma)."""

    def __init__(self, M, batch):
        n = M.k_max
        self.M, self.n = M, n
        self.Z = torch.randn(batch, 1, n, n, dtype=DT, device=M.hilbert.device, requires_grad=True)

    def parameters(self):
        return [self.Z]

    def elements(self):
        R = self.Z @ self.Z.mH
        return R / torch.diagonal(R, dim1=-2, dim2=-1).sum(-1).real[..., None, None].to(R.dtype)


def objective_for(pair, param):
    """sum_j F(x_j) per start, for the elements (B, m, n, n) of `param` (pair built for the batch B m)."""
    def objective():
        x = param.elements()
        B, m = x.shape[:2]
        return pair.terms(x.reshape(B * m, pair.n, pair.n)).reshape(B, m).sum(-1)
    return objective


def maximise(U, kind, m=None, starts=8, rounds=8, steps=150, lr=0.05, decay=0.6, sizes=None, lam=None,
             device='cpu', base=None, init=None):
    """
    Batched multistart maximisation of the Connes-Stormer functional for every unitary in U (nU, n, n), `starts`
    starts each (all in one batch). kind: 'partition' (m elements), 'basis' (m = n) or 'state' (sup of n F(rho) over
    states: the bound G* of the Lemma). Returns (best value per U over all starts and iterates, value at the final
    iterate per start (nU, starts), the parameter object).
    """
    nU, n = U.shape[0], U.shape[-1]
    B = nU * starts
    M = TypeIAlgebra([n], [n], precision='double', device=device)
    if kind == 'partition':
        param = PartitionParam(M, B, m, init=init)
    elif kind == 'basis':
        m = n
        param = BasisParam(M, B, base=base)
    elif kind == 'state':
        m = 1
        param = StateParam(M, B)
    else:
        raise ValueError(kind)
    Ub = U.to(device=device, dtype=DT).repeat_interleave(starts, 0).repeat_interleave(m, 0)
    pair = Pair(Ub, sizes=sizes, lam=lam, device=device)
    f = objective_for(pair, param)
    scale = float(n) if kind == 'state' else 1.0          # G(rho) = S(E_Q rho) - S(E_P rho) = n F(rho)
    best = opt.extremize(lambda: scale * f(), [param], maximize=True, rounds=rounds, steps=steps, lr=lr,
                         decay=decay)
    with torch.no_grad():
        final = scale * f()
    return best.reshape(nU, starts).max(1).values, final.reshape(nU, starts), param


# ----------------------------------------------------------------------------------------------------------------
# independent references (NumPy / SciPy / mpmath)
# ----------------------------------------------------------------------------------------------------------------
def eta_np(t):
    t = np.clip(np.asarray(t, dtype=float), 0.0, None)
    return -t * np.log(np.where(t > 0, t, 1.0))


def cs_numpy(U, xs):
    """sum_j [tau eta(E_Q x_j) - tau eta(E_P x_j)] for P = Delta_n, Q = U Delta_n U^*, tau = Tr/n (plain NumPy)."""
    U = np.asarray(U)
    n = U.shape[0]
    tot = 0.0
    for x in xs:
        dP = np.real(np.diag(x))
        dQ = np.real(np.diag(U.conj().T @ x @ U))
        tot += (eta_np(dQ).sum() - eta_np(dP).sum()) / n
    return tot


def unistochastic_entropy(U):
    """Choda's value h(U) = (1/n) sum_ij eta(|u_ij|^2): the functional at the minimal projections of Delta_n."""
    U = np.asarray(U)
    return float(eta_np(np.abs(U) ** 2).sum() / U.shape[0])


def h2(p):
    """binary entropy in nats"""
    return eta_np(p) + eta_np(1 - np.asarray(p, dtype=float))


def phi2(alpha, theta):
    """value of the basis (cos a, sin a), (-sin a, cos a) for U = rotation by theta: h(cos^2(a - theta)) - h(cos^2 a)"""
    return h2(np.cos(alpha - theta) ** 2) - h2(np.cos(alpha) ** 2)


def psi2(beta):
    """psi(b) = sin(2b) log tan^2(b) (continuous, period pi/2): phi2(a, t) = int_{a - t}^{a} psi."""
    beta = np.asarray(beta, dtype=float)
    s, t2 = np.sin(2 * beta), np.tan(beta) ** 2
    ok = (np.abs(s) > 0) & np.isfinite(t2) & (t2 > 0)
    return np.where(ok, s * np.log(np.where(ok, t2, 1.0)), 0.0)


def formula_n2(c, grid=4001):
    """
    H(Delta_2 | U Delta_2 U^*) for |u_11|^2 = c by the Proposition: max over alpha of phi2(alpha, theta), cos^2 theta
    = c, theta in [0, pi/2]. phi2 has period pi/2 in alpha. Returns (value, alpha* in (-pi/4, pi/4]).
    """
    from scipy.optimize import minimize_scalar
    theta = math.acos(math.sqrt(min(max(c, 0.0), 1.0)))
    a = np.linspace(-math.pi / 4, math.pi / 4, grid)
    v = phi2(a, theta)
    i = int(np.argmax(v))
    d = a[1] - a[0]
    r = minimize_scalar(lambda t: -float(phi2(t, theta)), bounds=(a[i] - d, a[i] + d), method='bounded',
                        options={'xatol': 1e-13})
    if -r.fun >= v[i]:
        val, al = -float(r.fun), float(r.x)
    else:
        val, al = float(v[i]), float(a[i])
    al = (al + math.pi / 4) % (math.pi / 2) - math.pi / 4
    return val, al


def kappa2():
    """kappa = max psi (attained in (-pi/4, 0)): H <= kappa theta and H / theta -> kappa as theta -> 0."""
    from scipy.optimize import minimize_scalar
    r = minimize_scalar(lambda b: -float(psi2(b)), bounds=(-math.pi / 4 + 1e-12, -1e-12), method='bounded',
                        options={'xatol': 1e-14})
    return -float(r.fun), float(r.x)


def rotation_with_phases(c, rng):
    """D1 R_theta D2 with random diagonal unitaries D1, D2 and cos^2 theta = c."""
    th = math.acos(math.sqrt(c))
    R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    D1, D2 = np.diag(np.exp(2j * np.pi * rng.random(2))), np.diag(np.exp(2j * np.pi * rng.random(2)))
    return D1 @ R @ D2


def fourier(n):
    j = np.arange(n)
    return np.exp(2j * np.pi * np.outer(j, j) / n) / math.sqrt(n)


def haar(nU, n, device='cpu'):
    """nU Haar-random unitaries of M_n from the library's sampler, (nU, n, n)."""
    M = TypeIAlgebra([n], [n], precision='double', device=device)
    return M.random_unitary_operator(nU).matrix[:, 0]


def certify(U, Z, digits=50):
    """
    The functional at the partition x_j = S^{-1/2} Z_j Z_j^* S^{-1/2} (Z: (m, n, r)) for P = Delta_n, Q = Ad(U') P,
    U' = U (U^* U)^{-1/2}, recomputed in mpmath: partition and unitary are exact at `digits` digits. Returns
    (value, |sum x_j - 1|, |U'^* U' - 1|) as floats.
    """
    import mpmath as mp
    with mp.workdps(digits):
        def mpm(a):
            return mp.matrix([[mp.mpc(complex(v)) for v in row] for row in np.asarray(a)])

        def herm_power(A, p):
            E, V = mp.eighe((A + A.transpose_conj()) / 2)
            return V * mp.diag([e ** p for e in E]) * V.transpose_conj()

        def eta_mp(t):
            return -t * mp.log(t) if t > 0 else mp.mpf(0)

        Um = mpm(U)
        Um = Um * herm_power(Um.transpose_conj() * Um, -0.5)
        n = Um.rows
        G = [mpm(z) * mpm(z).transpose_conj() for z in np.asarray(Z)]
        S = G[0]
        for g in G[1:]:
            S = S + g
        Si = herm_power(S, -0.5)
        xs = [Si * g * Si for g in G]
        tot = mp.mpf(0)
        for x in xs:
            y = Um.transpose_conj() * x * Um
            tot += sum(eta_mp(mp.re(y[i, i])) - eta_mp(mp.re(x[i, i])) for i in range(n)) / n
        res = xs[0]
        for x in xs[1:]:
            res = res + x
        res = res - mp.eye(n)
        err = max(abs(res[i, j]) for i in range(n) for j in range(n))
        un = Um.transpose_conj() * Um - mp.eye(n)
        uerr = max(abs(un[i, j]) for i in range(n) for j in range(n))
        return float(tot), float(err), float(uerr)


# ----------------------------------------------------------------------------------------------------------------
# the parts of the computation
# ----------------------------------------------------------------------------------------------------------------
def part_n2(cfg, dev, seed, log=print):
    rng = np.random.default_rng(seed)
    torch.manual_seed(seed)
    grid = np.linspace(0.0, 1.0, cfg['n2_grid'])
    curve = [formula_n2(c) for c in grid]
    rows = []
    Us = np.stack([rotation_with_phases(c, rng) for c in cfg['n2_table']])
    t0 = time.time()
    best, _, _ = maximise(torch.as_tensor(Us), 'partition', m=4, starts=cfg['n2_starts'], rounds=cfg['rounds'],
                          steps=cfg['steps'], device=dev)
    secs = time.time() - t0
    for c, U, b in zip(cfg['n2_table'], Us, best.tolist()):
        H, al = formula_n2(c)
        th = math.acos(math.sqrt(c))
        # the optimal basis of the Proposition, moved by the left phase D1 of U = D1 R D2 (independent check, NumPy):
        # the phases of the first column of U are those of D1 up to a common factor (R_11, R_21 >= 0)
        v = np.array([math.cos(al), math.sin(al)])
        w = np.array([-math.sin(al), math.cos(al)])
        D1 = np.diag(np.exp(1j * np.angle(U[:, 0]))) if c not in (0.0, 1.0) else np.eye(2)
        xs = [np.outer(D1 @ v, (D1 @ v).conj()), np.outer(D1 @ w, (D1 @ w).conj())]
        rows.append(dict(c=c, theta=th, choda=float(h2(c)), H=H, alpha=al, library=b, numpy_at_alpha=cs_numpy(U, xs),
                         kappa_theta=kappa2()[0] * th))
    dev_lib = max(abs(r['library'] - r['H']) for r in rows)
    dev_np = max(abs(r['numpy_at_alpha'] - r['H']) for r in rows)
    log(f"n=2: |library - formula| <= {dev_lib:.1e}, |numpy at alpha* - formula| <= {dev_np:.1e} "
        f"({len(rows)} values of c, complex U with random phases, {cfg['n2_starts']} starts, m = 4; {secs:.0f} s)")
    for r in rows:
        log(f"   c={r['c']:.3f}: H={r['H']:.10f}  alpha*={r['alpha']:+.6f}  library={r['library']:.10f}  "
            f"h(c)={r['choda']:.6f}")
    kap, bstar = kappa2()
    return dict(grid=grid.tolist(), H=[v for v, _ in curve], alpha=[a for _, a in curve],
                choda=h2(grid).tolist(), rows=rows, max_dev_library=dev_lib, max_dev_numpy=dev_np, kappa=kap,
                kappa_at=bstar, seconds=secs)


def directions_bases(x, count, rng):
    """Unitaries built from the rank-one directions of a partition x (m, n, n): `count` random choices of n
    directions (weighted by their trace), slightly perturbed and orthonormalised by the polar decomposition."""
    n = x.shape[-1]
    w, V = np.linalg.eigh(x)
    tr = w.sum(-1)
    keep = tr > 1e-8
    dirs = V[keep][:, :, -1]                          # principal eigenvector of every element
    p = tr[keep] / tr[keep].sum()
    out = []
    for _ in range(count):
        idx = rng.choice(len(dirs), size=n, replace=len(dirs) < n, p=p)
        A = dirs[idx].T + 1e-3 * (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
        Uu, _, Vh = np.linalg.svd(A)
        out.append(Uu @ Vh)
    return np.stack(out)


def structure(x, tol=1e-4):
    """Rank-one structure of a partition x (m, n, n): elements of non-negligible weight, whether they are all of rank
    one, the number of distinct directions (principal eigenvectors with |<d, d'>| > 1 - tol merged), and the largest
    overlap |<d, d'>| between distinct directions."""
    w, V = np.linalg.eigh(x)
    tr = w.sum(-1)
    keep = tr > 1e-7 * tr.max()
    rank_one = bool((w[keep][:, -2] <= 1e-6 * tr[keep]).all()) if x.shape[-1] > 1 else True
    clusters = []
    for d, t in zip(V[keep][:, :, -1], tr[keep]):
        for cl in clusters:
            if abs(np.vdot(cl[1], d)) > 1 - tol:
                cl[0] += t
                break
        else:
            clusters.append([t, d])
    ov = [abs(np.vdot(a[1], b[1])) for i, a in enumerate(clusters) for b in clusters[i + 1:]]
    return dict(elements=int(keep.sum()), rank_one=rank_one, directions=len(clusters),
                max_overlap=float(max(ov)) if ov else 0.0, weights=sorted(float(c[0]) for c in clusters))


def compare(n, nU, cfg, dev, seed, log=print, Us=None):
    """
    nU Haar unitaries of M_n: optimum over general partitions (m = n^2 elements) and over bases, Choda's value and G*.
    Where general partitions beat the best basis by more than cfg['flag'], the basis search is repeated with
    2 cfg['ext_starts'] more starts per U (Haar, and bases built from the directions of the best partition), and both
    explicit points (best partition, best basis) are recomputed in mpmath.
    """
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    U = haar(nU, n, dev) if Us is None else Us
    m, S, SB = n * n, cfg['starts_gen'], cfg['starts_basis']
    t0 = time.time()
    hg, fg, pg = maximise(U, 'partition', m=m, starts=S, rounds=cfg['rounds'], steps=cfg['steps'], device=dev)
    t1 = time.time()
    hb, fb, pb = maximise(U, 'basis', starts=SB, rounds=cfg['rounds'], steps=cfg['steps'], device=dev)
    t2 = time.time()
    gs, _, _ = maximise(U, 'state', starts=cfg['starts_state'], rounds=cfg['rounds'], steps=cfg['steps'], device=dev)
    t3 = time.time()
    with torch.no_grad():
        xall = pg.elements().reshape(nU, S, m, n, n).cpu().numpy()
        Zall = pg.Z.detach().reshape(nU, S, m, n, n).cpu().numpy()
        Wall = pb.matrix()[:, 0].reshape(nU, SB, n, n).cpu().numpy()
    hg, hb, gs = hg.cpu().numpy(), hb.cpu().numpy(), gs.cpu().numpy()
    fg, fb = fg.cpu().numpy(), fb.cpu().numpy()
    # G* >= G(rho_j) for the normalised elements of the best partition (some G(rho_j) >= its value, by the Lemma)
    Unp = U.cpu().numpy()
    for i in range(nU):
        x = xall[i, int(fg[i].argmax())]
        gs[i] = max(gs[i], max(n * cs_numpy(Unp[i], [xj / np.trace(xj).real]) for xj in x if np.trace(xj).real > 1e-12))
    first_flagged = [int(i) for i in np.nonzero(hg - hb > cfg['flag'])[0]]
    # explicit best points per U (final iterates)
    jg = fg.argmax(1)
    best_W = [Wall[i, int(fb[i].argmax())] for i in range(nU)]
    best_Wv = fb.max(1)
    hb_ext = hb.copy()
    K = cfg['ext_starts']
    if first_flagged:
        idx = torch.as_tensor(first_flagged)
        Uf = U[idx]
        hh, fh, ph = maximise(Uf, 'basis', starts=K, rounds=cfg['rounds'], steps=cfg['steps'], device=dev)
        base = np.concatenate([directions_bases(xall[i, jg[i]], K, rng) for i in first_flagged])
        hd, fd, pd = maximise(Uf, 'basis', starts=K, rounds=cfg['rounds'], steps=cfg['steps'], device=dev,
                              base=torch.as_tensor(base, dtype=DT, device=dev)[:, None])
        with torch.no_grad():
            Wh = ph.matrix()[:, 0].reshape(len(first_flagged), K, n, n).cpu().numpy()
            Wd = pd.matrix()[:, 0].reshape(len(first_flagged), K, n, n).cpu().numpy()
        fh, fd = fh.cpu().numpy(), fd.cpu().numpy()
        for q, i in enumerate(first_flagged):
            hb_ext[i] = max(hb[i], float(hh[q]), float(hd[q]))
            for fv, Wv in ((fh, Wh), (fd, Wd)):
                k = int(fv[q].argmax())
                if fv[q, k] > best_Wv[i]:
                    best_Wv[i], best_W[i] = float(fv[q, k]), Wv[q, k]
    t4 = time.time()
    flagged = [i for i in first_flagged if hg[i] - hb_ext[i] > cfg['flag']]
    # the U with the largest gap: a much larger basis search (cfg['stress_starts'] Haar starts and as many at bases
    # built from the directions of the partition)
    stress = None
    if flagged and cfg.get('stress_starts'):
        i = max(flagged, key=lambda i: hg[i] - hb_ext[i])
        K2 = cfg['stress_starts']
        hh, fh, ph = maximise(U[i:i + 1], 'basis', starts=K2, rounds=cfg['rounds'], steps=cfg['steps'], device=dev)
        base = torch.as_tensor(directions_bases(xall[i, jg[i]], K2, rng), dtype=DT, device=dev)[:, None]
        hd, fd, pd = maximise(U[i:i + 1], 'basis', starts=K2, rounds=cfg['rounds'], steps=cfg['steps'], device=dev,
                              base=base)
        with torch.no_grad():
            for fv, pv in ((fh, ph), (fd, pd)):
                k = int(fv[0].argmax())
                if float(fv[0, k]) > best_Wv[i]:
                    best_Wv[i], best_W[i] = float(fv[0, k]), pv.matrix()[k, 0].cpu().numpy()
        hb_ext[i] = max(hb_ext[i], float(hh[0]), float(hd[0]))
        stress = dict(index=i, starts=2 * K2, best_basis=float(max(hh[0], hd[0])), seconds=time.time() - t4)
        flagged = [i for i in flagged if hg[i] - hb_ext[i] > cfg['flag']]
    cases = []
    for i in flagged:
        Ui = U[i].cpu().numpy()
        cg, eg, ue = certify(Ui, Zall[i, jg[i]], cfg['digits'])
        cb, eb, _ = certify(Ui, best_W[i].T[:, :, None], cfg['digits'])
        st = structure(xall[i, jg[i]])
        cases.append(dict(index=i, general=float(fg[i, jg[i]]), basis=float(best_Wv[i]), certified_general=cg,
                          certified_basis=cb, certified_gap=cg - cb, partition_residual=eg, basis_residual=eb,
                          unitarity_residual=ue, U_real=Ui.real.tolist(), U_imag=Ui.imag.tolist(), **st))
    t5 = time.time()
    choda = np.array([unistochastic_entropy(u) for u in U.cpu().numpy()])
    H = np.maximum(hg, hb_ext)
    gap = hg - hb_ext
    out = dict(n=n, nU=nU, m=m, starts_gen=S, starts_basis=SB, ext_starts=2 * K, flag=cfg['flag'],
               H_general=hg.tolist(), H_bases=hb_ext.tolist(), H_bases_first=hb.tolist(), choda=choda.tolist(),
               G_star=gs.tolist(), first_flagged=first_flagged, flagged=flagged, cases=cases, stress=stress,
               n_flagged=len(flagged),
               certified_flagged=sum(c['certified_gap'] > cfg['flag'] and c['directions'] > n for c in cases),
               max_gap=float(gap.max()), min_gap=float(gap.min()),
               agree=float(np.abs(np.delete(gap, flagged)).max()) if len(flagged) < nU else 0.0,
               max_certified_gap=max((c['certified_gap'] for c in cases), default=0.0),
               directions=sorted({c['directions'] for c in cases}), all_rank_one=all(c['rank_one'] for c in cases),
               excess_min=float((H - choda).min()), excess_median=float(np.median(H - choda)),
               excess_max=float((H - choda).max()), naive_not_attained=int((gs - H > 1e-6).sum()),
               naive_gap_min=float((gs - H).min()),
               seconds=dict(general=t1 - t0, bases=t2 - t1, states=t3 - t2, extended=t4 - t3, certify=t5 - t4))
    log(f"n={n}, {nU} Haar U: general partitions beat the best basis by more than {cfg['flag']:g} for "
        f"{len(first_flagged)} U, after {2 * K} more basis starts for {len(flagged)} (certified in mpmath, with more "
        f"than n directions: {out['certified_flagged']}); largest certified gap {out['max_certified_gap']:.2e}; "
        f"otherwise min(H_gen - H_bases) = {out['min_gap']:.1e}; directions of the improving partitions: "
        f"{out['directions']}, "
        f"all of rank one: {out['all_rank_one']}; H - h_Choda in [{out['excess_min']:.4f}, {out['excess_max']:.4f}], "
        f"median {out['excess_median']:.4f}; G* > H + 1e-6 for {out['naive_not_attained']}/{nU} (min G* - H = "
        f"{out['naive_gap_min']:.3f}); {t1 - t0:.0f} + {t2 - t1:.0f} + {t3 - t2:.0f} + {t4 - t3:.0f} + "
        f"{t5 - t4:.0f} s")
    if stress:
        log(f"   stress test U#{stress['index']}: best basis over {stress['starts']} further starts "
            f"{stress['best_basis']:.12f}; {stress['seconds']:.0f} s")
    for c in cases:
        log(f"   U#{c['index']}: partition {c['certified_general']:.12f}, basis {c['certified_basis']:.12f} (mpmath), "
            f"gap {c['certified_gap']:.3e}; {c['elements']} elements, {c['directions']} directions, rank one "
            f"{c['rank_one']}, largest overlap {c['max_overlap']:.4f}")
    return out


def part_checks(cfg, dev, seed, log=print):
    torch.manual_seed(seed)
    out = {}
    # H(P|P) = 0
    v, _, _ = maximise(torch.eye(3, dtype=DT)[None], 'partition', m=9, starts=cfg['starts_gen'], rounds=2,
                       steps=cfg['steps'] // 2, device=dev)
    out['self'] = float(v.abs().max())
    # complementary MASAs: H = log n
    fr = []
    for n in (2, 3, 4):
        F = torch.as_tensor(fourier(n))[None]
        hg, _, _ = maximise(F, 'partition', m=n * n, starts=cfg['starts_gen'], rounds=cfg['rounds'],
                            steps=cfg['steps'], device=dev)
        fr.append(dict(n=n, H=float(hg[0]), log_n=math.log(n), choda=unistochastic_entropy(fourier(n))))
    out['fourier'] = fr
    out['fourier_max_dev'] = max(abs(r['H'] - r['log_n']) for r in fr)
    log(f"checks: max |H(P|P)| = {out['self']:.1e}; Fourier: max |H - log n| = {out['fourier_max_dev']:.1e} "
        f"(n = 2, 3, 4)")
    # non-abelian subalgebras against a rotated copy
    nonab = []
    for name, sizes, lam, n in (('M2+C in M3', [2, 1], [[1, 1]], 3), ('M2(x)1 in M4', [2], [[2]], 4)):
        U = haar(cfg['nonab_U'], n, dev)
        vals = {}
        for m in (n, 2 * n, n * n):
            hg, _, _ = maximise(U, 'partition', m=m, starts=cfg['starts_gen'], rounds=cfg['rounds'],
                                steps=cfg['steps'], sizes=sizes, lam=lam, device=dev)
            vals[m] = hg.cpu().numpy().tolist()
        gs, _, _ = maximise(U, 'state', starts=cfg['starts_state'], rounds=cfg['rounds'], steps=cfg['steps'],
                            sizes=sizes, lam=lam, device=dev)
        H = np.max(np.array(list(vals.values())), axis=0)
        spread = np.ptp(np.array(list(vals.values())), axis=0)
        nonab.append(dict(name=name, sizes=sizes, lam=lam, n=n, H_by_m={str(k): v for k, v in vals.items()},
                          H=H.tolist(), spread=spread.tolist(), G_star=gs.cpu().numpy().tolist()))
        log(f"   {name}: H = {np.round(H, 6).tolist()} (spread over m = n, 2n, n^2: {spread.max():.1e}); G* = "
            f"{np.round(gs.cpu().numpy(), 6).tolist()} (log 2 = {math.log(2):.6f})")
    out['nonabelian'] = nonab
    return out


CONFIGS = {
    # n2_table: values of c for the library check (n2_starts starts each); n3, n4, n5: number of Haar U; rounds x steps
    # (steps5 for n = 5) of optimize.extremize; starts per U for general partitions, bases and states; ext_starts:
    # the extended basis search (twice this number of starts) for the U flagged by a gap above `flag`
    'quick': dict(n2_grid=41, n2_table=[0.95, 0.8, 0.6], n2_starts=8, rounds=5, steps=80, n3=4, n4=2, n5=2,
                  starts_gen=4, starts_basis=4, starts_state=4, starts_gen5=4, steps5=80, ext_starts=8, flag=1e-8,
                  stress_starts=8, digits=30, nonab_U=1),
    'check': dict(n2_grid=101, n2_table=[0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5], n2_starts=4, rounds=2, steps=20,
                  n3=300, n4=150, n5=24, starts_gen=2, starts_basis=2, starts_state=2, starts_gen5=2, steps5=20,
                  ext_starts=2, flag=1e-8, stress_starts=2, digits=50, nonab_U=2),
    'full': dict(n2_grid=401, n2_table=[0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5], n2_starts=16, rounds=8, steps=150,
                 n3=300, n4=150, n5=24, starts_gen=8, starts_basis=8, starts_state=16, starts_gen5=16, steps5=250,
                 ext_starts=32, flag=1e-8, stress_starts=512, digits=50, nonab_U=2),
}


def run(mode, dev='cpu', seed=0, log=print):
    """All parts; returns the result dictionary (no files written)."""
    cfg = CONFIGS[mode]
    warnings.simplefilter('ignore', cost.CostWarning)
    t0 = time.time()
    res = dict(config=cfg)
    res['n2'] = part_n2(cfg, dev, seed, log)
    for n in (3, 4, 5):
        c = dict(cfg, starts_gen=cfg['starts_gen5'], starts_basis=cfg['starts_gen5'], steps=cfg['steps5']) \
            if n == 5 else cfg
        res[f'n{n}'] = compare(n, cfg[f'n{n}'], c, dev, seed + n, log)
    res['checks'] = part_checks(cfg, dev, seed + 7, log)
    res['seconds'] = time.time() - t0
    log(f"total {res['seconds']:.0f} s")
    return res


# ----------------------------------------------------------------------------------------------------------------
# output
# ----------------------------------------------------------------------------------------------------------------
def write_outputs(args, res):
    from common import save_json, write_tex, env_macro, sci
    cfg = res['config']
    save_json(args, 'connes_stormer', res)

    def f6(x):
        return f"{x:.6f}"

    tex = env_macro(args, 'CS')
    tex += "\\newcommand{\\CSTwoRows}{%\n"
    for r in res['n2']['rows']:
        tex += (f"{r['c']:g} & {f6(r['choda'])} & {f6(r['H'])} & {f6(r['library'])} & {r['alpha']:+.4f} & "
                f"{f6(r['kappa_theta'])} \\\\\n")
    tex += "}\n"
    n2 = res['n2']
    tex += (f"\\newcommand{{\\CSTwoDevLib}}{{{sci(max(n2['max_dev_library'], 1e-16), 1)}}}\n"
            f"\\newcommand{{\\CSTwoDevNumpy}}{{{sci(max(n2['max_dev_numpy'], 1e-16), 1)}}}\n"
            f"\\newcommand{{\\CSTwoStarts}}{{{cfg['n2_starts']}}}\n"
            f"\\newcommand{{\\CSKappa}}{{{n2['kappa']:.4f}}}\n")
    tex += "\\newcommand{\\CSRows}{%\n"
    for key in ('n3', 'n4', 'n5'):
        r = res[key]
        gap = f"${sci(r['max_certified_gap'], 1)}$" if r['cases'] else '--'
        tex += (f"{r['n']} & {r['nU']} & {r['starts_gen']}/{r['starts_basis']} & "
                f"{r['certified_flagged']} & {gap} & ${sci(max(r['agree'], 1e-16), 1)}$ & "
                f"{r['excess_min']:.4f} & {r['excess_median']:.4f} & {r['excess_max']:.4f} & "
                f"{r['naive_not_attained']} \\\\\n")
    tex += "}\n"
    tex += (f"\\newcommand{{\\CSSteps}}{{{cfg['rounds'] * cfg['steps']}}}\n"
            f"\\newcommand{{\\CSFiveSteps}}{{{cfg['rounds'] * cfg['steps5']}}}\n"
            f"\\newcommand{{\\CSExtStarts}}{{{2 * cfg['ext_starts']}}}\n"
            f"\\newcommand{{\\CSFlag}}{{{sci(cfg['flag'], 1)}}}\n"
            f"\\newcommand{{\\CSDigits}}{{{cfg['digits']}}}\n")
    for key, tag in (('n3', 'Three'), ('n4', 'Four'), ('n5', 'Five')):
        r = res[key]
        dirs = ', '.join(map(str, r['directions'])) if r['directions'] else '--'
        tex += (f"\\newcommand{{\\CS{tag}U}}{{{r['nU']}}}\n"
                f"\\newcommand{{\\CS{tag}Flagged}}{{{r['certified_flagged']}}}\n"
                f"\\newcommand{{\\CS{tag}Agree}}{{{sci(max(r['agree'], 1e-16), 1)}}}\n"
                f"\\newcommand{{\\CS{tag}Dirs}}{{{dirs}}}\n"
                f"\\newcommand{{\\CS{tag}NaiveMin}}{{{r['naive_gap_min']:.3f}}}\n")
        tex += f"\\newcommand{{\\CS{tag}StressStarts}}{{{r['stress']['starts'] if r['stress'] else 0}}}\n"
        if r['cases']:                       # every macro is defined, '--' without a case
            b = max(r['cases'], key=lambda c: c['certified_gap'])
            vals = dict(GapMax=sci(b['certified_gap'], 2),
                        GapMin=sci(min(c['certified_gap'] for c in r['cases']), 1),
                        General=f"{b['certified_general']:.10f}", Basis=f"{b['certified_basis']:.10f}",
                        Overlap=f"{b['max_overlap']:.3f}", Residual=sci(max(b['partition_residual'], 1e-99), 1))
        else:
            vals = {k: '\\text{--}' for k in ('GapMax', 'GapMin', 'General', 'Basis', 'Overlap', 'Residual')}
        for k, v in vals.items():
            tex += f"\\newcommand{{\\CS{tag}{k}}}{{{v}}}\n"
    ch = res['checks']
    tex += (f"\\newcommand{{\\CSSelf}}{{{sci(max(ch['self'], 1e-17), 1)}}}\n"
            f"\\newcommand{{\\CSFourierDev}}{{{sci(max(ch['fourier_max_dev'], 1e-16), 1)}}}\n")
    for r, tag in zip(ch['nonabelian'], ('A', 'B')):
        tex += (f"\\newcommand{{\\CSNonab{tag}}}{{{', '.join(f'{v:.4f}' for v in r['H'])}}}\n"
                f"\\newcommand{{\\CSNonab{tag}Spread}}{{{sci(max(max(r['spread']), 1e-16), 1)}}}\n"
                f"\\newcommand{{\\CSNonab{tag}G}}{{{', '.join(f'{v:.4f}' for v in r['G_star'])}}}\n")
    tex += f"\\newcommand{{\\CSSeconds}}{{{res['seconds']:.0f}}}\n"
    write_tex(args, 'connes_stormer', tex)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.5))
    g = np.array(n2['grid'])
    ax = axes[0]
    ax.plot(g, n2['H'], color='C0', lw=1.4, label='$H(\\Delta_2|U\\Delta_2U^*)$')
    ax.plot(g, n2['choda'], color='C1', lw=1.0, ls='--', label='$h(c)$ (Choda)')
    ax.plot([r['c'] for r in n2['rows']], [r['library'] for r in n2['rows']], 'o', ms=3.5, mfc='none', color='k',
            label='library, partitions')
    ax.set_xlabel('$c=|u_{11}|^2$')
    ax.set_ylabel('nats')
    ax.legend(frameon=False, fontsize=6.5, loc='lower center')
    ax.set_title('$n=2$', fontsize=9)
    ax = axes[1]
    th = np.arccos(np.sqrt(g))
    al = np.array(n2['alpha'])
    sel = (g > 0) & (g < 1) & (np.abs(g - 0.5) > 1e-9)
    ax.plot(g[sel], np.degrees(al[sel]), color='C0', lw=1.2, label='optimal $\\alpha^*$')
    ax.plot(g, np.degrees(np.minimum(th, np.pi / 2 - th)), color='0.5', lw=0.8, ls=':',
            label='$\\min(\\theta,\\frac{\\pi}{2}-\\theta)$')
    ax.axhline(0, color='0.8', lw=0.6)
    ax.set_xlabel('$c=|u_{11}|^2$')
    ax.set_ylabel('angle [degrees]')
    ax.legend(frameon=False, fontsize=6.5)
    ax.set_title('optimal basis, $n=2$', fontsize=9)
    ax = axes[2]
    for key, col in (('n3', 'C0'), ('n4', 'C2')):
        r = res[key]
        H = np.maximum(r['H_general'], r['H_bases'])
        ax.plot(r['choda'], H, '.', ms=2.5, color=col, label=f"$n={r['n']}$")
    lo = min(min(res[k]['choda']) for k in ('n3', 'n4'))
    hi = max(math.log(4), max(max(res[k]['H_general']) for k in ('n3', 'n4')))
    ax.plot([lo, hi], [lo, hi], color='0.6', lw=0.7, ls='--')
    ax.set_xlabel('$h(U)$ (Choda)')
    ax.set_ylabel('$H$')
    ax.legend(frameon=False, fontsize=6.5, loc='lower right', markerscale=2)
    ax.set_title('Haar $U$', fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{args.figdir}/connes_stormer.pdf", metadata={"CreationDate": None})


def main():
    from common import parse_args
    args = parse_args(__doc__)
    res = run(args.mode, args.device, args.seed, log=lambda s: print(s, flush=True))
    write_outputs(args, res)
    print("written:", args.out, args.figdir)


if __name__ == '__main__':
    main()
