#!/usr/bin/env python3
"""
Tomography in the algebra as training: learning the noise process of the XY register from measurement counts
(man/article, Section "Tomography in the algebra as training"; the register is that of Section "Noisy circuits with a
conserved charge").

A ring of n qubits evolves for a time T under the XY mixer H = (1/2) sum_<ij> (S^+_i S^-_j + S^-_i S^+_j) with
amplitude damping S^-_i (rate g_down), thermal excitation S^+_i (g_up) and dephasing S^z_i (g_phi), all sparse
SectorJumps. The true process Phi = exp(T L) is a channel of the algebra M = (+)_N M_binom(n, N). It is probed with
the basis states of all weights and k_N^2 Haar-random pure states in every sector N (an informationally complete set
of inputs for M); every output is measured in k_max + 3 Haar-random bases of M (POVM.from_unitary of block-diagonal
unitaries, combined with equal weights; every outcome also reveals the charge N, so POVM.sectors is a coarse-graining
of each setting). Counts are sampled with a fixed total number of shots per probe (1e2 ... 1e5), and the exact
probabilities serve as the limit of infinitely many shots. Four models are fitted by maximum likelihood
(learn.negative_log_likelihood; torch.optim.Adam with cosine decay, then torch.optim.LBFGS; several random starts, the
start with the smallest training loss is kept):

  rates     the noise model with known jump operators and four unknowns (J, g_down, g_up, g_phi);
  lindblad  learn.LindbladLayer: a free self-adjoint H in M, n-1 free jumps inside the sectors and n free jumps
            between neighbouring sectors N -> N +- 1 (contains the true generator);
  algebra   learn.ChannelLayer from M to M with all sector transitions and Kraus rank k_max^2 (every channel of M);
  dense     learn.ChannelLayer on M_{2^n} with Kraus rank 4^n (every channel of M_{2^n}), trained on the same probes
            embedded by interop.to_dense and measured in 2^n + 3 Haar-random bases of M_{2^n}: the experiment
            without knowledge of the symmetry, which has to learn that the process keeps the sectors.

Training evaluates the models through fast paths that are checked against the layers' own forward passes at the start
(exact exponential of the generator on M, of size sum_N binom(n, N)^2 instead of the padded (n+1) k_max^2; Kraus
operators applied to the pure probe vectors); evaluation uses the layers themselves. Reported per model and shot
budget: real parameters (learn.count_parameters), mean trace distance and infidelity of the predicted outputs of
held-out inputs, the error of the Choi matrix restricted to M, training time; for the rate model and the GKSL model
the learned rates (GKSL: least-squares projection of the learned generator onto -i[H, .] and the three dissipators),
with their spread over independent data sets and, for the rate model, the Cramer-Rao error bars from the Fisher
information. Writes man/article/generated/channel_learning.{json,tex} and man/article/figures/channel_learning.pdf.

    python scripts/channel_learning.py [--dims 3,4] [--adam 300] [--lbfgs 400] [--repeats 6] [--workers 4]
    python scripts/channel_learning.py --replot man/article/generated/channel_learning.json   # figure and macros only

About 15 minutes on 4 CPU cores (one fit per core, one thread each).
"""
import argparse
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import torch
from torch import nn

from common import environment                       # device, versions, matmul precision, commit

from torch_vn_algebra import SpinChain, TypeIAlgebra, dynamics, interop, learn
from torch_vn_algebra.measurement import POVM, fisher_information
from torch_vn_algebra.states import DensityMatrix

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'man' / 'article'
T = 1.0
J_TRUE = 1.0
RATES = {'down': 0.1, 'up': 0.02, 'phi': 0.05}
KINDS = ('down', 'up', 'phi')
SECTOR_JUMP = {'down': 'lower', 'up': 'raise', 'phi': 'sz'}
SHOTS = [100, 1000, 10000, 100000, None]               # per probe, all settings together; None: exact probabilities
MODELS = ('rates', 'lindblad', 'algebra', 'dense')
LR = {'rates': 0.05, 'lindblad': 0.02, 'algebra': 0.05, 'dense': 0.05}
H_SCALE = 0.1                                          # initial scale of the GKSL layer's Hamiltonian
COLORS = {'rates': '#1baf7a', 'lindblad': '#2a78d6', 'algebra': '#eda100', 'dense': '#eb6834'}
LABELS = {'rates': 'rates (4 parameters)', 'lindblad': r'GKSL in $\mathcal{M}$, jumps $N\to N, N\pm1$',
          'algebra': r'channel of $\mathcal{M}$', 'dense': r'channel of $M_{2^n}$ (dense)'}
GREY, INK = '#7d848c', '#0b0b0b'
OFFSET = {'rates': 0.7, 'lindblad': 0.85, 'algebra': 1.15, 'dense': 1.4}           # 'exact' markers side by side


# ----------------------------------------------------------------------
# the experiment
# ----------------------------------------------------------------------
def _hermitian_basis(alg, c):
    """Hermitian matrix units of sector c (E_aa, E_ab + E_ba, -i(E_ab - E_ba)) and the coefficients that recombine
    their images into the images of the matrix units E_ab (lindblad_evolve symmetrises its output, so the maps are
    applied to Hermitian inputs only)."""
    k = alg.k_factors[c]
    herm, recipe = [], []
    for a in range(k):
        for b in range(a, k):
            m = torch.zeros(alg.C, alg.k_max, alg.k_max, dtype=alg.dtype)
            if a == b:
                m[c, a, a] = 1
                recipe.append(((a, a), [(len(herm), 1.0)]))
                herm.append(m)
            else:
                m[c, a, b] = m[c, b, a] = 1
                m2 = torch.zeros_like(m)
                m2[c, a, b], m2[c, b, a] = -1j, 1j
                i = len(herm)
                herm += [m, m2]
                recipe.append(((a, b), [(i, 0.5), (i + 1, 0.5j)]))       # E_ab = (X + i Y)/2
                recipe.append(((b, a), [(i, 0.5), (i + 1, -0.5j)]))
    return torch.stack(herm), recipe


def _trace_norm_herm(X):
    return torch.linalg.eigvalsh(X).abs().sum(-1)


class Problem:
    """Register, true process, probes, measurements and exact outcome probabilities for one n (deterministic in
    the seed)."""

    def __init__(self, n, seed=0, test_pure=3, test_mixed=8):
        self.n = n
        g = torch.Generator().manual_seed(1000 * n + seed)
        torch.manual_seed(1000 * n + seed)                  # random_unitary_operator uses the global generator
        q = self.q = SpinChain(n, 'periodic', precision='double')
        alg = self.alg = q.algebra
        D = self.D = 2 ** n
        self.dense = TypeIAlgebra([D], [D], precision='double')
        self.H = q.xxz(J=1.0, Delta=0.0)
        self.ops = {k: [q.sector_jump(i, SECTOR_JUMP[k]) for i in range(n)] for k in KINDS}
        self.mask = alg._const('block_mask', torch.bool)                       # (C, k_max, k_max) active entries
        self.dim = int(self.mask.sum())                                        # sum_N binom(n, N)^2
        # training inputs: pure vectors (sector, vector): all basis states, then k_N^2 Haar-random states per sector
        vecs = [(c, torch.eye(k, dtype=alg.dtype)[a]) for c, k in enumerate(alg.k_factors) for a in range(k)]
        for c, k in enumerate(alg.k_factors):
            if k > 1:
                v = torch.randn(k * k, k, dtype=alg.dtype, generator=g)
                vecs += [(c, x / torch.linalg.vector_norm(x)) for x in v]
        self.train = self._pure_states(vecs)
        self.groups = []                                                       # per sector: probe indices, vectors
        for c, k in enumerate(alg.k_factors):
            idx = [i for i, (cc, _) in enumerate(vecs) if cc == c]
            self.groups.append((c, torch.tensor(idx), torch.stack([vecs[i][1] for i in idx])))
        pos = interop.chain_basis(q)
        start = np.cumsum([0] + list(alg.k_factors))
        self.dense_vecs = torch.zeros(len(vecs), D, dtype=alg.dtype)
        for i, (c, v) in enumerate(vecs):
            self.dense_vecs[i, pos[start[c]:start[c + 1]]] = v
        # held-out inputs: 3 pure states per sector with k_N > 1 and 8 random mixed states of M (all sectors)
        test = [(c, x / torch.linalg.vector_norm(x)) for c, k in enumerate(alg.k_factors) if k > 1
                for x in torch.randn(test_pure, k, dtype=alg.dtype, generator=g)]
        G = torch.randn(test_mixed, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype, generator=g) * self.mask
        R = G @ G.conj().transpose(-2, -1)
        R = R / torch.diagonal(R, dim1=-2, dim2=-1).sum((-2, -1))[:, None, None, None]
        self.test = DensityMatrix(alg, matrix=torch.cat([self._pure_states(test).matrix, R]), validate=False)
        self.n_train, self.n_test = len(vecs), self.test.matrix.shape[0]
        # measurements: k_max + 3 Haar-random bases of M; 2^n + 3 Haar-random bases of M_{2^n}
        self.settings_alg, self.settings_dense = alg.k_max + 3, D + 3
        U = alg.random_unitary_operator(self.settings_alg).matrix
        self.povm = POVM.combine([POVM.from_unitary(alg.operator(U[s:s + 1])) for s in range(self.settings_alg)])
        Ud = self.dense.random_unitary_operator(self.settings_dense).matrix
        self.povm_dense = POVM.combine([POVM.from_unitary(self.dense.operator(Ud[s:s + 1]))
                                        for s in range(self.settings_dense)])
        ones = torch.tensor([bin(x).count('1') for x in range(D)])              # n - N in QuTiP's convention
        self.block_mask_dense = ones[:, None] == ones[None, :]                  # entries inside the sectors
        self.fast_povm = MatrixPOVM(self.povm, self.mask)
        self.fast_povm_dense = MatrixPOVM(self.povm_dense, torch.ones(1, D, D, dtype=torch.bool))
        # generator matrices on M (unpadded) of the physical terms, true outputs and exact probabilities
        with torch.no_grad():
            self.units = torch.zeros(self.dim, *self.mask.shape, dtype=alg.dtype)    # matrix units of M
            flat = self.mask.reshape(-1).nonzero()[:, 0]
            self.units.reshape(self.dim, -1)[torch.arange(self.dim), flat] = 1
            self.G_phys = [self.superop(self.H, [])] + [self.superop(None, self.ops[k]) for k in KINDS]
            self.R_train = self.train.matrix[:, self.mask]                        # (probes, dim)
            self.out_train = self.evolve(self.train, J_TRUE, RATES)
            self.out_test = self.evolve(self.test, J_TRUE, RATES)
            self.p = self.povm.probabilities(self.out_train)
            self.p_dense = self.povm_dense.probabilities(self.to_dense_state(self.out_train))
            self.true_test = self.dense_matrix(self.out_test)
            self.choi_inputs = [_hermitian_basis(alg, c) for c in range(alg.C)]
            self.true_choi = self.choi(lambda rho: self.evolve(rho, J_TRUE, RATES), dense=False)

    def _pure_states(self, vecs):
        alg = self.alg
        mat = torch.zeros(len(vecs), alg.C, alg.k_max, alg.k_max, dtype=alg.dtype)
        for i, (c, v) in enumerate(vecs):
            k = len(v)
            mat[i, c, :k, :k] = v[:, None] * v.conj()[None, :]
        return DensityMatrix(alg, matrix=mat, validate=False)

    # ---- the true process
    def jumps(self, rates):
        return [L.scaled(torch.sqrt(torch.as_tensor(rates[k], dtype=torch.float64))) for k in KINDS
                for L in self.ops[k]]

    def evolve(self, rho, J, rates):
        """exp(T L)(rho) by dynamics.lindblad_evolve (exact exponential of the generator)."""
        H = self.alg.operator(J * self.H.matrix, is_self_adjoint=True)
        return dynamics.lindblad_evolve(rho, H, self.jumps(rates), [0.0, T], method='exact')[-1]

    # ---- generators as matrices on M (row-major entries of the blocks, sector after sector)
    def superop(self, H, jumps):
        """Matrix (dim x dim, complex) of -i[H, .] + sum_j D[L_j] on M, from dynamics.lindblad_rhs applied to the
        matrix units (the right-hand side is complex linear)."""
        Y = dynamics.lindblad_rhs(H, jumps, algebra=self.alg)(0.0, self.units)
        return Y[:, self.mask].transpose(0, 1)

    def propagate(self, G, R):
        """States exp(T G)(rho) for the rows R (probes x dim) of the input states."""
        out = R @ torch.linalg.matrix_exp(T * G).transpose(0, 1)
        mat = torch.zeros(R.shape[0], *self.mask.shape, dtype=out.dtype)
        mat[:, self.mask] = out
        return DensityMatrix(self.alg, matrix=mat, validate=False)

    # ---- dense representation (interop, computational basis of the n qubits)
    def dense_matrix(self, x):
        return interop.to_dense(x, layout=self.q)

    def to_dense_state(self, rho):
        return DensityMatrix(self.dense, matrix=self.dense_matrix(rho)[:, None], validate=False)

    # ---- Choi matrices of a map restricted to M
    def choi(self, fn, dense):
        """Normalised Choi matrices J_c = (1/k_c) sum_ab Phi(E_ab) (x) E_ab of the restriction to every input sector c
        (dense D k_c x D k_c matrices, outputs in the computational basis); fn acts on states of M (dense=False) or
        of M_{2^n} (dense=True)."""
        out = []
        for c, (herm, recipe) in enumerate(self.choi_inputs):
            k = self.alg.k_factors[c]
            X = DensityMatrix(self.alg, matrix=herm, validate=False)
            Y = fn(self.to_dense_state(X)).matrix[:, 0] if dense else self.dense_matrix(fn(X))
            J = torch.zeros(self.D, k, self.D, k, dtype=Y.dtype)
            for (a, b), terms in recipe:
                J[:, a, :, b] = sum(w * Y[i] for i, w in terms)
            out.append(J.reshape(self.D * k, self.D * k) / k)
        return out

    def dense_check(self):
        """Max deviation between the true outputs computed in M and an independent Lindblad evolution in M_{2^n}
        (H and the jumps converted by interop.to_dense and interop.channel_to_dense)."""
        Hd = self.dense.operator(self.dense_matrix(self.H)[:, None], is_self_adjoint=True)
        jumps = [self.dense.operator(interop.channel_to_dense(L, self.q, self.q)[:, 0][:, None])
                 for L in self.jumps(RATES)]
        out = dynamics.lindblad_evolve(self.to_dense_state(self.test), Hd, jumps, [0.0, T], method='exact')[-1]
        return float((out.matrix[:, 0] - self.true_test).abs().max())

    # ---- data
    def counts(self, shots, dense, seed):
        """Counts (probes, outcomes) with `shots` per probe (shots=None: the exact probabilities as weights)."""
        p = self.p_dense if dense else self.p
        if shots is None:
            return p.clone()
        gen = torch.Generator().manual_seed(seed)
        povm = self.povm_dense if dense else self.povm
        rho = self.to_dense_state(self.out_train) if dense else self.out_train
        return povm.sample(rho, shots, generator=gen).to(torch.float64)


# ----------------------------------------------------------------------
# models
# ----------------------------------------------------------------------
class RateModel(nn.Module):
    """The physical noise model with unknown coupling J and rates (log-parametrised): four real parameters."""

    def __init__(self, prob, gen):
        super().__init__()
        self.prob = prob
        u = torch.rand(4, dtype=torch.float64, generator=gen)
        init = torch.cat([0.5 + u[:1], math.log(1e-3) + u[1:] * (math.log(0.3) - math.log(1e-3))])
        self.theta = nn.Parameter(init)

    def rates(self):
        return {k: torch.exp(self.theta[1 + i]) for i, k in enumerate(KINDS)}

    def generator_matrix(self):
        r = self.rates()
        G = self.prob.G_phys
        return self.theta[0] * G[0] + sum(r[k] * G[1 + i] for i, k in enumerate(KINDS))

    def forward(self, rho):
        return self.prob.evolve(rho, self.theta[0], self.rates())


class MatrixPOVM:
    """A POVM of the module measurement as one matrix on the vectorised states (Tr rho E_j = sum_ab rho_ab (E_j)_ba
    over the active entries): the same probabilities and log-likelihood as the POVM, by one matrix product, for the
    training loop (learn.negative_log_likelihood calls its log_likelihood)."""

    def __init__(self, povm, mask):
        self.mask = mask
        self.A = povm.effects[0].transpose(-2, -1)[:, mask].transpose(0, 1).contiguous()     # (entries, outcomes)

    def probabilities(self, rho):
        return (rho.matrix[:, self.mask] @ self.A.to(rho.matrix.dtype)).real

    def log_likelihood(self, rho, counts, eps=1e-12):
        return (counts * torch.log(self.probabilities(rho).clamp(min=eps))).sum(dim=-1)


def tridiagonal(C):
    i = torch.arange(C)
    return (i[:, None] - i[None, :]).abs() == 1                                  # transitions[d, c]: N -> N +- 1


def build(model, prob, gen, h_scale=H_SCALE):
    """A fresh model with random parameters. The GKSL layer's Hamiltonian is rescaled by h_scale after construction
    (0.1: a start near the identity channel; 1: the library's default scale)."""
    n, alg = prob.n, prob.alg
    if model == 'rates':
        return RateModel(prob, gen)
    if model == 'lindblad':
        net = learn.LindbladLayer(alg, time=T, n_jumps=n - 1, n_inter=n, transitions=tridiagonal(alg.C),
                                  method='exact', generator=gen)
        with torch.no_grad():
            net.X.mul_(h_scale)
        return net
    if model == 'algebra':
        return learn.ChannelLayer(alg, kraus_rank=alg.k_max ** 2, algebra_out=alg, generator=gen)
    if model == 'dense':
        return learn.ChannelLayer(prob.dense, kraus_rank=prob.D ** 2, generator=gen)
    raise ValueError(model)


def train_outputs(model, net, prob):
    """Outputs of the training probes (a DensityMatrix of M, or of M_{2^n} for the dense model), computed by fast
    paths that agree with net(...) (see check_fast_paths)."""
    if model == 'rates':
        return prob.propagate(net.generator_matrix(), prob.R_train)
    if model == 'lindblad':
        return prob.propagate(prob.superop(*net.generator()), prob.R_train)
    if model == 'algebra':                    # Kraus operators applied to the pure probe vectors, sector by sector
        K = net.channel().kraus[0]                                             # (r, D, C, m, k)
        alg = prob.alg
        out = torch.zeros(prob.n_train, alg.C, alg.k_max, alg.k_max, dtype=K.dtype)
        for c, idx, V in prob.groups:
            W = torch.einsum('rdmk,pk->prdm', K[:, :, c, :, :V.shape[1]], V.to(K.dtype))
            out = out.index_add(0, idx, torch.einsum('prdm,prdn->pdmn', W, W.conj()))
        return DensityMatrix(alg, matrix=out, validate=False)
    K = net.channel().kraus[0, :, 0]                                          # (r, D, D)
    W = torch.einsum('rij,pj->pri', K, prob.dense_vecs.to(K.dtype))
    return DensityMatrix(prob.dense, matrix=torch.einsum('pri,prj->pij', W, W.conj())[:, None], validate=False)


def check_fast_paths(prob):
    """Max deviation between the fast training paths and the layers' forward passes, at a random parameter point."""
    dev = 0.0
    gen = torch.Generator().manual_seed(5)
    with torch.no_grad():
        for model in MODELS:
            net = build(model, prob, gen)
            fast = train_outputs(model, net, prob).matrix
            rho = prob.to_dense_state(prob.train) if model == 'dense' else prob.train
            dev = max(dev, float((fast - net(rho).matrix).abs().max()))
        for P, F, rho in ((prob.povm, prob.fast_povm, prob.out_train),
                          (prob.povm_dense, prob.fast_povm_dense, prob.to_dense_state(prob.out_train))):
            dev = max(dev, float((P.probabilities(rho) - F.probabilities(rho)).abs().max()))
    return dev


def fit(model, prob, counts, steps, seed, h_scale=H_SCALE):
    """Maximum likelihood: steps[0] steps of Adam (cosine learning-rate decay), then at most steps[1] iterations of
    L-BFGS (strong Wolfe line search). Loss: the negative log-likelihood per shot minus its value at the empirical
    frequencies, i.e. the mean Kullback-Leibler divergence of the predicted from the observed outcome distributions."""
    gen = torch.Generator().manual_seed(seed)
    net = build(model, prob, gen, h_scale)
    povm = prob.fast_povm_dense if model == 'dense' else prob.fast_povm
    N = counts.sum(-1)
    nll_min = -(counts * torch.log((counts / N[:, None]).clamp(min=1e-300))).sum()
    total = N.sum()

    def objective():
        return (learn.negative_log_likelihood(povm, train_outputs(model, net, prob), counts).sum() - nll_min) / total
    t0 = time.time()
    adam, lbfgs = steps
    opt = torch.optim.Adam(net.parameters(), lr=LR[model])
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, adam, eta_min=LR[model] / 100)
    history = []
    for step in range(adam):
        loss = objective()
        opt.zero_grad()
        loss.backward()
        opt.step()
        sched.step()
        if step % 25 == 0:
            history.append(float(loss.detach()))
    if lbfgs:
        opt = torch.optim.LBFGS(net.parameters(), lr=1.0, max_iter=lbfgs, history_size=50,
                                line_search_fn='strong_wolfe', tolerance_grad=1e-12, tolerance_change=1e-14)

        def closure():
            opt.zero_grad()
            loss = objective()
            loss.backward()
            history.append(float(loss.detach()))
            return loss
        opt.step(closure)
    with torch.no_grad():
        final = float(objective())
    history.append(final)
    return net, final, history[::max(1, len(history) // 60)] + [final], time.time() - t0


# ----------------------------------------------------------------------
# evaluation
# ----------------------------------------------------------------------
def evaluate(model, net, prob):
    """Held-out errors with the layers' own forward passes: mean trace distance and infidelity of the outputs, and
    the trace distance of the normalised Choi matrices of the restrictions to the input sectors (mean and max)."""
    dense = model == 'dense'
    with torch.no_grad():
        pred = net(prob.to_dense_state(prob.test)).matrix[:, 0] if dense else prob.dense_matrix(net(prob.test))
        true = prob.true_test
        td = 0.5 * _trace_norm_herm(pred - true)
        A = DensityMatrix(prob.dense, matrix=pred[:, None], validate=False)
        B = DensityMatrix(prob.dense, matrix=true[:, None], validate=False)
        inf = learn.infidelity(A, B).clamp(min=0)
        off = 0.5 * _trace_norm_herm(pred * ~prob.block_mask_dense)       # coherences between the sectors
        J = prob.choi(net, dense)
        choi = torch.stack([0.5 * _trace_norm_herm(a - b) for a, b in zip(J, prob.true_choi)])
    return dict(trace_distance=float(td.mean()), infidelity=float(inf.mean()), choi=float(choi.mean()),
                choi_max=float(choi.max()), off_sector=float(off.mean()))


def project_rates(prob, net):
    """Least-squares fit of the learned generator by J (-i[H_XY, .]) + sum_k g_k sum_i D[L_{k,i}] with real
    coefficients: the coefficients and the relative norm of the residual (the part of the learned generator outside
    the physical span)."""
    with torch.no_grad():
        G = prob.superop(*net.generator()).reshape(-1)
        A = torch.stack([b.reshape(-1) for b in prob.G_phys], dim=1)
        Ar, Gr = torch.cat([A.real, A.imag]), torch.cat([G.real, G.imag])
        coef = torch.linalg.lstsq(Ar, Gr[:, None]).solution[:, 0]
        res = float(torch.linalg.vector_norm(Ar @ coef - Gr) / torch.linalg.vector_norm(Gr))
    c = coef.tolist()
    return dict(J=c[0], down=c[1], up=c[2], phi=c[3], residual=res)


def rate_fisher(prob):
    """Fisher information per shot of (J, g_down, g_up, g_phi) at the true values, summed over the probes (equal
    shots per probe): measurement.fisher_information of the exact outcome probabilities."""
    def probs(theta):
        G = theta[0] * prob.G_phys[0] + sum(theta[1 + i] * prob.G_phys[1 + i] for i in range(3))
        return prob.fast_povm.probabilities(prob.propagate(G, prob.R_train))
    theta = torch.tensor([J_TRUE] + [RATES[k] for k in KINDS], dtype=torch.float64)
    return fisher_information(probs, [theta]).sum(0)                          # (4, 4)


# ----------------------------------------------------------------------
# tasks (run in worker processes)
# ----------------------------------------------------------------------
_PROBLEMS = {}


def _problem(n):
    if n not in _PROBLEMS:
        _PROBLEMS[n] = Problem(n)
    return _PROBLEMS[n]


def _worker_init():
    torch.set_num_threads(1)


def data_seed(n, shots, rep):
    return 7919 * n + 104729 * rep + (0 if shots is None else int(shots))


def run_task(task):
    """A task for a worker: ('fit', n, model, shots, repeat, starts, steps, h_scale) fits one data set with `starts`
    random starts and evaluates the best; ('check', n) and ('fisher', n) compute the consistency checks and the Fisher
    information of the rate model."""
    kind, n = task[0], task[1]
    prob = _problem(n)
    t0 = time.time()
    if kind == 'check':
        k = prob.alg.k_factors
        return dict(kind=kind, n=n, fast_path=check_fast_paths(prob), dense_check=prob.dense_check(),
                    probes=prob.n_train, test=prob.n_test, settings_alg=prob.settings_alg,
                    settings_dense=prob.settings_dense, dim_alg=prob.dim, dim_dense=prob.D ** 2,
                    set_alg=prob.dim ** 2 - prob.dim, set_dense=prob.D ** 4 - prob.D ** 2,
                    identifiable_dense=prob.dim * prob.D ** 2 - prob.dim, k_factors=list(k),
                    outcomes_alg=len(prob.povm), outcomes_dense=len(prob.povm_dense), seconds=time.time() - t0)
    if kind == 'fisher':
        F = rate_fisher(prob)
        return dict(kind=kind, n=n, fisher=F.tolist(), seconds=time.time() - t0)
    _, n, model, shots, rep, starts, steps, h_scale = task
    seed = data_seed(n, shots, rep)
    counts = prob.counts(shots, model == 'dense', seed + (1 if model == 'dense' else 0))
    best, losses = None, []
    for s in range(starts):
        net, loss, hist, _ = fit(model, prob, counts, steps, 31 * seed + 17 * MODELS.index(model) + s, h_scale)
        losses.append(loss)
        if best is None or loss < best[1]:
            best = (net, loss, hist)
    seconds = time.time() - t0
    net, loss, hist = best
    res = dict(kind=kind, n=n, model=model, shots=shots, repeat=rep, h_scale=h_scale, loss=loss, losses=losses,
               history=hist, seconds=seconds, starts=starts, params=learn.count_parameters(net),
               **evaluate(model, net, prob))
    if model == 'rates':
        res['rates'] = dict(J=net.theta[0].item(), **{k: v.item() for k, v in net.rates().items()})
    elif model == 'lindblad':
        res['rates'] = project_rates(prob, net)
    return res


# ----------------------------------------------------------------------
# output
# ----------------------------------------------------------------------
def shots_key(shots):
    return 'exact' if shots is None else str(int(shots))


def figure(main, spread, nmax, repeat_shots, repeats, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.edgecolor': GREY, 'axes.labelcolor': INK, 'xtick.color': GREY,
                         'ytick.color': GREY, 'text.color': INK})
    finite = [s for s in SHOTS if s is not None]
    xs = finite + [10 * finite[-1]]                                             # 'exact' drawn one decade right
    floor = 1.5e-4
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw={'width_ratios': [1.25, 1]})
    ns = sorted({r['n'] for r in main})
    for model in MODELS:
        for n in ns:
            ys = [next(r['trace_distance'] for r in main if r['n'] == n and r['model'] == model and r['shots'] == s)
                  for s in SHOTS]
            solid = n == nmax
            ax.plot(xs[:-1], ys[:-1], '-' if solid else '--', color=COLORS[model], lw=1.5 if solid else 0.9,
                    marker='o' if solid else None, ms=3.5, alpha=1 if solid else 0.6,
                    label=LABELS[model] if solid else None)
            if solid:                                          # exact data: the optimisation floor (clipped)
                ax.plot(xs[-1] * OFFSET[model], max(ys[-1], floor), marker='o' if ys[-1] >= floor else 'v', color=COLORS[model],
                        ms=4, ls='none')
    ref = np.array(finite, dtype=float)
    ax.plot(ref, 2.5e-2 * (ref / ref[0]) ** -0.5, ':', color=GREY, lw=0.9)
    ax.text(ref[2], 2.5e-2 * (ref[2] / ref[0]) ** -0.5 * 1.3, r'$\propto 1/\sqrt{\mathrm{shots}}$', color=GREY,
            fontsize=7)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 1)
    ax.set_xticks(xs)
    ax.set_xticklabels([r'$10^{%d}$' % round(math.log10(s)) for s in finite] + ['exact'])
    ax.axvline(math.sqrt(finite[-1] * xs[-1]), color='#e6e5e1', lw=0.8)
    ax.set_xlabel('shots per probe state')
    ax.set_ylabel('held-out trace distance')
    ax.set_title(f'(a) prediction error, $n={nmax}$ (dashed: $n={min(ns)}$)', fontsize=8, loc='left')
    ax.spines[['top', 'right']].set_visible(False)
    # (b) learned parameters over independent data sets
    sp = spread[nmax]
    keys = ('J',) + KINDS
    true = dict(J=J_TRUE, **RATES)
    names = {'J': '$J$', 'down': r'$\gamma_\downarrow$', 'up': r'$\gamma_\uparrow$', 'phi': r'$\gamma_\varphi$'}
    for i, k in enumerate(keys):
        cr = sp['cramer_rao'][k] / true[k]
        bx.add_patch(plt.Rectangle((i - 0.3, 1 - cr), 0.6, 2 * cr, color='#e6e5e1', lw=0))
    bx.axhline(1, color=GREY, lw=0.8)
    for model, dx, filled in (('rates', -0.12, True), ('lindblad', 0.12, False)):
        m = [sp[model][k]['mean'] / true[k] for k in keys]
        e = [sp[model][k]['std'] / true[k] for k in keys]
        bx.errorbar(np.arange(len(keys)) + dx, m, yerr=e, fmt='o', color=COLORS[model], ms=4, capsize=2, lw=1,
                    mfc=COLORS[model] if filled else 'white', label=LABELS[model] if model == 'rates'
                    else 'GKSL layer, projected')
    bx.add_patch(plt.Rectangle((0, 0), 0, 0, color='#e6e5e1', label='Cramér–Rao $\\pm\\sigma$ (rate model)'))
    bx.set_xticks(range(len(keys)))
    bx.set_xticklabels([names[k] for k in keys])
    bx.set_xlim(-0.6, len(keys) - 0.4)
    bx.set_ylabel('learned / true')
    bx.set_title(f'(b) parameters, $n={nmax}$, $10^{{{round(math.log10(repeat_shots))}}}$ shots, '
                 f'{repeats} data sets', fontsize=8, loc='left')
    bx.legend(frameon=False, fontsize=6.5, loc='lower left')
    bx.spines[['top', 'right']].set_visible(False)
    handles = [plt.Line2D([], [], color=COLORS[m], marker='o', ms=3.5, lw=1.5, label=LABELS[m]) for m in MODELS]
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, fontsize=7)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(path, metadata={'CreationDate': None})
    plt.close(fig)


def sci(x, digits=1):
    if x == 0:
        return '0'
    m, e = f'{x:.{digits}e}'.split('e')
    e = int(e)
    return m if e == 0 else f'{m}{{\\times}}10^{{{e}}}'


WORD = {100: 'Two', 1000: 'Three', 10000: 'Four', 100000: 'Five', None: 'Exact'}
MACRO = {'rates': 'Rates', 'lindblad': 'Lindblad', 'algebra': 'Algebra', 'dense': 'Dense'}
NAME = {'rates': 'rates', 'lindblad': 'GKSL layer', 'algebra': r'channel of $\M$', 'dense': r'channel of $M_{2^n}$'}


def shots_to_reach(main, n, model, target):
    """Shots per probe at which the held-out trace distance reaches `target` (log-log interpolation of the curve;
    None if it does not within the budgets)."""
    finite = [s for s in SHOTS if s is not None]
    ys = [next(r['trace_distance'] for r in main if r['n'] == n and r['model'] == model and r['shots'] == s)
          for s in finite]
    if ys[0] <= target:
        return float(finite[0])
    for (s0, y0), (s1, y1) in zip(zip(finite, ys), zip(finite[1:], ys[1:])):
        if y1 <= target < y0:
            t = (math.log(target) - math.log(y0)) / (math.log(y1) - math.log(y0))
            return float(math.exp(math.log(s0) + t * (math.log(s1) - math.log(s0))))
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dims', default='3,4')
    ap.add_argument('--adam', type=int, default=300, help='Adam steps per start')
    ap.add_argument('--lbfgs', type=int, default=400, help='L-BFGS iterations per start')
    ap.add_argument('--starts', default='rates:3,lindblad:3,algebra:2,dense:2')
    ap.add_argument('--repeats', type=int, default=6, help='independent data sets for the spread of the rates')
    ap.add_argument('--repeat-shots', type=int, default=10000)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--out', default=str(OUT / 'generated'))
    ap.add_argument('--figdir', default=str(OUT / 'figures'))
    ap.add_argument('--replot', default=None, help='only redraw the figure and rewrite the macros from this JSON')
    args = ap.parse_args()
    if args.replot:
        payload = json.loads(Path(args.replot).read_text())
        for key in ('checks', 'fisher_per_shot', 'spread'):                  # JSON object keys are strings
            payload[key] = {int(n): v for n, v in payload[key].items()}
        report(payload, args.out, args.figdir)
        return
    t_start = time.time()
    dims = [int(x) for x in args.dims.split(',')]
    nmax = max(dims)
    starts = {kv.split(':')[0]: int(kv.split(':')[1]) for kv in args.starts.split(',')}
    steps = (args.adam, args.lbfgs)
    tasks = [('check', n) for n in dims] + [('fisher', n) for n in dims]
    tasks += [('fit', n, m, s, 0, starts[m], steps, H_SCALE) for n in dims for m in MODELS for s in SHOTS]
    tasks += [('fit', n, m, args.repeat_shots, rep, 2, steps, H_SCALE) for n in dims for m in ('rates', 'lindblad')
              for rep in range(1, args.repeats)]
    tasks += [('fit', nmax, 'lindblad', None, 0, 3, steps, 1.0)]                 # the library's default H scale
    cost = {'dense': 4, 'algebra': 2, 'lindblad': 2, 'rates': 1}
    order = sorted(range(len(tasks)), key=lambda i: -(tasks[i][1] ** 2) * (cost[tasks[i][2]] * tasks[i][5]
                                                                            if tasks[i][0] == 'fit' else 1))
    results = [None] * len(tasks)
    with ProcessPoolExecutor(max_workers=args.workers, initializer=_worker_init) as pool:
        futures = {pool.submit(run_task, tasks[i]): i for i in order}
        done = 0
        for fut in as_completed(futures):
            i = futures[fut]
            results[i] = fut.result()
            done += 1
            r = results[i]
            if r['kind'] == 'fit':
                print(f"[{done}/{len(tasks)}] n={r['n']} {r['model']:8s} shots={shots_key(r['shots']):>6s} "
                      f"rep={r['repeat']} h={r['h_scale']}: loss {r['loss']:.2e}, trace distance "
                      f"{r['trace_distance']:.2e}, Choi {r['choi']:.2e}, {r['seconds']:.0f} s", flush=True)
            else:
                print(f"[{done}/{len(tasks)}] {r['kind']} n={r['n']} ({r['seconds']:.0f} s)", flush=True)
    payload = dict(environment=environment('cpu'), runtime_seconds=time.time() - t_start, workers=args.workers,
                   setup=dict(T=T, J=J_TRUE, rates=RATES, shots=[shots_key(s) for s in SHOTS], dims=dims,
                              adam=args.adam, lbfgs=args.lbfgs, starts=starts, lr=LR, h_scale=H_SCALE,
                              repeats=args.repeats, repeat_shots=args.repeat_shots),
                   checks={r['n']: r for r in results if r['kind'] == 'check'},
                   fisher_per_shot={r['n']: r['fisher'] for r in results if r['kind'] == 'fisher'},
                   fits=[r for r in results if r['kind'] == 'fit'])
    payload['spread'] = spread(payload)
    gen = Path(args.out)
    gen.mkdir(parents=True, exist_ok=True)
    (gen / 'channel_learning.json').write_text(json.dumps(payload, indent=1))
    report(payload, args.out, args.figdir)
    print(f"wrote figures/channel_learning.pdf and generated/channel_learning.{{json,tex}}; "
          f"{payload['runtime_seconds']:.0f} s")


def spread(payload):
    """Mean and standard deviation of the learned rates over the independent data sets, and the Cramer-Rao standard
    deviations of the rate model (inverse Fisher information) for the same number of shots."""
    setup, out = payload['setup'], {}
    for n in setup['dims']:
        cov = np.linalg.inv(np.array(payload['fisher_per_shot'][n])) / setup['repeat_shots']
        out[n] = {}
        for m in ('rates', 'lindblad'):
            rs = [r for r in payload['fits'] if r['n'] == n and r['model'] == m
                  and r['shots'] == setup['repeat_shots'] and r['h_scale'] == H_SCALE]
            out[n][m] = {k: dict(mean=float(np.mean([r['rates'][k] for r in rs])),
                                 std=float(np.std([r['rates'][k] for r in rs], ddof=1)), samples=len(rs))
                         for k in ('J',) + KINDS}
            if m == 'lindblad':
                out[n][m]['residual'] = float(np.mean([r['rates']['residual'] for r in rs]))
        out[n]['cramer_rao'] = {k: float(math.sqrt(cov[i, i])) for i, k in enumerate(('J',) + KINDS)}
    return out


def report(payload, outdir, figdir):
    """Figure and TeX macros from the results (also from a saved JSON: --replot)."""
    setup = payload['setup']
    dims = setup['dims']
    nmax = max(dims)
    starts = setup['starts']
    checks = payload['checks']
    runtime = payload['runtime_seconds']
    sp_all = payload['spread']
    fits = payload['fits']
    main_fits = [r for r in fits if r['repeat'] == 0 and r['h_scale'] == H_SCALE]
    init_study = [r for r in fits if r['h_scale'] != H_SCALE][0]
    start_spread = max((max(r['losses']) - min(r['losses'])) / min(r['losses']) for r in fits
                       if r['shots'] is not None and r['h_scale'] == H_SCALE)
    Path(figdir).mkdir(parents=True, exist_ok=True)
    figure(main_fits, sp_all, nmax, setup['repeat_shots'], setup['repeats'], Path(figdir) / 'channel_learning.pdf')
    gen = Path(outdir)
    # ---- TeX
    get = lambda n, m, s: next(r for r in main_fits if r['n'] == n and r['model'] == m and r['shots'] == s)
    rows = []
    for n in dims:
        for m in MODELS:
            tds = ' & '.join(f'${sci(get(n, m, s)["trace_distance"])}$' for s in SHOTS)
            secs = sum(get(n, m, s)['seconds'] for s in SHOTS) / len(SHOTS) / starts[m]
            rows.append(f'{n} & {NAME[m]} & {get(n, m, None)["params"]} & {tds} & '
                        f'${sci(get(n, m, 10000)["choi"])}$ & {secs:.0f}\\\\')
        if n != dims[-1]:
            rows.append('\\midrule')
    c = checks[nmax]
    tex = ['% generated by scripts/channel_learning.py -- do not edit by hand',
           '\\newcommand{\\CLRows}{%', *rows, '}',
           f'\\newcommand{{\\CLT}}{{{T:g}}}',
           f'\\newcommand{{\\CLGdown}}{{{RATES["down"]:g}}}',
           f'\\newcommand{{\\CLGup}}{{{RATES["up"]:g}}}',
           f'\\newcommand{{\\CLGphi}}{{{RATES["phi"]:g}}}',
           f'\\newcommand{{\\CLnMin}}{{{min(dims)}}}',
           f'\\newcommand{{\\CLnMax}}{{{nmax}}}',
           f'\\newcommand{{\\CLProbes}}{{{c["probes"]}}}',
           f'\\newcommand{{\\CLProbesMin}}{{{checks[min(dims)]["probes"]}}}',
           f'\\newcommand{{\\CLTest}}{{{c["test"]}}}',
           f'\\newcommand{{\\CLSettingsAlg}}{{{c["settings_alg"]}}}',
           f'\\newcommand{{\\CLSettingsDense}}{{{c["settings_dense"]}}}',
           f'\\newcommand{{\\CLDimAlg}}{{{c["dim_alg"]}}}',
           f'\\newcommand{{\\CLDimDense}}{{{c["dim_dense"]}}}',
           f'\\newcommand{{\\CLSetAlg}}{{{c["set_alg"]}}}',
           f'\\newcommand{{\\CLSetDense}}{{{c["set_dense"]}}}',
           f'\\newcommand{{\\CLIdentDense}}{{{c["identifiable_dense"]}}}',
           f'\\newcommand{{\\CLDenseCheck}}{{{sci(max(r["dense_check"] for r in checks.values()))}}}',
           f'\\newcommand{{\\CLFastPath}}{{{sci(max(r["fast_path"] for r in checks.values()))}}}',
           f'\\newcommand{{\\CLAdam}}{{{setup['adam']}}}',
           f'\\newcommand{{\\CLLbfgs}}{{{setup['lbfgs']}}}',
           f'\\newcommand{{\\CLStartsRates}}{{{starts["rates"]}}}',
           f'\\newcommand{{\\CLStartsLindblad}}{{{starts["lindblad"]}}}',
           f'\\newcommand{{\\CLStartsAlgebra}}{{{starts["algebra"]}}}',
           f'\\newcommand{{\\CLStartsDense}}{{{starts["dense"]}}}',
           f'\\newcommand{{\\CLRepeats}}{{{setup['repeats']}}}',
           f'\\newcommand{{\\CLRuntime}}{{{runtime / 60:.0f}}}',
           f'\\newcommand{{\\CLWorkers}}{{{payload['workers']}}}',
           f'\\newcommand{{\\CLHScale}}{{{H_SCALE:g}}}',
           f'\\newcommand{{\\CLInitConverged}}{{{sum(l < 1e-6 for l in init_study["losses"])}}}',
           f'\\newcommand{{\\CLInitStarts}}{{{len(init_study["losses"])}}}',
           f'\\newcommand{{\\CLInitLoss}}{{{sci(max(init_study["losses"]))}}}',
           f'\\newcommand{{\\CLStartSpread}}{{{sci(start_spread)}}}']
    for m in MODELS:
        M = MACRO[m]
        tex.append(f'\\newcommand{{\\CLPar{M}}}{{{get(nmax, m, None)["params"]}}}')
        tex.append(f'\\newcommand{{\\CLParMin{M}}}{{{get(min(dims), m, None)["params"]}}}')
        for s in SHOTS:
            r = get(nmax, m, s)
            tex.append(f'\\newcommand{{\\CLTD{M}{WORD[s]}}}{{{sci(r["trace_distance"])}}}')
            tex.append(f'\\newcommand{{\\CLInf{M}{WORD[s]}}}{{{sci(r["infidelity"])}}}')
            tex.append(f'\\newcommand{{\\CLChoi{M}{WORD[s]}}}{{{sci(r["choi"])}}}')
            tex.append(f'\\newcommand{{\\CLTDMin{M}{WORD[s]}}}{{{sci(get(min(dims), m, s)["trace_distance"])}}}')
        tex.append(f'\\newcommand{{\\CLOff{M}}}{{{sci(get(nmax, m, 10000)["off_sector"])}}}')
        tex.append(f'\\newcommand{{\\CLSec{M}}}{{{get(nmax, m, 10000)["seconds"] / starts[m]:.0f}}}')
    finite = [x for x in SHOTS if x is not None]
    for m in MODELS:                                   # log-log slope of the error over the finite shot budgets
        slope = -np.polyfit(np.log10(finite), np.log10([get(nmax, m, x)['trace_distance'] for x in finite]), 1)[0]
        tex.append(f'\\newcommand{{\\CLSlope{MACRO[m]}}}{{{slope:.2f}}}')
    ratio = get(nmax, 'algebra', None)['params'] / get(nmax, 'lindblad', None)['params']
    tex.append(f'\\newcommand{{\\CLParRatio}}{{{ratio:.0f}}}')
    target = get(nmax, 'algebra', 1000)['trace_distance']
    need = {m: shots_to_reach(main_fits, nmax, m, target) for m in MODELS}
    tex.append(f'\\newcommand{{\\CLTarget}}{{{sci(target)}}}')
    for m in MODELS:
        v = need[m]
        tex.append(f'\\newcommand{{\\CLNeed{MACRO[m]}}}{{{"--" if v is None else sci(v)}}}')
    ratio = None if need['dense'] is None else need['dense'] / need['algebra']
    tex.append(f'\\newcommand{{\\CLNeedRatio}}{{{"--" if ratio is None else f"{ratio:.0f}"}}}')
    sp = sp_all[nmax]
    for k, K in (('J', 'J'), ('down', 'Down'), ('up', 'Up'), ('phi', 'Phi')):
        for m, M in (('rates', 'Rates'), ('lindblad', 'Lindblad')):
            tex.append(f'\\newcommand{{\\CL{M}{K}}}{{{sp[m][k]["mean"]:.4f}}}')
            tex.append(f'\\newcommand{{\\CL{M}{K}Sd}}{{{sp[m][k]["std"]:.4f}}}')
        tex.append(f'\\newcommand{{\\CLCR{K}}}{{{sp["cramer_rao"][k]:.4f}}}')
    tex.append(f'\\newcommand{{\\CLRepeatShots}}{{{sci(setup['repeat_shots'], 0)}}}')
    tex.append(f'\\newcommand{{\\CLResidual}}{{{sci(sp["lindblad"]["residual"])}}}')
    (gen / 'channel_learning.tex').write_text('\n'.join(tex) + '\n')




if __name__ == '__main__':
    main()
