#!/usr/bin/env python3
"""
Measured effect of the design decisions of torch_vn_algebra described in man/article, Section "Architecture": each row
times the same result computed with and without one mechanism (the library path against the naive path), on the
device of the run (CPU here; the GPU columns of the performance section are separate).

  known spectrum      lambda_max of an operator built from its eigenvalues: no decomposition vs. eigvalsh
  eigenvalues only    eigenvalues() (eigvalsh) vs. eigh() with eigenvectors
  spectral cache      exp(-i t A) for 50 times t: one eigendecomposition (cache_spectrum) vs. 50
  size groups         eigvalsh of 64 sectors of size 8: one batched call vs. a Python loop over the sectors
  trace of a product  Tr(rho A) in O(k^2) vs. forming rho A (blocks of size 256)
  sparse jumps        Lindblad evolution of 8 qubits with S^-_i, S^z_i: SectorJump vs. dense jumps between sectors

Writes man/article/generated/architecture.{json,tex}.

    python scripts/architecture_bench.py [--repeats 5]
"""
import argparse
import json
import time
from pathlib import Path

import torch

from common import environment                       # device, versions, matmul precision, commit

from torch_vn_algebra import TypeIAlgebra, SpinChain, dynamics, states

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'man' / 'article' / 'generated'


def best(fn, repeats):
    fn()                                                       # warm-up
    ts = []
    for _ in range(repeats):
        t = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t)
    return min(ts)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--repeats', type=int, default=5)
    args = ap.parse_args()
    torch.manual_seed(0)
    R = args.repeats
    rows = []

    alg = TypeIAlgebra(n_factors=[32] * 8, k_factors=[32] * 8, precision='double')
    B = 256

    def sampler(d):
        lam = 0.2 + torch.rand(B, d, dtype=torch.float64)
        return lam
    Xk = alg.operator_from_eigenvalues(sampler, batch_size=B, force_self_adjoint=True, force_positive=True)
    Xk.matrix
    Xu = alg.operator(Xk.matrix.clone(), is_self_adjoint=True)

    def lam(Y):
        Y._lambda_max = Y._lambda_min = None                   # drop the cached value, keep everything else
        return Y.lambda_max
    rows.append(('known spectrum', r'$\lambda_{\max}$, 256 operators in $M_{32}^{\oplus8}$',
                 best(lambda: lam(Xk), R), best(lambda: lam(Xu), R)))
    M = Xk.matrix.clone()

    X = alg.operator(M.clone(), is_self_adjoint=True)
    rows.append(('eigenvalues only', 'spectra of the same operators',
                 best(lambda: X.eigenvalues(), R), best(lambda: X.eigh(), R)))

    ts = torch.linspace(0.1, 5.0, 50, dtype=torch.float64)

    def props(cached):
        Y = alg.operator(M.clone(), is_self_adjoint=True)
        if cached:
            Y.cache_spectrum()
        for t in ts.tolist():
            Y.apply_function(lambda w, t=t: torch.exp(-1j * t * w)).matrix
    rows.append(('spectral cache', r'$e^{-itA}$ for 50 values of $t$', best(lambda: props(True), 2),
                 best(lambda: props(False), 2)))

    alg64 = TypeIAlgebra(n_factors=[8] * 64, k_factors=[8] * 64, precision='double')
    M64 = alg64.operator_from_eigenvalues(lambda d: torch.rand(B, d, dtype=torch.float64), batch_size=B,
                                          force_self_adjoint=True).matrix.clone()
    Y64 = alg64.operator(M64, is_self_adjoint=True)

    def loop():
        return [torch.linalg.eigvalsh(M64[:, c]) for c in range(64)]
    rows.append(('size groups', r'spectra, 256 operators in $M_8^{\oplus64}$', best(lambda: Y64.eigenvalues(), R),
                 best(loop, R)))

    alg2 = TypeIAlgebra(n_factors=[256] * 2, k_factors=[256] * 2, precision='double')
    rho = states.random_density_matrix(alg2, batch_size=32)
    rho.matrix
    X2 = states.random_density_matrix(alg2, batch_size=32)
    X2.matrix
    rows.append(('trace of a product', r'$\Tr(\rho A)$, 32 pairs in $M_{256}^{\oplus2}$',
                 best(lambda: states.trace_of_product(rho, X2), R),
                 best(lambda: torch.diagonal(rho.matrix @ X2.matrix, dim1=-2, dim2=-1).sum((-1, -2)), R)))

    n = 8
    q = SpinChain(n, 'periodic', precision='double')
    H = q.xxz(J=1.0, Delta=0.0)
    rho0 = q.basis_state('11110000')
    sparse = [q.sector_jump(i, 'lower') for i in range(n)] + [q.sector_jump(i, 'sz') for i in range(n)]
    dense = [q.lowering(i) for i in range(n)] + [q.sz(i) for i in range(n)]
    rates = [0.05] * (2 * n)

    def evolve(js):
        return dynamics.lindblad_evolve(rho0, H, js, [0.0, 1.0], rates=rates, substeps=10, method='kraus')[-1].matrix
    diff = float((evolve(sparse) - evolve(dense)).abs().max())
    rows.append(('sparse jumps', r'Lindblad, 8 qubits, $S^-_i$ and $S^z_i$', best(lambda: evolve(sparse), 2),
                 best(lambda: evolve(dense), 2)))
    kd = sum(J.kraus.numel() * J.kraus.element_size() for J in dense[:n])
    ks = sum(J.val.numel() * J.val.element_size() + J.dst.numel() * 8 for J in sparse[:n])

    OUT.mkdir(parents=True, exist_ok=True)
    data = dict(environment=environment('cpu'), device='cpu', torch=torch.__version__, threads=torch.get_num_threads(),
                rows=[dict(name=a, what=b, library=c, naive=d, speedup=d / c) for a, b, c, d in rows],
                jump_bytes_dense=kd, jump_bytes_sparse=ks, jump_max_diff=diff)
    (OUT / 'architecture.json').write_text(json.dumps(data, indent=1))
    lines = [f"{a} & {b} & {1e3 * c:.2f} & {1e3 * d:.2f} & {d / c:.0f}\\\\" if d / c >= 10 else
             f"{a} & {b} & {1e3 * c:.2f} & {1e3 * d:.2f} & {d / c:.1f}\\\\" for a, b, c, d in rows]
    tex = ['% generated by scripts/architecture_bench.py -- do not edit by hand',
           '\\newcommand{\\ArchRows}{%', *lines, '}',
           f'\\newcommand{{\\ArchThreads}}{{{torch.get_num_threads()}}}',
           f'\\newcommand{{\\ArchJumpDense}}{{{kd / 2 ** 20:.1f}}}',
           f'\\newcommand{{\\ArchJumpSparse}}{{{ks / 2 ** 20:.2f}}}',
           ]
    m, e = f'{diff:.0e}'.split('e')
    tex.append(f'\\newcommand{{\\ArchJumpDiff}}{{{m}\\times10^{{{int(e)}}}}}')
    (OUT / 'architecture.tex').write_text('\n'.join(tex) + '\n')
    for a, b, c, d in rows:
        print(f'{a:20s} library {1e3 * c:9.2f} ms   naive {1e3 * d:9.2f} ms   x{d / c:.1f}')
    print(f'jump storage: dense {kd / 2 ** 20:.1f} MiB, sparse {ks / 2 ** 20:.2f} MiB; max difference {diff:.1e}')


if __name__ == '__main__':
    main()
