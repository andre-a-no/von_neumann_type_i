#!/usr/bin/env python3
"""
Section "Number field and precision" of paper: what real / complex and single / double cost and
what they buy, for the same computations in the four configurations

    float32 (real, single)   float64 (real, double)   complex64 (single)   complex128 (double)

  * resources: bytes per stored entry, time per sample of library operations;
  * accuracy:  unitarity of Haar samples, SVD square root, entropy of an almost pure state,
               round-off floor of an inequality that holds with equality;
  * modelling: real symmetric versus complex Hermitian random matrices have different level
               statistics (GOE <r> = 0.5307 versus GUE <r> = 0.5996), i.e. the field changes the model.

Writes generated/numerics.tex, generated/numerics.json and figures/numerics.pdf.
"""
import math
import statistics
import time

import torch

from common import parse_args, save_json, write_tex, env_macro, sci, sync

from torch_vn_algebra import TypeIAlgebra, channels, cost, dynamics, optimize as opt, states

args = parse_args(__doc__)
dev = args.device
FULL = args.mode in ('full', 'check')            # check: the sizes of full, minimal repetitions
CHECK = args.mode == 'check'
CONFIGS = [('float32', False, 'single'), ('float64', False, 'double'),
           ('complex64', True, 'single'), ('complex128', True, 'double')]
K, C = (32, 16) if FULL else (16, 8)
TARGET = 2 ** 22 if FULL else 2 ** 18
B = max(1, TARGET // (C * K * K))
REPEATS = 1 if CHECK else 5 if FULL else 3
cost.set_limits(warn_seconds=float('inf'), probe_flops=float('inf'))


def timeit(fn):
    fn()
    sync(dev)
    ts = []
    for _ in range(REPEATS):
        t0 = time.perf_counter()
        fn()
        sync(dev)
        ts.append(time.perf_counter() - t0)
    return statistics.median(ts) / B


rows = []
for name, cplx, prec in CONFIGS:
    torch.manual_seed(args.seed)
    alg = TypeIAlgebra([K] * C, [K] * C, complex_valued=cplx, precision=prec, device=dev)
    gen = lambda: alg.operator_from_eigenvalues(lambda d: 0.1 + torch.rand(B, d, device=dev), batch_size=B,
                                                force_positive=True, force_self_adjoint=True)
    X = gen()
    X.matrix
    rho = states.random_density_matrix(alg, B)
    Phi = channels.random_channel(alg, 4, B)
    t = {
        'generate': timeit(lambda: gen().matrix),
        'eigh f(A)': timeit(lambda: alg.operator(X.matrix, is_self_adjoint=True).apply_function(torch.sqrt).matrix),
        'SVD |A|': timeit(lambda: alg.operator(X.matrix).abs().matrix),
        'channel': timeit(lambda: Phi(rho).matrix),
        'expm': timeit(lambda: alg.operator(X.matrix).expm(-0.1).matrix),
    }
    Xp, Yp, Up = opt.PositiveParam(alg, B, delta=torch.full((B,), 0.5, device=dev)), opt.PositiveParam(alg, B), \
        opt.UnitaryParam(alg, B)

    def step():
        x, y = Xp.operator(), Yp.operator()
        z = (x @ Up.operator() @ y).trace.abs() - (x @ y).trace.real
        z.sum().backward()
    t['opt. step'] = timeit(step)

    # accuracy ------------------------------------------------------------------------------
    n = 64
    a1 = TypeIAlgebra([n], [n], complex_valued=cplx, precision=prec, device=dev)
    U = a1.random_unitary(n, batch_size=256)
    eye = torch.eye(n, dtype=U.dtype, device=dev)
    unitarity = (U @ U.conj().transpose(-2, -1) - eye).abs().amax().item()
    P = a1.operator_from_eigenvalues(lambda d: torch.rand(256, d, device=dev), batch_size=256, force_positive=True)
    R = P.sqrt().matrix
    sqrt_err = (torch.linalg.matrix_norm(R @ R - P.matrix) / torch.linalg.matrix_norm(P.matrix)).max().item()
    eps = 1e-7                                     # rho = (1 - eps) |0><0| + eps 1/n
    lam = torch.full((n,), eps / n, dtype=torch.float64)
    lam[0] += 1 - eps
    S_exact = -(lam * torch.log(lam)).sum().item()
    rho_np = a1.operator_from_eigenvalues(lambda d: lam.clone(), batch_size=64, force_positive=True)
    S = states.von_neumann_entropy(states.DensityMatrix(a1, matrix=rho_np.matrix, validate=False))
    ent_err = ((S.double() - S_exact).abs().max() / S_exact).item()
    # Exp. 3 with commuting X, Y: z = Tr(Y|X|Y) - Tr|YXY| = 0 exactly, so |z| is pure round-off
    a2 = TypeIAlgebra([K] * C, [K] * C, complex_valued=cplx, precision=prec, device=dev)
    V = a2.random_unitary_operator(512).matrix
    dx = torch.randn(512, C, K, device=dev).to(V.dtype)
    dy = (0.1 + torch.rand(512, C, K, device=dev)).to(V.dtype)
    Xc = a2.operator(V @ torch.diag_embed(dx) @ V.conj().transpose(-2, -1), is_self_adjoint=True)
    Yc = a2.operator(V @ torch.diag_embed(dy) @ V.conj().transpose(-2, -1), is_self_adjoint=True, is_positive=True)
    z = (Yc @ Xc.abs() @ Yc).trace.real - (Yc @ Xc @ Yc).abs().trace.real
    floor = z.abs().max().item()

    bytes_entry = torch.empty((), dtype=alg.dtype).element_size()
    row = dict(config=name, bytes_per_entry=bytes_entry, times=t, unitarity=unitarity, sqrt_error=sqrt_err,
               entropy_rel_error=ent_err, inequality_floor=floor,
               operator_MB=alg.bytes_per_operator(B) / 2 ** 20)
    rows.append(row)
    print(f"{name:10s} {bytes_entry:2d} B/entry | " + " ".join(f"{k} {1e6 * v:8.1f}us" for k, v in t.items())
          + f" | unitarity {unitarity:.1e} sqrt {sqrt_err:.1e} entropy {ent_err:.1e} z-floor {floor:.1e}")

# modelling: GOE versus GUE ----------------------------------------------------------------------
NM, SM = (200, 400) if FULL else (100, 100)
r_vals = {}
for name, cplx in (('real symmetric (GOE)', False), ('complex Hermitian (GUE)', True)):
    alg = TypeIAlgebra([NM], [NM], complex_valued=cplx, precision='double', device=dev)
    G = torch.randn(SM, NM, NM, dtype=alg.dtype, device=dev)
    H = (G + G.conj().transpose(-2, -1)) / 2
    w = torch.linalg.eigvalsh(H)
    s = torch.diff(w[:, NM // 4: 3 * NM // 4], dim=-1)
    r = (torch.minimum(s[:, :-1], s[:, 1:]) / torch.maximum(s[:, :-1], s[:, 1:])).mean(-1)
    r_vals[name] = (r.mean().item(), (r.std() / math.sqrt(SM)).item())
    print(f"{name}: <r> = {r_vals[name][0]:.4f} +- {r_vals[name][1]:.4f}")

save_json(args, 'numerics', dict(rows=rows, K=K, C=C, batch=B, r=r_vals))
ref = rows[0]['times']
tex = env_macro(args, 'Num')
tex += f"\\newcommand{{\\NumK}}{{{K}}}\n\\newcommand{{\\NumC}}{{{C}}}\n\\newcommand{{\\NumB}}{{{B}}}\n"
tex += (f"\\newcommand{{\\NumRGOE}}{{{r_vals['real symmetric (GOE)'][0]:.4f}}}\n"
        f"\\newcommand{{\\NumRGUE}}{{{r_vals['complex Hermitian (GUE)'][0]:.4f}}}\n")
tex += "\\newcommand{\\NumTimeRows}{%\n"
for r in rows:
    tex += (f"\\texttt{{{r['config']}}} & {r['bytes_per_entry']} & {r['operator_MB']:.3g} & "
            + " & ".join(f"{r['times'][k] / ref[k]:.2g}" for k in ref) + " \\\\\n")
tex += "}\n\\newcommand{\\NumAccRows}{%\n"
for r in rows:
    tex += (f"\\texttt{{{r['config']}}} & ${sci(r['unitarity'])}$ & ${sci(r['sqrt_error'])}$ & "
            f"${sci(r['entropy_rel_error'])}$ & ${sci(max(r['inequality_floor'], 1e-300))}$ \\\\\n")
tex += "}\n"
tex += "\\newcommand{\\NumTimeHeader}{" + " & ".join(ref.keys()) + "}\n"
tex += "\\newcommand{\\NumRefTimes}{" + ", ".join(f"{k}: {1e6 * v:.3g}\\,$\\mu$s" for k, v in ref.items()) + "}\n"
write_tex(args, 'numerics', tex)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(7, 3))
ops = list(ref.keys())
w = 0.2
for j, r in enumerate(rows):
    ax.bar([i + j * w for i in range(len(ops))], [r['times'][k] / ref[k] for k in ops], w, label=r['config'])
ax.set_xticks([i + 1.5 * w for i in range(len(ops))])
ax.set_xticklabels(ops, fontsize=8)
ax.set_ylabel('time relative to float32')
ax.set_yscale('log', base=2)
ax.axhline(1, color='k', lw=0.6)
ax.legend(frameon=False, fontsize=7, ncol=4)
ax.set_title(f'$k={K}$, $C={C}$, batch {B}', fontsize=9)
fig.tight_layout()
fig.savefig(f"{args.figdir}/numerics.pdf", metadata={"CreationDate": None})
print("written:", args.out, args.figdir)
