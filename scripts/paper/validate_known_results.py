#!/usr/bin/env python3
"""
Section "Validation against known results" of paper/v2.

Every check compares a Monte Carlo estimate obtained with the library to an exact formula:

  1. Haar moments of U(n) and O(n) (Weingarten calculus)
  2. E tau((A U B U*)^2) = -1/(N^2 - 1) for traceless involutions A, B and Haar U (Weingarten)
  3. Nearest-neighbour spacings of COE / CUE / CSE eigenphases vs. the Wigner surmise
  4. Page's formula for the mean entanglement entropy of random pure states (tensor products,
     partial trace)
  5. Mean purity of induced random states, E Tr rho^2 = (k + r) / (k r + 1)
  6. Spontaneous emission (Lindblad): exp(tL) and RK4 vs the analytic solution

Writes generated/validation.tex (table rows + macros), generated/validation.json and
figures/spacings.pdf, figures/page.pdf.

    python scripts/paper/validate_known_results.py --mode full --device cuda
"""
import math

import numpy as np
import torch

from common import parse_args, save_json, write_tex, env_macro, sci, fmt

from torch_vn_algebra import TypeIAlgebra, DensityMatrix, dynamics, states, tensor_product, partial_trace

args = parse_args(__doc__)
dev = args.device
FULL = args.mode == 'full'
rows, record = [], {}


def add(name, params, theory, est, err, approx=False):
    """approx=True: `theory` is an approximation (no z-score)."""
    z = (est - theory) / err if err > 0 and not approx else float('nan')
    rows.append((name, params, theory, est, err, z))
    record.setdefault(name, []).append(dict(params=params, theory=theory, estimate=est, stderr=err, z=z))
    print(f"{name:28s} {params:22s} theory {theory: .6g}  MC {est: .6g} +- {err:.2g}  (z = {z:+.2f})")


def mean_err(x):
    x = x.double()
    return x.mean().item(), (x.std() / math.sqrt(x.numel())).item()


def chunks(total, size):
    while total > 0:
        yield min(size, total)
        total -= size


# ----------------------------------------------------------------------------------------------
# 1. Haar moments
# ----------------------------------------------------------------------------------------------
S_HAAR = 400_000 if FULL else 40_000
for cplx in (True, False):
    for n in (2, 4, 8, 16, 32):
        alg = TypeIAlgebra([n], [n], complex_valued=cplx, precision='double', device=dev)
        vals2, vals4 = [], []
        for b in chunks(S_HAAR, 50_000):
            u = alg.random_unitary(n, batch_size=b)[:, 0, 0].abs()
            vals2.append(u ** 2)
            vals4.append(u ** 4)
        v2, v4 = torch.cat(vals2), torch.cat(vals4)
        group = 'U' if cplx else 'O'
        t4 = 2 / (n * (n + 1)) if cplx else 3 / (n * (n + 2))
        add(f"E|U_11|^4, {group}(n)", f"n = {n}", t4, *mean_err(v4))
        if n in (2, 32):
            add(f"E|U_11|^2, {group}(n)", f"n = {n}", 1 / n, *mean_err(v2))

# ----------------------------------------------------------------------------------------------
# 2. Weingarten: second-order mixed moment
# ----------------------------------------------------------------------------------------------
S_W = 400_000 if FULL else 50_000
for N in (4, 16, 64):
    alg = TypeIAlgebra([N], [N], precision='double', device=dev)
    A = torch.diag(torch.tensor([1.0] * (N // 2) + [-1.0] * (N // 2), dtype=torch.complex128, device=dev))
    vals = []
    for b in chunks(S_W if N < 64 else S_W // 4, 20_000):
        U = alg.random_unitary(N, batch_size=b)
        M = A @ U @ A @ U.conj().transpose(-2, -1)
        vals.append((M @ M).diagonal(dim1=-2, dim2=-1).sum(-1).real / N)
    add("E tau((AUBU*)^2)", f"N = {N}", -1 / (N ** 2 - 1), *mean_err(torch.cat(vals)))

# ----------------------------------------------------------------------------------------------
# 3. Circular ensembles: spacing distribution
# ----------------------------------------------------------------------------------------------
SURMISE = {
    1: (lambda s: math.pi / 2 * s * np.exp(-math.pi * s ** 2 / 4), 4 / math.pi),
    2: (lambda s: 32 / math.pi ** 2 * s ** 2 * np.exp(-4 * s ** 2 / math.pi), 3 * math.pi / 8),
    4: (lambda s: 2 ** 18 / (3 ** 6 * math.pi ** 3) * s ** 4 * np.exp(-64 * s ** 2 / (9 * math.pi)),
        45 * math.pi / 128),
}
N_CE = 64 if FULL else 32
S_CE = 20_000 if FULL else 2_000
alg_ce = TypeIAlgebra([1], [1], precision='double', device=dev)
spacings = {}
for beta, measure in ((1, 'coe'), (2, 'haar'), (4, 'cse')):
    out = []
    for b in chunks(S_CE, 2_000):
        S = alg_ce.random_unitary(N_CE, measure=measure, batch_size=b)
        theta = torch.sort(torch.angle(torch.linalg.eigvals(S)), dim=-1)[0]
        if beta == 4:                       # Kramers degeneracy: every eigenvalue twice
            theta = theta[:, ::2]
        n_eff = theta.shape[-1]
        gaps = torch.diff(theta, dim=-1, append=theta[:, :1] + 2 * math.pi)
        out.append(gaps * n_eff / (2 * math.pi))
    s = torch.cat(out).flatten().cpu()
    spacings[beta] = s.numpy()
    p, s2 = SURMISE[beta]
    hist, edges = np.histogram(spacings[beta], bins=60, range=(0, 3), density=True)
    centers = 0.5 * (edges[1:] + edges[:-1])
    l1 = float(np.sum(np.abs(hist - p(centers))) * (edges[1] - edges[0]))
    m, e = mean_err(s ** 2)
    record.setdefault('spacing_l1', {})[beta] = l1
    add(f"<s^2>, beta = {beta}", f"N = {N_CE} (surmise)", s2, m, e, approx=True)

# ----------------------------------------------------------------------------------------------
# 4. Page's formula
# ----------------------------------------------------------------------------------------------
S_PAGE = 50_000 if FULL else 5_000
n_page = 16 if FULL else 8
page_rows = []
for m in [mm for mm in (2, 3, 4, 6, 8, 12, 16) if mm <= n_page]:
    A = TypeIAlgebra([m], [m], precision='double', device=dev)
    Bs = TypeIAlgebra([n_page], [n_page], precision='double', device=dev)
    AB = tensor_product(A, Bs)
    vals = []
    for b in chunks(S_PAGE, 5_000):
        psi = states.random_density_matrix(AB, batch_size=b, rank=1)
        vals.append(states.von_neumann_entropy(partial_trace(psi, keep=1)))
    est, err = mean_err(torch.cat(vals))
    page = sum(1.0 / k for k in range(n_page + 1, m * n_page + 1)) - (m - 1) / (2 * n_page)
    page_rows.append((m, page, est, err))
    if m in (2, n_page // 2, n_page):
        add("Page entropy <S_A>", f"m = {m}, n = {n_page}", page, est, err)

# ----------------------------------------------------------------------------------------------
# 5. Purity of induced random states
# ----------------------------------------------------------------------------------------------
S_PUR = 100_000 if FULL else 10_000
for k, r in ((4, 2), (4, 4), (16, 4), (16, 16)):
    alg = TypeIAlgebra([k], [k], precision='double', device=dev)
    vals = torch.cat([states.purity(states.random_density_matrix(alg, batch_size=b, rank=r))
                      for b in chunks(S_PUR, 20_000)])
    add("E Tr rho^2 (induced)", f"k = {k}, r = {r}", (k + r) / (k * r + 1), *mean_err(vals))

# ----------------------------------------------------------------------------------------------
# 6. Spontaneous emission (deterministic, error only)
# ----------------------------------------------------------------------------------------------
alg = TypeIAlgebra([2], [2], precision='double', device=dev)
gamma = 0.3
H = alg.from_blocks([torch.diag(torch.tensor([0.5, -0.5]))])
L = alg.from_blocks([torch.tensor([[0.0, 0.0], [1.0, 0.0]])])
psi = torch.tensor([0.6, 0.8], dtype=torch.complex128)
rho0 = DensityMatrix(alg, matrix=torch.outer(psi, psi.conj())[None, None].to(dev))
ts = [0.0, 1.0, 2.0, 5.0, 10.0]
rk = dynamics.lindblad_evolve(rho0, H, [L], ts, rates=[gamma], substeps=100)
err_exp, err_rk = 0.0, 0.0
for t, r in zip(ts, rk):
    exact_pe, exact_coh = 0.36 * math.exp(-gamma * t), 0.48 * math.exp(-gamma * t / 2)
    ex = dynamics.lindblad_channel(alg, H, [L], t=t, rates=[gamma])(rho0).matrix[0, 0]
    err_exp = max(err_exp, abs(ex[0, 0].real.item() - exact_pe), abs(ex[0, 1].abs().item() - exact_coh))
    m = r.matrix[0, 0]
    err_rk = max(err_rk, abs(m[0, 0].real.item() - exact_pe), abs(m[0, 1].abs().item() - exact_coh))
record['spontaneous_emission'] = dict(max_error_exp_tL=err_exp, max_error_rk4=err_rk)
print(f"spontaneous emission: max error exp(tL) {err_exp:.1e}, RK4 {err_rk:.1e}")

# ----------------------------------------------------------------------------------------------
# 7. The three trace functionals (identities that hold exactly)
# ----------------------------------------------------------------------------------------------
alg = TypeIAlgebra([2, 3, 5], [2, 3, 5], precision='double', device=dev)
Bt = 2000
def rand_op():
    m = torch.zeros(Bt, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype, device=dev)
    for c, k in enumerate(alg.k_factors):
        m[:, c, :k, :k] = torch.randn(Bt, k, k, dtype=alg.dtype, device=dev)
    return alg.operator(m)
A, Bop = rand_op(), rand_op()
U = alg.random_unitary_operator(Bt)
one = alg.identity(1)
k = torch.tensor(alg.k_factors, dtype=torch.float64)
trace_checks = []
for name, f in (('Tr_blunt', lambda X: X.Tr_blunt()), ('Tr_norm', lambda X: X.Tr_norm()),
                ('tau_vN', lambda X: X.tau_vN())):
    norm1 = {'Tr_blunt': float(k.sum()), 'Tr_norm': float(alg.C), 'tau_vN': 1.0}[name]
    def val(X):
        return f(X)

    def cmp(a, b):
        return (a - b).abs().max().item()
    unit = abs(val(one).real.item() - norm1)
    tracial = cmp(val(A @ Bop), val(Bop @ A))
    unitary = cmp(val(U @ A @ U.adjoint()), val(A))
    positive = (val(A.adjoint() @ A).real / (A.frobenius_norm() ** 2 + 1e-300)).min().item()
    trace_checks.append(dict(trace=name, unit=unit, tracial=tracial, unitary_invariance=unitary,
                             min_positivity_ratio=positive))
    print(f"{name:9s}: f(1) error {unit:.1e}, |f(AB)-f(BA)| {tracial:.1e}, |f(UAU*)-f(A)| {unitary:.1e}, "
          f"min f(A*A)/||A||_F^2 = {positive:.3f} > 0")
rel = (A.Tr_norm() - alg.C * A.tau_vN()).abs().max().item()
record['trace_functionals'] = dict(checks=trace_checks, Tr_norm_equals_C_tau=rel)

# ----------------------------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------------------------
save_json(args, 'validation', record)
tex = env_macro(args, 'Val')
tex += f"\\newcommand{{\\ValEmissionExp}}{{{sci(err_exp)}}}\n\\newcommand{{\\ValEmissionRK}}{{{sci(err_rk)}}}\n"
tex += f"\\newcommand{{\\ValNce}}{{{N_CE}}}\n\\newcommand{{\\ValSce}}{{{S_CE}}}\n"
tex += f"\\newcommand{{\\ValSHaar}}{{{S_HAAR}}}\n\\newcommand{{\\ValSPage}}{{{S_PAGE}}}\n"
tex += f"\\newcommand{{\\ValNPage}}{{{n_page}}}\n"
for beta in (1, 2, 4):
    name = {1: 'One', 2: 'Two', 4: 'Four'}[beta]
    tex += f"\\newcommand{{\\ValLone{name}}}{{{record['spacing_l1'][beta]:.3f}}}\n"
tex += "\\newcommand{\\ValTraceRows}{%\n"
for tc in trace_checks:
    nm = {'Tr_blunt': '$\\Tr_{\\mathrm{blunt}}$', 'Tr_norm': '$\\Tr_{\\mathrm{norm}}$',
          'tau_vN': '$\\tau_{\\mathrm{vN}}$'}[tc['trace']]
    tex += (f"{nm} & ${sci(max(tc['unit'], 1e-300)) if tc['unit'] else '0'}$ & ${sci(tc['tracial'])}$ "
            f"& ${sci(tc['unitary_invariance'])}$ & {tc['min_positivity_ratio']:.3f} \\\\\n")
tex += "}\n"
tex += "\\newcommand{\\ValidationRows}{%\n"
labels = {
    'E|U_11|^4, U(n)': r'$\mathbb E|U_{11}|^4$, $U(n)$',
    'E|U_11|^2, U(n)': r'$\mathbb E|U_{11}|^2$, $U(n)$',
    'E|U_11|^4, O(n)': r'$\mathbb E|U_{11}|^4$, $O(n)$',
    'E|U_11|^2, O(n)': r'$\mathbb E|U_{11}|^2$, $O(n)$',
    'E tau((AUBU*)^2)': r'$\mathbb E\,\tau((AUBU^*)^2)$',
    '<s^2>, beta = 1': r'$\langle s^2\rangle$, COE',
    '<s^2>, beta = 2': r'$\langle s^2\rangle$, CUE',
    '<s^2>, beta = 4': r'$\langle s^2\rangle$, CSE',
    'Page entropy <S_A>': r'$\langle S(\rho_A)\rangle$',
    'E Tr rho^2 (induced)': r'$\mathbb E\Tr\rho^2$',
}
for name, params, theory, est, err, z in rows:
    p = params.replace('(surmise)', '').strip()
    p = p.replace(' = ', '=')
    tex += (f"{labels[name]} & ${p}$ & {fmt(theory, 6)} & {fmt(est, 6)} & ${sci(err)}$ "
            f"& {'--' if math.isnan(z) else f'{z:+.1f}'} \\\\\n")
tex += "}\n"
write_tex(args, 'validation', tex)

# figures
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(10, 3), sharey=True)
xs = np.linspace(0, 3, 300)
for ax, beta, name in zip(axes, (1, 2, 4), ('COE', 'CUE', 'CSE')):
    ax.hist(spacings[beta], bins=60, range=(0, 3), density=True, color='#9ab', label='library')
    ax.plot(xs, SURMISE[beta][0](xs), 'k-', lw=1.2, label='Wigner surmise')
    ax.set_title(f'{name} ($\\beta={beta}$), $N={N_CE}$')
    ax.set_xlabel('$s$')
axes[0].set_ylabel('$P(s)$')
axes[0].legend(frameon=False)
fig.tight_layout()
fig.savefig(f"{args.figdir}/spacings.pdf")

fig, ax = plt.subplots(figsize=(4.5, 3.2))
ms = [r[0] for r in page_rows]
ax.plot(ms, [r[1] for r in page_rows], 'k-', label="Page's formula")
ax.errorbar(ms, [r[2] for r in page_rows], yerr=[2 * r[3] for r in page_rows], fmt='o', ms=4, label='library')
ax.set_xlabel('$m = \\dim H_A$')
ax.set_ylabel('$\\langle S(\\rho_A)\\rangle$')
ax.set_title(f'$n = \\dim H_B = {n_page}$')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(f"{args.figdir}/page.pdf")
print("written:", args.out, args.figdir)
