#!/usr/bin/env python3
"""
Section "Symmetric tensor products and spin chains" of paper/v2: exact checks of SpinChain.

  1. XX chain (open): the spectrum of every sector N equals the sums of N distinct single-particle
     energies cos(pi k / (L + 1)) (Jordan-Wigner, free fermions)
  2. Heisenberg ring: ground-state energies against the tabulated exact values; SU(2) nesting of
     sector spectra
  3. Domain-wall melting in the XX chain: <S^z_i(t)> from many-body evolution in the sector
     N = L/2 against the single-particle correlation matrix
  4. Entanglement entropy of the XX ground state against Peschel's formula
  5. Random-field Heisenberg chain: mean level-spacing ratio <r> versus disorder W, from GOE
     (0.5307) to Poisson (2 ln 2 - 1 = 0.3863); disorder realisations form the batch

Writes generated/chains.tex, generated/chains.json, figures/chain_*.pdf.
"""
import itertools
import math
import time

import numpy as np
import torch

from common import parse_args, save_json, write_tex, env_macro, sci, sync

import warnings

from torch_vn_algebra import SpinChain, DensityMatrix, dynamics, krylov, cost

warnings.simplefilter('ignore', cost.CostWarning)

args = parse_args(__doc__)
dev = args.device
FULL = args.mode in ('full', 'check')            # check: the sizes of full, minimal repetitions
CHECK = args.mode == 'check'
rec = {}

# exact ground-state energies of the S=1/2 Heisenberg ring, H = sum S_i . S_{i+1}
HEISENBERG_RING = {4: -2.0, 6: -2.802775638, 8: -3.651093409, 10: -4.515446354, 12: -5.387390917,
                   14: -6.263549534, 16: -7.142296361}
BETHE = 0.25 - math.log(2)                                     # infinite-chain energy per site


def hop(L):
    h = torch.zeros(L, L, dtype=torch.complex128)
    for i in range(L - 1):
        h[i, i + 1] = h[i + 1, i] = 0.5
    return h


# 1. free fermions -----------------------------------------------------------------------------
L1 = 12 if FULL else 10
ch = SpinChain(L1, device=dev)
eps = [math.cos(math.pi * k / (L1 + 1)) for k in range(1, L1 + 1)]
err_ff = 0.0
for N in range(L1 + 1):
    w = torch.linalg.eigvalsh(ch.xxz(J=1.0, Delta=0.0, sector=N).matrix[0, 0])
    exact = sorted(sum(c) for c in itertools.combinations(eps, N)) if N else [0.0]
    err_ff = max(err_ff, (w.cpu() - torch.tensor(exact, dtype=torch.float64)).abs().max().item())
rec['free_fermions'] = dict(L=L1, max_error=err_ff, levels=2 ** L1)
print(f"XX chain L={L1}: all {2 ** L1} levels match free fermions to {err_ff:.1e}")

# 2. Heisenberg ring ---------------------------------------------------------------------------
ring = []
for L, e_exact in HEISENBERG_RING.items():
    if L > (16 if FULL else 12):
        continue
    t0 = time.time()
    c = SpinChain(L, 'periodic', complex_valued=False, device=dev)
    e0 = torch.linalg.eigvalsh(c.heisenberg(sector=L // 2).matrix[0, 0])[0].item()
    sync(dev)
    ring.append(dict(L=L, E0=e0, exact=e_exact, error=abs(e0 - e_exact), seconds=time.time() - t0,
                     dim=c.sector_dim(L // 2)))
    print(f"Heisenberg ring L={L:2d}: E0 = {e0:.9f} (exact {e_exact}), dim {c.sector_dim(L // 2)}, "
          f"{time.time() - t0:.1f} s")
c = SpinChain(10, 'periodic', complex_valued=False, device=dev)
E = [torch.sort(torch.linalg.eigvalsh(c.heisenberg(sector=N).matrix[0, 0]))[0] for N in range(11)]
su2 = max(max((E[N + 1] - e).abs().min().item() for e in E[N]) for N in range(5))
rec['heisenberg_ring'] = ring
rec['su2_nesting_error'] = su2

# 3. domain wall --------------------------------------------------------------------------------
L3 = 16 if FULL else 12
N3 = L3 // 2
c3 = SpinChain(L3, device=dev)
H = c3.xxz(J=1.0, Delta=0.0, sector=N3)
psi0 = c3.vector_in_sector('1' * N3 + '0' * (L3 - N3))
times = np.linspace(0, 0.6 * L3, 25)
t0 = time.time()
out = dynamics.schrodinger(psi0[None, None], H, times)
bits = c3.bits(N3).to(torch.float64)
mag = torch.stack([(p[0, 0].abs() ** 2).double() @ bits for p in out]).cpu() - 0.5      # <S^z_i(t)>
C0 = torch.diag(torch.tensor([1.0] * N3 + [0.0] * (L3 - N3), dtype=torch.complex128))
mag_ff = torch.stack([(torch.linalg.matrix_exp(-1j * t * hop(L3)).conj() @ C0
                       @ torch.linalg.matrix_exp(-1j * t * hop(L3)).T).diagonal().real for t in times]) - 0.5
err_dw = (mag - mag_ff).abs().max().item()
rec['domain_wall'] = dict(L=L3, sector_dim=c3.sector_dim(N3), max_error=err_dw, seconds=time.time() - t0)
print(f"domain wall L={L3} (sector dim {c3.sector_dim(N3)}): max error {err_dw:.1e}")

# 4. Peschel ------------------------------------------------------------------------------------
L4 = 14 if FULL else 12
N4 = L4 // 2
c4 = SpinChain(L4, device=dev)
w, V = c4.xxz(J=1.0, Delta=0.0, sector=N4).eigh()[0]
psi = V[0, :, 0]
rho = DensityMatrix(c4.sector_algebra(N4), matrix=torch.outer(psi, psi.conj())[None, None])
e, U = torch.linalg.eigh(hop(L4).real)
Cm = U[:, :N4] @ U[:, :N4].T
S_lib, S_pes = [], []
for nA in range(1, L4):
    S_lib.append(c4.entanglement_entropy(rho, nA, sector=N4).item())
    nu = torch.linalg.eigvalsh(Cm[:nA, :nA]).clamp(1e-14, 1 - 1e-14)
    S_pes.append((-(nu * torch.log(nu) + (1 - nu) * torch.log(1 - nu)).sum()).item())
err_pes = max(abs(a - b) for a, b in zip(S_lib, S_pes))
rec['peschel'] = dict(L=L4, max_error=err_pes, S=S_lib)
print(f"entanglement L={L4}: max deviation from Peschel {err_pes:.1e}")

# 5. level statistics ---------------------------------------------------------------------------
SIZES = (10, 12, 14) if FULL else (8, 10)
REAL = 20 if CHECK else 400 if FULL else 40
WS = [0.5, 12] if CHECK else [0.5, 1, 2, 3, 4, 6, 8, 12]
mbl = {}
g = torch.Generator().manual_seed(args.seed)
for L in SIZES:
    c = SpinChain(L, 'periodic', complex_valued=False, device=dev)
    rs = []
    t0 = time.time()
    for W in WS:
        vals = []
        left = REAL
        while left > 0:
            b = min(left, 100 if L <= 12 else 20)
            Hd = c.random_field_heisenberg(W, b, sector=L // 2, generator=g)
            Ew = torch.linalg.eigvalsh(Hd.matrix[:, 0])
            k = Ew.shape[-1]
            vals.append(SpinChain.level_spacing_ratio(Ew[:, k // 4: 3 * k // 4]).cpu())
            left -= b
        v = torch.cat(vals)
        rs.append((v.mean().item(), (v.std() / math.sqrt(len(v))).item()))
    mbl[L] = rs
    print(f"L={L}: <r>(W) = " + ", ".join(f"{m:.3f}" for m, _ in rs) + f"  ({time.time() - t0:.0f} s)")
rec['level_statistics'] = dict(W=WS, realisations=REAL, r={str(L): v for L, v in mbl.items()})

# 6. Lanczos beyond dense matrices ---------------------------------------------------------------
kry_rows = []
for L in ((16, 20, 24) if FULL else (12, 16, 18)):
    c = SpinChain(L, complex_valued=False, device=dev)
    t0 = time.time()
    E, _, res = krylov.ground_state(c.xxz_sparse(J=1.0, Delta=0.0, sector=L // 2))
    sync(dev)
    eps_k = sorted(math.cos(math.pi * k / (L + 1)) for k in range(1, L + 1))
    exact = sum(eps_k[:L // 2])
    kry_rows.append(dict(model=f'XX open, $L={L}$', dim=c.sector_dim(L // 2), E=E.item(), exact=exact,
                         error=abs(E.item() - exact), seconds=time.time() - t0))
c = SpinChain(16, 'periodic', complex_valued=False, device=dev)
t0 = time.time()
E, _, _ = krylov.ground_state(c.xxz_sparse(sector=8))
kry_rows.append(dict(model='Heisenberg ring, $L=16$', dim=c.sector_dim(8), E=E.item(), exact=-7.142296361,
                     error=abs(E.item() + 7.142296361), seconds=time.time() - t0, digits=9))
for r in kry_rows:
    print(f"Lanczos {r['model']}: dim {r['dim']}, E0 {r['E']:.10f}, error {r['error']:.1e}, {r['seconds']:.1f} s")

# 7. Krylov evolution beyond dense matrices --------------------------------------------------------
L7 = 22 if FULL else 16
N7 = L7 // 2
c7 = SpinChain(L7, device=dev)
t0 = time.time()
out7 = krylov.evolve(c7.xxz_sparse(J=1.0, Delta=0.0, sector=N7), c7.vector_in_sector('1' * N7 + '0' * N7)[None],
                     [0.0, 2.0, 4.0])
bits7 = c7.bits(N7).to(torch.float64)
C07 = torch.diag(torch.tensor([1.0] * N7 + [0.0] * N7, dtype=torch.complex128))
err_kdw = 0.0
for t, p in zip([0.0, 2.0, 4.0], out7):
    U = torch.linalg.matrix_exp(-1j * t * hop(L7))
    err_kdw = max(err_kdw, ((p[0].abs() ** 2).double().cpu() @ bits7.cpu()
                            - (U.conj() @ C07 @ U.T).diagonal().real).abs().max().item())
rec['krylov_domain_wall'] = dict(L=L7, dim=c7.sector_dim(N7), error=err_kdw, seconds=time.time() - t0)
print(f"Krylov domain wall L={L7} (dim {c7.sector_dim(N7)}): error {err_kdw:.1e}, {time.time() - t0:.0f} s")
rec['lanczos'] = kry_rows

# 8. momentum sectors ----------------------------------------------------------------------------
L8 = 14 if FULL else 12
c8 = SpinChain(L8, 'periodic', device=dev)
t0 = time.time()
wk = torch.sort(torch.cat([w[0] for w in c8.xxz_momentum(1.0, 1.0, sector=L8 // 2).eigenvalues()]))[0]
t_mom = time.time() - t0
t0 = time.time()
wd = torch.linalg.eigvalsh(c8.xxz(1.0, 1.0, sector=L8 // 2).matrix[0, 0])
t_dense = time.time() - t0
err_mom = (wk - wd).abs().max().item()
dims8 = c8.momentum_dims(L8 // 2)
rec['momentum'] = dict(L=L8, dims=dims8, error=err_mom, seconds_momentum=t_mom, seconds_dense=t_dense)
print(f"momentum sectors L={L8}: blocks {dims8}, spectra agree to {err_mom:.1e} "
      f"({t_mom:.2f} s vs dense {t_dense:.2f} s)")

# 9. particle loss between sectors -----------------------------------------------------------------
c9 = SpinChain(6, 'periodic', device=dev)
g = 0.4
ts9 = [0.0, 0.5, 1.0, 2.0]
rhos = dynamics.lindblad_evolve(c9.basis_state('110111'), c9.xxz(1.0, 0.7), [c9.lowering(i) for i in range(6)],
                                ts9, rates=[g] * 6, substeps=80)
err_loss = 0.0
for t, r in zip(ts9, rhos):
    q = math.exp(-g * t)
    binom = torch.tensor([math.comb(5, n) * q ** n * (1 - q) ** (5 - n) for n in range(6)] + [0.0],
                         dtype=torch.float64)
    err_loss = max(err_loss, (r.sector_probabilities()[0].cpu() - binom).abs().max().item())
rec['loss'] = dict(error=err_loss)
print(f"particle loss: sector distribution vs binomial death process, max error {err_loss:.1e}")

# output ---------------------------------------------------------------------------------------
save_json(args, 'chains', rec)
tex = env_macro(args, 'Chain')
tex += (f"\\newcommand{{\\ChainFFL}}{{{L1}}}\n\\newcommand{{\\ChainFFErr}}{{{sci(err_ff)}}}\n"
        f"\\newcommand{{\\ChainSUErr}}{{{sci(max(su2, 1e-16))}}}\n"
        f"\\newcommand{{\\ChainDWL}}{{{L3}}}\n\\newcommand{{\\ChainDWDim}}{{{c3.sector_dim(N3)}}}\n"
        f"\\newcommand{{\\ChainDWErr}}{{{sci(err_dw)}}}\n"
        f"\\newcommand{{\\ChainPesL}}{{{L4}}}\n\\newcommand{{\\ChainPesErr}}{{{sci(err_pes)}}}\n"
        f"\\newcommand{{\\ChainReal}}{{{REAL}}}\n"
        f"\\newcommand{{\\ChainKDWL}}{{{L7}}}\n\\newcommand{{\\ChainKDWDim}}{{{c7.sector_dim(N7)}}}\n"
        f"\\newcommand{{\\ChainKDWErr}}{{{sci(max(err_kdw, 1e-16))}}}\n"
        f"\\newcommand{{\\ChainMomL}}{{{L8}}}\n\\newcommand{{\\ChainMomErr}}{{{sci(max(err_mom, 1e-16))}}}\n"
        f"\\newcommand{{\\ChainMomDims}}{{{', '.join(map(str, dims8))}}}\n"
        f"\\newcommand{{\\ChainLossErr}}{{{sci(max(err_loss, 1e-16))}}}\n")
tex += "\\newcommand{\\ChainKrylovRows}{%\n"
def diff_cell(err, digits):
    """An exact value given to `digits` decimals cannot resolve differences below half its last digit."""
    if digits and err < 0.5 * 10.0 ** -digits:
        return f"$<{sci(0.5 * 10.0 ** -digits, 1)}$"
    return f"${sci(max(err, 1e-16))}$"


for r in kry_rows:
    d = r.get('digits')
    tex += (f"{r['model']} & {r['dim']} & {r['E']:.10f} & {r['exact']:.{d or 10}f} & {diff_cell(r['error'], d)} "
            f"& {r['seconds']:.2g} \\\\\n")
tex += "}\n"
tex += "\\newcommand{\\ChainRingRows}{%\n"
for r in ring:
    tex += (f"{r['L']} & {r['dim']} & {r['E0']:.9f} & {r['exact']:.9f} & {diff_cell(r['error'], 9)} "
            f"& {r['E0'] / r['L']:.5f} \\\\\n")
tex += "}\n"
write_tex(args, 'chains', tex)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
ax = axes[0]
im = ax.imshow(mag.numpy(), aspect='auto', origin='lower', cmap='RdBu_r', vmin=-0.5, vmax=0.5,
               extent=[-0.5, L3 - 0.5, times[0], times[-1]])
ax.set_xlabel('site $i$')
ax.set_ylabel('$t$')
ax.set_title(f'domain wall, $L={L3}$: $\\langle S^z_i(t)\\rangle$', fontsize=9)
fig.colorbar(im, ax=ax)
ax = axes[1]
ax.plot(range(1, L4), S_pes, 'k-', label='Peschel (free fermions)')
ax.plot(range(1, L4), S_lib, 'o', ms=4, label='library (reduced state)')
ax.set_xlabel('$n$ (left block)')
ax.set_ylabel('$S(\\rho_{1..n})$')
ax.set_title(f'XX ground state, $L={L4}$', fontsize=9)
ax.legend(frameon=False, fontsize=7)
ax = axes[2]
for L, rs in mbl.items():
    ax.errorbar(WS, [m for m, _ in rs], yerr=[2 * e for _, e in rs], fmt='o-', ms=3, label=f'$L={L}$')
ax.axhline(0.5307, color='k', ls='--', lw=0.8)
ax.axhline(2 * math.log(2) - 1, color='k', ls=':', lw=0.8)
ax.text(WS[-1], 0.5307, 'GOE', ha='right', va='bottom', fontsize=7)
ax.text(WS[-1], 0.3863, 'Poisson', ha='right', va='bottom', fontsize=7)
ax.set_xscale('log')
ax.set_xlabel('disorder $W$')
ax.set_ylabel('$\\langle r\\rangle$')
ax.set_title(f'random-field Heisenberg, {REAL} realisations', fontsize=9)
ax.legend(frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(f"{args.figdir}/chains.pdf", metadata={"CreationDate": None})
print("written:", args.out, args.figdir)
