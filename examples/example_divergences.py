#!/usr/bin/env python3
"""
Example: quantum Renyi divergences, data processing and recovery on M = M_3 (+) M_2 (+) M_4.

1. Petz and sandwiched Renyi divergences of random states for several alpha, with the limits (Umegaki, D_max,
   -log F at alpha = 1/2) and the ordering ~D_alpha <= D_alpha.
2. The sector decomposition D(rho||sigma) = D(p||q) + sum_c p_c D(rho_c/p_c||sigma_c/q_c).
3. Data processing under a random channel and under the conditional expectation onto the centre.
4. Petz recovery: R(Phi(sigma)) = sigma, and the recoverability bound of Junge et al. (2018) with the universal
   recovery map.
5. Noncommutative L^p norms for the three traces, and the state-weighted (Kosaki) norms, which do not depend on the
   trace.
"""
import math

import torch

from torch_vn_algebra import TypeIAlgebra, states, channels, divergences as dv

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'


def main():
    torch.manual_seed(0)
    B = 1000
    alg = TypeIAlgebra([3, 2, 4], [3, 2, 4], precision='double', device=DEVICE)
    rho, sigma = states.random_density_matrix(alg, B), states.random_density_matrix(alg, B)

    print("1. Renyi divergences (mean over", B, "random pairs of states)")
    print(f"{'alpha':>7s} {'Petz D_a':>10s} {'sandwiched':>11s} {'max(~D - D)':>12s}")
    for a in (0.25, 0.5, 0.9, 1.5, 2.0, 5.0):
        P, S = dv.petz_renyi_divergence(rho, sigma, a), dv.sandwiched_renyi_divergence(rho, sigma, a)
        print(f"{a:7.2f} {P.mean().item():10.4f} {S.mean().item():11.4f} {(S - P).max().item():12.1e}")
    D = dv.umegaki_relative_entropy(rho, sigma)
    print(f"  Umegaki D = {D.mean().item():.4f}; D_max = {dv.max_relative_entropy(rho, sigma).mean().item():.4f}; "
          f"|~D_(1/2) + log F| <= "
          f"{(dv.sandwiched_renyi_divergence(rho, sigma, 0.5) + torch.log(states.fidelity(rho, sigma))).abs().max().item():.1e}")
    # one value of alpha per sample: the whole curve alpha -> ~D_alpha of one pair in one call
    alphas = torch.linspace(0.5, 4.0, 8, dtype=torch.float64, device=DEVICE)
    one = lambda X: states.DensityMatrix(alg, matrix=X.matrix[:1], validate=False)
    curve = dv.sandwiched_renyi_divergence(one(rho), one(sigma), alphas)
    print("  ~D_alpha of one pair for alpha = 0.5, 1 (Umegaki), ..., 4:", " ".join(f"{v:.3f}" for v in curve.tolist()),
          "(nondecreasing)")

    print("2. Sector decomposition of the relative entropy")
    sc = dv.sector_contributions(rho, sigma, kind='umegaki')
    classical = (sc.p * torch.log(sc.p / sc.q)).sum(-1)
    quantum = (sc.p * sc.conditional).sum(-1)
    print(f"  mean D(p||q) = {classical.mean().item():.4f}, mean sum_c p_c D(rho_c||sigma_c) = "
          f"{quantum.mean().item():.4f}, max |sum - D| = {(classical + quantum - D).abs().max().item():.1e}")

    print("3. Data processing: largest D(Phi rho||Phi sigma) - D(rho||sigma) over the samples (<= 0 by the theorems)")
    for name, Phi in (('random channel', channels.random_channel(alg, 2, B)),
                      ('centre expectation', channels.center_expectation(alg))):
        r, s = Phi(rho), Phi(sigma)
        worst_s = max((dv.sandwiched_renyi_divergence(r, s, a) - dv.sandwiched_renyi_divergence(rho, sigma, a))
                      .max().item() for a in (0.5, 1.5, 3.0))
        worst_p = max((dv.petz_renyi_divergence(r, s, a) - dv.petz_renyi_divergence(rho, sigma, a)).max().item()
                      for a in (0.3, 1.5, 2.0))
        print(f"  {name:20s} sandwiched (alpha = 0.5, 1.5, 3): {worst_s:+.2e}   Petz (alpha = 0.3, 1.5, 2): "
              f"{worst_p:+.2e}")

    print("4. Petz recovery map")
    Phi = channels.random_channel(alg, 2, B)
    R = dv.PetzRecovery(sigma, Phi)
    print(f"  max |R(Phi(sigma)) - sigma| = {(R(Phi(sigma)).matrix - sigma.matrix).abs().max().item():.1e}")
    small = slice(0, 50)                                     # the universal map integrates ~200 rotated maps
    r50, s50 = (states.DensityMatrix(alg, matrix=X.matrix[small], validate=False) for X in (rho, sigma))
    Phi50 = channels.Channel(alg, Phi.kraus[small], trace_preserving=True)
    gap = dv.umegaki_relative_entropy(r50, s50) - dv.umegaki_relative_entropy(Phi50(r50), Phi50(s50))
    bound = -torch.log(states.fidelity(r50, dv.universal_recovery(s50, Phi50, Phi50(r50))))
    print(f"  D(rho||sigma) - D(Phi rho||Phi sigma) >= -log F(rho, R_univ(Phi rho)) on {gap.shape[0]} pairs: "
          f"min(gap - bound) = {(gap - bound).min().item():.2e}, mean gap {gap.mean().item():.3f}, "
          f"mean bound {bound.mean().item():.3f}")

    print("5. L^p norms of A = rho - sigma (mean) and Kosaki norms of a random x")
    A = rho - sigma
    for trace in ('blunt', 'norm', 'tau_vN'):
        vals = " ".join(f"{dv.lp_norm(A, p, trace).mean().item():.4f}" for p in (1, 2, 4, math.inf))
        print(f"  {trace:7s} p = 1, 2, 4, inf: {vals}")
    x = alg.operator(torch.randn(B, alg.C, alg.k_max, alg.k_max, dtype=alg.dtype, device=DEVICE))
    k_blunt = dv.kosaki_norm(x, rho, 2.0)
    d_norm = alg.operator(rho.density('norm'), is_self_adjoint=True, is_positive=True)
    k_norm = dv.kosaki_norm(x, d_norm, 2.0, trace='norm')
    print(f"  ||x||_(2,rho) from the Tr_blunt and the Tr_norm densities agree to "
          f"{(k_blunt - k_norm).abs().max().item():.1e}")


if __name__ == '__main__':
    main()
