"""Violation envelope of Gardner's inequality |phi(X)| <= phi(|X|) for phi = Tr(W .) on M_2 (item 1 of sec:tc)."""
import torch
from torch_vn_algebra import TypeIAlgebra, optimize as opt

torch.manual_seed(0)                                                  # reproducible output
alg = TypeIAlgebra([2], [2], precision='double')
grid, starts = torch.tensor([0.1, 0.3, 0.5, 0.7, 0.9], dtype=torch.float64), 16
B = len(grid) * starts                                                # one batch: contrasts x starts

W = opt.SpectrumParam(alg, B, delta=grid.repeat_interleave(starts))   # Delta(W) exact
P = opt.PositiveParam(alg, B, delta_range=(0.0, 1.0))                 # |X|, any rank
U = opt.UnitaryParam(alg, B)                                          # X = U P


def objective():
    w = W.eigenvalues() / W.eigenvalues().sum((-2, -1), keepdim=True)
    phi = lambda op: (torch.diagonal(op.matrix, dim1=-2, dim2=-1) * w).sum((-2, -1))
    p = phi(P.operator()).real
    return (phi(U.operator() @ P.operator()).abs() - p) / p


V = opt.extremize(objective, [W, P, U]).reshape(-1, starts).max(1).values
kappa = (1 - grid ** 2) ** -0.5                                       # Kantorovich constant
for d, v, k in zip(grid.tolist(), V.tolist(), kappa.tolist()):
    print(f"delta={d:.1f}: V={v:.6f}   kappa-1={k - 1:.6f}")
