"""Random operators with a prescribed spectrum in M = M_2 (+) M_3 (+) M_4."""
import torch
from torch_vn_algebra import TypeIAlgebra

torch.manual_seed(0)                               # reproducible output
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
alg = TypeIAlgebra(n_factors=[2, 3, 4], k_factors=[2, 3, 4], device=dev)
B = 10_000                                     # Monte Carlo samples, processed as one batch


def contrast_sampler(delta):
    """Spectra in [lo, 1] with lambda_max = 1, lambda_min = lo, i.e. Michelson contrast delta."""
    lo = (1 - delta) / (1 + delta)

    def sampler(dim):
        lam = lo + (1 - lo) * torch.rand(B, dim, device=dev)
        lam[:, 0], lam[:, -1] = 1.0, lo
        return lam
    return sampler


# A_c = U_c diag(lambda) U_c^*, with independent Haar unitaries U_c in every block and sample
X = alg.operator_from_eigenvalues(contrast_sampler(0.6), batch_size=B,
                                  force_positive=True, force_self_adjoint=True)
print(X.matrix.shape)                          # (10000, 3, 4, 4): batch, sector, padded block
print(X.michelson_contrast[:4])                # 0.6 for every sample
print(X.Tr_blunt().real.mean(), X.Tr_norm().real.mean(), X.tau_vN().real.mean())

# functional calculus through the spectral theorem, block by block
R = X.sqrt()
print((R @ R - X).frobenius_norm().max())      # ~1e-6 in single precision
print(X.apply_function(torch.log).tau_vN().real.mean())   # tau_vN(log X)
