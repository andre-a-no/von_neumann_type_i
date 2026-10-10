#!/usr/bin/env python3
"""
Validation runs reported in Section 4.1 of the paper:

1. Haar measure: E|U_11|^2 = 1/n and E[U_11^2] = 0 for U ~ Haar(U(n)),
   using the library sampler TypeIAlgebra.random_unitary.
2. Power iteration: number of iterations to reach 1e-8 for a small and a large spectral gap.
3. Square root via SVD: relative error ||sqrt(A)^2 - A||_F / ||A||_F for random positive A.

Usage:  python scripts/validation.py [--device cuda] [--samples 10000]
"""
import argparse
import os
import sys

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from torch_vn_algebra import TypeIAlgebra


def test_haar(device, num_samples):
    print(f"=== Haar measure on U(n), {num_samples} samples ===")
    print("n\t1/n\t\tmean |U11|^2\trel. error\t|mean U11^2|\tstd. error of mean")
    alg = TypeIAlgebra([1], [1], complex_valued=True, device=device)
    for n in [2, 4, 8, 16, 32]:
        u = alg.random_unitary(n, measure='haar', batch_size=num_samples)[:, 0, 0]
        p = (u.abs() ** 2).double()
        avg = p.mean().item()
        sem = p.std().item() / np.sqrt(num_samples)
        rel = abs(avg - 1 / n) * n
        print(f"{n}\t{1/n:.5f}\t\t{avg:.5f}\t\t{rel:.2e}\t{(u ** 2).mean().abs().item():.2e}\t{sem * n:.2e} (rel.)")
    print()


def test_power_iteration_gap(tol=1e-8):
    print("=== Power iteration vs spectral gap ===")
    rng = np.random.default_rng(0)
    for eigvals in ([1.0, 0.99, 0.98], [2.0, 1.0, 0.5]):
        dim = len(eigvals)
        q, _ = np.linalg.qr(rng.standard_normal((dim, dim)))
        mat = q @ np.diag(eigvals) @ q.T
        v = rng.standard_normal(dim)
        v /= np.linalg.norm(v)
        lambda_old = 0.0
        for i in range(100000):
            w = mat @ v
            lambda_est = v @ w
            if abs(lambda_est - lambda_old) < tol:
                break
            lambda_old = lambda_est
            v = w / np.linalg.norm(w)
        print(f"spectrum {eigvals}: gap {eigvals[0] - eigvals[1]:.2f}, iterations to {tol:g}: {i + 1}, "
              f"error {abs(lambda_est - eigvals[0]):.1e}")
    print()


def test_svd_sqrt(device, dim=100, num_matrices=1000):
    print(f"=== SVD square root, {num_matrices} positive {dim}x{dim} matrices, eigenvalues U[0,1] ===")
    alg = TypeIAlgebra([dim], [dim], complex_valued=False, device=device)
    A = alg.operator_from_eigenvalues(lambda d: torch.rand(num_matrices, d, device=device),
                                      batch_size=num_matrices, force_positive=True, force_self_adjoint=True)
    for name, dtype in (('float32 (library default)', None), ('float64', torch.float64)):
        if dtype is None:
            mats, roots = A.matrix[:, 0], A.sqrt().matrix[:, 0]
        else:
            mats = A.matrix[:, 0].to(dtype)
            _, S, Vh = torch.linalg.svd(mats)
            roots = Vh.mT @ torch.diag_embed(S.sqrt()) @ Vh
        err = torch.linalg.matrix_norm(roots @ roots - mats) / torch.linalg.matrix_norm(mats)
        print(f"{name}: mean rel. error {err.mean().item():.2e}, max {err.max().item():.2e}")
    print()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--device', default='cuda' if torch.cuda.is_available() else 'cpu')
    parser.add_argument('--samples', type=int, default=10000)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    torch.manual_seed(args.seed)
    test_haar(args.device, args.samples)
    test_power_iteration_gap()
    test_svd_sqrt(args.device)
