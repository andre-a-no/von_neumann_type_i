#!/usr/bin/env python3
"""
Validation tests for Haar measure, power iteration, and SVD-based sqrt.
"""
import torch
import numpy as np
from structures.Algebra import TypeIAlgebra

# ------------------------------------------------------------
# 1. Haar measure test (independent, не использует TypeIAlgebra)
# ------------------------------------------------------------
def test_haar():
    print("=== Haar measure test (independent) ===")
    ns = [2, 4, 8, 16, 32]
    num_samples = 10000
    print("n\tExpected\tMean |U11|^2\tRelative error")
    for n in ns:
        total = 0.0
        for _ in range(num_samples):
            Z = torch.randn(n, n, dtype=torch.float32)
            Q, R = torch.linalg.qr(Z)
            d = torch.diag(R).sign()
            Q = Q * d.unsqueeze(0)
            total += (Q[0, 0].abs() ** 2).item()
        avg = total / num_samples
        expected = 1.0 / n
        rel_err = abs(avg - expected) / expected
        print(f"{n}\t{expected:.5f}\t{avg:.5f}\t{rel_err:.2e}")
    print()

# ------------------------------------------------------------
# 2. Power iteration with small gap (можно оставить как есть)
# ------------------------------------------------------------
def test_power_iteration_small_gap():
    print("=== Power iteration with small gap ===")
    dim = 3
    eigvals = [1.0, 0.99, 0.98]
    alg = TypeIAlgebra(n_factors=[dim], k_factors=[dim], batch_size=1,
                       complex_valued=False, device='cpu')
    def sampler(d):
        return torch.tensor([eigvals], dtype=torch.float32)
    A = alg.operator_from_eigenvalues(lambda d: sampler(d)[0, :d], batch_size=1,
                                      force_self_adjoint=True)
    mat = A.matrix[0, 0].cpu().numpy()
    v = np.random.randn(dim)
    v = v / np.linalg.norm(v)
    tol = 1e-8
    lambda_old = 0.0
    for i in range(10000):
        v_new = mat @ v
        lambda_est = np.dot(v, v_new) / np.dot(v, v)
        v_new = v_new / np.linalg.norm(v_new)
        if abs(lambda_est - lambda_old) < tol:
            print(f"Spectral gap = {eigvals[0]-eigvals[1]:.4f}, iterations to converge: {i+1}")
            break
        lambda_old = lambda_est
        v = v_new
    else:
        print("Did not converge within 10000 iterations")
    print()

# ------------------------------------------------------------
# 3. SVD‑based sqrt с двойной точностью
# ------------------------------------------------------------
def test_svd_sqrt():
    print("=== SVD sqrt accuracy (float64) ===")
    dim = 100
    num_matrices = 1000
    batch_size = num_matrices
    alg = TypeIAlgebra(n_factors=[dim], k_factors=[dim], batch_size=batch_size,
                       complex_valued=False, device='cuda')
    def sampler(d):
        # генерируем float32, как требует библиотека
        return torch.rand(batch_size, d, device='cuda', dtype=torch.float32) + 1e-8
    ops = alg.operator_from_eigenvalues(sampler, batch_size=batch_size, force_positive=True)
    mats = ops.matrix[:, 0].double()  # переводим в float64 для высокой точности
    errors = []
    for i in range(num_matrices):
        A = mats[i]
        U, S, V = torch.linalg.svd(A)
        sqrtA = U @ torch.diag_embed(torch.sqrt(S)) @ U.T
        diff = sqrtA @ sqrtA - A
        rel_err = torch.norm(diff).item() / (torch.norm(A).item() + 1e-12)
        errors.append(rel_err)
    print(f"Number of matrices tested: {len(errors)}")
    print(f"Mean relative error: {np.mean(errors):.2e}")
    print(f"Max relative error:  {np.max(errors):.2e}")

if __name__ == '__main__':
    test_haar()
    test_power_iteration_small_gap()
    test_svd_sqrt()