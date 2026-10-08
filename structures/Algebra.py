"""
Type I von Neumann algebra M = ⊕_{c=1}^{C} M_{n_c}(ℂ) with subspace dimensions k_c.

Storage shape: (batch, C, k_max, k_max)
- batch: parallel Monte Carlo samples
- C: number of factors (channels)
- k_max = max(k_c): maximum subspace dimension across channels

Each channel c has an active block of size k_c × k_c in the top-left corner.
The remaining (k_max - k_c) dimensions are padding.
"""

import torch
import numpy as np
from typing import Optional, List, Tuple, Callable, Union
from .Hilbert_space import HilbertSpace


class TypeIAlgebra:
    """
    Type I von Neumann algebra M = ⊕_{c=1}^{C} M_{n_c}(ℂ).
    
    Operator storage shape: (batch, C, k_max, k_max)
    """

    def __init__(
        self,
        n_factors: List[int],
        k_factors: List[int],
        hilbert: Optional[HilbertSpace] = None,
        batch_size: int = 1,
        complex_valued: bool = True,
        device: Optional[torch.device] = None
    ):
        C = len(n_factors)
        assert len(k_factors) == C, "n_factors and k_factors must have same length"

        for n, k in zip(n_factors, k_factors):
            assert k <= n, f"k={k} exceeds n={n}"
            assert n > 0 and k >= 0

        self.n_factors = n_factors
        self.k_factors = k_factors
        self.C = C
        self.k_max = max(k_factors) if k_factors else 0
        self.total_subspace_dim = sum(k_factors)

        if hilbert is None:
            self.hilbert = HilbertSpace(
                n=self.total_subspace_dim,
                k=self.total_subspace_dim,
                batch_size=batch_size,
                n_channels=C,
                complex_valued=complex_valued,
                device=device
            )
        else:
            self.hilbert = hilbert

    # ========================================================================
    # Random Unitary / Orthogonal / SU(n)
    # ========================================================================

    def random_unitary(self, n: int, measure: str = 'haar') -> torch.Tensor:
        device = self.hilbert.device
        complex_ = self.hilbert.complex_valued
        if measure == 'haar':
            if complex_:
                Z = torch.randn(n, n, dtype=torch.complex64, device=device)
                Z = Z + 1j * torch.randn(n, n, device=device)
            else:
                Z = torch.randn(n, n, dtype=torch.float32, device=device)
            Q, R = torch.linalg.qr(Z)
            # Исправлено: для комплексных используем .sgn(), для вещественных .sign()
            if complex_:
                d = torch.diag(R).sgn()
            else:
                d = torch.diag(R).sign()
            Q = Q * d.unsqueeze(0)
            return Q
        elif measure == 'haar_su':
            Q = self.random_unitary(n, measure='haar')
            det = torch.linalg.det(Q)
            Q[..., 0, :] = Q[..., 0, :] / det
            return Q
        # Остальные меры (coe, cse, diag) оставить без изменений
        elif measure == 'coe':
            Z = torch.randn(n, n, dtype=torch.float32, device=device)
            Z = (Z + Z.T) / 2
            Q, _ = torch.linalg.qr(Z)
            return Q
        elif measure == 'cse':
            assert n % 2 == 0, "CSE requires even dimension"
            half = n // 2
            A = torch.randn(half, half, dtype=torch.float32, device=device)
            B = torch.randn(half, half, dtype=torch.float32, device=device)
            Cmat = torch.randn(half, half, dtype=torch.float32, device=device)
            D = torch.randn(half, half, dtype=torch.float32, device=device)
            Z = torch.cat([torch.cat([A, B], dim=-1), torch.cat([Cmat, D], dim=-1)], dim=-2)
            Q, _ = torch.linalg.qr(Z)
            return Q
        elif measure == 'diag':
            if not complex_:
                raise ValueError("Diagonal measure requires complex-valued matrices")
            phases = torch.exp(2j * torch.pi * torch.rand(n, device=device))
            return torch.diag(phases)
        else:
            raise ValueError(f"Unknown measure: {measure}")

    # ========================================================================
    # Operator Class
    # ========================================================================

    class Operator:
        """
        Operator in the von Neumann algebra.
        
        Storage: (batch, C, k_max, k_max)
        """

        def __init__(
            self,
            algebra: 'TypeIAlgebra',
            generator: Optional[Callable[[], torch.Tensor]] = None,
            matrix: Optional[torch.Tensor] = None,
            is_self_adjoint: Optional[bool] = None,
            is_normal: Optional[bool] = None,
            is_positive: Optional[bool] = None,
            is_invertible: Optional[bool] = None,
            is_projection: Optional[bool] = None,
            eigenvalues: Optional[List[torch.Tensor]] = None
        ):
            self.algebra = algebra
            self._generator = generator
            self._matrix = matrix
            self._is_materialized = matrix is not None

            self._is_self_adjoint = is_self_adjoint
            self._is_self_adjoint_set = is_self_adjoint is not None
            self._is_normal = is_normal
            self._is_normal_set = is_normal is not None
            self._is_positive = is_positive
            self._is_positive_set = is_positive is not None
            self._is_invertible = is_invertible
            self._is_invertible_set = is_invertible is not None
            self._is_projection = is_projection
            self._is_projection_set = is_projection is not None
            self._eigenvalues = eigenvalues

            self._lambda_max = None
            self._lambda_min = None
            self._trace = None
            self._inverse = None
            self._abs = None

            self._update_memory()

        def _update_memory(self):
            if self._matrix is not None:
                self._bytes_held = self._matrix.element_size() * self._matrix.numel()
            else:
                self._bytes_held = 0

        @property
        def matrix(self) -> torch.Tensor:
            if not self._is_materialized:
                self._matrix = self._generator()
                self._is_materialized = True
                self._validate_properties()
                self._update_memory()
            return self._matrix

        def _validate_properties(self):
            mat = self._matrix
            if self._is_self_adjoint_set and self._is_self_adjoint:
                diff = mat - mat.conj().transpose(-2, -1)
                if torch.norm(diff) > 1e-6:
                    raise ValueError("Operator marked self-adjoint but not Hermitian")
            if self._is_normal_set and self._is_normal:
                mat_d = mat.conj().transpose(-2, -1)
                left = mat @ mat_d
                right = mat_d @ mat
                if torch.norm(left - right) > 1e-6:
                    raise ValueError("Operator marked normal but not normal")
            if self._is_positive_set and self._is_positive:
                lambda_min = self.lambda_min
                if torch.any(lambda_min < -1e-12).item():
                    raise ValueError("Operator marked positive but has negative eigenvalues")
            if self._is_invertible_set and self._is_invertible:
                lambda_min = self.lambda_min
                if torch.any(lambda_min < 1e-12).item():
                    raise ValueError("Operator marked invertible but has zero eigenvalue")
            if self._is_projection_set and self._is_projection:
                mat2 = mat @ mat
                if torch.norm(mat2 - mat) > 1e-6:
                    raise ValueError("Operator marked projection but P^2 != P")

        @property
        def shape(self) -> Tuple[int, ...]:
            return self.matrix.shape

        # ----- type checks -----
        @property
        def is_self_adjoint(self) -> bool:
            if self._is_self_adjoint is not None:
                return self._is_self_adjoint
            mat = self.matrix
            diff = mat - mat.conj().transpose(-2, -1)
            self._is_self_adjoint = torch.norm(diff) < 1e-6
            return self._is_self_adjoint

        @is_self_adjoint.setter
        def is_self_adjoint(self, value: bool):
            self._is_self_adjoint = value
            self._is_self_adjoint_set = True

        @property
        def is_normal(self) -> bool:
            if self._is_normal is not None:
                return self._is_normal
            mat = self.matrix
            mat_d = mat.conj().transpose(-2, -1)
            left = mat @ mat_d
            right = mat_d @ mat
            diff = left - right
            self._is_normal = torch.norm(diff) < 1e-6
            return self._is_normal

        @is_normal.setter
        def is_normal(self, value: bool):
            self._is_normal = value
            self._is_normal_set = True

        @property
        def is_positive(self) -> bool:
            if self._is_positive is not None:
                return self._is_positive
            if not self.is_self_adjoint:
                self._is_positive = False
                return False
            lambda_min = self.lambda_min
            self._is_positive = torch.all(lambda_min >= -1e-12).item()
            return self._is_positive

        @is_positive.setter
        def is_positive(self, value: bool):
            self._is_positive = value
            self._is_positive_set = True

        @property
        def is_invertible(self) -> bool:
            if self._is_invertible is not None:
                return self._is_invertible
            if not self.is_self_adjoint:
                self._is_invertible = False
                return False
            lambda_min = self.lambda_min
            self._is_invertible = torch.all(lambda_min > 1e-12).item()
            return self._is_invertible

        @is_invertible.setter
        def is_invertible(self, value: bool):
            self._is_invertible = value
            self._is_invertible_set = True

        @property
        def is_projection(self) -> bool:
            if self._is_projection is not None:
                return self._is_projection
            if not self.is_self_adjoint:
                self._is_projection = False
                return False
            mat = self.matrix
            mat2 = mat @ mat
            diff = mat2 - mat
            self._is_projection = (torch.norm(diff) < 1e-6).item()
            return self._is_projection

        @is_projection.setter
        def is_projection(self, value: bool):
            self._is_projection = value
            self._is_projection_set = True

        # ----- power iteration (legacy, kept for operator_norm) -----
        @staticmethod
        def _power_iteration(mat: torch.Tensor, num_iters: int = 100, tol: float = 1e-8) -> torch.Tensor:
            batch, C, k_max, _ = mat.shape
            device = mat.device
            dtype = mat.dtype
            v = torch.randn(batch, C, k_max, 1, dtype=torch.float32, device=device)
            if torch.is_complex(mat):
                v = v.to(torch.complex64)
                v = v + 1j * torch.randn(batch, C, k_max, 1, device=device)
            else:
                v = v.to(dtype)
            norm_sq = (v.conj() * v).real.sum(dim=(2, 3), keepdim=True)
            v = v / torch.sqrt(norm_sq + 1e-12)
            lambda_old = torch.zeros(batch, C, dtype=dtype, device=device)
            for _ in range(num_iters):
                v_new = torch.matmul(mat, v)
                vh_v_new = (v.conj() * v_new).real.sum(dim=(2, 3))
                vh_v = (v.conj() * v).real.sum(dim=(2, 3))
                lambda_est = vh_v_new / (vh_v + 1e-12)
                norm_sq = (v_new.conj() * v_new).real.sum(dim=(2, 3), keepdim=True)
                v_new = v_new / torch.sqrt(norm_sq + 1e-12)
                if torch.max(torch.abs(lambda_est - lambda_old)) < tol:
                    break
                lambda_old = lambda_est
                v = v_new
            return lambda_est

        # ----- helper for single‑block power iteration (used by lambda_max/lambda_min for large sizes) -----
        @staticmethod
        def _power_iteration_single(block: torch.Tensor, num_iters: int = 5000, tol: float = 1e-12) -> torch.Tensor:
            batch, dim, _ = block.shape
            device = block.device
            dtype = block.dtype
            # double precision for stability
            if dtype in (torch.float32, torch.complex64):
                block = block.to(torch.float64 if not torch.is_complex(block) else torch.complex128)
                dtype = block.dtype
            v = torch.randn(batch, dim, 1, dtype=torch.float64, device=device)
            if torch.is_complex(block):
                v = v.to(torch.complex128)
                v = v + 1j * torch.randn(batch, dim, 1, device=device)
            else:
                v = v.to(dtype)
            norm_sq = (v.conj() * v).real.sum(dim=(1, 2), keepdim=True)
            v = v / torch.sqrt(norm_sq + 1e-12)
            lambda_old = torch.zeros(batch, dtype=dtype, device=device)
            for _ in range(num_iters):
                v_new = torch.matmul(block, v)
                vh_v_new = (v.conj() * v_new).real.sum(dim=(1, 2))
                vh_v = (v.conj() * v).real.sum(dim=(1, 2))
                lambda_est = vh_v_new / (vh_v + 1e-12)
                norm_sq = (v_new.conj() * v_new).real.sum(dim=(1, 2), keepdim=True)
                v_new = v_new / torch.sqrt(norm_sq + 1e-12)
                if torch.max(torch.abs(lambda_est - lambda_old)) < tol:
                    break
                lambda_old = lambda_est
                v = v_new
            return lambda_est.to(block.dtype)

        # ----- lambda_max and lambda_min (global extremes) -----
        @property
        def lambda_max(self) -> torch.Tensor:
            if self._lambda_max is not None:
                return self._lambda_max
            if not self.is_self_adjoint:
                raise RuntimeError("lambda_max requires self-adjoint operator")
            A = self.matrix  # (batch, C, k_max, k_max)
            batch, C, _, _ = A.shape
            max_k = max(self.algebra.k_factors) if self.algebra.k_factors else 0
            if max_k <= 64:
                # Точная диагонализация для малых размеров
                all_eigvals = []
                for c in range(C):
                    k_c = self.algebra.k_factors[c]
                    if k_c == 0:
                        continue
                    block = A[:, c, :k_c, :k_c]
                    # Не приводим к real, работаем с комплексными эрмитовыми матрицами
                    eigvals = torch.linalg.eigvalsh(block)   # (batch, k_c)
                    all_eigvals.append(eigvals)
                all_eigvals = torch.cat(all_eigvals, dim=-1)   # (batch, total_dim)
                self._lambda_max = all_eigvals.max(dim=1)[0]
                self._lambda_min = all_eigvals.min(dim=1)[0]
                return self._lambda_max
            else:
                # Степенной метод для больших размеров (ваш алгоритм с μ_global и ν_global)
                # 1. Вычисляем μ_c для каждого канала
                mu_c = torch.zeros(batch, C, dtype=A.real.dtype, device=A.device)
                for c in range(C):
                    k_c = self.algebra.k_factors[c]
                    if k_c == 0:
                        mu_c[:, c] = 0.0
                        continue
                    block = A[:, c, :k_c, :k_c]
                    # Не обрезаем мнимую часть
                    mu_c[:, c] = self._power_iteration_single(block)   # (batch,)
                # 2. Глобальное μ: максимальное по модулю (при равенстве предпочитаем положительное)
                abs_mu = torch.abs(mu_c)
                _, idx = torch.max(abs_mu, dim=1)
                mu_global = torch.gather(mu_c, 1, idx.unsqueeze(1)).squeeze(1)
                # 3. Для всех каналов строим B_c = μ_global * I_c - A_c и находим ν_c (модуль)
                nu_c = torch.zeros(batch, C, dtype=A.real.dtype, device=A.device)
                for c in range(C):
                    k_c = self.algebra.k_factors[c]
                    if k_c == 0:
                        nu_c[:, c] = 0.0
                        continue
                    block = A[:, c, :k_c, :k_c]
                    mu_exp = mu_global.unsqueeze(-1).unsqueeze(-1)   # (batch,1,1)
                    I = torch.eye(k_c, dtype=block.dtype, device=block.device).unsqueeze(0)
                    B = mu_exp * I - block
                    nu_c[:, c] = self._power_iteration_single(B).abs()
                # 4. Глобальное ν = максимум по каналам
                nu_global, _ = torch.max(nu_c, dim=1)
                # 5. Восстанавливаем λ_max и λ_min
                is_pos = mu_global > 0
                lmax = torch.where(is_pos, mu_global, mu_global + nu_global)
                lmin = torch.where(is_pos, mu_global - nu_global, mu_global)
                # Защита от отрицательных (для положительных операторов)
                lmin = torch.clamp(lmin, min=0.0)
                self._lambda_max = lmax
                self._lambda_min = lmin
                return lmax

        @property
        def lambda_min(self) -> torch.Tensor:
            if self._lambda_min is not None:
                return self._lambda_min
            _ = self.lambda_max
            return self._lambda_min

        # ----- other methods (trace, norms, etc.) -----
        @property
        def trace(self) -> torch.Tensor:
            if self._trace is not None:
                return self._trace
            diag = torch.diagonal(self.matrix, dim1=-2, dim2=-1)
            self._trace = torch.sum(diag, dim=(-2, -1))
            return self._trace

        @trace.setter
        def trace(self, value: torch.Tensor):
            self._trace = value

        def trace_norm(self) -> torch.Tensor:
            U, S, Vh = torch.linalg.svd(self.matrix)
            return torch.sum(S, dim=(-2, -1))

        def frobenius_norm(self) -> torch.Tensor:
            mat = self.matrix
            total_sq = 0.0
            for c, k_c in enumerate(self.algebra.k_factors):
                if k_c == 0:
                    continue
                block = mat[:, c, :k_c, :k_c]
                total_sq += torch.sum(block.conj() * block, dim=(-2, -1)).real
            return torch.sqrt(total_sq)

        def operator_norm(self) -> torch.Tensor:
            mat = self.matrix
            AHA = mat.conj().transpose(-2, -1) @ mat
            lambda_max_AHA = self._power_iteration(AHA)   # shape (batch, C)
            return torch.sqrt(torch.max(lambda_max_AHA, dim=1)[0])

        def inverse(self, tol: float = 1e-12) -> 'TypeIAlgebra.Operator':
            def inv_generator():
                U, S, Vh = torch.linalg.svd(self.matrix)
                S_inv = torch.where(S > tol, 1.0 / S, torch.zeros_like(S))
                S_inv = S_inv.to(dtype=U.dtype)
                Vh_conj = Vh.conj().transpose(-2, -1)
                U_conj = U.conj().transpose(-2, -1)
                return Vh_conj @ torch.diag_embed(S_inv) @ U_conj
            op = TypeIAlgebra.Operator(self.algebra, generator=inv_generator)
            op._is_self_adjoint = self.is_self_adjoint
            op._is_normal = self.is_normal
            op._is_positive = self.is_positive
            op._is_invertible = self.is_invertible
            return op

        @property
        def inv(self) -> 'TypeIAlgebra.Operator':
            if self._inverse is None:
                self._inverse = self.inverse()
            return self._inverse

        def abs(self) -> 'TypeIAlgebra.Operator':
            def abs_generator():
                U, S, Vh = torch.linalg.svd(self.matrix)
                S = S.to(dtype=U.dtype)
                Vh_conj = Vh.conj().transpose(-2, -1)
                return Vh_conj @ torch.diag_embed(S) @ Vh
            op = TypeIAlgebra.Operator(self.algebra, generator=abs_generator)
            op._is_self_adjoint = True
            op._is_positive = True
            op._is_normal = True
            return op

        def sqrt(self) -> 'TypeIAlgebra.Operator':
            if not self.is_positive:
                raise RuntimeError("sqrt requires positive operator")
            def sqrt_generator():
                U, S, Vh = torch.linalg.svd(self.matrix)
                S = S.to(dtype=U.dtype)
                Vh_conj = Vh.conj().transpose(-2, -1)
                return Vh_conj @ torch.diag_embed(torch.sqrt(S)) @ Vh
            op = TypeIAlgebra.Operator(self.algebra, generator=sqrt_generator)
            op._is_self_adjoint = True
            op._is_positive = True
            op._is_normal = True
            return op

        def trace_a_log_a(self, tol: float = 1e-12) -> torch.Tensor:
            if self._is_positive is not None and not self._is_positive:
                raise RuntimeError("trace_a_log_a requires positive operator")
            U, S, Vh = torch.linalg.svd(self.matrix)
            S_safe = torch.clamp(S, min=tol)
            S_log_S = S_safe * torch.log(S_safe)
            return torch.sum(S_log_S, dim=(-2, -1))

        def entropy(self) -> torch.Tensor:
            return -self.trace_a_log_a()

        @property
        def michelson_contrast(self) -> float:
            if not self.is_positive:
                raise RuntimeError("Michelson contrast requires a positive operator.")
            mat = self.matrix
            batch, C, _, _ = mat.shape
            device = mat.device
            dtype = mat.dtype
            contrasts = []
            for b in range(batch):
                for c in range(C):
                    k_c = self.algebra.k_factors[c]
                    if k_c == 0:
                        continue
                    block = mat[b, c, :k_c, :k_c]
                    block_4d = block.reshape(1, 1, k_c, k_c)
                    lambda_max = self._power_iteration(block_4d).item()
                    lambda_max = lambda_max.real if isinstance(lambda_max, complex) else lambda_max
                    I = torch.eye(k_c, dtype=dtype, device=device)
                    B_mat = lambda_max * I - block
                    B_4d = B_mat.reshape(1, 1, k_c, k_c)
                    lambda_max_B = self._power_iteration(B_4d).item()
                    lambda_max_B = lambda_max_B.real if isinstance(lambda_max_B, complex) else lambda_max_B
                    lambda_min = lambda_max - lambda_max_B
                    denom = lambda_max + lambda_min
                    if denom.real < 1e-12:
                        contrast = 0.0
                    else:
                        contrast = ((lambda_max - lambda_min) / denom).real
                    contrasts.append(contrast)
            if not contrasts:
                return 0.0
            return sum(contrasts) / len(contrasts)

        def __add__(self, other: 'TypeIAlgebra.Operator') -> 'TypeIAlgebra.Operator':
            assert self.algebra is other.algebra
            return TypeIAlgebra.Operator(self.algebra, generator=lambda: self.matrix + other.matrix)

        def __matmul__(self, other: 'TypeIAlgebra.Operator') -> 'TypeIAlgebra.Operator':
            assert self.algebra is other.algebra
            return TypeIAlgebra.Operator(self.algebra, generator=lambda: torch.matmul(self.matrix, other.matrix))

        def __mul__(self, scalar: Union[float, complex, torch.Tensor]) -> 'TypeIAlgebra.Operator':
            return TypeIAlgebra.Operator(self.algebra, generator=lambda: self.matrix * scalar)

        __rmul__ = __mul__

        def adjoint(self) -> 'TypeIAlgebra.Operator':
            return TypeIAlgebra.Operator(self.algebra, generator=lambda: self.matrix.conj().transpose(-2, -1))

        def Tr_blunt(self) -> torch.Tensor:
            return self.trace

        def Tr_norm(self, basis=None) -> torch.Tensor:
            mat = self.matrix
            total = torch.zeros(mat.shape[0], dtype=mat.real.dtype, device=mat.device)
            for c, k_c in enumerate(self.algebra.k_factors):
                if k_c > 0:
                    block = mat[:, c, :k_c, :k_c]
                    tr_c = torch.diagonal(block, dim1=-2, dim2=-1).sum(-1).real
                    total += tr_c / k_c
            return total

        def tau_vN(self) -> torch.Tensor:
            mat = self.matrix
            B = mat.shape[0]
            total = torch.zeros(B, dtype=mat.real.dtype, device=mat.device)
            for c, k_c in enumerate(self.algebra.k_factors):
                if k_c > 0:
                    block = mat[:, c, :k_c, :k_c]
                    tr_c = torch.diagonal(block, dim1=-2, dim2=-1).sum(-1).real
                    total = total + tr_c / k_c
            return total / self.algebra.C

        def __repr__(self) -> str:
            status = "mat" if self._is_materialized else "lazy"
            mem = self._bytes_held / 1024
            return f"Operator({self.shape}, {status}, {mem:.1f}KB)"

    # ========================================================================
    # Factory Methods
    # ========================================================================

    def operator_from_eigenvalues(
        self,
        eigenvalue_sampler: Callable[[int], torch.Tensor],
        batch_size: int = 1,
        unitary_measure: str = 'haar',
        force_self_adjoint: Optional[bool] = None,
        force_positive: Optional[bool] = None,
        force_normal: Optional[bool] = None,
        force_invertible: Optional[bool] = None,
        force_projection: Optional[bool] = None,
    ) -> Operator:
        all_eig = []
        inferred_self_adjoint = True
        inferred_positive = True
        inferred_normal = True
        inferred_invertible = True
        inferred_projection = True

        for c in range(self.C):
            k_c = self.k_factors[c]
            eig = eigenvalue_sampler(k_c)
            if isinstance(eig, np.ndarray):
                eig = torch.from_numpy(eig)
            eig = eig.to(device=self.hilbert.device, dtype=self.hilbert.dtype)
            all_eig.append(eig)

            if torch.is_complex(eig):
                if torch.any(eig.imag != 0):
                    inferred_self_adjoint = False
                    inferred_normal = False
                if torch.any(eig.real < -1e-12):
                    inferred_positive = False
                if torch.any(torch.abs(eig) < 1e-12):
                    inferred_invertible = False
                proj_ok = torch.all((torch.abs(eig.real) < 1e-12) | (torch.abs(eig.real - 1) < 1e-12))
                if not proj_ok:
                    inferred_projection = False
            else:
                if torch.any(eig < -1e-12):
                    inferred_positive = False
                if torch.any(torch.abs(eig) < 1e-12):
                    inferred_invertible = False
                proj_ok = torch.all((torch.abs(eig) < 1e-12) | (torch.abs(eig - 1) < 1e-12))
                if not proj_ok:
                    inferred_projection = False

        if force_self_adjoint is not None:
            if force_self_adjoint and not inferred_self_adjoint:
                for c, eig in enumerate(all_eig):
                    if torch.is_complex(eig) and torch.any(eig.imag != 0):
                        raise ValueError(f"force_self_adjoint=True but eigenvalues for channel {c} have non-zero imaginary part")
            inferred_self_adjoint = force_self_adjoint

        if force_positive is not None:
            if force_positive and not inferred_positive:
                for c, eig in enumerate(all_eig):
                    if torch.is_complex(eig):
                        if torch.any(eig.real < -1e-12):
                            bad_indices = torch.where(eig.real < -1e-12)[0].tolist()
                            raise ValueError(
                                f"force_positive=True but eigenvalues for channel {c} are negative (real part): "
                                f"indices {bad_indices}"
                            )
                    else:
                        if torch.any(eig < -1e-12):
                            bad_indices = torch.where(eig < -1e-12)[0].tolist()
                            raise ValueError(
                                f"force_positive=True but eigenvalues for channel {c} are negative: "
                                f"indices {bad_indices}"
                            )
            inferred_positive = force_positive

        if force_normal is not None:
            inferred_normal = force_normal

        if force_invertible is not None:
            if force_invertible and not inferred_invertible:
                for c, eig in enumerate(all_eig):
                    if torch.any(torch.abs(eig) < 1e-12):
                        raise ValueError(f"force_invertible=True but eigenvalues for channel {c} contain zeros")
            inferred_invertible = force_invertible

        if force_projection is not None:
            if force_projection and not inferred_projection:
                for c, eig in enumerate(all_eig):
                    if torch.is_complex(eig):
                        bad_mask = (torch.abs(eig.real) > 1e-12) & (torch.abs(eig.real - 1) > 1e-12)
                    else:
                        bad_mask = (torch.abs(eig) > 1e-12) & (torch.abs(eig - 1) > 1e-12)
                    if torch.any(bad_mask):
                        raise ValueError(f"force_projection=True but eigenvalues for channel {c} are not 0 or 1")
            inferred_projection = force_projection

        def generator():
            batch_matrices = []
            for c in range(self.C):
                k_c = self.k_factors[c]
                eig = all_eig[c]                     # (batch, k_c)
                # Generate batch of unitary matrices
                U_batch = torch.stack([self.random_unitary(k_c, measure=unitary_measure) for _ in range(batch_size)])
                diag_eig = torch.diag_embed(eig)      # (batch, k_c, k_c)
                A_c = U_batch @ diag_eig @ U_batch.conj().transpose(-2, -1)  # (batch, k_c, k_c)
                A_full = torch.zeros(batch_size, self.k_max, self.k_max, dtype=A_c.dtype, device=self.hilbert.device)
                A_full[:, :k_c, :k_c] = A_c
                batch_matrices.append(A_full)
            matrix = torch.stack(batch_matrices, dim=1)   # (batch, C, k_max, k_max)
            return matrix

        op = self.Operator(self, generator=generator)
        op._eigenvalues = all_eig
        op._is_self_adjoint = inferred_self_adjoint
        op._is_positive = inferred_positive
        op._is_normal = inferred_normal
        op._is_invertible = inferred_invertible
        op._is_projection = inferred_projection
        return op

    def identity(self, batch_size: int = 1) -> Operator:
        def generator():
            matrices = []
            for k_c in self.k_factors:
                mat = torch.zeros(self.k_max, self.k_max, dtype=self.hilbert.dtype, device=self.hilbert.device)
                mat[:k_c, :k_c] = torch.eye(k_c, dtype=self.hilbert.dtype, device=self.hilbert.device)
                matrices.append(mat)
            return torch.stack(matrices, dim=0).unsqueeze(0).expand(batch_size, -1, -1, -1)
        op = self.Operator(self, generator=generator)
        op._is_self_adjoint = True
        op._is_normal = True
        op._is_positive = True
        op._is_invertible = False
        op._is_projection = True
        return op

    def zero(self, batch_size: int = 1) -> Operator:
        shape = (batch_size, self.C, self.k_max, self.k_max)
        def generator():
            return torch.zeros(shape, dtype=self.hilbert.dtype, device=self.hilbert.device)
        op = self.Operator(self, generator=generator)
        op._is_self_adjoint = True
        op._is_normal = True
        op._is_positive = True
        op._is_invertible = False
        op._is_projection = True
        return op

    def central(self, scalars: List[float], batch_size: int = 1) -> Operator:
        assert len(scalars) == self.C
        def generator():
            matrices = []
            for c, (k_c, scalar) in enumerate(zip(self.k_factors, scalars)):
                mat = torch.zeros(self.k_max, self.k_max, dtype=self.hilbert.dtype, device=self.hilbert.device)
                mat[:k_c, :k_c] = scalar * torch.eye(k_c, dtype=self.hilbert.dtype, device=self.hilbert.device)
                matrices.append(mat)
            return torch.stack(matrices, dim=0).unsqueeze(0).expand(batch_size, -1, -1, -1)
        op = self.Operator(self, generator=generator)
        op._is_self_adjoint = True
        op._is_normal = True
        op._is_positive = all(s >= 0 for s in scalars)
        op._is_invertible = all(s != 0 for s in scalars) and all(k_c == self.k_max for k_c in self.k_factors)
        op._is_projection = all(s == 0 or s == 1 for s in scalars)
        return op

    def __repr__(self) -> str:
        return f"TypeIAlgebra(C={self.C}, n={self.n_factors}, k={self.k_factors})"