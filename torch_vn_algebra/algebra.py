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
from .hilbert_space import HilbertSpace


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
        # Blocks up to this size are diagonalised exactly (torch.linalg.eigvalsh) in
        # lambda_max / lambda_min; larger ones fall back to shifted power iteration.
        self.exact_eig_max_dim = 256

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

    UNITARY_MEASURES = ('haar', 'haar_su', 'coe', 'cse', 'diag')

    def random_unitary(self, n: int, measure: str = 'haar',
                       batch_size: Optional[int] = None) -> torch.Tensor:
        """
        Random n x n unitary (orthogonal if the algebra is real).

        measure:
            'haar'    - Haar measure on U(n) (complex) or O(n) (real)
            'haar_su' - Haar measure on SU(n) (complex) or SO(n) (real)
            'coe'     - Circular Orthogonal Ensemble, S = W^T W (complex only)
            'cse'     - Circular Symplectic Ensemble, S = W^R W with
                        W^R = J W^T J^T (complex only, n even)
            'diag'    - diagonal matrix of i.i.d. uniform phases (complex only)

        Returns a tensor of shape (n, n), or (batch_size, n, n) if batch_size is given.
        """
        shape = () if batch_size is None else (batch_size,)
        device = self.hilbert.device
        complex_ = self.hilbert.complex_valued
        if measure in ('coe', 'cse', 'diag') and not complex_:
            raise ValueError(f"measure '{measure}' requires a complex-valued algebra")

        if measure == 'haar':
            # Mezzadri's recipe: QR of a Ginibre matrix with the phases of diag(R) fixed.
            # torch.randn with a complex dtype already draws circular N(0, 1) entries.
            dtype = torch.complex64 if complex_ else torch.float32
            Z = torch.randn(*shape, n, n, dtype=dtype, device=device)
            Q, R = torch.linalg.qr(Z)
            d = torch.diagonal(R, dim1=-2, dim2=-1)
            d = d.sgn() if complex_ else d.sign()
            return Q * d.unsqueeze(-2)
        elif measure == 'haar_su':
            Q = self.random_unitary(n, 'haar', batch_size)
            det = torch.linalg.det(Q)
            if complex_:
                # Multiply by det^{-1/n}: a Haar U(n) matrix times a phase is Haar on SU(n).
                phase = torch.exp(-1j * torch.angle(det) / n).to(Q.dtype)
                return Q * phase[..., None, None]
            # Real case: flip the first column when det = -1 to land in SO(n).
            Q = Q.clone()
            Q[..., :, 0] = Q[..., :, 0] * det.sign().unsqueeze(-1)
            return Q
        elif measure == 'coe':
            W = self.random_unitary(n, 'haar', batch_size)
            return W.transpose(-2, -1) @ W
        elif measure == 'cse':
            if n % 2 != 0:
                raise ValueError("CSE requires even dimension")
            W = self.random_unitary(n, 'haar', batch_size)
            half = n // 2
            I = torch.eye(half, dtype=W.dtype, device=device)
            O = torch.zeros_like(I)
            J = torch.cat([torch.cat([O, I], dim=-1), torch.cat([-I, O], dim=-1)], dim=-2)
            W_dual = J @ W.transpose(-2, -1) @ J.transpose(-2, -1)
            return W_dual @ W
        elif measure == 'diag':
            phases = torch.exp(2j * torch.pi * torch.rand(*shape, n, device=device))
            return torch.diag_embed(phases.to(torch.complex64))
        else:
            raise ValueError(f"Unknown measure: {measure}. Choose from {self.UNITARY_MEASURES}")

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

        def _tol(self) -> float:
            """Absolute tolerance for property checks, scaled to the precision and size of the entries."""
            mat = self._matrix
            eps = 1e-5 if mat.dtype in (torch.float32, torch.complex64) else 1e-10
            scale = mat.abs().max().item() if mat.numel() else 0.0
            return eps * max(1.0, scale)

        def _max_abs(self, t: torch.Tensor) -> float:
            return t.abs().max().item() if t.numel() else 0.0

        def _validate_properties(self):
            mat = self._matrix
            tol = self._tol()
            if self._is_self_adjoint_set and self._is_self_adjoint:
                if self._max_abs(mat - mat.conj().transpose(-2, -1)) > tol:
                    raise ValueError("Operator marked self-adjoint but not Hermitian")
            if self._is_normal_set and self._is_normal:
                mat_d = mat.conj().transpose(-2, -1)
                if self._max_abs(mat @ mat_d - mat_d @ mat) > tol * max(1.0, self._max_abs(mat)):
                    raise ValueError("Operator marked normal but not normal")
            if self._is_positive_set and self._is_positive:
                if torch.any(self.lambda_min < -tol).item():
                    raise ValueError("Operator marked positive but has negative eigenvalues")
            if self._is_invertible_set and self._is_invertible:
                if torch.any(self._block_svdvals_min() <= tol).item():
                    raise ValueError("Operator marked invertible but has zero eigenvalue")
            if self._is_projection_set and self._is_projection:
                if self._max_abs(mat @ mat - mat) > tol:
                    raise ValueError("Operator marked projection but P^2 != P")

        def _block_svdvals_min(self) -> torch.Tensor:
            """Smallest singular value over the active blocks, shape (batch,)."""
            mat = self.matrix
            mins = []
            for c, k_c in enumerate(self.algebra.k_factors):
                if k_c > 0:
                    mins.append(torch.linalg.svdvals(mat[:, c, :k_c, :k_c]).min(dim=-1)[0])
            return torch.stack(mins, dim=-1).min(dim=-1)[0]

        @property
        def shape(self) -> Tuple[int, ...]:
            return self.matrix.shape

        # ----- type checks -----
        @property
        def is_self_adjoint(self) -> bool:
            if self._is_self_adjoint is not None:
                return self._is_self_adjoint
            mat = self.matrix
            self._is_self_adjoint = self._max_abs(mat - mat.conj().transpose(-2, -1)) <= self._tol()
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
            diff = mat @ mat_d - mat_d @ mat
            self._is_normal = self._max_abs(diff) <= self._tol() * max(1.0, self._max_abs(mat))
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
            self.matrix
            self._is_positive = torch.all(self.lambda_min >= -self._tol()).item()
            return self._is_positive

        @is_positive.setter
        def is_positive(self, value: bool):
            self._is_positive = value
            self._is_positive_set = True

        @property
        def is_invertible(self) -> bool:
            if self._is_invertible is not None:
                return self._is_invertible
            self.matrix
            self._is_invertible = torch.all(self._block_svdvals_min() > self._tol()).item()
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
            self._is_projection = self._max_abs(mat @ mat - mat) <= self._tol()
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
            if max_k <= self.algebra.exact_eig_max_dim:
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
        def michelson_contrast(self) -> torch.Tensor:
            """
            Michelson contrast Delta(A) = (lambda_max - lambda_min) / (lambda_max + lambda_min)
            of a positive operator, with the extreme eigenvalues taken over the whole
            algebra (all channels). Returns a tensor of shape (batch,); Delta = 0 for A = 0.
            """
            if not self.is_positive:
                raise RuntimeError("Michelson contrast requires a positive operator.")
            lmax, lmin = self.lambda_max, self.lambda_min
            denom = lmax + lmin
            safe = torch.where(denom > 1e-12, denom, torch.ones_like(denom))
            return torch.where(denom > 1e-12, (lmax - lmin) / safe, torch.zeros_like(denom))

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
            if not self._is_materialized:   # do not trigger materialisation
                return f"Operator(C={self.algebra.C}, k_max={self.algebra.k_max}, lazy)"
            mem = self._bytes_held / 1024
            return f"Operator({tuple(self._matrix.shape)}, materialized, {mem:.1f}KB)"

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
                if k_c == 0:
                    batch_matrices.append(torch.zeros(batch_size, self.k_max, self.k_max,
                                                      dtype=self.hilbert.dtype, device=self.hilbert.device))
                    continue
                U_batch = self.random_unitary(k_c, measure=unitary_measure, batch_size=batch_size)
                diag_eig = torch.diag_embed(eig)      # (batch, k_c, k_c)
                A_c = U_batch @ diag_eig @ U_batch.conj().transpose(-2, -1)  # (batch, k_c, k_c)
                if not torch.is_complex(eig) or not torch.any(eig.imag != 0):
                    # real spectrum: remove round-off so that A_c is exactly Hermitian
                    A_c = 0.5 * (A_c + A_c.conj().transpose(-2, -1))
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

    def random_unitary_operator(self, batch_size: int = 1, measure: str = 'haar') -> Operator:
        """
        Block-diagonal unitary U = (+)_c U_c in the algebra, with independent
        U_c drawn from `measure` (see random_unitary) for every channel and sample.
        """
        def generator():
            mat = torch.zeros(batch_size, self.C, self.k_max, self.k_max,
                              dtype=self.hilbert.dtype, device=self.hilbert.device)
            for c, k_c in enumerate(self.k_factors):
                if k_c > 0:
                    mat[:, c, :k_c, :k_c] = self.random_unitary(k_c, measure=measure, batch_size=batch_size)
            return mat
        op = self.Operator(self, generator=generator)
        op._is_normal = True
        op._is_invertible = True
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
        op._is_invertible = True
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
        op._is_invertible = all(s != 0 for s, k_c in zip(scalars, self.k_factors) if k_c > 0)
        op._is_projection = all(s == 0 or s == 1 for s in scalars)
        return op

    def __repr__(self) -> str:
        return f"TypeIAlgebra(C={self.C}, n={self.n_factors}, k={self.k_factors})"